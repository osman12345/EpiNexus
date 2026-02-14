"""
Super-Enhancer Detection Module

Implements ROSE (Rank Ordering of Super-Enhancers) algorithm:
1. Stitch nearby enhancer peaks
2. Rank by signal intensity
3. Identify inflection point
4. Classify super-enhancers vs typical enhancers

Reference: Lovén et al., Cell 2013
"""

import numpy as np
import pandas as pd
from typing import Tuple, Dict, Any
from dataclasses import dataclass


@dataclass
class SuperEnhancerResult:
    """Result of super-enhancer analysis."""

    all_enhancers: pd.DataFrame
    super_enhancers: pd.DataFrame
    typical_enhancers: pd.DataFrame
    cutoff_value: float
    cutoff_rank: int
    statistics: Dict[str, Any]


class SuperEnhancerDetector:
    """
    Detect super-enhancers using the ROSE algorithm.

    Super-enhancers are large clusters of enhancers that:
    - Drive cell identity gene expression
    - Are enriched for disease-associated variants
    - Show exceptionally high levels of H3K27ac, BRD4, Mediator
    """

    def __init__(self, stitch_distance: int = 12500, tss_exclusion: int = 2500, min_enhancer_size: int = 500):
        """
        Initialize detector.

        Args:
            stitch_distance: Max distance to stitch enhancers (default 12.5kb)
            tss_exclusion: Distance from TSS to exclude (promoters)
            min_enhancer_size: Minimum enhancer size in bp
        """
        self.stitch_distance = stitch_distance
        self.tss_exclusion = tss_exclusion
        self.min_enhancer_size = min_enhancer_size

    def detect(
        self, peaks: pd.DataFrame, signal_col: str = "signal", tss_regions: pd.DataFrame = None
    ) -> SuperEnhancerResult:
        """
        Run super-enhancer detection.

        Args:
            peaks: DataFrame with chr, start, end, and signal columns
            signal_col: Column name for signal intensity
            tss_regions: Optional TSS regions to exclude (promoters)

        Returns:
            SuperEnhancerResult with classified enhancers
        """

        # Step 1: Filter and prepare peaks
        peaks = self._prepare_peaks(peaks, signal_col)

        # Step 2: Exclude TSS-proximal regions (optional)
        if tss_regions is not None:
            peaks = self._exclude_tss(peaks, tss_regions)

        # Step 3: Stitch nearby enhancers
        stitched = self._stitch_enhancers(peaks, signal_col)

        # Step 4: Rank by signal
        stitched = self._rank_enhancers(stitched, signal_col)

        # Step 5: Find inflection point (hockey stick)
        cutoff_value, cutoff_rank = self._find_cutoff(stitched, signal_col)

        # Step 6: Classify enhancers
        stitched["is_super_enhancer"] = stitched[signal_col] >= cutoff_value
        stitched["enhancer_class"] = stitched["is_super_enhancer"].map(
            {True: "Super-enhancer", False: "Typical enhancer"}
        )

        # Split results
        super_enhancers = stitched[stitched["is_super_enhancer"]].copy()
        typical_enhancers = stitched[~stitched["is_super_enhancer"]].copy()

        # Calculate statistics
        stats = self._calculate_statistics(stitched, super_enhancers, typical_enhancers, signal_col)

        return SuperEnhancerResult(
            all_enhancers=stitched,
            super_enhancers=super_enhancers,
            typical_enhancers=typical_enhancers,
            cutoff_value=cutoff_value,
            cutoff_rank=cutoff_rank,
            statistics=stats,
        )

    def _prepare_peaks(self, peaks: pd.DataFrame, signal_col: str) -> pd.DataFrame:
        """Prepare and validate peak data."""

        required_cols = ["chr", "start", "end", signal_col]
        for col in required_cols:
            if col not in peaks.columns:
                # Try common alternatives
                if col == "chr" and "chrom" in peaks.columns:
                    peaks = peaks.rename(columns={"chrom": "chr"})
                elif col == signal_col and "score" in peaks.columns:
                    peaks = peaks.rename(columns={"score": signal_col})
                else:
                    raise ValueError(f"Missing required column: {col}")

        # Filter by size
        peaks = peaks.copy()
        peaks["size"] = peaks["end"] - peaks["start"]
        peaks = peaks[peaks["size"] >= self.min_enhancer_size]

        # Sort by chromosome and position
        peaks = peaks.sort_values(["chr", "start"]).reset_index(drop=True)

        return peaks

    def _exclude_tss(self, peaks: pd.DataFrame, tss: pd.DataFrame) -> pd.DataFrame:
        """Exclude peaks near TSS (promoter regions).

        Uses interval-tree overlap detection (O(n log n)) for performance.
        """
        from .genomic_utils import exclude_overlapping

        # Expand TSS regions
        tss = tss.copy()
        tss["start"] = tss["start"] - self.tss_exclusion
        tss["end"] = tss["end"] + self.tss_exclusion

        return exclude_overlapping(
            peaks,
            tss,
            chrom_col="chr",
            start_col="start",
            end_col="end",
        )

    def _stitch_enhancers(self, peaks: pd.DataFrame, signal_col: str) -> pd.DataFrame:
        """Stitch nearby enhancers within stitch_distance."""

        if len(peaks) == 0:
            return peaks

        stitched_regions = []
        current_region = None

        for _, peak in peaks.iterrows():
            if current_region is None:
                current_region = {
                    "chr": peak["chr"],
                    "start": peak["start"],
                    "end": peak["end"],
                    "signal": peak[signal_col],
                    "constituent_count": 1,
                    "constituent_peaks": [peak.name],
                }
            elif peak["chr"] == current_region["chr"] and peak["start"] - current_region["end"] <= self.stitch_distance:
                # Stitch: extend region and sum signal
                current_region["end"] = max(current_region["end"], peak["end"])
                current_region["signal"] += peak[signal_col]
                current_region["constituent_count"] += 1
                current_region["constituent_peaks"].append(peak.name)
            else:
                # Save current region and start new one
                stitched_regions.append(current_region)
                current_region = {
                    "chr": peak["chr"],
                    "start": peak["start"],
                    "end": peak["end"],
                    "signal": peak[signal_col],
                    "constituent_count": 1,
                    "constituent_peaks": [peak.name],
                }

        # Don't forget last region
        if current_region is not None:
            stitched_regions.append(current_region)

        # Convert to DataFrame
        stitched_df = pd.DataFrame(stitched_regions)
        stitched_df["size"] = stitched_df["end"] - stitched_df["start"]

        # Rename signal column back
        stitched_df = stitched_df.rename(columns={"signal": signal_col})

        return stitched_df

    def _rank_enhancers(self, enhancers: pd.DataFrame, signal_col: str) -> pd.DataFrame:
        """Rank enhancers by signal intensity."""

        enhancers = enhancers.sort_values(signal_col, ascending=True).reset_index(drop=True)
        enhancers["rank"] = range(1, len(enhancers) + 1)
        enhancers["rank_normalized"] = enhancers["rank"] / len(enhancers)

        return enhancers

    def _find_cutoff(self, enhancers: pd.DataFrame, signal_col: str) -> Tuple[float, int]:
        """
        Find the inflection point using the tangent line method.

        Draws a line from the origin to the max point, finds the enhancer
        with maximum distance from this line.
        """

        if len(enhancers) < 10:
            # Not enough data for meaningful cutoff
            return enhancers[signal_col].median(), len(enhancers) // 2

        # Normalize data for geometric calculation
        x = enhancers["rank"].values
        y = enhancers[signal_col].values

        # Normalize to 0-1 range
        x_norm = (x - x.min()) / (x.max() - x.min())
        y_norm = (y - y.min()) / (y.max() - y.min() + 1e-10)

        # Line from origin (0,0) to max point (1, 1)
        # For normalized data, the line is y = x

        # Calculate perpendicular distance from each point to line y=x
        # Distance = |y - x| / sqrt(2)
        distances = np.abs(y_norm - x_norm) / np.sqrt(2)

        # Find point with maximum distance
        cutoff_idx = np.argmax(distances)
        cutoff_value = enhancers.iloc[cutoff_idx][signal_col]
        cutoff_rank = enhancers.iloc[cutoff_idx]["rank"]

        return cutoff_value, int(cutoff_rank)

    def _calculate_statistics(
        self, all_enh: pd.DataFrame, super_enh: pd.DataFrame, typical_enh: pd.DataFrame, signal_col: str
    ) -> Dict[str, Any]:
        """Calculate summary statistics."""

        return {
            "total_enhancers": len(all_enh),
            "super_enhancer_count": len(super_enh),
            "typical_enhancer_count": len(typical_enh),
            "super_enhancer_fraction": len(super_enh) / len(all_enh) if len(all_enh) > 0 else 0,
            "mean_se_size": super_enh["size"].mean() if len(super_enh) > 0 else 0,
            "mean_te_size": typical_enh["size"].mean() if len(typical_enh) > 0 else 0,
            "mean_se_signal": super_enh[signal_col].mean() if len(super_enh) > 0 else 0,
            "mean_te_signal": typical_enh[signal_col].mean() if len(typical_enh) > 0 else 0,
            "se_signal_fraction": super_enh[signal_col].sum() / all_enh[signal_col].sum() if len(all_enh) > 0 else 0,
            "median_constituents_se": super_enh["constituent_count"].median() if len(super_enh) > 0 else 0,
            "median_constituents_te": typical_enh["constituent_count"].median() if len(typical_enh) > 0 else 0,
        }


def generate_hockey_stick_data(result: SuperEnhancerResult, signal_col: str = "signal") -> pd.DataFrame:
    """Generate data for hockey stick plot visualization."""

    plot_data = result.all_enhancers[["rank", signal_col, "is_super_enhancer"]].copy()
    plot_data["cutoff"] = result.cutoff_value

    return plot_data
