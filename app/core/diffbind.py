"""
DiffBind Analysis - Python Implementation

This module provides a DiffBind-compatible interface using pure Python.
No R dependencies required.

Uses:
- PyDESeq2 for differential analysis (equivalent to DESeq2)
- pysam/pybedtools for read counting
- scipy for statistics fallback
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
import pandas as pd

# Import from our Python differential module
from app.core.differential import (
    DifferentialAnalyzer,
    DifferentialConfig,
    DifferentialResults as DiffResults,
    Sample,
)

logger = logging.getLogger(__name__)


@dataclass
class DiffBindConfig:
    """Configuration for DiffBind-style analysis (Python implementation)."""
    comparison_name: str
    group1: str
    group2: str
    histone_mark: str
    genome: str = "hg38"

    # Analysis parameters
    fdr_threshold: float = 0.05
    lfc_threshold: float = 1.0
    min_overlap: int = 2
    summit_size: int = 250

    # Normalization: RLE (DESeq2-style), TMM, or median_ratio
    normalize_method: str = "RLE"

    # Gene annotation
    tss_upstream: int = 3000
    tss_downstream: int = 3000

    # Output options
    output_dir: Optional[str] = None


@dataclass
class DiffBindResults:
    """Results from DiffBind-style analysis."""
    comparison_name: str
    total_peaks: int
    significant_peaks: int
    gained_peaks: int
    lost_peaks: int

    # File paths
    all_peaks_file: Optional[str] = None
    significant_peaks_file: Optional[str] = None
    genes_file: Optional[str] = None

    # DataFrames
    all_peaks: Optional[pd.DataFrame] = None
    significant_df: Optional[pd.DataFrame] = None

    # Summary statistics
    summary: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if not self.summary:
            self.summary = {
                "comparison": self.comparison_name,
                "total_peaks": self.total_peaks,
                "significant_peaks": self.significant_peaks,
                "gained": self.gained_peaks,
                "lost": self.lost_peaks
            }


class DiffBindRunner:
    """
    Python-native DiffBind-style analysis runner.

    Provides equivalent functionality to R DiffBind without R dependencies.
    Uses PyDESeq2 for differential testing.
    """

    def __init__(self):
        """Initialize the DiffBind runner."""
        self.analyzer = DifferentialAnalyzer()

    def create_sample_sheet(
        self,
        samples: List[Dict],
        output_path: str
    ) -> str:
        """
        Create sample sheet CSV.

        Args:
            samples: List of sample dictionaries
            output_path: Path to save CSV

        Returns:
            Path to created sample sheet
        """
        required_cols = ["SampleID", "Condition", "Factor", "Replicate"]

        df = pd.DataFrame(samples)

        # Validate
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            raise ValueError(f"Sample sheet missing columns: {missing}")

        df.to_csv(output_path, index=False)
        logger.info(f"Created sample sheet: {output_path}")

        return output_path

    def run(
        self,
        samples: List[Dict],
        config: DiffBindConfig,
        output_dir: Optional[str] = None
    ) -> DiffBindResults:
        """
        Run differential peak analysis.

        Args:
            samples: List of sample dictionaries with keys:
                - SampleID: Unique identifier
                - Condition: Group label (Treatment/Control)
                - Factor: Histone mark name
                - Replicate: Replicate number
                - bamReads: Path to BAM file (optional)
                - Peaks: Path to peak file
            config: Analysis configuration
            output_dir: Output directory (overrides config.output_dir)

        Returns:
            DiffBindResults object
        """
        if output_dir:
            config.output_dir = output_dir

        logger.info(f"Running DiffBind analysis: {config.comparison_name}")
        logger.info(f"Comparing {config.group1} vs {config.group2}")

        # Convert to Sample objects
        sample_objects = []
        for s in samples:
            sample_objects.append(Sample(
                sample_id=s.get('SampleID', s.get('sample_id')),
                condition=s.get('Condition', s.get('condition')),
                histone_mark=s.get('Factor', s.get('histone_mark', config.histone_mark)),
                replicate=s.get('Replicate', s.get('replicate', 1)),
                bam_file=s.get('bamReads', s.get('bam_file')),
                peak_file=s.get('Peaks', s.get('peak_file'))
            ))

        # Create differential config
        diff_config = DifferentialConfig(
            comparison_name=config.comparison_name,
            group1=config.group1,
            group2=config.group2,
            min_overlap=config.min_overlap,
            summit_extend=config.summit_size,
            normalize_method=config.normalize_method,
            fdr_threshold=config.fdr_threshold,
            lfc_threshold=config.lfc_threshold,
            output_dir=config.output_dir
        )

        # Run analysis
        results = self.analyzer.run(sample_objects, diff_config)

        # Convert to DiffBindResults format
        out_dir = config.output_dir

        return DiffBindResults(
            comparison_name=results.comparison_name,
            total_peaks=results.total_peaks,
            significant_peaks=results.significant_peaks,
            gained_peaks=results.gained_peaks,
            lost_peaks=results.lost_peaks,
            all_peaks_file=f"{out_dir}/all_peaks.csv" if out_dir else None,
            significant_peaks_file=f"{out_dir}/significant_peaks.csv" if out_dir else None,
            all_peaks=results.all_peaks,
            significant_df=results.significant,
            summary=results.to_dict()
        )

    def run_from_samplesheet(
        self,
        samplesheet_path: str,
        config: DiffBindConfig
    ) -> DiffBindResults:
        """
        Run analysis from an existing sample sheet.

        Args:
            samplesheet_path: Path to sample sheet CSV
            config: Analysis configuration

        Returns:
            DiffBindResults object
        """
        samples = pd.read_csv(samplesheet_path).to_dict("records")
        return self.run(samples, config)


# Convenience function for direct usage
def run_diffbind(
    samples: List[Dict],
    comparison_name: str,
    group1: str,
    group2: str,
    histone_mark: str = "H3K27ac",
    output_dir: Optional[str] = None,
    **kwargs
) -> DiffBindResults:
    """
    Run DiffBind-style differential peak analysis.

    This is a pure Python implementation - no R required.

    Args:
        samples: List of sample dictionaries
        comparison_name: Name for the comparison
        group1: Treatment/test group name
        group2: Control/reference group name
        histone_mark: Histone mark being analyzed
        output_dir: Directory for output files
        **kwargs: Additional config options (fdr_threshold, lfc_threshold, etc.)

    Returns:
        DiffBindResults object with differential peaks
    """
    config = DiffBindConfig(
        comparison_name=comparison_name,
        group1=group1,
        group2=group2,
        histone_mark=histone_mark,
        output_dir=output_dir,
        **kwargs
    )

    runner = DiffBindRunner()
    return runner.run(samples, config)
