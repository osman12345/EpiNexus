"""
Quality Control Module for Histone Mark Analysis.

Provides QC metrics for CUT&Tag and ChIP-seq data:
- Mapping statistics
- Fragment size distribution
- FRiP (Fraction of Reads in Peaks)
- Library complexity
- Spike-in calibration
"""

import subprocess
import logging
from pathlib import Path
from typing import Dict, Optional, List, Tuple
from dataclasses import dataclass
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class QCMetrics:
    """QC metrics for a sample."""

    sample_name: str

    # Mapping metrics
    total_reads: int = 0
    mapped_reads: int = 0
    mapping_rate: float = 0.0

    # Duplication
    unique_reads: int = 0
    duplicate_reads: int = 0
    duplicate_rate: float = 0.0

    # Peaks
    peak_count: int = 0
    reads_in_peaks: int = 0
    frip_score: float = 0.0

    # Fragment size
    median_fragment_size: int = 0
    fragment_size_std: float = 0.0

    # Spike-in (if applicable)
    spikein_reads: Optional[int] = None
    scale_factor: Optional[float] = None

    # Quality flags
    pass_qc: bool = True
    qc_warnings: List[str] = None

    def __post_init__(self):
        if self.qc_warnings is None:
            self.qc_warnings = []

    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            "sample_name": self.sample_name,
            "total_reads": self.total_reads,
            "mapped_reads": self.mapped_reads,
            "mapping_rate": self.mapping_rate,
            "unique_reads": self.unique_reads,
            "duplicate_rate": self.duplicate_rate,
            "peak_count": self.peak_count,
            "frip_score": self.frip_score,
            "median_fragment_size": self.median_fragment_size,
            "spikein_reads": self.spikein_reads,
            "scale_factor": self.scale_factor,
            "pass_qc": self.pass_qc,
            "qc_warnings": self.qc_warnings,
        }


class QCAnalyzer:
    """Quality control analyzer for histone mark data."""

    # QC thresholds
    MIN_MAPPING_RATE = 0.70
    MAX_DUPLICATE_RATE = 0.30
    MIN_FRIP = 0.10
    MIN_PEAK_COUNT = 1000
    OPTIMAL_FRAGMENT_SIZE = (150, 300)

    def __init__(self, samtools_path: str = "samtools", bedtools_path: str = "bedtools"):
        """
        Initialize QC analyzer.

        Args:
            samtools_path: Path to samtools executable
            bedtools_path: Path to bedtools executable
        """
        self.samtools = samtools_path
        self.bedtools = bedtools_path

    def analyze_bam(self, bam_path: str, sample_name: str) -> QCMetrics:
        """
        Analyze BAM file for mapping statistics.

        Args:
            bam_path: Path to BAM file
            sample_name: Sample identifier

        Returns:
            QCMetrics with mapping statistics
        """
        metrics = QCMetrics(sample_name=sample_name)

        try:
            # Get flagstat
            result = subprocess.run([self.samtools, "flagstat", bam_path], capture_output=True, text=True, check=True)

            # Parse flagstat output
            for line in result.stdout.split("\n"):
                if "in total" in line:
                    metrics.total_reads = int(line.split()[0])
                elif "mapped (" in line and "primary" not in line:
                    metrics.mapped_reads = int(line.split()[0])
                elif "duplicates" in line:
                    metrics.duplicate_reads = int(line.split()[0])

            # Calculate rates
            if metrics.total_reads > 0:
                metrics.mapping_rate = metrics.mapped_reads / metrics.total_reads
                metrics.unique_reads = metrics.mapped_reads - metrics.duplicate_reads
                if metrics.mapped_reads > 0:
                    metrics.duplicate_rate = metrics.duplicate_reads / metrics.mapped_reads

        except subprocess.CalledProcessError as e:
            logger.error(f"Error running samtools flagstat: {e}")
            metrics.qc_warnings.append("Failed to get mapping statistics")

        return metrics

    def calculate_fragment_sizes(self, bam_path: str, sample_size: int = 100000) -> Tuple[int, float, List[int]]:
        """
        Calculate fragment size distribution from BAM file.

        Args:
            bam_path: Path to BAM file
            sample_size: Number of reads to sample

        Returns:
            Tuple of (median, std, distribution)
        """
        try:
            # Use samtools to extract fragment sizes
            result = subprocess.run(
                [self.samtools, "view", "-f", "66", bam_path],  # Properly paired, first in pair
                capture_output=True,
                text=True,
                check=True,
            )

            sizes = []
            for i, line in enumerate(result.stdout.split("\n")):
                if i >= sample_size:
                    break
                if line:
                    fields = line.split("\t")
                    if len(fields) > 8:
                        tlen = abs(int(fields[8]))
                        if 0 < tlen < 700:  # Filter reasonable sizes (CUT&Tag typically <700bp)
                            sizes.append(tlen)

            if sizes:
                median = int(np.median(sizes))
                std = float(np.std(sizes))
                return median, std, sizes

        except subprocess.CalledProcessError as e:
            logger.error(f"Error calculating fragment sizes: {e}")

        return 0, 0.0, []

    def calculate_frip(self, bam_path: str, peak_path: str) -> Tuple[int, float]:
        """
        Calculate Fraction of Reads in Peaks (FRiP).

        Uses ``samtools view -c -L`` for accurate read counting rather than
        line-counting bedtools intersect output, which can miscount with
        multi-overlap scenarios.

        Args:
            bam_path: Path to BAM file
            peak_path: Path to peak BED file

        Returns:
            Tuple of (reads_in_peaks, frip_score)
        """
        try:
            # Count reads overlapping peaks using samtools view -c -L
            # This accurately counts reads (not lines) that overlap the BED regions
            result = subprocess.run(
                [self.samtools, "view", "-c", "-F", "4", "-L", peak_path, bam_path],
                capture_output=True,
                text=True,
                check=True,
            )
            reads_in_peaks = int(result.stdout.strip()) if result.stdout.strip() else 0

            # Get total mapped reads
            total_result = subprocess.run(
                [self.samtools, "view", "-c", "-F", "4", bam_path], capture_output=True, text=True, check=True
            )
            total_mapped = int(total_result.stdout.strip())

            frip = reads_in_peaks / total_mapped if total_mapped > 0 else 0.0

            return reads_in_peaks, frip

        except subprocess.CalledProcessError as e:
            logger.error(f"Error calculating FRiP: {e}")
            return 0, 0.0

    def count_peaks(self, peak_path: str) -> int:
        """Count number of peaks in BED file."""
        try:
            with open(peak_path) as f:
                return sum(1 for line in f if line.strip() and not line.startswith("#"))
        except Exception as e:
            logger.error(f"Error counting peaks: {e}")
            return 0

    def run_full_qc(
        self,
        sample_name: str,
        bam_path: Optional[str] = None,
        peak_path: Optional[str] = None,
        spikein_bam: Optional[str] = None,
    ) -> QCMetrics:
        """
        Run full QC analysis on a sample.

        Args:
            sample_name: Sample identifier
            bam_path: Path to BAM file
            peak_path: Path to peak file
            spikein_bam: Path to spike-in aligned BAM (optional)

        Returns:
            Complete QCMetrics
        """
        metrics = QCMetrics(sample_name=sample_name)

        # BAM QC
        if bam_path and Path(bam_path).exists():
            bam_metrics = self.analyze_bam(bam_path, sample_name)
            metrics.total_reads = bam_metrics.total_reads
            metrics.mapped_reads = bam_metrics.mapped_reads
            metrics.mapping_rate = bam_metrics.mapping_rate
            metrics.duplicate_reads = bam_metrics.duplicate_reads
            metrics.duplicate_rate = bam_metrics.duplicate_rate
            metrics.unique_reads = bam_metrics.unique_reads

            # Fragment sizes
            median, std, _ = self.calculate_fragment_sizes(bam_path)
            metrics.median_fragment_size = median
            metrics.fragment_size_std = std

        # Peak QC
        if peak_path and Path(peak_path).exists():
            metrics.peak_count = self.count_peaks(peak_path)

            if bam_path and Path(bam_path).exists():
                reads_in_peaks, frip = self.calculate_frip(bam_path, peak_path)
                metrics.reads_in_peaks = reads_in_peaks
                metrics.frip_score = frip

        # Spike-in
        if spikein_bam and Path(spikein_bam).exists():
            spikein_metrics = self.analyze_bam(spikein_bam, sample_name)
            metrics.spikein_reads = spikein_metrics.mapped_reads
            if metrics.spikein_reads and metrics.spikein_reads > 0:
                metrics.scale_factor = 10000 / metrics.spikein_reads

        # Check QC thresholds
        metrics = self._check_qc_thresholds(metrics)

        return metrics

    def _check_qc_thresholds(self, metrics: QCMetrics) -> QCMetrics:
        """Check QC metrics against thresholds."""
        warnings = []

        if metrics.mapping_rate < self.MIN_MAPPING_RATE:
            warnings.append(f"Low mapping rate: {metrics.mapping_rate:.1%} (threshold: {self.MIN_MAPPING_RATE:.1%})")

        if metrics.duplicate_rate > self.MAX_DUPLICATE_RATE:
            warnings.append(
                f"High duplication: {metrics.duplicate_rate:.1%} (threshold: {self.MAX_DUPLICATE_RATE:.1%})"
            )

        if metrics.frip_score == 0 and metrics.peak_count and metrics.peak_count > 0:
            warnings.append("FRiP score is 0.0 — no reads overlap peaks; possible failed sample")

        elif metrics.frip_score < self.MIN_FRIP and metrics.frip_score > 0:
            warnings.append(f"Low FRiP: {metrics.frip_score:.2%} (threshold: {self.MIN_FRIP:.1%})")

        if metrics.peak_count < self.MIN_PEAK_COUNT and metrics.peak_count > 0:
            warnings.append(f"Low peak count: {metrics.peak_count} (threshold: {self.MIN_PEAK_COUNT})")

        if metrics.median_fragment_size > 0:
            if metrics.median_fragment_size < self.OPTIMAL_FRAGMENT_SIZE[0]:
                warnings.append(f"Small fragments: {metrics.median_fragment_size}bp")
            elif metrics.median_fragment_size > self.OPTIMAL_FRAGMENT_SIZE[1]:
                warnings.append(f"Large fragments: {metrics.median_fragment_size}bp")

        metrics.qc_warnings = warnings
        metrics.pass_qc = len(warnings) == 0

        return metrics

    def generate_qc_report(self, samples: List[QCMetrics], output_path: str) -> pd.DataFrame:
        """
        Generate a QC report for multiple samples.

        Args:
            samples: List of QCMetrics objects
            output_path: Path to save CSV report

        Returns:
            DataFrame with QC metrics
        """
        data = [s.to_dict() for s in samples]
        df = pd.DataFrame(data)

        # Add summary columns
        df["qc_status"] = df["pass_qc"].map({True: "PASS", False: "WARN"})

        # Save report
        df.to_csv(output_path, index=False)
        logger.info(f"QC report saved to: {output_path}")

        return df
