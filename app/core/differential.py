"""
Differential Peak Analysis Module (Pure Python)

Replaces DiffBind with a Python-native implementation using:
- PyDESeq2 for differential analysis
- pybedtools for peak processing
- pyranges for genomic operations

Provides equivalent functionality to DiffBind:
- Consensus peak generation
- Read counting in peaks
- Normalization (RLE/TMM)
- Differential analysis with DESeq2
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
import numpy as np
import pandas as pd
from collections import defaultdict

logger = logging.getLogger(__name__)

# Optional imports with fallbacks
try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False
    logger.warning("pysam not available - BAM reading disabled")

try:
    import pyranges as pr
    PYRANGES_AVAILABLE = True
except ImportError:
    PYRANGES_AVAILABLE = False
    logger.warning("pyranges not available - using pandas fallback")

try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    PYDESEQ2_AVAILABLE = True
except ImportError:
    PYDESEQ2_AVAILABLE = False
    logger.warning("pydeseq2 not available - using scipy fallback")


@dataclass
class Sample:
    """Sample metadata."""
    sample_id: str
    condition: str
    histone_mark: str
    replicate: int = 1
    bam_file: Optional[str] = None
    peak_file: Optional[str] = None

    def __post_init__(self):
        self.sample_id = str(self.sample_id)


@dataclass
class DifferentialConfig:
    """Configuration for differential analysis."""
    comparison_name: str
    group1: str  # Treatment
    group2: str  # Control

    # Peak processing
    min_overlap: int = 2  # Minimum samples with peak
    summit_extend: int = 250  # Extend summit in each direction

    # Normalization
    normalize_method: str = "RLE"  # RLE, TMM, or median_ratio

    # Significance thresholds
    fdr_threshold: float = 0.1
    lfc_threshold: float = 0.5

    # Output
    output_dir: Optional[str] = None


@dataclass
class DifferentialResults:
    """Results from differential analysis."""
    comparison_name: str
    total_peaks: int
    significant_peaks: int
    gained_peaks: int
    lost_peaks: int

    # DataFrames
    all_peaks: pd.DataFrame = field(default_factory=pd.DataFrame)
    significant: pd.DataFrame = field(default_factory=pd.DataFrame)

    # File paths (if saved)
    output_dir: Optional[str] = None

    def to_dict(self) -> Dict:
        return {
            "comparison_name": self.comparison_name,
            "total_peaks": self.total_peaks,
            "significant_peaks": self.significant_peaks,
            "gained_peaks": self.gained_peaks,
            "lost_peaks": self.lost_peaks
        }


class DifferentialAnalyzer:
    """
    Pure Python differential peak analysis.

    Workflow:
    1. Load peak files and create consensus peaks
    2. Count reads in consensus peaks for each sample
    3. Normalize counts (RLE normalization)
    4. Run differential analysis (PyDESeq2)
    5. Annotate results with fold changes and p-values
    """

    def __init__(self):
        """Initialize the analyzer."""
        self._check_dependencies()

    def _check_dependencies(self):
        """Check for required dependencies."""
        if not PYSAM_AVAILABLE:
            logger.warning("pysam not installed - install with: pip install pysam")
        if not PYDESEQ2_AVAILABLE:
            logger.warning("pydeseq2 not installed - install with: pip install pydeseq2")

    def load_peaks(self, peak_file: str) -> pd.DataFrame:
        """
        Load peaks from BED/narrowPeak file.

        Args:
            peak_file: Path to peak file

        Returns:
            DataFrame with columns: chrom, start, end, name, score, summit
        """
        # Detect format by extension and column count
        with open(peak_file) as f:
            first_line = f.readline().strip()
            n_cols = len(first_line.split('\t'))

        if n_cols >= 10:  # narrowPeak format
            cols = ['chrom', 'start', 'end', 'name', 'score', 'strand',
                   'signalValue', 'pValue', 'qValue', 'summit']
            df = pd.read_csv(peak_file, sep='\t', header=None, names=cols[:n_cols])
        else:  # BED format
            cols = ['chrom', 'start', 'end', 'name', 'score', 'strand']
            df = pd.read_csv(peak_file, sep='\t', header=None, names=cols[:n_cols])

        # Ensure required columns
        df['chrom'] = df['chrom'].astype(str)
        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)

        # Add summit if not present (use peak center)
        if 'summit' not in df.columns:
            df['summit'] = (df['start'] + df['end']) // 2 - df['start']

        return df

    def create_consensus_peaks(
        self,
        samples: List[Sample],
        min_overlap: int = 2,
        summit_extend: int = 250
    ) -> pd.DataFrame:
        """
        Create consensus peak set from multiple samples.

        Args:
            samples: List of Sample objects with peak_file paths
            min_overlap: Minimum number of samples with overlapping peak
            summit_extend: Base pairs to extend from summit

        Returns:
            DataFrame with consensus peaks
        """
        all_peaks = []

        for sample in samples:
            if not sample.peak_file or not Path(sample.peak_file).exists():
                logger.warning(f"Peak file not found for {sample.sample_id}")
                continue

            peaks = self.load_peaks(sample.peak_file)
            peaks['sample'] = sample.sample_id
            all_peaks.append(peaks)

        if not all_peaks:
            raise ValueError("No valid peak files found")

        combined = pd.concat(all_peaks, ignore_index=True)

        # Use pyranges if available for efficient overlap detection
        if PYRANGES_AVAILABLE:
            consensus = self._merge_peaks_pyranges(combined, min_overlap, summit_extend)
        else:
            consensus = self._merge_peaks_pandas(combined, min_overlap, summit_extend)

        logger.info(f"Created {len(consensus)} consensus peaks from {len(samples)} samples")
        return consensus

    def _merge_peaks_pyranges(
        self,
        peaks: pd.DataFrame,
        min_overlap: int,
        summit_extend: int
    ) -> pd.DataFrame:
        """Merge peaks using pyranges."""
        # Create summit-centered regions
        peaks = peaks.copy()
        peaks['summit_pos'] = peaks['start'] + peaks['summit']
        peaks['Start'] = peaks['summit_pos'] - summit_extend
        peaks['End'] = peaks['summit_pos'] + summit_extend
        peaks['Chromosome'] = peaks['chrom']

        # Ensure no negative starts
        peaks['Start'] = peaks['Start'].clip(lower=0)

        gr = pr.PyRanges(peaks[['Chromosome', 'Start', 'End', 'sample']])

        # Cluster overlapping peaks
        clustered = gr.cluster()
        cluster_df = clustered.df

        # Count samples per cluster
        cluster_counts = cluster_df.groupby('Cluster').agg({
            'Chromosome': 'first',
            'Start': 'min',
            'End': 'max',
            'sample': 'nunique'
        }).reset_index()

        cluster_counts.columns = ['cluster', 'chrom', 'start', 'end', 'n_samples']

        # Filter by minimum overlap
        consensus = cluster_counts[cluster_counts['n_samples'] >= min_overlap].copy()
        consensus['peak_id'] = [f"peak_{i:06d}" for i in range(len(consensus))]

        return consensus[['peak_id', 'chrom', 'start', 'end', 'n_samples']]

    def _merge_peaks_pandas(
        self,
        peaks: pd.DataFrame,
        min_overlap: int,
        summit_extend: int
    ) -> pd.DataFrame:
        """Merge peaks using pandas (fallback)."""
        # Sort peaks
        peaks = peaks.copy()
        peaks['summit_pos'] = peaks['start'] + peaks['summit']
        peaks['new_start'] = (peaks['summit_pos'] - summit_extend).clip(lower=0)
        peaks['new_end'] = peaks['summit_pos'] + summit_extend

        peaks = peaks.sort_values(['chrom', 'new_start'])

        # Simple greedy merge
        merged = []
        current_chrom = None
        current_start = None
        current_end = None
        current_samples = set()

        for _, row in peaks.iterrows():
            if current_chrom is None:
                current_chrom = row['chrom']
                current_start = row['new_start']
                current_end = row['new_end']
                current_samples = {row['sample']}
            elif row['chrom'] == current_chrom and row['new_start'] <= current_end:
                # Overlapping - extend
                current_end = max(current_end, row['new_end'])
                current_samples.add(row['sample'])
            else:
                # Save current and start new
                if len(current_samples) >= min_overlap:
                    merged.append({
                        'chrom': current_chrom,
                        'start': current_start,
                        'end': current_end,
                        'n_samples': len(current_samples)
                    })
                current_chrom = row['chrom']
                current_start = row['new_start']
                current_end = row['new_end']
                current_samples = {row['sample']}

        # Don't forget last region
        if current_chrom and len(current_samples) >= min_overlap:
            merged.append({
                'chrom': current_chrom,
                'start': current_start,
                'end': current_end,
                'n_samples': len(current_samples)
            })

        consensus = pd.DataFrame(merged)
        consensus['peak_id'] = [f"peak_{i:06d}" for i in range(len(consensus))]

        return consensus[['peak_id', 'chrom', 'start', 'end', 'n_samples']]

    def count_reads_in_peaks(
        self,
        consensus_peaks: pd.DataFrame,
        samples: List[Sample]
    ) -> pd.DataFrame:
        """
        Count reads in consensus peaks for each sample.

        Args:
            consensus_peaks: DataFrame with consensus peaks
            samples: List of samples with bam_file paths

        Returns:
            Count matrix (peaks x samples)
        """
        if not PYSAM_AVAILABLE:
            raise ImportError("pysam required for read counting")

        counts = {}

        for sample in samples:
            if not sample.bam_file or not Path(sample.bam_file).exists():
                logger.warning(f"BAM file not found for {sample.sample_id}")
                continue

            sample_counts = []

            with pysam.AlignmentFile(sample.bam_file, "rb") as bam:
                for _, peak in consensus_peaks.iterrows():
                    try:
                        count = bam.count(
                            peak['chrom'],
                            int(peak['start']),
                            int(peak['end'])
                        )
                    except ValueError:
                        count = 0
                    sample_counts.append(count)

            counts[sample.sample_id] = sample_counts
            logger.info(f"Counted reads for {sample.sample_id}")

        count_df = pd.DataFrame(counts, index=consensus_peaks['peak_id'])
        return count_df

    def normalize_counts(
        self,
        counts: pd.DataFrame,
        method: str = "RLE"
    ) -> Tuple[pd.DataFrame, pd.Series]:
        """
        Normalize count matrix.

        Args:
            counts: Raw count matrix (peaks x samples)
            method: Normalization method (RLE, TMM, median_ratio)

        Returns:
            Tuple of (normalized counts, size factors)
        """
        if method == "RLE":
            size_factors = self._rle_normalization(counts)
        elif method == "TMM":
            size_factors = self._tmm_normalization(counts)
        else:  # median_ratio
            size_factors = self._median_ratio_normalization(counts)

        normalized = counts.div(size_factors, axis=1)

        return normalized, size_factors

    def _rle_normalization(self, counts: pd.DataFrame) -> pd.Series:
        """RLE (DESeq2-style) normalization."""
        # Calculate geometric mean per peak
        log_counts = np.log(counts.replace(0, np.nan))
        geo_means = log_counts.mean(axis=1)

        # Calculate size factors
        size_factors = {}
        for sample in counts.columns:
            log_ratios = log_counts[sample] - geo_means
            size_factors[sample] = np.exp(np.nanmedian(log_ratios))

        return pd.Series(size_factors)

    def _tmm_normalization(self, counts: pd.DataFrame) -> pd.Series:
        """TMM (edgeR-style) normalization."""
        # Use first sample as reference
        ref = counts.iloc[:, 0]

        size_factors = {}
        for sample in counts.columns:
            obs = counts[sample]

            # M values (log ratios)
            with np.errstate(divide='ignore', invalid='ignore'):
                m = np.log2(obs / ref)
                a = 0.5 * np.log2(obs * ref)

            # Trim extremes
            mask = np.isfinite(m) & np.isfinite(a)
            m_trimmed = m[mask]

            if len(m_trimmed) > 0:
                # Weighted mean
                size_factors[sample] = 2 ** np.median(m_trimmed)
            else:
                size_factors[sample] = 1.0

        return pd.Series(size_factors)

    def _median_ratio_normalization(self, counts: pd.DataFrame) -> pd.Series:
        """Simple median ratio normalization."""
        median_counts = counts.median(axis=1)

        size_factors = {}
        for sample in counts.columns:
            with np.errstate(divide='ignore', invalid='ignore'):
                ratios = counts[sample] / median_counts
            size_factors[sample] = np.nanmedian(ratios)

        return pd.Series(size_factors)

    def run_deseq2(
        self,
        counts: pd.DataFrame,
        samples: List[Sample],
        group1: str,
        group2: str
    ) -> pd.DataFrame:
        """
        Run differential analysis using PyDESeq2.

        Args:
            counts: Count matrix (peaks x samples)
            samples: Sample metadata
            group1: Treatment group name
            group2: Control group name

        Returns:
            DataFrame with differential analysis results
        """
        if not PYDESEQ2_AVAILABLE:
            logger.warning("PyDESeq2 not available, using fallback")
            return self._run_differential_fallback(counts, samples, group1, group2)

        # Prepare metadata
        sample_ids = [s.sample_id for s in samples if s.sample_id in counts.columns]
        metadata = pd.DataFrame({
            'sample': sample_ids,
            'condition': [s.condition for s in samples if s.sample_id in counts.columns]
        }).set_index('sample')

        # Subset counts to matching samples
        counts_subset = counts[sample_ids].T  # PyDESeq2 expects samples x genes

        # Create DESeq dataset
        dds = DeseqDataSet(
            counts=counts_subset,
            metadata=metadata,
            design_factors="condition"
        )

        # Run DESeq2
        dds.deseq2()

        # Extract results for the comparison
        stat_res = DeseqStats(dds, contrast=["condition", group1, group2])
        stat_res.summary()

        results = stat_res.results_df.copy()
        results['peak_id'] = results.index
        results = results.reset_index(drop=True)

        # Rename columns to match expected output
        results = results.rename(columns={
            'log2FoldChange': 'log2FC',
            'pvalue': 'pvalue',
            'padj': 'FDR',
            'baseMean': 'baseMean'
        })

        # Add direction
        results['direction'] = np.where(results['log2FC'] > 0, 'Gained', 'Lost')

        return results

    def _run_differential_fallback(
        self,
        counts: pd.DataFrame,
        samples: List[Sample],
        group1: str,
        group2: str
    ) -> pd.DataFrame:
        """Fallback differential analysis using scipy."""
        from scipy import stats

        # Get sample groups
        group1_samples = [s.sample_id for s in samples if s.condition == group1]
        group2_samples = [s.sample_id for s in samples if s.condition == group2]

        results = []

        for peak_id in counts.index:
            g1_counts = counts.loc[peak_id, group1_samples].values
            g2_counts = counts.loc[peak_id, group2_samples].values

            # Calculate fold change (add pseudocount)
            g1_mean = np.mean(g1_counts) + 1
            g2_mean = np.mean(g2_counts) + 1
            log2fc = np.log2(g1_mean / g2_mean)

            # Mann-Whitney U test
            try:
                _, pval = stats.mannwhitneyu(g1_counts, g2_counts, alternative='two-sided')
            except (ValueError, TypeError) as e:
                # Can fail if all values are identical or arrays are empty
                logger.debug(f"Mann-Whitney test failed for {peak_id}: {e}")
                pval = 1.0

            results.append({
                'peak_id': peak_id,
                'log2FC': log2fc,
                'pvalue': pval,
                'baseMean': (g1_mean + g2_mean) / 2
            })

        results_df = pd.DataFrame(results)

        # FDR correction
        try:
            from scipy.stats import false_discovery_control
            results_df['FDR'] = false_discovery_control(results_df['pvalue'])
        except (ImportError, AttributeError):
            # false_discovery_control not available in older scipy versions
            # Fall back to manual Benjamini-Hochberg
            logger.info("Using manual Benjamini-Hochberg FDR correction")
            pvals = results_df['pvalue'].values
            n = len(pvals)
            ranks = np.argsort(np.argsort(pvals)) + 1
            results_df['FDR'] = np.minimum(1, pvals * n / ranks)

        results_df['direction'] = np.where(results_df['log2FC'] > 0, 'Gained', 'Lost')

        return results_df

    def run(
        self,
        samples: List[Sample],
        config: DifferentialConfig
    ) -> DifferentialResults:
        """
        Run complete differential analysis pipeline.

        Args:
            samples: List of Sample objects
            config: Analysis configuration

        Returns:
            DifferentialResults object
        """
        logger.info(f"Starting differential analysis: {config.comparison_name}")
        logger.info(f"Comparing {config.group1} vs {config.group2}")

        # Step 1: Create consensus peaks
        logger.info("Creating consensus peaks...")
        consensus = self.create_consensus_peaks(
            samples,
            min_overlap=config.min_overlap,
            summit_extend=config.summit_extend
        )

        # Step 2: Count reads
        logger.info("Counting reads in peaks...")
        counts = self.count_reads_in_peaks(consensus, samples)

        # Step 3: Normalize
        logger.info(f"Normalizing with {config.normalize_method}...")
        normalized, size_factors = self.normalize_counts(counts, config.normalize_method)

        # Step 4: Differential analysis
        logger.info("Running differential analysis...")
        results = self.run_deseq2(counts, samples, config.group1, config.group2)

        # Step 5: Merge with peak coordinates
        results = results.merge(consensus, on='peak_id', how='left')

        # Step 6: Filter significant
        results['significant'] = (
            (results['FDR'] < config.fdr_threshold) &
            (np.abs(results['log2FC']) > config.lfc_threshold)
        )

        significant = results[results['significant']].copy()

        # Calculate summary stats
        gained = (significant['direction'] == 'Gained').sum()
        lost = (significant['direction'] == 'Lost').sum()

        # Save results if output_dir specified
        if config.output_dir:
            output_path = Path(config.output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            results.to_csv(output_path / "all_peaks.csv", index=False)
            significant.to_csv(output_path / "significant_peaks.csv", index=False)

            # Save size factors
            size_factors.to_csv(output_path / "size_factors.csv")

            logger.info(f"Results saved to {output_path}")

        return DifferentialResults(
            comparison_name=config.comparison_name,
            total_peaks=len(results),
            significant_peaks=len(significant),
            gained_peaks=int(gained),
            lost_peaks=int(lost),
            all_peaks=results,
            significant=significant,
            output_dir=config.output_dir
        )


# Convenience function
def run_differential_analysis(
    samples: List[Dict],
    comparison_name: str,
    group1: str,
    group2: str,
    output_dir: str = None,
    **kwargs
) -> DifferentialResults:
    """
    Convenience function to run differential analysis.

    Args:
        samples: List of sample dictionaries
        comparison_name: Name for the comparison
        group1: Treatment group
        group2: Control group
        output_dir: Output directory
        **kwargs: Additional config options

    Returns:
        DifferentialResults
    """
    sample_objects = [Sample(**s) for s in samples]

    config = DifferentialConfig(
        comparison_name=comparison_name,
        group1=group1,
        group2=group2,
        output_dir=output_dir,
        **kwargs
    )

    analyzer = DifferentialAnalyzer()
    return analyzer.run(sample_objects, config)
