"""
ATAC-seq Analysis Module

Provides:
- ATAC-seq QC metrics (TSS enrichment, fragment size)
- Peak calling integration
- Differential accessibility
- Footprinting analysis
- Integration with histone marks
"""

import numpy as np
import pandas as pd
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass


@dataclass
class ATACQCMetrics:
    """ATAC-seq quality control metrics."""
    total_reads: int
    mapped_reads: int
    mapping_rate: float
    mitochondrial_rate: float
    duplicate_rate: float
    tss_enrichment: float
    frip: float
    nfr_ratio: float  # Nucleosome-free region ratio
    peak_count: int
    median_fragment_size: int


@dataclass
class FragmentSizeDistribution:
    """Fragment size distribution data."""
    sizes: np.ndarray
    counts: np.ndarray
    nfr_fraction: float  # <100bp
    mono_fraction: float  # 180-247bp
    di_fraction: float  # 315-473bp
    tri_fraction: float  # 558-615bp


class ATACSeqAnalyzer:
    """
    Analyze ATAC-seq data for chromatin accessibility.

    ATAC-seq (Assay for Transposase-Accessible Chromatin) measures
    open chromatin regions across the genome.
    """

    def __init__(self):
        self.nfr_range = (0, 100)  # Nucleosome-free
        self.mono_range = (180, 247)  # Mononucleosome
        self.di_range = (315, 473)  # Dinucleosome

    def calculate_qc_metrics(
        self,
        bam_stats: Dict[str, Any],
        peaks: pd.DataFrame,
        tss_scores: np.ndarray = None
    ) -> ATACQCMetrics:
        """Calculate comprehensive QC metrics."""

        total = bam_stats.get('total_reads', 0)
        mapped = bam_stats.get('mapped_reads', 0)
        mito = bam_stats.get('mitochondrial_reads', 0)
        dups = bam_stats.get('duplicates', 0)

        return ATACQCMetrics(
            total_reads=total,
            mapped_reads=mapped,
            mapping_rate=mapped / total if total > 0 else 0,
            mitochondrial_rate=mito / mapped if mapped > 0 else 0,
            duplicate_rate=dups / mapped if mapped > 0 else 0,
            tss_enrichment=tss_scores.mean() if tss_scores is not None else 0,
            frip=bam_stats.get('frip', 0),
            nfr_ratio=bam_stats.get('nfr_ratio', 0),
            peak_count=len(peaks),
            median_fragment_size=bam_stats.get('median_fragment', 0)
        )

    def analyze_fragment_sizes(
        self,
        fragment_sizes: np.ndarray
    ) -> FragmentSizeDistribution:
        """Analyze fragment size distribution."""

        # Calculate histogram
        bins = np.arange(0, 1000, 10)
        counts, _ = np.histogram(fragment_sizes, bins=bins)
        sizes = bins[:-1] + 5  # Bin centers

        # Calculate fractions
        total = len(fragment_sizes)

        nfr_mask = fragment_sizes < self.nfr_range[1]
        mono_mask = (fragment_sizes >= self.mono_range[0]) & (fragment_sizes <= self.mono_range[1])
        di_mask = (fragment_sizes >= self.di_range[0]) & (fragment_sizes <= self.di_range[1])
        tri_mask = (fragment_sizes >= 558) & (fragment_sizes <= 615)

        return FragmentSizeDistribution(
            sizes=sizes,
            counts=counts,
            nfr_fraction=nfr_mask.sum() / total,
            mono_fraction=mono_mask.sum() / total,
            di_fraction=di_mask.sum() / total,
            tri_fraction=tri_mask.sum() / total
        )

    def calculate_tss_enrichment(
        self,
        coverage_at_tss: np.ndarray,
        window_size: int = 2000
    ) -> Tuple[np.ndarray, float]:
        """
        Calculate TSS enrichment score.

        TSS enrichment measures signal around transcription start sites
        compared to flanking regions.
        """

        # Average across all TSS
        avg_signal = coverage_at_tss.mean(axis=0)

        # Normalize to flanking regions
        flank_mean = np.concatenate([
            avg_signal[:window_size//4],
            avg_signal[-window_size//4:]
        ]).mean()

        if flank_mean > 0:
            normalized = avg_signal / flank_mean
            enrichment_score = normalized[window_size//2]  # Score at TSS
        else:
            normalized = avg_signal
            enrichment_score = 0

        return normalized, enrichment_score

    def differential_accessibility(
        self,
        count_matrix: pd.DataFrame,
        condition1: List[str],
        condition2: List[str],
        method: str = "deseq2"
    ) -> pd.DataFrame:
        """
        Identify differentially accessible regions.

        Args:
            count_matrix: Peaks x Samples count matrix
            condition1: Sample names for condition 1
            condition2: Sample names for condition 2
            method: Statistical method ('deseq2', 'edger', 'wilcoxon')

        Returns:
            DataFrame with differential accessibility results
        """

        # Get counts for each condition
        counts1 = count_matrix[condition1].values
        counts2 = count_matrix[condition2].values

        # Calculate fold changes
        mean1 = counts1.mean(axis=1) + 1  # Pseudocount
        mean2 = counts2.mean(axis=1) + 1
        log2fc = np.log2(mean2 / mean1)

        # Simple statistical test (for demo - real implementation would use DESeq2)
        from scipy import stats

        pvalues = []
        for i in range(len(count_matrix)):
            if counts1[i].std() > 0 or counts2[i].std() > 0:
                _, pval = stats.mannwhitneyu(counts1[i], counts2[i], alternative='two-sided')
            else:
                pval = 1.0
            pvalues.append(pval)

        # FDR correction
        pvalues = np.array(pvalues)
        fdr = self._benjamini_hochberg(pvalues)

        results = pd.DataFrame({
            'peak_id': count_matrix.index,
            'log2FoldChange': log2fc,
            'pvalue': pvalues,
            'padj': fdr,
            'baseMean': (mean1 + mean2) / 2
        })

        return results.sort_values('padj')

    def _benjamini_hochberg(self, pvalues: np.ndarray) -> np.ndarray:
        """Apply Benjamini-Hochberg FDR correction."""
        n = len(pvalues)
        ranked = np.argsort(pvalues)
        fdr = np.zeros(n)

        for i, rank in enumerate(ranked):
            fdr[rank] = pvalues[rank] * n / (i + 1)

        # Ensure monotonicity
        fdr = np.minimum.accumulate(fdr[::-1])[::-1]
        fdr = np.minimum(fdr, 1.0)

        return fdr

    def footprint_analysis(
        self,
        coverage: np.ndarray,
        motif_positions: np.ndarray,
        window: int = 100
    ) -> Dict[str, np.ndarray]:
        """
        Perform TF footprinting analysis.

        TF footprints appear as local dips in ATAC-seq signal
        where the TF protects DNA from transposase.
        """

        # Aggregate signal around motif positions
        n_motifs = len(motif_positions)
        aggregated = np.zeros((n_motifs, 2 * window))

        for i, pos in enumerate(motif_positions):
            start = max(0, pos - window)
            end = min(len(coverage), pos + window)

            # Extract and center
            signal = coverage[start:end]
            if len(signal) == 2 * window:
                aggregated[i] = signal

        # Calculate mean footprint
        mean_footprint = aggregated.mean(axis=0)

        # Calculate footprint depth (dip at center)
        flank_signal = np.concatenate([mean_footprint[:20], mean_footprint[-20:]]).mean()
        center_signal = mean_footprint[window-10:window+10].mean()
        footprint_depth = (flank_signal - center_signal) / flank_signal if flank_signal > 0 else 0

        return {
            'positions': np.arange(-window, window),
            'mean_footprint': mean_footprint,
            'individual_footprints': aggregated,
            'footprint_depth': footprint_depth
        }

    def integrate_with_histones(
        self,
        atac_peaks: pd.DataFrame,
        histone_peaks: Dict[str, pd.DataFrame],
        overlap_threshold: int = 1
    ) -> pd.DataFrame:
        """
        Integrate ATAC-seq with histone modification data.

        Classify accessible regions by their histone modification patterns.
        """

        results = atac_peaks.copy()

        # Check overlap with each histone mark
        for mark, peaks in histone_peaks.items():
            overlaps = self._count_overlaps(atac_peaks, peaks, overlap_threshold)
            results[f'{mark}_overlap'] = overlaps > 0
            results[f'{mark}_count'] = overlaps

        # Classify regions
        def classify_region(row):
            has_k4me3 = row.get('H3K4me3_overlap', False)
            has_k27ac = row.get('H3K27ac_overlap', False)
            has_k4me1 = row.get('H3K4me1_overlap', False)
            has_k27me3 = row.get('H3K27me3_overlap', False)

            if has_k4me3 and has_k27ac:
                return 'Active Promoter'
            elif has_k4me1 and has_k27ac:
                return 'Active Enhancer'
            elif has_k4me1 and not has_k27ac:
                return 'Poised Enhancer'
            elif has_k4me3 and has_k27me3:
                return 'Bivalent'
            elif has_k27me3:
                return 'Repressed'
            else:
                return 'Other Accessible'

        results['chromatin_state'] = results.apply(classify_region, axis=1)

        return results

    def _count_overlaps(
        self,
        peaks_a: pd.DataFrame,
        peaks_b: pd.DataFrame,
        min_overlap: int = 1
    ) -> np.ndarray:
        """Count overlapping peaks (simplified version)."""

        overlaps = np.zeros(len(peaks_a))

        for i, (_, peak_a) in enumerate(peaks_a.iterrows()):
            for _, peak_b in peaks_b.iterrows():
                if peak_a['chr'] == peak_b['chr']:
                    overlap = min(peak_a['end'], peak_b['end']) - max(peak_a['start'], peak_b['start'])
                    if overlap >= min_overlap:
                        overlaps[i] += 1

        return overlaps


def generate_demo_atac_data() -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Generate demo ATAC-seq data for testing."""
    np.random.seed(42)

    n_peaks = 50000

    # Generate peaks
    starts = np.random.randint(1000000, 200000000, n_peaks)
    widths = np.random.randint(200, 1000, n_peaks)

    peaks = pd.DataFrame({
        'chr': np.random.choice([f'chr{i}' for i in range(1, 23)], n_peaks),
        'start': starts,
        'end': starts + widths,
        'score': np.random.lognormal(3, 1, n_peaks),
        'summit': starts + widths // 2
    })

    # Generate QC stats
    bam_stats = {
        'total_reads': 80000000,
        'mapped_reads': 76000000,
        'mitochondrial_reads': 3800000,
        'duplicates': 7600000,
        'frip': 0.35,
        'nfr_ratio': 0.45,
        'median_fragment': 150
    }

    return peaks, bam_stats
