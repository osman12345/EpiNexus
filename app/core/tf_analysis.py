"""
Transcription Factor ChIP-seq Analysis Module

Specialized analysis for TF ChIP-seq data:
- Motif enrichment at binding sites
- Target gene identification
- Co-binding and co-localization analysis
- Binding site annotation
- Differential binding analysis
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from enum import Enum


class TFCategory(Enum):
    """Categories of transcription factors."""
    PIONEER = "pioneer"  # FOXA1, PU.1
    CHROMATIN_MODIFIER = "chromatin_modifier"  # BRD4, CBP
    SEQUENCE_SPECIFIC = "sequence_specific"  # MYC, P53
    ARCHITECTURAL = "architectural"  # CTCF, Cohesin
    SIGNALING = "signaling"  # STAT, SMAD
    HOUSEKEEPING = "housekeeping"  # SP1, NRF1


@dataclass
class TFBindingSite:
    """A transcription factor binding site."""
    chr: str
    start: int
    end: int
    peak_id: str
    tf_name: str
    signal: float
    summit: int
    motif_score: float = 0.0
    nearest_gene: str = ""
    distance_to_tss: int = 0
    genomic_annotation: str = ""


class TFAnalysisEngine:
    """
    Engine for transcription factor ChIP-seq analysis.

    Provides specialized analysis for TF binding data including
    motif analysis, target identification, and co-binding.
    """

    def __init__(
        self,
        tf_name: str,
        genome: str = "hg38",
        promoter_distance: int = 2000
    ):
        self.tf_name = tf_name
        self.genome = genome
        self.promoter_distance = promoter_distance

    def analyze_peaks(
        self,
        peaks: pd.DataFrame,
        genes: pd.DataFrame = None
    ) -> Dict:
        """
        Comprehensive TF peak analysis.

        Args:
            peaks: DataFrame with chr, start, end, signal columns
            genes: Optional gene annotation DataFrame

        Returns:
            Dictionary with analysis results
        """
        results = {
            'peak_summary': self._summarize_peaks(peaks),
            'genomic_distribution': self._analyze_genomic_distribution(peaks),
            'signal_analysis': self._analyze_signal(peaks),
        }

        if genes is not None:
            results['target_genes'] = self._identify_target_genes(peaks, genes)
            results['promoter_binding'] = self._analyze_promoter_binding(peaks, genes)

        return results

    def _summarize_peaks(self, peaks: pd.DataFrame) -> Dict:
        """Summarize peak characteristics."""
        peak_widths = peaks['end'] - peaks['start']

        return {
            'total_peaks': len(peaks),
            'mean_width': float(peak_widths.mean()),
            'median_width': float(peak_widths.median()),
            'min_width': int(peak_widths.min()),
            'max_width': int(peak_widths.max()),
            'total_coverage_bp': int(peak_widths.sum()),
            'chromosomes': peaks['chr'].nunique(),
            'peaks_per_chr': peaks.groupby('chr').size().to_dict()
        }

    def _analyze_genomic_distribution(self, peaks: pd.DataFrame) -> Dict:
        """Analyze genomic distribution of binding sites."""
        # Simulated distribution - in production would use actual annotation
        np.random.seed(42)
        n = len(peaks)

        # Typical TF distribution
        distribution = {
            'Promoter (<=1kb)': int(n * np.random.uniform(0.15, 0.35)),
            'Promoter (1-2kb)': int(n * np.random.uniform(0.05, 0.15)),
            '5\' UTR': int(n * np.random.uniform(0.02, 0.05)),
            '3\' UTR': int(n * np.random.uniform(0.03, 0.08)),
            'Exon': int(n * np.random.uniform(0.03, 0.08)),
            'Intron': int(n * np.random.uniform(0.25, 0.40)),
            'Intergenic': int(n * np.random.uniform(0.15, 0.30)),
        }

        # Normalize to sum to n
        total = sum(distribution.values())
        if total != n:
            distribution['Intergenic'] += (n - total)

        # Calculate percentages
        distribution_pct = {
            k: round(v / n * 100, 1) for k, v in distribution.items()
        }

        return {
            'counts': distribution,
            'percentages': distribution_pct
        }

    def _analyze_signal(self, peaks: pd.DataFrame) -> Dict:
        """Analyze peak signal characteristics."""
        if 'signal' not in peaks.columns:
            peaks['signal'] = np.random.lognormal(3, 1, len(peaks))

        signals = peaks['signal']

        return {
            'mean_signal': float(signals.mean()),
            'median_signal': float(signals.median()),
            'std_signal': float(signals.std()),
            'signal_range': [float(signals.min()), float(signals.max())],
            'high_confidence_peaks': int((signals > signals.quantile(0.75)).sum()),
            'low_confidence_peaks': int((signals < signals.quantile(0.25)).sum())
        }

    def _identify_target_genes(
        self,
        peaks: pd.DataFrame,
        genes: pd.DataFrame,
        max_distance: int = 100000
    ) -> pd.DataFrame:
        """Identify potential target genes for TF binding sites."""
        targets = []

        for _, peak in peaks.iterrows():
            peak_chr = peak['chr']
            peak_center = (peak['start'] + peak['end']) // 2

            # Find nearby genes
            nearby = genes[
                (genes['chr'] == peak_chr) &
                (np.abs(genes['tss'] - peak_center) <= max_distance)
            ]

            for _, gene in nearby.iterrows():
                distance = peak_center - gene['tss']

                targets.append({
                    'peak_id': peak.get('peak_id', f"peak_{peak.name}"),
                    'gene_id': gene.get('gene_id', ''),
                    'gene_symbol': gene.get('gene_symbol', ''),
                    'distance_to_tss': distance,
                    'abs_distance': abs(distance),
                    'peak_signal': peak.get('signal', 1.0),
                    'is_promoter': abs(distance) <= self.promoter_distance,
                    'binding_location': 'upstream' if distance < 0 else 'downstream'
                })

        if not targets:
            return pd.DataFrame()

        df = pd.DataFrame(targets)
        return df.sort_values('abs_distance')

    def _analyze_promoter_binding(
        self,
        peaks: pd.DataFrame,
        genes: pd.DataFrame
    ) -> Dict:
        """Analyze TF binding at promoters."""
        targets = self._identify_target_genes(peaks, genes, max_distance=self.promoter_distance)

        if len(targets) == 0:
            return {
                'promoter_bound_genes': 0,
                'fraction_genes_bound': 0.0,
                'fraction_peaks_at_promoters': 0.0
            }

        promoter_targets = targets[targets['is_promoter']]

        return {
            'promoter_bound_genes': promoter_targets['gene_symbol'].nunique(),
            'fraction_genes_bound': promoter_targets['gene_symbol'].nunique() / len(genes),
            'fraction_peaks_at_promoters': len(promoter_targets) / len(peaks),
            'top_promoter_targets': promoter_targets.groupby('gene_symbol')['peak_signal'].max().nlargest(20).to_dict()
        }


class TFCoBindingAnalyzer:
    """Analyze co-binding between multiple TFs."""

    def __init__(self, overlap_threshold: int = 100):
        self.overlap_threshold = overlap_threshold

    def analyze_cobinding(
        self,
        tf1_peaks: pd.DataFrame,
        tf2_peaks: pd.DataFrame,
        tf1_name: str = "TF1",
        tf2_name: str = "TF2"
    ) -> Dict:
        """
        Analyze co-localization between two TFs.

        Args:
            tf1_peaks: Peaks for first TF
            tf2_peaks: Peaks for second TF
            tf1_name: Name of first TF
            tf2_name: Name of second TF

        Returns:
            Dictionary with co-binding statistics
        """
        # Find overlapping peaks
        overlaps = self._find_overlaps(tf1_peaks, tf2_peaks)

        n1 = len(tf1_peaks)
        n2 = len(tf2_peaks)
        n_overlap = len(overlaps)

        # Calculate enrichment
        genome_size = 3e9
        tf1_coverage = (tf1_peaks['end'] - tf1_peaks['start']).sum()
        tf2_coverage = (tf2_peaks['end'] - tf2_peaks['start']).sum()

        expected_overlap = (tf1_coverage * tf2_coverage) / genome_size
        fold_enrichment = n_overlap / expected_overlap if expected_overlap > 0 else 0

        return {
            f'{tf1_name}_peaks': n1,
            f'{tf2_name}_peaks': n2,
            'overlapping_peaks': n_overlap,
            f'{tf1_name}_only': n1 - n_overlap,
            f'{tf2_name}_only': n2 - n_overlap,
            f'{tf1_name}_overlap_fraction': n_overlap / n1 if n1 > 0 else 0,
            f'{tf2_name}_overlap_fraction': n_overlap / n2 if n2 > 0 else 0,
            'fold_enrichment': fold_enrichment,
            'jaccard_index': n_overlap / (n1 + n2 - n_overlap) if (n1 + n2 - n_overlap) > 0 else 0,
            'overlap_peaks': overlaps
        }

    def _find_overlaps(
        self,
        peaks1: pd.DataFrame,
        peaks2: pd.DataFrame
    ) -> pd.DataFrame:
        """Find overlapping peaks between two sets."""
        overlaps = []

        for _, p1 in peaks1.iterrows():
            chr1 = p1['chr']
            start1 = p1['start'] - self.overlap_threshold
            end1 = p1['end'] + self.overlap_threshold

            # Find overlapping peaks in set 2
            matching = peaks2[
                (peaks2['chr'] == chr1) &
                (peaks2['start'] < end1) &
                (peaks2['end'] > start1)
            ]

            if len(matching) > 0:
                for _, p2 in matching.iterrows():
                    overlaps.append({
                        'chr': chr1,
                        'start': max(p1['start'], p2['start']),
                        'end': min(p1['end'], p2['end']),
                        'peak1_id': p1.get('peak_id', ''),
                        'peak2_id': p2.get('peak_id', ''),
                        'peak1_signal': p1.get('signal', 1.0),
                        'peak2_signal': p2.get('signal', 1.0)
                    })

        return pd.DataFrame(overlaps) if overlaps else pd.DataFrame()


class MotifEnrichmentAnalyzer:
    """Analyze motif enrichment at TF binding sites."""

    KNOWN_TF_MOTIFS = {
        'MYC': {'consensus': 'CACGTG', 'name': 'E-box'},
        'P53': {'consensus': 'RRRCWWGYYY', 'name': 'p53 RE'},
        'CTCF': {'consensus': 'CCGCGNGGNGGCAG', 'name': 'CTCF motif'},
        'STAT3': {'consensus': 'TTCNNNGAA', 'name': 'STAT binding'},
        'NFkB': {'consensus': 'GGGRNWYYCC', 'name': 'NFkB motif'},
        'AP1': {'consensus': 'TGASTCA', 'name': 'AP-1 site'},
        'SP1': {'consensus': 'GGGCGG', 'name': 'GC box'},
        'GATA': {'consensus': 'WGATAR', 'name': 'GATA motif'},
        'ETS': {'consensus': 'GGAA', 'name': 'ETS motif'},
        'FOX': {'consensus': 'TRTTKRY', 'name': 'Forkhead box'}
    }

    def analyze_motif_enrichment(
        self,
        peaks: pd.DataFrame,
        tf_name: str,
        background: pd.DataFrame = None
    ) -> Dict:
        """
        Analyze motif enrichment at binding sites.

        In production, this would use HOMER or MEME-ChIP.
        Here we simulate results for demonstration.
        """
        np.random.seed(42)
        n_peaks = len(peaks)

        # Get expected motif for this TF
        expected_motif = self.KNOWN_TF_MOTIFS.get(
            tf_name.upper(),
            {'consensus': 'NNNNNNN', 'name': 'Unknown'}
        )

        # Simulate motif finding results
        results = {
            'primary_motif': {
                'consensus': expected_motif['consensus'],
                'name': expected_motif['name'],
                'peaks_with_motif': int(n_peaks * np.random.uniform(0.4, 0.8)),
                'enrichment_pvalue': 10 ** -np.random.uniform(20, 100),
                'centrality_score': np.random.uniform(0.7, 0.95)
            },
            'secondary_motifs': [],
            'de_novo_motifs': []
        }

        # Add secondary known motifs
        other_tfs = [tf for tf in self.KNOWN_TF_MOTIFS.keys() if tf != tf_name.upper()]
        np.random.shuffle(other_tfs)

        for tf in other_tfs[:5]:
            motif = self.KNOWN_TF_MOTIFS[tf]
            results['secondary_motifs'].append({
                'tf_name': tf,
                'consensus': motif['consensus'],
                'name': motif['name'],
                'peaks_with_motif': int(n_peaks * np.random.uniform(0.05, 0.3)),
                'enrichment_pvalue': 10 ** -np.random.uniform(2, 30),
                'possible_cobinder': np.random.random() > 0.5
            })

        # De novo motifs
        for i in range(3):
            results['de_novo_motifs'].append({
                'motif_id': f'denovo_{i+1}',
                'consensus': ''.join(np.random.choice(['A', 'C', 'G', 'T', 'N'], 8)),
                'peaks_with_motif': int(n_peaks * np.random.uniform(0.1, 0.4)),
                'enrichment_pvalue': 10 ** -np.random.uniform(5, 50),
                'best_match': np.random.choice(other_tfs) if other_tfs else 'Unknown'
            })

        return results


def generate_demo_tf_peaks(tf_name: str = "MYC", n_peaks: int = 5000) -> pd.DataFrame:
    """Generate demo TF binding site data."""
    np.random.seed(42)

    # TF-specific peak characteristics
    tf_params = {
        'MYC': {'width_mean': 200, 'width_std': 50, 'signal_mean': 4},
        'P53': {'width_mean': 300, 'width_std': 100, 'signal_mean': 3},
        'CTCF': {'width_mean': 150, 'width_std': 30, 'signal_mean': 5},
        'STAT3': {'width_mean': 250, 'width_std': 80, 'signal_mean': 3.5},
        'NFkB': {'width_mean': 350, 'width_std': 120, 'signal_mean': 3}
    }

    params = tf_params.get(tf_name.upper(), {'width_mean': 250, 'width_std': 75, 'signal_mean': 3.5})

    chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX']
    chr_weights = np.array([8, 7, 6, 5, 5, 5, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 3])
    chr_weights = chr_weights / chr_weights.sum()

    peaks = pd.DataFrame({
        'chr': np.random.choice(chromosomes, n_peaks, p=chr_weights),
        'start': np.random.randint(1000000, 200000000, n_peaks),
        'peak_id': [f'{tf_name}_peak_{i}' for i in range(n_peaks)],
        'signal': np.random.lognormal(params['signal_mean'], 1, n_peaks),
        'qvalue': 10 ** -np.random.uniform(2, 20, n_peaks)
    })

    widths = np.maximum(50, np.random.normal(params['width_mean'], params['width_std'], n_peaks).astype(int))
    peaks['end'] = peaks['start'] + widths
    peaks['summit'] = peaks['start'] + widths // 2

    return peaks[['chr', 'start', 'end', 'peak_id', 'signal', 'summit', 'qvalue']]


def generate_demo_genes() -> pd.DataFrame:
    """Generate demo gene annotation data."""
    np.random.seed(42)

    # Known genes
    known_genes = [
        {'gene_id': 'ENSG00000136997', 'gene_symbol': 'MYC', 'chr': 'chr8', 'tss': 127736231, 'strand': '+'},
        {'gene_id': 'ENSG00000141510', 'gene_symbol': 'TP53', 'chr': 'chr17', 'tss': 7687538, 'strand': '-'},
        {'gene_id': 'ENSG00000012048', 'gene_symbol': 'BRCA1', 'chr': 'chr17', 'tss': 43125364, 'strand': '-'},
        {'gene_id': 'ENSG00000146648', 'gene_symbol': 'EGFR', 'chr': 'chr7', 'tss': 55191822, 'strand': '+'},
        {'gene_id': 'ENSG00000133703', 'gene_symbol': 'KRAS', 'chr': 'chr12', 'tss': 25245384, 'strand': '-'},
        {'gene_id': 'ENSG00000171862', 'gene_symbol': 'PTEN', 'chr': 'chr10', 'tss': 87863113, 'strand': '+'},
        {'gene_id': 'ENSG00000157764', 'gene_symbol': 'BRAF', 'chr': 'chr7', 'tss': 140924764, 'strand': '-'},
        {'gene_id': 'ENSG00000181019', 'gene_symbol': 'NQO1', 'chr': 'chr16', 'tss': 69710984, 'strand': '+'},
    ]

    genes = known_genes.copy()

    # Add random genes
    for i in range(200):
        genes.append({
            'gene_id': f'ENSG{i:011d}',
            'gene_symbol': f'GENE{i}',
            'chr': f'chr{np.random.randint(1, 23)}',
            'tss': np.random.randint(1000000, 200000000),
            'strand': np.random.choice(['+', '-'])
        })

    return pd.DataFrame(genes)
