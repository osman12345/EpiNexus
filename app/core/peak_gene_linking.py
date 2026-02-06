"""
Peak-to-Gene Linking Module

Implements Activity-by-Contact (ABC) model and other methods:
- ABC score calculation
- Distance-based linking
- Correlation-based linking
- Hi-C integration
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass


@dataclass
class EnhancerGeneLink:
    """A predicted enhancer-gene regulatory link."""
    enhancer_id: str
    enhancer_chr: str
    enhancer_start: int
    enhancer_end: int
    gene_id: str
    gene_symbol: str
    gene_tss: int
    distance: int
    activity: float
    contact: float
    abc_score: float
    method: str


class PeakGeneLinkingEngine:
    """
    Link enhancers/peaks to target genes using multiple methods.

    The Activity-by-Contact (ABC) model predicts that an enhancer's effect
    on a gene is proportional to the enhancer's activity multiplied by
    the frequency of 3D contact between them.

    ABC Score = (Activity × Contact) / Σ(Activity × Contact) for all enhancers

    Reference: Fulco et al., Nature Genetics 2019
    """

    def __init__(
        self,
        max_distance: int = 1000000,  # 1 Mb default
        activity_column: str = 'signal',
        use_hic: bool = False
    ):
        self.max_distance = max_distance
        self.activity_column = activity_column
        self.use_hic = use_hic

    def link_peaks_to_genes(
        self,
        peaks: pd.DataFrame,
        genes: pd.DataFrame,
        hic_matrix: np.ndarray = None,
        method: str = 'abc'
    ) -> pd.DataFrame:
        """
        Link peaks to genes using specified method.

        Args:
            peaks: DataFrame with chr, start, end, signal columns
            genes: DataFrame with chr, tss, gene_id, gene_symbol columns
            hic_matrix: Optional Hi-C contact matrix
            method: 'abc', 'distance', 'nearest', or 'correlation'

        Returns:
            DataFrame of enhancer-gene links
        """

        if method == 'abc':
            return self._abc_linking(peaks, genes, hic_matrix)
        elif method == 'distance':
            return self._distance_linking(peaks, genes)
        elif method == 'nearest':
            return self._nearest_gene_linking(peaks, genes)
        elif method == 'correlation':
            return self._correlation_linking(peaks, genes)
        else:
            raise ValueError(f"Unknown method: {method}")

    def _abc_linking(
        self,
        peaks: pd.DataFrame,
        genes: pd.DataFrame,
        hic_matrix: np.ndarray = None
    ) -> pd.DataFrame:
        """
        Activity-by-Contact (ABC) model linking.

        ABC Score = (Activity × Contact) / Σ(Activity × Contact)
        """

        links = []

        for _, gene in genes.iterrows():
            gene_chr = gene['chr']
            gene_tss = gene['tss']

            # Find peaks within max_distance on same chromosome
            nearby_peaks = peaks[
                (peaks['chr'] == gene_chr) &
                (np.abs(peaks['start'] - gene_tss) <= self.max_distance)
            ]

            if len(nearby_peaks) == 0:
                continue

            # Calculate activity and contact for each peak
            activities = []
            contacts = []

            for _, peak in nearby_peaks.iterrows():
                # Activity = enhancer signal (H3K27ac, ATAC, etc.)
                activity = peak.get(self.activity_column, 1.0)

                # Contact = inverse distance (or Hi-C if available)
                peak_center = (peak['start'] + peak['end']) // 2
                distance = abs(peak_center - gene_tss)

                if hic_matrix is not None and self.use_hic:
                    # Use Hi-C contact frequency (simplified)
                    contact = self._get_hic_contact(peak_center, gene_tss, hic_matrix)
                else:
                    # Power-law decay with distance
                    contact = 1.0 / (1.0 + (distance / 5000) ** 0.87)

                activities.append(activity)
                contacts.append(contact)

            activities = np.array(activities)
            contacts = np.array(contacts)

            # Calculate ABC scores
            abc_products = activities * contacts
            abc_sum = abc_products.sum()

            if abc_sum > 0:
                abc_scores = abc_products / abc_sum
            else:
                abc_scores = np.zeros_like(abc_products)

            # Create links
            for i, (_, peak) in enumerate(nearby_peaks.iterrows()):
                peak_center = (peak['start'] + peak['end']) // 2
                distance = abs(peak_center - gene_tss)

                links.append({
                    'enhancer_id': peak.get('peak_id', f"peak_{peak.name}"),
                    'enhancer_chr': peak['chr'],
                    'enhancer_start': peak['start'],
                    'enhancer_end': peak['end'],
                    'gene_id': gene.get('gene_id', ''),
                    'gene_symbol': gene.get('gene_symbol', ''),
                    'gene_tss': gene_tss,
                    'distance': distance,
                    'activity': activities[i],
                    'contact': contacts[i],
                    'abc_score': abc_scores[i],
                    'method': 'ABC'
                })

        return pd.DataFrame(links)

    def _distance_linking(
        self,
        peaks: pd.DataFrame,
        genes: pd.DataFrame,
        distance_threshold: int = 50000
    ) -> pd.DataFrame:
        """Simple distance-based linking within threshold."""

        links = []

        for _, peak in peaks.iterrows():
            peak_chr = peak['chr']
            peak_center = (peak['start'] + peak['end']) // 2

            # Find genes within distance
            nearby_genes = genes[
                (genes['chr'] == peak_chr) &
                (np.abs(genes['tss'] - peak_center) <= distance_threshold)
            ]

            for _, gene in nearby_genes.iterrows():
                distance = abs(peak_center - gene['tss'])

                links.append({
                    'enhancer_id': peak.get('peak_id', f"peak_{peak.name}"),
                    'enhancer_chr': peak['chr'],
                    'enhancer_start': peak['start'],
                    'enhancer_end': peak['end'],
                    'gene_id': gene.get('gene_id', ''),
                    'gene_symbol': gene.get('gene_symbol', ''),
                    'gene_tss': gene['tss'],
                    'distance': distance,
                    'activity': peak.get(self.activity_column, 1.0),
                    'contact': 1.0 / (1.0 + distance / 10000),
                    'abc_score': 0.0,  # Not calculated for distance method
                    'method': 'Distance'
                })

        return pd.DataFrame(links)

    def _nearest_gene_linking(
        self,
        peaks: pd.DataFrame,
        genes: pd.DataFrame
    ) -> pd.DataFrame:
        """Link each peak to its nearest gene."""

        links = []

        for _, peak in peaks.iterrows():
            peak_chr = peak['chr']
            peak_center = (peak['start'] + peak['end']) // 2

            # Find genes on same chromosome
            chr_genes = genes[genes['chr'] == peak_chr]

            if len(chr_genes) == 0:
                continue

            # Find nearest gene
            distances = np.abs(chr_genes['tss'].values - peak_center)
            nearest_idx = np.argmin(distances)
            gene = chr_genes.iloc[nearest_idx]
            distance = distances[nearest_idx]

            links.append({
                'enhancer_id': peak.get('peak_id', f"peak_{peak.name}"),
                'enhancer_chr': peak['chr'],
                'enhancer_start': peak['start'],
                'enhancer_end': peak['end'],
                'gene_id': gene.get('gene_id', ''),
                'gene_symbol': gene.get('gene_symbol', ''),
                'gene_tss': gene['tss'],
                'distance': distance,
                'activity': peak.get(self.activity_column, 1.0),
                'contact': 1.0,
                'abc_score': 0.0,
                'method': 'Nearest'
            })

        return pd.DataFrame(links)

    def _correlation_linking(
        self,
        peaks: pd.DataFrame,
        genes: pd.DataFrame,
        expression_matrix: pd.DataFrame = None,
        signal_matrix: pd.DataFrame = None
    ) -> pd.DataFrame:
        """
        Link based on correlation between enhancer activity and gene expression.

        Requires matched samples with both ChIP-seq and RNA-seq data.
        """

        # Placeholder - would need expression and signal matrices
        return self._distance_linking(peaks, genes)

    def _get_hic_contact(
        self,
        pos1: int,
        pos2: int,
        hic_matrix: np.ndarray,
        resolution: int = 10000
    ) -> float:
        """Get Hi-C contact frequency between two positions."""

        bin1 = pos1 // resolution
        bin2 = pos2 // resolution

        if bin1 < hic_matrix.shape[0] and bin2 < hic_matrix.shape[1]:
            return hic_matrix[bin1, bin2]
        return 0.0

    def filter_links(
        self,
        links: pd.DataFrame,
        min_abc_score: float = 0.02,
        max_distance: int = None
    ) -> pd.DataFrame:
        """Filter links by ABC score and distance."""

        filtered = links.copy()

        if min_abc_score > 0:
            filtered = filtered[filtered['abc_score'] >= min_abc_score]

        if max_distance:
            filtered = filtered[filtered['distance'] <= max_distance]

        return filtered.sort_values('abc_score', ascending=False)

    def summarize_links(self, links: pd.DataFrame) -> Dict:
        """Generate summary statistics for links."""

        return {
            'total_links': len(links),
            'unique_enhancers': links['enhancer_id'].nunique(),
            'unique_genes': links['gene_symbol'].nunique(),
            'mean_abc_score': links['abc_score'].mean(),
            'median_distance': links['distance'].median(),
            'links_per_gene': links.groupby('gene_symbol').size().mean(),
            'links_per_enhancer': links.groupby('enhancer_id').size().mean()
        }


def generate_demo_genes() -> pd.DataFrame:
    """Generate demo gene annotation data."""
    np.random.seed(42)

    genes = [
        {'gene_id': 'ENSG00000136997', 'gene_symbol': 'MYC', 'chr': 'chr8', 'tss': 127736231},
        {'gene_id': 'ENSG00000141510', 'gene_symbol': 'TP53', 'chr': 'chr17', 'tss': 7687538},
        {'gene_id': 'ENSG00000012048', 'gene_symbol': 'BRCA1', 'chr': 'chr17', 'tss': 43125364},
        {'gene_id': 'ENSG00000146648', 'gene_symbol': 'EGFR', 'chr': 'chr7', 'tss': 55191822},
        {'gene_id': 'ENSG00000133703', 'gene_symbol': 'KRAS', 'chr': 'chr12', 'tss': 25245384},
    ]

    # Add more random genes
    for i in range(50):
        genes.append({
            'gene_id': f'ENSG{i:011d}',
            'gene_symbol': f'GENE{i}',
            'chr': f'chr{np.random.randint(1, 23)}',
            'tss': np.random.randint(1000000, 200000000)
        })

    return pd.DataFrame(genes)
