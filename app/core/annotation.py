"""
Peak Annotation Module (Pure Python)

Replaces ChIPseeker with a Python-native implementation using:
- pyranges for genomic interval operations
- GTF/GFF parsing for gene models

Provides:
- Peak to gene annotation (nearest gene, TSS distance)
- Genomic feature classification (promoter, exon, intron, intergenic)
- Regulatory region overlap detection
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)

try:
    import pyranges as pr
    PYRANGES_AVAILABLE = True
except ImportError:
    PYRANGES_AVAILABLE = False


@dataclass
class AnnotationConfig:
    """Configuration for peak annotation."""
    # TSS region definition
    tss_upstream: int = 3000
    tss_downstream: int = 3000

    # Feature priorities (lower = higher priority)
    feature_priority: Dict[str, int] = None

    def __post_init__(self):
        if self.feature_priority is None:
            self.feature_priority = {
                'Promoter': 1,
                '5UTR': 2,
                '3UTR': 3,
                'Exon': 4,
                'Intron': 5,
                'Downstream': 6,
                'Intergenic': 7
            }


class GeneAnnotation:
    """
    Gene annotation from GTF/GFF files.

    Provides genomic features for peak annotation:
    - Gene locations
    - Transcript structures
    - Promoter regions
    - Exons and introns
    """

    # Pre-defined URLs for common genomes
    GTF_URLS = {
        "mm10": "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz",
        "mm39": "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.annotation.gtf.gz",
        "hg38": "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz",
        "hg19": "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh37_mapping/gencode.v44lift37.annotation.gtf.gz",
    }

    def __init__(self, gtf_file: Optional[str] = None, genome: str = "mm10"):
        """
        Initialize gene annotation.

        Args:
            gtf_file: Path to GTF file (optional)
            genome: Genome name (used for default GTF if file not provided)
        """
        self.genome = genome
        self.genes = None
        self.transcripts = None
        self.exons = None
        self.gene_info = None

        if gtf_file:
            self.load_gtf(gtf_file)

    def load_gtf(self, gtf_file: str):
        """
        Load gene annotation from GTF file.

        Args:
            gtf_file: Path to GTF/GFF file
        """
        logger.info(f"Loading GTF: {gtf_file}")

        # Parse GTF using pandas
        df = self._parse_gtf(gtf_file)

        # Extract genes
        self.genes = df[df['feature'] == 'gene'].copy()
        self.genes = self._standardize_columns(self.genes)

        # Extract transcripts
        self.transcripts = df[df['feature'] == 'transcript'].copy()

        # Extract exons
        self.exons = df[df['feature'] == 'exon'].copy()

        # Create gene info lookup
        self.gene_info = self.genes.set_index('gene_id')[
            ['gene_name', 'chrom', 'start', 'end', 'strand']
        ].to_dict('index')

        logger.info(f"Loaded {len(self.genes)} genes, {len(self.exons)} exons")

    def _parse_gtf(self, gtf_file: str) -> pd.DataFrame:
        """Parse GTF file into DataFrame."""
        records = []

        # Handle gzipped files
        import gzip
        opener = gzip.open if gtf_file.endswith('.gz') else open

        with opener(gtf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                chrom, source, feature, start, end, score, strand, frame, attributes = fields

                # Parse attributes
                attrs = self._parse_attributes(attributes)

                records.append({
                    'chrom': chrom,
                    'source': source,
                    'feature': feature,
                    'start': int(start),
                    'end': int(end),
                    'strand': strand,
                    **attrs
                })

        return pd.DataFrame(records)

    def _parse_attributes(self, attr_string: str) -> Dict[str, str]:
        """Parse GTF attribute string."""
        attrs = {}
        for item in attr_string.strip().split(';'):
            item = item.strip()
            if not item:
                continue

            # Handle both GTF and GFF formats
            if '=' in item:  # GFF3
                key, value = item.split('=', 1)
            elif ' ' in item:  # GTF
                parts = item.split(' ', 1)
                key = parts[0]
                value = parts[1].strip('"') if len(parts) > 1 else ''
            else:
                continue

            attrs[key] = value.strip('"')

        return attrs

    def _standardize_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Standardize column names across GTF versions."""
        df = df.copy()

        # Standardize gene_id
        if 'gene_id' not in df.columns:
            for col in ['ID', 'Name', 'gene']:
                if col in df.columns:
                    df['gene_id'] = df[col]
                    break

        # Standardize gene_name
        if 'gene_name' not in df.columns:
            for col in ['gene_symbol', 'Name', 'gene_id']:
                if col in df.columns:
                    df['gene_name'] = df[col]
                    break

        return df

    def get_promoter_regions(
        self,
        upstream: int = 3000,
        downstream: int = 3000
    ) -> pd.DataFrame:
        """
        Get promoter regions for all genes.

        Args:
            upstream: Base pairs upstream of TSS
            downstream: Base pairs downstream of TSS

        Returns:
            DataFrame with promoter regions
        """
        if self.genes is None:
            raise ValueError("Gene annotation not loaded")

        promoters = self.genes.copy()

        # Calculate TSS based on strand
        promoters['tss'] = np.where(
            promoters['strand'] == '+',
            promoters['start'],
            promoters['end']
        )

        # Define promoter region
        promoters['prom_start'] = np.where(
            promoters['strand'] == '+',
            promoters['tss'] - upstream,
            promoters['tss'] - downstream
        ).clip(min=0)

        promoters['prom_end'] = np.where(
            promoters['strand'] == '+',
            promoters['tss'] + downstream,
            promoters['tss'] + upstream
        )

        return promoters[['gene_id', 'gene_name', 'chrom', 'prom_start', 'prom_end', 'strand', 'tss']]


class PeakAnnotator:
    """
    Annotate peaks with genomic features.

    Provides ChIPseeker-like functionality:
    - Nearest gene assignment
    - TSS distance calculation
    - Genomic feature classification
    """

    def __init__(
        self,
        gene_annotation: GeneAnnotation,
        config: AnnotationConfig = None
    ):
        """
        Initialize peak annotator.

        Args:
            gene_annotation: Loaded GeneAnnotation object
            config: Annotation configuration
        """
        self.genes = gene_annotation
        self.config = config or AnnotationConfig()

        # Pre-compute promoter regions
        self.promoters = self.genes.get_promoter_regions(
            upstream=self.config.tss_upstream,
            downstream=self.config.tss_downstream
        )

    def annotate_peaks(
        self,
        peaks: pd.DataFrame,
        chrom_col: str = 'chrom',
        start_col: str = 'start',
        end_col: str = 'end'
    ) -> pd.DataFrame:
        """
        Annotate peaks with nearest gene and genomic features.

        Args:
            peaks: DataFrame with peak coordinates
            chrom_col: Chromosome column name
            start_col: Start position column name
            end_col: End position column name

        Returns:
            Annotated peaks DataFrame
        """
        peaks = peaks.copy()

        # Calculate peak centers
        peaks['peak_center'] = (peaks[start_col] + peaks[end_col]) // 2

        # Find nearest gene for each peak
        annotations = []

        for idx, peak in peaks.iterrows():
            ann = self._annotate_single_peak(
                peak[chrom_col],
                peak[start_col],
                peak[end_col],
                peak['peak_center']
            )
            ann['peak_idx'] = idx
            annotations.append(ann)

        ann_df = pd.DataFrame(annotations)

        # Merge with original peaks
        result = peaks.merge(ann_df, left_index=True, right_on='peak_idx', how='left')
        result = result.drop(columns=['peak_idx'], errors='ignore')

        return result

    def _annotate_single_peak(
        self,
        chrom: str,
        start: int,
        end: int,
        center: int
    ) -> Dict[str, Any]:
        """Annotate a single peak."""
        result = {
            'gene_id': None,
            'gene_name': None,
            'distance_to_tss': None,
            'annotation': 'Intergenic',
            'gene_strand': None
        }

        # Get genes on same chromosome
        chrom_genes = self.genes.genes[self.genes.genes['chrom'] == chrom]

        if len(chrom_genes) == 0:
            return result

        # Get promoters on same chromosome
        chrom_promoters = self.promoters[self.promoters['chrom'] == chrom]

        # Check promoter overlap first
        promoter_hits = chrom_promoters[
            (chrom_promoters['prom_start'] <= end) &
            (chrom_promoters['prom_end'] >= start)
        ]

        if len(promoter_hits) > 0:
            # Take closest promoter
            promoter_hits = promoter_hits.copy()
            promoter_hits['dist'] = np.abs(promoter_hits['tss'] - center)
            closest = promoter_hits.loc[promoter_hits['dist'].idxmin()]

            result['gene_id'] = closest['gene_id']
            result['gene_name'] = closest['gene_name']
            result['distance_to_tss'] = int(center - closest['tss'])
            result['gene_strand'] = closest['strand']
            result['annotation'] = 'Promoter'
            return result

        # Find nearest gene
        chrom_genes = chrom_genes.copy()
        chrom_genes['tss'] = np.where(
            chrom_genes['strand'] == '+',
            chrom_genes['start'],
            chrom_genes['end']
        )
        chrom_genes['dist'] = np.abs(chrom_genes['tss'] - center)
        closest_gene = chrom_genes.loc[chrom_genes['dist'].idxmin()]

        result['gene_id'] = closest_gene.get('gene_id')
        result['gene_name'] = closest_gene.get('gene_name')
        result['distance_to_tss'] = int(center - closest_gene['tss'])
        result['gene_strand'] = closest_gene['strand']

        # Determine annotation type
        gene_start = closest_gene['start']
        gene_end = closest_gene['end']

        if start >= gene_start and end <= gene_end:
            # Inside gene body
            result['annotation'] = 'Intron'  # Simplified - would need exon data for precise
        elif center > gene_end:
            if closest_gene['strand'] == '+':
                result['annotation'] = 'Downstream'
            else:
                result['annotation'] = 'Upstream'
        else:
            if closest_gene['strand'] == '+':
                result['annotation'] = 'Upstream'
            else:
                result['annotation'] = 'Downstream'

        return result

    def annotate_to_genes(
        self,
        peaks: pd.DataFrame,
        summarize: bool = True
    ) -> pd.DataFrame:
        """
        Annotate peaks and optionally create gene-level summary.

        Args:
            peaks: DataFrame with peaks (must include diff analysis columns)
            summarize: If True, return gene-level summary

        Returns:
            Annotated peaks or gene summary
        """
        # First annotate peaks
        annotated = self.annotate_peaks(peaks)

        if not summarize:
            return annotated

        # Create gene-level summary
        gene_summary = annotated.groupby(['gene_id', 'gene_name']).agg({
            'peak_id': 'count',
            'annotation': lambda x: ','.join(x.unique()),
            'distance_to_tss': 'min',
            'FDR': 'min' if 'FDR' in annotated.columns else lambda x: None,
            'log2FC': 'mean' if 'log2FC' in annotated.columns else lambda x: None,
            'direction': lambda x: ','.join(x.unique()) if 'direction' in annotated.columns else None
        }).reset_index()

        gene_summary.columns = [
            'gene_id', 'gene_name', 'n_peaks', 'annotations',
            'min_distance_to_tss', 'min_FDR', 'mean_log2FC', 'directions'
        ]

        return gene_summary.sort_values('min_FDR' if 'min_FDR' in gene_summary.columns else 'n_peaks')


def annotate_differential_peaks(
    peaks_df: pd.DataFrame,
    gtf_file: str,
    tss_upstream: int = 3000,
    tss_downstream: int = 3000
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Convenience function to annotate differential peaks.

    Args:
        peaks_df: DataFrame from differential analysis
        gtf_file: Path to GTF annotation file
        tss_upstream: TSS upstream distance
        tss_downstream: TSS downstream distance

    Returns:
        Tuple of (annotated_peaks, gene_summary)
    """
    # Load annotation
    genes = GeneAnnotation(gtf_file)

    # Configure annotator
    config = AnnotationConfig(
        tss_upstream=tss_upstream,
        tss_downstream=tss_downstream
    )

    annotator = PeakAnnotator(genes, config)

    # Annotate peaks
    annotated_peaks = annotator.annotate_peaks(peaks_df)

    # Create gene summary
    gene_summary = annotator.annotate_to_genes(peaks_df, summarize=True)

    return annotated_peaks, gene_summary


class AnnotationDatabase:
    """
    Pre-built annotation databases for common genomes.

    Provides easy access to gene annotations without needing
    to download and parse GTF files each time.
    """

    CACHE_DIR = Path.home() / ".histone_analyzer" / "annotations"

    def __init__(self):
        """Initialize annotation database."""
        self.CACHE_DIR.mkdir(parents=True, exist_ok=True)

    def get_annotation(self, genome: str) -> GeneAnnotation:
        """
        Get gene annotation for a genome.

        Downloads and caches GTF file if not present.

        Args:
            genome: Genome name (mm10, mm39, hg38, hg19)

        Returns:
            GeneAnnotation object
        """
        cache_dir = self.CACHE_DIR / f"{genome}_annotation"

        # Support new parquet cache; fall back to legacy .pkl if present
        if cache_dir.is_dir() and (cache_dir / "genes.parquet").exists():
            logger.info(f"Loading cached annotation for {genome} (parquet)")
            return self._load_cached(cache_dir)

        # Legacy pickle cache – load then migrate to parquet
        legacy_pkl = self.CACHE_DIR / f"{genome}_genes.pkl"
        if legacy_pkl.exists():
            logger.info(f"Migrating legacy pickle cache for {genome} to parquet")
            annotation = self._load_cached_pickle(legacy_pkl)
            self._save_cached(annotation, cache_dir)
            legacy_pkl.unlink()
            return annotation

        # Download and parse GTF
        if genome in GeneAnnotation.GTF_URLS:
            gtf_url = GeneAnnotation.GTF_URLS[genome]
            logger.info(f"Downloading GTF for {genome}...")
            gtf_file = self._download_gtf(gtf_url, genome)

            annotation = GeneAnnotation(gtf_file, genome)

            # Cache for future use
            self._save_cached(annotation, cache_dir)

            return annotation
        else:
            raise ValueError(f"No annotation available for genome: {genome}")

    def _download_gtf(self, url: str, genome: str) -> str:
        """Download GTF file."""
        import urllib.request

        local_file = self.CACHE_DIR / f"{genome}.gtf.gz"

        if not local_file.exists():
            logger.info(f"Downloading from {url}...")
            urllib.request.urlretrieve(url, local_file)

        return str(local_file)

    def _save_cached(self, annotation: GeneAnnotation, cache_dir: Path):
        """Save annotation DataFrames to parquet files.

        Parquet is safer than pickle (no arbitrary code execution on load),
        faster to read, and produces smaller files with columnar compression.
        """
        cache_dir.mkdir(parents=True, exist_ok=True)

        for name in ("genes", "transcripts", "exons", "gene_info"):
            df = getattr(annotation, name, None)
            if df is not None:
                df.to_parquet(cache_dir / f"{name}.parquet", index=True)

        # Store genome name as a small metadata file
        (cache_dir / "genome.txt").write_text(annotation.genome)

    def _load_cached(self, cache_dir: Path) -> GeneAnnotation:
        """Load annotation from parquet cache directory."""
        annotation = GeneAnnotation.__new__(GeneAnnotation)

        for name in ("genes", "transcripts", "exons", "gene_info"):
            parquet_path = cache_dir / f"{name}.parquet"
            if parquet_path.exists():
                setattr(annotation, name, pd.read_parquet(parquet_path))
            else:
                setattr(annotation, name, None)

        genome_file = cache_dir / "genome.txt"
        annotation.genome = genome_file.read_text().strip() if genome_file.exists() else "unknown"

        return annotation

    @staticmethod
    def _load_cached_pickle(cache_file: Path) -> "GeneAnnotation":
        """Load legacy pickle cache (for one-time migration to parquet)."""
        import pickle
        with open(cache_file, 'rb') as f:
            data = pickle.load(f)

        annotation = GeneAnnotation.__new__(GeneAnnotation)
        annotation.genes = data['genes']
        annotation.transcripts = data['transcripts']
        annotation.exons = data['exons']
        annotation.gene_info = data['gene_info']
        annotation.genome = data['genome']

        return annotation
