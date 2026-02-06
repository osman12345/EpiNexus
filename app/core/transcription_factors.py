"""
Transcription Factor Analysis Module

Comprehensive TF analysis integrated with histone marks:
- Motif enrichment in differential peaks
- TF binding site scanning (JASPAR database)
- TF-target gene prediction
- TF + histone mark co-occupancy analysis

Pure Python implementation using:
- gimmemotifs or MOODS for motif scanning
- JASPAR database via pyjaspar or local PWMs
- Custom integration algorithms
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Set
from dataclasses import dataclass, field
import numpy as np
import pandas as pd
from collections import defaultdict
import re

logger = logging.getLogger(__name__)

# Optional imports
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    logger.warning("Biopython not available - install with: pip install biopython")

try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class MotifMatch:
    """A TF motif match in a sequence."""
    tf_name: str
    motif_id: str
    chrom: str
    start: int
    end: int
    strand: str
    score: float
    p_value: float
    sequence: str


@dataclass
class TFBindingSite:
    """A predicted TF binding site."""
    tf_name: str
    chrom: str
    start: int
    end: int
    strand: str
    score: float
    peak_id: Optional[str] = None
    nearest_gene: Optional[str] = None
    distance_to_tss: Optional[int] = None


@dataclass
class MotifEnrichmentResult:
    """Result of motif enrichment analysis."""
    tf_name: str
    motif_id: str
    target_count: int
    target_percent: float
    background_count: int
    background_percent: float
    fold_enrichment: float
    p_value: float
    adjusted_p_value: float
    target_sequences: List[str] = field(default_factory=list)


@dataclass
class TFTargetPrediction:
    """Predicted TF-target gene relationship."""
    tf_name: str
    target_gene: str
    binding_score: float
    distance_to_tss: int
    peak_id: str
    histone_context: Optional[Dict[str, str]] = None  # e.g., {"H3K27ac": "gained"}
    confidence: str = "medium"  # low, medium, high


# =============================================================================
# Position Weight Matrix (PWM) Handling
# =============================================================================

class PositionWeightMatrix:
    """
    Position Weight Matrix for TF motif representation.

    Supports:
    - JASPAR format
    - MEME format
    - Custom matrix format
    """

    def __init__(
        self,
        matrix: np.ndarray,
        name: str,
        motif_id: str = None,
        alphabet: str = "ACGT"
    ):
        """
        Initialize PWM.

        Args:
            matrix: 4 x n numpy array (rows: A, C, G, T; cols: positions)
            name: TF name
            motif_id: Motif identifier (e.g., JASPAR ID)
            alphabet: Nucleotide alphabet
        """
        self.matrix = matrix
        self.name = name
        self.motif_id = motif_id or name
        self.alphabet = alphabet
        self.length = matrix.shape[1]

        # Convert to log-odds (PWM) if needed
        self._pwm = self._to_log_odds()

    def _to_log_odds(self, pseudocount: float = 0.01) -> np.ndarray:
        """Convert frequency matrix to log-odds PWM."""
        # Normalize columns to sum to 1
        col_sums = self.matrix.sum(axis=0, keepdims=True)
        freq = (self.matrix + pseudocount) / (col_sums + 4 * pseudocount)

        # Log-odds against background (0.25 each)
        return np.log2(freq / 0.25)

    def score_sequence(self, sequence: str) -> float:
        """
        Score a sequence against this PWM.

        Args:
            sequence: DNA sequence (must be same length as motif)

        Returns:
            Log-odds score
        """
        if len(sequence) != self.length:
            return float('-inf')

        sequence = sequence.upper()
        base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

        score = 0.0
        for i, base in enumerate(sequence):
            if base in base_to_idx:
                score += self._pwm[base_to_idx[base], i]
            else:
                score += 0  # N or other ambiguous bases

        return score

    def max_score(self) -> float:
        """Get maximum possible score."""
        return self._pwm.max(axis=0).sum()

    def min_score(self) -> float:
        """Get minimum possible score."""
        return self._pwm.min(axis=0).sum()

    @classmethod
    def from_jaspar(cls, jaspar_entry: str) -> 'PositionWeightMatrix':
        """
        Parse PWM from JASPAR format.

        Format:
        >MA0001.1 AGL3
        A [ 0 3 79 ... ]
        C [ 94 75 4 ... ]
        G [ 1 0 3 ... ]
        T [ 2 19 11 ... ]
        """
        lines = jaspar_entry.strip().split('\n')

        # Parse header
        header = lines[0]
        if header.startswith('>'):
            parts = header[1:].split()
            motif_id = parts[0] if parts else "unknown"
            name = parts[1] if len(parts) > 1 else motif_id
        else:
            motif_id = "unknown"
            name = "unknown"

        # Parse matrix
        matrix_lines = [l for l in lines[1:] if l.strip()]
        matrix = []

        for line in matrix_lines:
            # Extract numbers from line like "A [ 1 2 3 ]"
            numbers = re.findall(r'[\d.]+', line)
            if numbers:
                matrix.append([float(n) for n in numbers])

        return cls(np.array(matrix), name, motif_id)

    def to_consensus(self) -> str:
        """Get consensus sequence."""
        bases = ['A', 'C', 'G', 'T']
        consensus = ""
        for i in range(self.length):
            max_idx = self.matrix[:, i].argmax()
            consensus += bases[max_idx]
        return consensus


# =============================================================================
# JASPAR Database
# =============================================================================

class JASPARDatabase:
    """
    Interface to JASPAR transcription factor motif database.

    Provides access to TF motifs without requiring external tools.
    """

    # Common TF motifs (subset of JASPAR CORE vertebrates)
    # Format: {motif_id: (name, consensus, matrix)}
    BUILTIN_MOTIFS = {
        # AP-1 family
        "MA0476.1": ("FOS", "TGAGTCA", np.array([
            [0.1, 0.0, 0.8, 0.0, 0.1, 0.0, 0.8],  # A
            [0.1, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0],  # C
            [0.1, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0],  # G
            [0.7, 0.1, 0.1, 1.0, 0.9, 0.1, 0.2],  # T
        ])),
        "MA0099.3": ("JUN", "TGAGTCA", np.array([
            [0.1, 0.0, 0.8, 0.0, 0.1, 0.0, 0.8],
            [0.1, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0],
            [0.1, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0],
            [0.7, 0.1, 0.1, 1.0, 0.9, 0.1, 0.2],
        ])),
        # NF-kB
        "MA0061.1": ("NF-kappaB", "GGGAATTTCC", np.array([
            [0.0, 0.0, 0.0, 0.9, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.9, 0.9],
            [1.0, 1.0, 1.0, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.9, 0.9, 0.1, 0.1],
        ])),
        # GATA factors
        "MA0035.4": ("GATA1", "AGATAA", np.array([
            [0.8, 0.0, 0.9, 0.0, 0.9, 0.9],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.1, 0.9, 0.1, 0.0, 0.0, 0.0],
            [0.1, 0.1, 0.0, 1.0, 0.1, 0.1],
        ])),
        "MA0036.3": ("GATA2", "AGATAA", np.array([
            [0.8, 0.0, 0.9, 0.0, 0.9, 0.9],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.1, 0.9, 0.1, 0.0, 0.0, 0.0],
            [0.1, 0.1, 0.0, 1.0, 0.1, 0.1],
        ])),
        # ETS factors
        "MA0098.3": ("ETS1", "GGAA", np.array([
            [0.0, 0.0, 0.9, 0.9],
            [0.0, 0.0, 0.0, 0.0],
            [1.0, 1.0, 0.1, 0.1],
            [0.0, 0.0, 0.0, 0.0],
        ])),
        # MEF2
        "MA0052.4": ("MEF2A", "CTAAAAATAG", np.array([
            [0.0, 0.0, 0.9, 0.9, 0.9, 0.9, 0.9, 0.0, 0.9, 0.0],
            [0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9],
            [0.1, 0.9, 0.1, 0.1, 0.1, 0.1, 0.1, 0.9, 0.1, 0.1],
        ])),
        # SP1
        "MA0079.4": ("SP1", "GGGCGG", np.array([
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.9, 0.0, 0.0],
            [1.0, 1.0, 1.0, 0.1, 1.0, 1.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ])),
        # CREB/ATF
        "MA0018.3": ("CREB1", "TGACGTCA", np.array([
            [0.0, 0.0, 0.9, 0.0, 0.0, 0.0, 0.0, 0.9],
            [0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.9, 0.0],
            [0.0, 0.9, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0],
            [1.0, 0.1, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1],
        ])),
        # NRF2
        "MA0150.2": ("NRF2", "TGACNNNGC", np.array([
            [0.0, 0.0, 0.9, 0.0, 0.25, 0.25, 0.25, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.9, 0.25, 0.25, 0.25, 0.0, 0.9],
            [0.0, 0.9, 0.0, 0.0, 0.25, 0.25, 0.25, 0.9, 0.0],
            [1.0, 0.1, 0.1, 0.1, 0.25, 0.25, 0.25, 0.1, 0.1],
        ])),
        # p53
        "MA0106.3": ("TP53", "CATGCCCGGGCATG", np.array([
            [0.0, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0],
            [0.9, 0.0, 0.0, 0.0, 0.9, 0.9, 0.9, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0, 0.9, 0.9, 0.9, 0.0, 0.0, 0.0, 0.9],
            [0.1, 0.1, 0.9, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.9, 0.1],
        ])),
    }

    def __init__(self, custom_motifs_file: Optional[str] = None):
        """
        Initialize JASPAR database.

        Args:
            custom_motifs_file: Path to custom motifs in JASPAR format
        """
        self.motifs: Dict[str, PositionWeightMatrix] = {}

        # Load built-in motifs
        self._load_builtin_motifs()

        # Load custom motifs if provided
        if custom_motifs_file:
            self.load_motifs_file(custom_motifs_file)

    def _load_builtin_motifs(self):
        """Load built-in common TF motifs."""
        for motif_id, (name, consensus, matrix) in self.BUILTIN_MOTIFS.items():
            self.motifs[motif_id] = PositionWeightMatrix(matrix, name, motif_id)

        logger.info(f"Loaded {len(self.motifs)} built-in TF motifs")

    def load_motifs_file(self, filepath: str):
        """Load motifs from JASPAR format file."""
        with open(filepath) as f:
            content = f.read()

        # Split by motif entries
        entries = content.split('>')[1:]  # Skip empty first element

        for entry in entries:
            entry = '>' + entry
            try:
                pwm = PositionWeightMatrix.from_jaspar(entry)
                self.motifs[pwm.motif_id] = pwm
            except Exception as e:
                logger.warning(f"Failed to parse motif: {e}")

        logger.info(f"Loaded {len(entries)} motifs from {filepath}")

    def get_motif(self, motif_id: str) -> Optional[PositionWeightMatrix]:
        """Get a specific motif by ID."""
        return self.motifs.get(motif_id)

    def search_by_name(self, name: str) -> List[PositionWeightMatrix]:
        """Search motifs by TF name (case-insensitive)."""
        name_lower = name.lower()
        return [m for m in self.motifs.values() if name_lower in m.name.lower()]

    def list_motifs(self) -> List[Tuple[str, str]]:
        """List all available motifs."""
        return [(m.motif_id, m.name) for m in self.motifs.values()]


# =============================================================================
# Motif Scanner
# =============================================================================

class MotifScanner:
    """
    Scan sequences for TF binding motifs.

    Provides efficient motif scanning with configurable thresholds.
    """

    def __init__(
        self,
        motif_db: JASPARDatabase,
        score_threshold: float = 0.8,  # Fraction of max score
        p_value_threshold: float = 1e-4
    ):
        """
        Initialize motif scanner.

        Args:
            motif_db: JASPAR database instance
            score_threshold: Minimum score as fraction of max possible
            p_value_threshold: P-value cutoff for reporting matches
        """
        self.db = motif_db
        self.score_threshold = score_threshold
        self.p_value_threshold = p_value_threshold

    def scan_sequence(
        self,
        sequence: str,
        motif_ids: Optional[List[str]] = None
    ) -> List[MotifMatch]:
        """
        Scan a sequence for all motifs.

        Args:
            sequence: DNA sequence to scan
            motif_ids: Specific motifs to scan for (None = all)

        Returns:
            List of MotifMatch objects
        """
        matches = []
        sequence = sequence.upper()

        motifs_to_scan = (
            [self.db.motifs[mid] for mid in motif_ids if mid in self.db.motifs]
            if motif_ids else list(self.db.motifs.values())
        )

        for motif in motifs_to_scan:
            matches.extend(self._scan_single_motif(sequence, motif))

        return matches

    def _scan_single_motif(
        self,
        sequence: str,
        motif: PositionWeightMatrix
    ) -> List[MotifMatch]:
        """Scan for a single motif."""
        matches = []

        max_score = motif.max_score()
        min_acceptable = max_score * self.score_threshold

        # Scan forward strand
        for i in range(len(sequence) - motif.length + 1):
            subseq = sequence[i:i + motif.length]
            score = motif.score_sequence(subseq)

            if score >= min_acceptable:
                # Estimate p-value (simplified)
                p_value = 10 ** (-(score / max_score) * 4)  # Rough approximation

                if p_value <= self.p_value_threshold:
                    matches.append(MotifMatch(
                        tf_name=motif.name,
                        motif_id=motif.motif_id,
                        chrom="",  # Set by caller
                        start=i,
                        end=i + motif.length,
                        strand="+",
                        score=score,
                        p_value=p_value,
                        sequence=subseq
                    ))

        # Scan reverse strand
        rev_comp = self._reverse_complement(sequence)
        for i in range(len(rev_comp) - motif.length + 1):
            subseq = rev_comp[i:i + motif.length]
            score = motif.score_sequence(subseq)

            if score >= min_acceptable:
                p_value = 10 ** (-(score / max_score) * 4)

                if p_value <= self.p_value_threshold:
                    # Adjust position for reverse strand
                    real_start = len(sequence) - i - motif.length
                    matches.append(MotifMatch(
                        tf_name=motif.name,
                        motif_id=motif.motif_id,
                        chrom="",
                        start=real_start,
                        end=real_start + motif.length,
                        strand="-",
                        score=score,
                        p_value=p_value,
                        sequence=subseq
                    ))

        return matches

    def _reverse_complement(self, seq: str) -> str:
        """Get reverse complement of sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(b, 'N') for b in reversed(seq.upper()))


# =============================================================================
# Motif Enrichment Analysis
# =============================================================================

class MotifEnrichmentAnalyzer:
    """
    Analyze motif enrichment in peak sets.

    Compares motif occurrence in target peaks vs background.
    """

    def __init__(
        self,
        motif_db: JASPARDatabase,
        genome_fasta: Optional[str] = None
    ):
        """
        Initialize enrichment analyzer.

        Args:
            motif_db: JASPAR database
            genome_fasta: Path to reference genome FASTA
        """
        self.db = motif_db
        self.scanner = MotifScanner(motif_db)
        self.genome_fasta = genome_fasta
        self._genome = None

    def _load_genome(self):
        """Load reference genome."""
        if self._genome is None and self.genome_fasta:
            if BIOPYTHON_AVAILABLE:
                self._genome = SeqIO.to_dict(SeqIO.parse(self.genome_fasta, "fasta"))
            else:
                # Simple FASTA parser fallback
                self._genome = self._parse_fasta_simple(self.genome_fasta)

    def _parse_fasta_simple(self, fasta_path: str) -> Dict[str, str]:
        """Simple FASTA parser without Biopython."""
        sequences = {}
        current_id = None
        current_seq = []

        with open(fasta_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        sequences[current_id] = ''.join(current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)

        if current_id:
            sequences[current_id] = ''.join(current_seq)

        return sequences

    def get_peak_sequences(
        self,
        peaks: pd.DataFrame,
        chrom_col: str = 'chrom',
        start_col: str = 'start',
        end_col: str = 'end'
    ) -> List[str]:
        """Extract sequences for peaks from reference genome."""
        self._load_genome()

        sequences = []
        for _, peak in peaks.iterrows():
            chrom = str(peak[chrom_col])
            start = int(peak[start_col])
            end = int(peak[end_col])

            if self._genome and chrom in self._genome:
                if BIOPYTHON_AVAILABLE:
                    seq = str(self._genome[chrom].seq[start:end])
                else:
                    seq = self._genome[chrom][start:end]
                sequences.append(seq.upper())
            else:
                sequences.append("")

        return sequences

    def analyze_enrichment(
        self,
        target_sequences: List[str],
        background_sequences: List[str],
        motif_ids: Optional[List[str]] = None
    ) -> List[MotifEnrichmentResult]:
        """
        Analyze motif enrichment in target vs background.

        Args:
            target_sequences: Sequences from target peaks
            background_sequences: Background sequences
            motif_ids: Specific motifs to test (None = all)

        Returns:
            List of enrichment results
        """
        results = []

        motifs = (
            [self.db.motifs[mid] for mid in motif_ids if mid in self.db.motifs]
            if motif_ids else list(self.db.motifs.values())
        )

        for motif in motifs:
            result = self._test_single_motif(motif, target_sequences, background_sequences)
            results.append(result)

        # Adjust p-values for multiple testing (Benjamini-Hochberg)
        results = self._adjust_pvalues(results)

        # Sort by adjusted p-value
        results.sort(key=lambda x: x.adjusted_p_value)

        return results

    def _test_single_motif(
        self,
        motif: PositionWeightMatrix,
        target_seqs: List[str],
        bg_seqs: List[str]
    ) -> MotifEnrichmentResult:
        """Test enrichment for a single motif."""
        # Count sequences with motif
        target_with_motif = 0
        target_matches = []

        for seq in target_seqs:
            if seq:
                matches = self.scanner._scan_single_motif(seq, motif)
                if matches:
                    target_with_motif += 1
                    target_matches.append(seq)

        bg_with_motif = 0
        for seq in bg_seqs:
            if seq:
                matches = self.scanner._scan_single_motif(seq, motif)
                if matches:
                    bg_with_motif += 1

        # Calculate percentages
        target_pct = (target_with_motif / len(target_seqs) * 100) if target_seqs else 0
        bg_pct = (bg_with_motif / len(bg_seqs) * 100) if bg_seqs else 0

        # Fold enrichment
        fold_enrichment = (target_pct / bg_pct) if bg_pct > 0 else float('inf')

        # Fisher's exact test
        from scipy import stats
        contingency = [
            [target_with_motif, len(target_seqs) - target_with_motif],
            [bg_with_motif, len(bg_seqs) - bg_with_motif]
        ]
        _, p_value = stats.fisher_exact(contingency, alternative='greater')

        return MotifEnrichmentResult(
            tf_name=motif.name,
            motif_id=motif.motif_id,
            target_count=target_with_motif,
            target_percent=target_pct,
            background_count=bg_with_motif,
            background_percent=bg_pct,
            fold_enrichment=fold_enrichment,
            p_value=p_value,
            adjusted_p_value=p_value,  # Will be adjusted later
            target_sequences=target_matches[:10]  # Keep top 10 examples
        )

    def _adjust_pvalues(
        self,
        results: List[MotifEnrichmentResult]
    ) -> List[MotifEnrichmentResult]:
        """Benjamini-Hochberg p-value adjustment."""
        n = len(results)
        if n == 0:
            return results

        # Sort by p-value
        sorted_results = sorted(results, key=lambda x: x.p_value)

        # Adjust
        for i, result in enumerate(sorted_results):
            result.adjusted_p_value = min(1.0, result.p_value * n / (i + 1))

        # Enforce monotonicity
        for i in range(n - 2, -1, -1):
            sorted_results[i].adjusted_p_value = min(
                sorted_results[i].adjusted_p_value,
                sorted_results[i + 1].adjusted_p_value
            )

        return results


# =============================================================================
# TF-Target Gene Prediction
# =============================================================================

class TFTargetPredictor:
    """
    Predict TF-target gene relationships.

    Integrates:
    - TF binding site locations
    - Distance to TSS
    - Histone mark context
    - Expression data (optional)
    """

    def __init__(
        self,
        max_distance: int = 100000,  # Max distance from TSS
        promoter_distance: int = 3000
    ):
        """
        Initialize predictor.

        Args:
            max_distance: Maximum distance to consider as potential target
            promoter_distance: Distance defining promoter region
        """
        self.max_distance = max_distance
        self.promoter_distance = promoter_distance

    def predict_targets(
        self,
        tf_peaks: pd.DataFrame,
        gene_annotation: pd.DataFrame,
        histone_data: Optional[Dict[str, pd.DataFrame]] = None,
        expression_data: Optional[pd.DataFrame] = None
    ) -> List[TFTargetPrediction]:
        """
        Predict TF target genes.

        Args:
            tf_peaks: DataFrame with TF binding peaks
            gene_annotation: DataFrame with gene info (gene_name, chrom, tss)
            histone_data: Dict of histone mark DataFrames keyed by mark name
            expression_data: Optional expression/DE data

        Returns:
            List of TF-target predictions
        """
        predictions = []

        for _, peak in tf_peaks.iterrows():
            peak_chrom = peak.get('chrom', peak.get('Chromosome'))
            peak_center = (peak.get('start', 0) + peak.get('end', 0)) // 2
            peak_id = peak.get('peak_id', f"{peak_chrom}:{peak.get('start')}-{peak.get('end')}")
            tf_name = peak.get('tf_name', 'Unknown_TF')

            # Find nearby genes
            nearby_genes = gene_annotation[
                (gene_annotation['chrom'] == peak_chrom) &
                (abs(gene_annotation['tss'] - peak_center) <= self.max_distance)
            ]

            for _, gene in nearby_genes.iterrows():
                distance = peak_center - gene['tss']

                # Get histone context if available
                histone_context = None
                if histone_data:
                    histone_context = self._get_histone_context(
                        peak_chrom, peak.get('start'), peak.get('end'),
                        histone_data
                    )

                # Calculate confidence
                confidence = self._calculate_confidence(
                    distance, histone_context, expression_data, gene['gene_name']
                )

                predictions.append(TFTargetPrediction(
                    tf_name=tf_name,
                    target_gene=gene['gene_name'],
                    binding_score=peak.get('score', 0),
                    distance_to_tss=distance,
                    peak_id=peak_id,
                    histone_context=histone_context,
                    confidence=confidence
                ))

        return predictions

    def _get_histone_context(
        self,
        chrom: str,
        start: int,
        end: int,
        histone_data: Dict[str, pd.DataFrame]
    ) -> Dict[str, str]:
        """Get histone mark status at TF binding site."""
        context = {}

        for mark, data in histone_data.items():
            # Check for overlapping significant peaks
            overlapping = data[
                (data['chrom'] == chrom) &
                (data['start'] <= end) &
                (data['end'] >= start) &
                (data.get('significant', True))
            ]

            if len(overlapping) > 0:
                direction = overlapping['direction'].iloc[0] if 'direction' in overlapping.columns else 'present'
                context[mark] = direction

        return context if context else None

    def _calculate_confidence(
        self,
        distance: int,
        histone_context: Optional[Dict],
        expression_data: Optional[pd.DataFrame],
        gene_name: str
    ) -> str:
        """Calculate prediction confidence."""
        score = 0

        # Distance score
        if abs(distance) <= self.promoter_distance:
            score += 3  # Promoter binding
        elif abs(distance) <= 10000:
            score += 2  # Proximal
        else:
            score += 1  # Distal

        # Histone context score
        if histone_context:
            if 'H3K27ac' in histone_context:
                score += 2  # Active mark present
            if 'H3K4me1' in histone_context:
                score += 1  # Enhancer mark

        # Expression concordance
        if expression_data is not None and gene_name in expression_data.index:
            if expression_data.loc[gene_name, 'significant']:
                score += 2

        # Convert to confidence level
        if score >= 6:
            return "high"
        elif score >= 3:
            return "medium"
        else:
            return "low"


# =============================================================================
# TF + Histone Integration
# =============================================================================

class TFHistoneIntegrator:
    """
    Integrate TF binding with histone modifications.

    Analyzes:
    - Co-occupancy of TFs and histone marks
    - TF binding at differential histone mark regions
    - Regulatory logic (TF + marks → gene expression)
    """

    def __init__(self):
        """Initialize integrator."""
        pass

    def analyze_cooccupancy(
        self,
        tf_peaks: pd.DataFrame,
        histone_peaks: Dict[str, pd.DataFrame],
        overlap_fraction: float = 0.5
    ) -> pd.DataFrame:
        """
        Analyze co-occupancy of TF and histone marks.

        Args:
            tf_peaks: TF binding peaks
            histone_peaks: Dict of histone mark peaks
            overlap_fraction: Minimum overlap fraction

        Returns:
            DataFrame with co-occupancy statistics
        """
        results = []

        for mark_name, mark_peaks in histone_peaks.items():
            # Count overlaps
            overlaps = self._count_overlaps(tf_peaks, mark_peaks, overlap_fraction)

            total_tf = len(tf_peaks)
            total_mark = len(mark_peaks)

            results.append({
                'histone_mark': mark_name,
                'tf_peaks_total': total_tf,
                'mark_peaks_total': total_mark,
                'overlapping': overlaps,
                'tf_overlap_pct': (overlaps / total_tf * 100) if total_tf > 0 else 0,
                'mark_overlap_pct': (overlaps / total_mark * 100) if total_mark > 0 else 0
            })

        return pd.DataFrame(results)

    def _count_overlaps(
        self,
        peaks1: pd.DataFrame,
        peaks2: pd.DataFrame,
        min_overlap: float = 0.5
    ) -> int:
        """Count overlapping peaks."""
        count = 0

        # Simple O(n*m) overlap - could be optimized with interval trees
        for _, p1 in peaks1.iterrows():
            chrom1 = p1.get('chrom', p1.get('Chromosome'))
            start1 = p1.get('start', p1.get('Start', 0))
            end1 = p1.get('end', p1.get('End', 0))
            len1 = end1 - start1

            for _, p2 in peaks2.iterrows():
                chrom2 = p2.get('chrom', p2.get('Chromosome'))
                if chrom1 != chrom2:
                    continue

                start2 = p2.get('start', p2.get('Start', 0))
                end2 = p2.get('end', p2.get('End', 0))

                # Calculate overlap
                overlap_start = max(start1, start2)
                overlap_end = min(end1, end2)

                if overlap_end > overlap_start:
                    overlap_len = overlap_end - overlap_start
                    if overlap_len / len1 >= min_overlap:
                        count += 1
                        break

        return count

    def find_tf_at_differential_marks(
        self,
        tf_peaks: pd.DataFrame,
        diff_histone_peaks: pd.DataFrame,
        direction: str = 'both'  # 'gained', 'lost', or 'both'
    ) -> pd.DataFrame:
        """
        Find TF binding sites at differential histone mark regions.

        Args:
            tf_peaks: TF binding peaks
            diff_histone_peaks: Differential histone peaks with 'direction' column
            direction: Filter by direction of histone change

        Returns:
            TF peaks overlapping differential histone regions
        """
        if direction != 'both':
            diff_peaks = diff_histone_peaks[
                diff_histone_peaks['direction'].str.lower() == direction.lower()
            ]
        else:
            diff_peaks = diff_histone_peaks

        # Find overlapping TF peaks
        overlapping_tf = []

        for _, tf_peak in tf_peaks.iterrows():
            tf_chrom = tf_peak.get('chrom')
            tf_start = tf_peak.get('start', 0)
            tf_end = tf_peak.get('end', 0)

            for _, h_peak in diff_peaks.iterrows():
                h_chrom = h_peak.get('chrom')
                if tf_chrom != h_chrom:
                    continue

                h_start = h_peak.get('start', 0)
                h_end = h_peak.get('end', 0)

                # Check overlap
                if tf_start < h_end and tf_end > h_start:
                    result = tf_peak.to_dict()
                    result['histone_peak_id'] = h_peak.get('peak_id')
                    result['histone_direction'] = h_peak.get('direction')
                    result['histone_log2FC'] = h_peak.get('log2FC', h_peak.get('Fold'))
                    overlapping_tf.append(result)
                    break

        return pd.DataFrame(overlapping_tf)


# =============================================================================
# Convenience Functions
# =============================================================================

def run_tf_motif_analysis(
    peaks_df: pd.DataFrame,
    genome_fasta: str,
    background_peaks: Optional[pd.DataFrame] = None,
    custom_motifs: Optional[str] = None
) -> Tuple[List[MotifEnrichmentResult], pd.DataFrame]:
    """
    Convenience function to run complete TF motif analysis.

    Args:
        peaks_df: Target peaks DataFrame
        genome_fasta: Path to reference genome FASTA
        background_peaks: Background peaks (optional, will generate if None)
        custom_motifs: Path to custom motifs file

    Returns:
        Tuple of (enrichment results, peaks with motif annotations)
    """
    # Initialize database
    db = JASPARDatabase(custom_motifs)

    # Initialize analyzer
    analyzer = MotifEnrichmentAnalyzer(db, genome_fasta)

    # Get sequences
    target_seqs = analyzer.get_peak_sequences(peaks_df)

    # Get or generate background
    if background_peaks is not None:
        bg_seqs = analyzer.get_peak_sequences(background_peaks)
    else:
        # Simple background: shuffle target sequences
        import random
        bg_seqs = [''.join(random.sample(s, len(s))) for s in target_seqs if s]

    # Run enrichment
    enrichment_results = analyzer.analyze_enrichment(target_seqs, bg_seqs)

    # Scan all peaks for motifs
    scanner = MotifScanner(db)
    peaks_annotated = peaks_df.copy()

    motif_annotations = []
    for seq in target_seqs:
        if seq:
            matches = scanner.scan_sequence(seq)
            tf_names = list(set(m.tf_name for m in matches))
            motif_annotations.append(';'.join(tf_names) if tf_names else '')
        else:
            motif_annotations.append('')

    peaks_annotated['tf_motifs'] = motif_annotations

    return enrichment_results, peaks_annotated
