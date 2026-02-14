"""
Shared genomic utilities for EpiNexus.

Provides high-performance interval overlap detection using NCLS (Nested
Containment List) or sorted-endpoint sweepline, replacing O(n*m) nested
loops throughout the codebase.

Also contains shared helpers for peak file parsing, column name detection,
and chromosome sorting.
"""

import logging
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Try to import NCLS for fast interval overlap; fall back to sweep-line
# ---------------------------------------------------------------------------
try:
    from ncls import NCLS

    _HAS_NCLS = True
except ImportError:
    _HAS_NCLS = False
    logger.info("ncls not installed – using sweep-line overlap algorithm")

try:
    import pyranges as pr

    _HAS_PYRANGES = True
except ImportError:
    _HAS_PYRANGES = False


# ============================================================================
# Core overlap functions
# ============================================================================


def _build_ncls_index(
    starts: np.ndarray, ends: np.ndarray
) -> "NCLS":
    """Build an NCLS index from start/end arrays."""
    ids = np.arange(len(starts), dtype=np.int64)
    return NCLS(
        starts.astype(np.int64),
        ends.astype(np.int64),
        ids,
    )


def find_overlaps(
    query_df: pd.DataFrame,
    subject_df: pd.DataFrame,
    chrom_col: str = "chr",
    start_col: str = "start",
    end_col: str = "end",
    min_overlap_bp: int = 1,
    min_overlap_frac: float = 0.0,
    report: str = "all",
) -> pd.DataFrame:
    """Find overlapping intervals between two DataFrames.

    Replaces all O(n*m) nested-loop overlap detection in the codebase with
    an O(n log n) algorithm (NCLS when available, sweep-line otherwise).

    Parameters
    ----------
    query_df : pd.DataFrame
        Query intervals (the "left" set).
    subject_df : pd.DataFrame
        Subject intervals (the "right" set to search against).
    chrom_col : str
        Column name for chromosome in both DataFrames.
    start_col, end_col : str
        Column names for interval boundaries.
    min_overlap_bp : int
        Minimum overlap in base pairs (default 1).
    min_overlap_frac : float
        Minimum overlap as a fraction of the *query* interval length
        (0.0–1.0, default 0.0 = any overlap).
    report : str
        "all" – return all overlapping pairs.
        "first" – return only the first hit per query.
        "count" – return a count of overlaps per query.

    Returns
    -------
    pd.DataFrame
        If report="all" or "first": DataFrame with columns
            [query_idx, subject_idx, overlap_bp].
        If report="count": DataFrame with columns
            [query_idx, count].
    """
    if query_df.empty or subject_df.empty:
        if report == "count":
            return pd.DataFrame({"query_idx": pd.RangeIndex(len(query_df)), "count": 0})
        return pd.DataFrame(columns=["query_idx", "subject_idx", "overlap_bp"])

    results: List[Tuple[int, int, int]] = []

    # Group by chromosome for both sets
    query_groups = query_df.groupby(chrom_col)
    subject_groups = {name: grp for name, grp in subject_df.groupby(chrom_col)}

    for chrom, q_grp in query_groups:
        if chrom not in subject_groups:
            continue
        s_grp = subject_groups[chrom]

        q_starts = q_grp[start_col].values
        q_ends = q_grp[end_col].values
        q_indices = q_grp.index.values

        s_starts = s_grp[start_col].values
        s_ends = s_grp[end_col].values
        s_indices = s_grp.index.values

        if _HAS_NCLS:
            ncls = _build_ncls_index(s_starts, s_ends)
            for i in range(len(q_starts)):
                qs, qe = int(q_starts[i]), int(q_ends[i])
                q_len = qe - qs
                # NCLS.find_overlap returns an iterator of (start, end, id) tuples
                for hit in ncls.find_overlap(qs, qe):
                    _s_start, _s_end, s_local_idx = hit
                    ovlp = min(qe, int(_s_end)) - max(qs, int(_s_start))
                    if ovlp < min_overlap_bp:
                        continue
                    if min_overlap_frac > 0 and q_len > 0 and ovlp / q_len < min_overlap_frac:
                        continue
                    results.append((int(q_indices[i]), int(s_indices[int(s_local_idx)]), ovlp))
                    if report == "first":
                        break
        else:
            # Sweep-line fallback – still O(n log n) via sorting
            results.extend(
                _sweepline_overlaps(
                    q_starts, q_ends, q_indices,
                    s_starts, s_ends, s_indices,
                    min_overlap_bp, min_overlap_frac, report,
                )
            )

    if report == "count":
        if results:
            df = pd.DataFrame(results, columns=["query_idx", "subject_idx", "overlap_bp"])
            counts = df.groupby("query_idx").size().reset_index(name="count")
        else:
            counts = pd.DataFrame({"query_idx": pd.Series(dtype=int), "count": pd.Series(dtype=int)})
        # Fill in zeros for queries with no overlaps
        all_q = pd.DataFrame({"query_idx": query_df.index})
        return all_q.merge(counts, on="query_idx", how="left").fillna(0).astype({"count": int})

    if not results:
        return pd.DataFrame(columns=["query_idx", "subject_idx", "overlap_bp"])

    return pd.DataFrame(results, columns=["query_idx", "subject_idx", "overlap_bp"])


def _sweepline_overlaps(
    q_starts, q_ends, q_indices,
    s_starts, s_ends, s_indices,
    min_overlap_bp, min_overlap_frac, report,
) -> List[Tuple[int, int, int]]:
    """Sweep-line overlap detection (fallback when NCLS unavailable)."""
    results = []
    # Sort subject by start coordinate
    s_order = np.argsort(s_starts)
    s_starts_sorted = s_starts[s_order]
    s_ends_sorted = s_ends[s_order]
    s_indices_sorted = s_indices[s_order]

    for i in range(len(q_starts)):
        qs, qe = int(q_starts[i]), int(q_ends[i])
        q_len = qe - qs

        # Binary search for first possible overlapping subject
        lo = np.searchsorted(s_starts_sorted, qs - (s_ends_sorted.max() - s_starts_sorted.min() if len(s_starts_sorted) > 0 else 0), side="left")
        lo = max(0, lo)

        for j in range(lo, len(s_starts_sorted)):
            ss = int(s_starts_sorted[j])
            if ss >= qe:
                break
            se = int(s_ends_sorted[j])
            ovlp = min(qe, se) - max(qs, ss)
            if ovlp < min_overlap_bp:
                continue
            if min_overlap_frac > 0 and q_len > 0 and ovlp / q_len < min_overlap_frac:
                continue
            results.append((int(q_indices[i]), int(s_indices_sorted[j]), ovlp))
            if report == "first":
                break
    return results


# ============================================================================
# Convenience wrappers replacing specific O(n²) patterns
# ============================================================================


def count_overlaps(
    query_df: pd.DataFrame,
    subject_df: pd.DataFrame,
    chrom_col: str = "chr",
    start_col: str = "start",
    end_col: str = "end",
    min_overlap_bp: int = 1,
) -> np.ndarray:
    """Count overlaps per query interval.

    Direct replacement for atacseq._count_overlaps() and similar patterns.

    Returns
    -------
    np.ndarray
        Array of length len(query_df) with overlap counts.
    """
    counts_df = find_overlaps(
        query_df, subject_df,
        chrom_col=chrom_col, start_col=start_col, end_col=end_col,
        min_overlap_bp=min_overlap_bp,
        report="count",
    )
    result = np.zeros(len(query_df), dtype=int)
    for _, row in counts_df.iterrows():
        idx = query_df.index.get_loc(row["query_idx"])
        result[idx] = row["count"]
    return result


def count_overlaps_with_fraction(
    query_df: pd.DataFrame,
    subject_df: pd.DataFrame,
    min_overlap_frac: float = 0.5,
    chrom_col: str = "chr",
    start_col: str = "start",
    end_col: str = "end",
) -> int:
    """Count how many query intervals have at least one subject overlap.

    Direct replacement for transcription_factors._count_overlaps().
    """
    hits = find_overlaps(
        query_df, subject_df,
        chrom_col=chrom_col, start_col=start_col, end_col=end_col,
        min_overlap_frac=min_overlap_frac,
        report="first",
    )
    return len(hits)


def find_first_overlaps(
    query_df: pd.DataFrame,
    subject_df: pd.DataFrame,
    chrom_col: str = "chr",
    start_col: str = "start",
    end_col: str = "end",
) -> pd.DataFrame:
    """Find the first subject overlap for each query.

    Direct replacement for transcription_factors.find_tf_at_differential_marks().

    Returns
    -------
    pd.DataFrame
        Rows from query_df that have at least one overlap, with the index of
        the matching subject row in column ``subject_idx``.
    """
    return find_overlaps(
        query_df, subject_df,
        chrom_col=chrom_col, start_col=start_col, end_col=end_col,
        report="first",
    )


def exclude_overlapping(
    peaks: pd.DataFrame,
    exclusion_regions: pd.DataFrame,
    chrom_col: str = "chr",
    start_col: str = "start",
    end_col: str = "end",
) -> pd.DataFrame:
    """Remove peaks that overlap any exclusion region.

    Direct replacement for super_enhancers._exclude_tss().
    """
    hits = find_overlaps(
        peaks, exclusion_regions,
        chrom_col=chrom_col, start_col=start_col, end_col=end_col,
        report="first",
    )
    overlapping_indices = set(hits["query_idx"].values)
    keep_mask = ~peaks.index.isin(overlapping_indices)
    return peaks[keep_mask].reset_index(drop=True)


# ============================================================================
# Peak file parsing utilities
# ============================================================================

# Standard column name mappings (used across many pages/modules)
CHROM_COLS = ["chr", "chrom", "chromosome", "seqnames", "#chr"]
START_COLS = ["start", "chromStart", "peak_start"]
END_COLS = ["end", "chromEnd", "peak_end"]
SIGNAL_COLS = ["signal", "signalValue", "score", "fold_enrichment", "enrichment"]
FC_COLS = ["log2FC", "log2FoldChange", "log2fc", "logFC", "Fold"]
PVAL_COLS = ["pvalue", "pval", "p.value", "PValue", "P.Value"]
FDR_COLS = ["FDR", "padj", "p.adj", "adj.P.Val", "q_value", "qvalue"]


def detect_column(df: pd.DataFrame, candidates: List[str], required: bool = False) -> Optional[str]:
    """Find the first matching column name from a list of candidates.

    Parameters
    ----------
    df : pd.DataFrame
    candidates : list of str
        Column names to search for (case-insensitive).
    required : bool
        If True, raise ValueError when not found.

    Returns
    -------
    str or None
    """
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    if required:
        raise ValueError(
            f"Could not find any of {candidates} in columns: {list(df.columns)}"
        )
    return None


def standardize_peak_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename common peak column variants to standardized names.

    Produces columns: chr, start, end (and optionally signal, log2FC, pvalue, FDR).
    """
    mapping = {}
    for std_name, candidates in [
        ("chr", CHROM_COLS),
        ("start", START_COLS),
        ("end", END_COLS),
        ("signal", SIGNAL_COLS),
        ("log2FC", FC_COLS),
        ("pvalue", PVAL_COLS),
        ("FDR", FDR_COLS),
    ]:
        col = detect_column(df, candidates)
        if col and col != std_name:
            mapping[col] = std_name
    return df.rename(columns=mapping)


def load_peak_file(filepath_or_buffer, sep: str = "\t") -> pd.DataFrame:
    """Load a BED/narrowPeak/broadPeak/CSV file into a standardized DataFrame.

    Handles:
    - BED (3-6+ columns, no header)
    - narrowPeak / broadPeak (ENCODE format)
    - CSV/TSV with headers

    Returns a DataFrame with at least: chr, start, end
    """
    import io

    if hasattr(filepath_or_buffer, "read"):
        content = filepath_or_buffer.read()
        if isinstance(content, bytes):
            content = content.decode("utf-8")
        buf = io.StringIO(content)
    else:
        buf = str(filepath_or_buffer)

    # Sniff header: if first line starts with '#' or has alpha characters in column 2
    try:
        peek = pd.read_csv(buf, sep=sep, nrows=2, header=None)
        if hasattr(buf, "seek"):
            buf.seek(0)

        first_val = str(peek.iloc[0, 1]) if peek.shape[1] > 1 else ""
        has_header = not first_val.replace(".", "").replace("-", "").isdigit()
    except Exception:
        has_header = True
        if hasattr(buf, "seek"):
            buf.seek(0)

    df = pd.read_csv(buf, sep=sep, header=0 if has_header else None, comment="#")

    if not has_header:
        # Assign BED-style column names
        bed_cols = ["chr", "start", "end", "name", "score", "strand",
                    "signalValue", "pValue", "qValue", "peak"]
        df.columns = bed_cols[: len(df.columns)]

    df = standardize_peak_columns(df)
    return df


# ============================================================================
# Chromosome utilities
# ============================================================================

_CHROM_ORDER = {f"chr{i}": i for i in range(1, 23)}
_CHROM_ORDER.update({"chrX": 23, "chrY": 24, "chrM": 25, "chrMT": 25})


def sort_chromosomes(chroms: List[str]) -> List[str]:
    """Sort chromosome names in natural order (1,2,...,22,X,Y,M)."""
    def _sort_key(c: str) -> Tuple[int, str]:
        c_stripped = c.replace("chr", "") if c.startswith("chr") else c
        if c in _CHROM_ORDER:
            return (_CHROM_ORDER[c], c)
        try:
            return (int(c_stripped), c)
        except ValueError:
            return (100, c)
    return sorted(chroms, key=_sort_key)


def filter_standard_chroms(df: pd.DataFrame, chrom_col: str = "chr") -> pd.DataFrame:
    """Filter to standard chromosomes (1-22, X, Y), removing random/Un/hap."""
    standard = set(_CHROM_ORDER.keys()) | {str(i) for i in range(1, 23)} | {"X", "Y"}
    return df[df[chrom_col].isin(standard)].copy()
