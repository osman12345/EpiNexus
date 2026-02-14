"""
Unit tests for genomic_utils shared module.

Tests cover:
- Interval overlap detection (NCLS and sweep-line fallback)
- Convenience wrappers (count_overlaps, exclude_overlapping, etc.)
- Peak file parsing and column standardization
- Chromosome utilities
"""

import io

import numpy as np
import pandas as pd
import pytest

from app.core.genomic_utils import (
    find_overlaps,
    count_overlaps,
    count_overlaps_with_fraction,
    find_first_overlaps,
    exclude_overlapping,
    detect_column,
    standardize_peak_columns,
    load_peak_file,
    sort_chromosomes,
    filter_standard_chroms,
    CHROM_COLS,
    _sweepline_overlaps,
)


# ============================================================================
# Core overlap detection
# ============================================================================


class TestFindOverlaps:
    """Tests for the main find_overlaps function."""

    def test_basic_overlap(self):
        """Two intervals on the same chromosome that overlap."""
        query = pd.DataFrame({"chr": ["chr1"], "start": [100], "end": [300]})
        subject = pd.DataFrame({"chr": ["chr1"], "start": [200], "end": [400]})
        result = find_overlaps(query, subject)
        assert len(result) == 1
        assert result.iloc[0]["overlap_bp"] == 100  # min(300,400) - max(100,200)

    def test_no_overlap(self):
        """Intervals on the same chromosome that do not overlap."""
        query = pd.DataFrame({"chr": ["chr1"], "start": [100], "end": [200]})
        subject = pd.DataFrame({"chr": ["chr1"], "start": [300], "end": [400]})
        result = find_overlaps(query, subject)
        assert len(result) == 0

    def test_different_chromosomes(self):
        """Intervals on different chromosomes never overlap."""
        query = pd.DataFrame({"chr": ["chr1"], "start": [100], "end": [300]})
        subject = pd.DataFrame({"chr": ["chr2"], "start": [100], "end": [300]})
        result = find_overlaps(query, subject)
        assert len(result) == 0

    def test_multiple_overlaps(self, overlapping_peaks):
        """Use the shared fixture with known overlaps."""
        query, subject = overlapping_peaks
        result = find_overlaps(query, subject)
        # chr1:100-300 overlaps chr1:250-350 (50 bp)
        assert len(result) >= 1

    def test_min_overlap_bp_filter(self):
        """Filter by minimum overlap base pairs."""
        query = pd.DataFrame({"chr": ["chr1"], "start": [100], "end": [300]})
        subject = pd.DataFrame({"chr": ["chr1"], "start": [290], "end": [400]})
        # Overlap is 10 bp
        result_low = find_overlaps(query, subject, min_overlap_bp=1)
        result_high = find_overlaps(query, subject, min_overlap_bp=50)
        assert len(result_low) == 1
        assert len(result_high) == 0

    def test_min_overlap_frac_filter(self):
        """Filter by minimum overlap fraction."""
        query = pd.DataFrame({"chr": ["chr1"], "start": [100], "end": [200]})  # length 100
        subject = pd.DataFrame({"chr": ["chr1"], "start": [180], "end": [300]})  # overlap 20
        # Fraction = 20/100 = 0.2
        result_low = find_overlaps(query, subject, min_overlap_frac=0.1)
        result_high = find_overlaps(query, subject, min_overlap_frac=0.5)
        assert len(result_low) == 1
        assert len(result_high) == 0

    def test_report_count(self):
        """report='count' returns per-query overlap counts."""
        query = pd.DataFrame(
            {
                "chr": ["chr1", "chr1", "chr2"],
                "start": [100, 500, 100],
                "end": [300, 700, 300],
            }
        )
        subject = pd.DataFrame(
            {
                "chr": ["chr1", "chr1"],
                "start": [200, 600],
                "end": [400, 800],
            }
        )
        result = find_overlaps(query, subject, report="count")
        assert "count" in result.columns
        assert len(result) == 3
        # First query overlaps subject 0, second overlaps subject 1, third has 0
        counts = result.set_index("query_idx")["count"]
        assert counts.loc[0] >= 1
        assert counts.loc[2] == 0

    def test_report_first(self):
        """report='first' returns at most one hit per query."""
        query = pd.DataFrame({"chr": ["chr1"], "start": [100], "end": [500]})
        subject = pd.DataFrame(
            {
                "chr": ["chr1", "chr1"],
                "start": [150, 300],
                "end": [250, 400],
            }
        )
        result = find_overlaps(query, subject, report="first")
        assert len(result) == 1

    def test_empty_query(self):
        """Empty query returns empty result."""
        query = pd.DataFrame(columns=["chr", "start", "end"])
        subject = pd.DataFrame({"chr": ["chr1"], "start": [100], "end": [200]})
        result = find_overlaps(query, subject)
        assert len(result) == 0

    def test_empty_subject(self):
        """Empty subject returns empty result."""
        query = pd.DataFrame({"chr": ["chr1"], "start": [100], "end": [200]})
        subject = pd.DataFrame(columns=["chr", "start", "end"])
        result = find_overlaps(query, subject)
        assert len(result) == 0

    def test_empty_count_report(self):
        """Empty subject with report='count' returns zeros for all queries."""
        query = pd.DataFrame({"chr": ["chr1", "chr2"], "start": [100, 200], "end": [300, 400]})
        subject = pd.DataFrame(columns=["chr", "start", "end"])
        result = find_overlaps(query, subject, report="count")
        assert len(result) == 2
        assert result["count"].sum() == 0


class TestSweeplineOverlaps:
    """Direct tests for the sweep-line fallback algorithm."""

    def test_basic_sweepline(self):
        """Verify sweep-line finds overlaps correctly."""
        q_starts = np.array([100, 500])
        q_ends = np.array([300, 700])
        q_indices = np.array([0, 1])
        s_starts = np.array([200, 600])
        s_ends = np.array([400, 800])
        s_indices = np.array([0, 1])

        results = _sweepline_overlaps(
            q_starts,
            q_ends,
            q_indices,
            s_starts,
            s_ends,
            s_indices,
            min_overlap_bp=1,
            min_overlap_frac=0.0,
            report="all",
        )
        assert len(results) == 2

    def test_sweepline_no_overlap(self):
        """Sweep-line with non-overlapping intervals."""
        q_starts = np.array([100])
        q_ends = np.array([200])
        q_indices = np.array([0])
        s_starts = np.array([300])
        s_ends = np.array([400])
        s_indices = np.array([0])

        results = _sweepline_overlaps(
            q_starts,
            q_ends,
            q_indices,
            s_starts,
            s_ends,
            s_indices,
            min_overlap_bp=1,
            min_overlap_frac=0.0,
            report="all",
        )
        assert len(results) == 0


# ============================================================================
# Convenience wrappers
# ============================================================================


class TestCountOverlaps:
    """Tests for count_overlaps wrapper."""

    def test_returns_array(self, sample_peaks):
        """Return type is a numpy array with correct length."""
        subject = pd.DataFrame(
            {
                "chr": ["chr1", "chr2"],
                "start": [1500, 2500],
                "end": [2500, 3500],
            }
        )
        result = count_overlaps(sample_peaks, subject)
        assert isinstance(result, np.ndarray)
        assert len(result) == len(sample_peaks)

    def test_known_counts(self):
        """Verify correct counts for known overlap pattern."""
        query = pd.DataFrame(
            {
                "chr": ["chr1", "chr1", "chr2"],
                "start": [100, 500, 100],
                "end": [300, 700, 300],
            }
        )
        subject = pd.DataFrame(
            {
                "chr": ["chr1", "chr1"],
                "start": [200, 250],
                "end": [400, 350],
            }
        )
        result = count_overlaps(query, subject)
        # First query (100-300) overlaps both subjects, second (500-700) overlaps none
        assert result[0] == 2
        assert result[1] == 0
        assert result[2] == 0


class TestCountOverlapsWithFraction:
    """Tests for count_overlaps_with_fraction."""

    def test_returns_int(self):
        """Return type is an integer."""
        query = pd.DataFrame({"chr": ["chr1"], "start": [100], "end": [200]})
        subject = pd.DataFrame({"chr": ["chr1"], "start": [100], "end": [200]})
        result = count_overlaps_with_fraction(query, subject, min_overlap_frac=0.5)
        assert isinstance(result, int)

    def test_fraction_filtering(self):
        """Fractional overlap threshold works."""
        query = pd.DataFrame({"chr": ["chr1"], "start": [100], "end": [200]})  # 100bp
        subject = pd.DataFrame({"chr": ["chr1"], "start": [190], "end": [300]})  # 10bp overlap = 10%
        assert count_overlaps_with_fraction(query, subject, min_overlap_frac=0.05) == 1
        assert count_overlaps_with_fraction(query, subject, min_overlap_frac=0.5) == 0


class TestFindFirstOverlaps:
    """Tests for find_first_overlaps."""

    def test_returns_first_only(self):
        """Only the first overlap per query is returned."""
        query = pd.DataFrame({"chr": ["chr1"], "start": [100], "end": [500]})
        subject = pd.DataFrame(
            {
                "chr": ["chr1", "chr1"],
                "start": [150, 300],
                "end": [250, 400],
            }
        )
        result = find_first_overlaps(query, subject)
        assert len(result) == 1


class TestExcludeOverlapping:
    """Tests for exclude_overlapping."""

    def test_removes_overlapping(self):
        """Peaks overlapping exclusion regions are removed."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1", "chr1", "chr2"],
                "start": [100, 500, 100],
                "end": [200, 600, 200],
            }
        )
        exclusion = pd.DataFrame(
            {
                "chr": ["chr1"],
                "start": [150],
                "end": [250],
            }
        )
        result = exclude_overlapping(peaks, exclusion)
        # First peak overlaps exclusion; peaks 2 and 3 should remain
        assert len(result) == 2

    def test_no_exclusion(self):
        """Empty exclusion keeps all peaks."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1", "chr2"],
                "start": [100, 200],
                "end": [200, 300],
            }
        )
        exclusion = pd.DataFrame(columns=["chr", "start", "end"])
        result = exclude_overlapping(peaks, exclusion)
        assert len(result) == len(peaks)


# ============================================================================
# Column detection and standardization
# ============================================================================


class TestDetectColumn:
    """Tests for detect_column."""

    def test_finds_standard_name(self):
        """Finds column when an exact match exists."""
        df = pd.DataFrame({"chr": [1], "start": [100], "end": [200]})
        assert detect_column(df, CHROM_COLS) == "chr"

    def test_case_insensitive(self):
        """Detection is case-insensitive."""
        df = pd.DataFrame({"CHR": [1], "START": [100], "END": [200]})
        assert detect_column(df, CHROM_COLS) == "CHR"

    def test_returns_none_when_missing(self):
        """Returns None when no candidate matches."""
        df = pd.DataFrame({"x": [1], "y": [2]})
        assert detect_column(df, CHROM_COLS) is None

    def test_required_raises(self):
        """Raises ValueError when required=True and not found."""
        df = pd.DataFrame({"x": [1]})
        with pytest.raises(ValueError, match="Could not find"):
            detect_column(df, CHROM_COLS, required=True)


class TestStandardizePeakColumns:
    """Tests for standardize_peak_columns."""

    def test_renames_variants(self):
        """Non-standard column names get renamed."""
        df = pd.DataFrame(
            {
                "chrom": ["chr1"],
                "chromStart": [100],
                "chromEnd": [200],
                "signalValue": [5.0],
            }
        )
        result = standardize_peak_columns(df)
        assert "chr" in result.columns
        assert "start" in result.columns
        assert "end" in result.columns
        assert "signal" in result.columns

    def test_already_standard(self):
        """Standard names are left unchanged."""
        df = pd.DataFrame({"chr": ["chr1"], "start": [100], "end": [200]})
        result = standardize_peak_columns(df)
        assert list(result.columns) == ["chr", "start", "end"]


# ============================================================================
# Peak file loading
# ============================================================================


class TestLoadPeakFile:
    """Tests for load_peak_file."""

    def test_load_bed(self, sample_bed_file):
        """Load a standard BED file."""
        df = load_peak_file(sample_bed_file)
        assert len(df) == 5
        assert "chr" in df.columns
        assert "start" in df.columns
        assert "end" in df.columns

    def test_load_narrowpeak(self, sample_narrowpeak_file):
        """Load a narrowPeak file."""
        df = load_peak_file(sample_narrowpeak_file)
        assert len(df) == 3
        assert "chr" in df.columns

    def test_load_csv_with_header(self, temp_dir):
        """Load a CSV with column headers."""
        csv_path = temp_dir / "peaks.csv"
        csv_path.write_text("chrom,chromStart,chromEnd,score\nchr1,100,200,5\nchr2,300,400,10\n")
        df = load_peak_file(csv_path, sep=",")
        assert len(df) == 2
        assert "chr" in df.columns
        assert "start" in df.columns

    def test_load_from_buffer(self):
        """Load peak data from a file-like object."""
        content = "chr1\t100\t200\tpeak1\t50\nchr2\t300\t400\tpeak2\t60\n"
        buf = io.StringIO(content)
        df = load_peak_file(buf)
        assert len(df) == 2


# ============================================================================
# Chromosome utilities
# ============================================================================


class TestSortChromosomes:
    """Tests for sort_chromosomes."""

    def test_natural_order(self):
        """Chromosomes sort in natural genomic order."""
        chroms = ["chr3", "chr1", "chrX", "chr22", "chr2", "chrY"]
        result = sort_chromosomes(chroms)
        assert result == ["chr1", "chr2", "chr3", "chr22", "chrX", "chrY"]

    def test_unknown_chromosomes(self):
        """Unknown contigs sort to the end."""
        chroms = ["chr1", "chrUn_random", "chr2"]
        result = sort_chromosomes(chroms)
        assert result[-1] == "chrUn_random"


class TestFilterStandardChroms:
    """Tests for filter_standard_chroms."""

    def test_removes_nonstandard(self):
        """Random and haplotype contigs are removed."""
        df = pd.DataFrame(
            {
                "chr": ["chr1", "chr2", "chrUn_gl000220", "chr1_random", "chrX"],
                "start": [100, 200, 300, 400, 500],
                "end": [200, 300, 400, 500, 600],
            }
        )
        result = filter_standard_chroms(df)
        assert set(result["chr"]) == {"chr1", "chr2", "chrX"}

    def test_empty_after_filter(self):
        """If all chroms are non-standard, result is empty."""
        df = pd.DataFrame(
            {
                "chr": ["chrUn", "chrGL"],
                "start": [100, 200],
                "end": [200, 300],
            }
        )
        result = filter_standard_chroms(df)
        assert len(result) == 0
