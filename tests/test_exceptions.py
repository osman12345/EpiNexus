"""
Unit tests for custom exception classes and validation helpers.
"""

import pandas as pd
import pytest

from app.core.exceptions import (
    EpiNexusError,
    FileFormatError,
    PeakFileFormatError,
    GTFParseError,
    ExpressionFileFormatError,
    ValidationError,
    MissingColumnError,
    InsufficientSamplesError,
    EmptyDataError,
    InvalidParameterError,
    AnalysisError,
    DifferentialAnalysisError,
    NormalizationError,
    PipelineError,
    AlignmentError,
    PeakCallingError,
    validate_dataframe,
    validate_numeric_param,
)


# ============================================================================
# Exception hierarchy
# ============================================================================


class TestExceptionHierarchy:
    """Verify the exception inheritance chain."""

    def test_base_exception(self):
        with pytest.raises(EpiNexusError):
            raise EpiNexusError("base error")

    def test_file_format_is_epinexus(self):
        with pytest.raises(EpiNexusError):
            raise FileFormatError("bad file")

    def test_peak_file_is_file_format(self):
        with pytest.raises(FileFormatError):
            raise PeakFileFormatError("bad peak file")

    def test_gtf_is_file_format(self):
        with pytest.raises(FileFormatError):
            raise GTFParseError("bad GTF")

    def test_validation_is_epinexus(self):
        with pytest.raises(EpiNexusError):
            raise ValidationError("invalid")

    def test_analysis_is_epinexus(self):
        with pytest.raises(EpiNexusError):
            raise AnalysisError("failed")

    def test_differential_is_analysis(self):
        with pytest.raises(AnalysisError):
            raise DifferentialAnalysisError("diff failed")

    def test_normalization_is_analysis(self):
        with pytest.raises(AnalysisError):
            raise NormalizationError("norm failed")

    def test_pipeline_is_epinexus(self):
        with pytest.raises(EpiNexusError):
            raise PipelineError("pipe error")

    def test_alignment_is_pipeline(self):
        with pytest.raises(PipelineError):
            raise AlignmentError("align failed")

    def test_peak_calling_is_pipeline(self):
        with pytest.raises(PipelineError):
            raise PeakCallingError("peaks failed")


# ============================================================================
# Rich exception attributes
# ============================================================================


class TestMissingColumnError:
    """Tests for MissingColumnError with custom attributes."""

    def test_message_includes_column(self):
        err = MissingColumnError("log2FC", "peaks_df")
        assert "log2FC" in str(err)
        assert "peaks_df" in str(err)

    def test_available_columns_shown(self):
        err = MissingColumnError("log2FC", available=["chr", "start", "end"])
        assert "chr" in str(err)
        assert err.column == "log2FC"
        assert err.available == ["chr", "start", "end"]

    def test_is_validation_error(self):
        with pytest.raises(ValidationError):
            raise MissingColumnError("score")


class TestInsufficientSamplesError:
    """Tests for InsufficientSamplesError attributes."""

    def test_message(self):
        err = InsufficientSamplesError(required=3, actual=1, context="DiffBind")
        assert "3" in str(err)
        assert "1" in str(err)
        assert "DiffBind" in str(err)
        assert err.required == 3
        assert err.actual == 1


class TestEmptyDataError:
    """Tests for EmptyDataError."""

    def test_message(self):
        err = EmptyDataError("peak matrix")
        assert "peak matrix" in str(err)
        assert err.data_name == "peak matrix"


class TestInvalidParameterError:
    """Tests for InvalidParameterError attributes."""

    def test_with_range(self):
        err = InvalidParameterError("fdr", -0.5, ">= 0")
        assert "fdr" in str(err)
        assert "-0.5" in str(err)
        assert ">= 0" in str(err)
        assert err.param == "fdr"
        assert err.value == -0.5

    def test_without_range(self):
        err = InvalidParameterError("method", "BAD")
        assert "BAD" in str(err)


# ============================================================================
# Validation helpers
# ============================================================================


class TestValidateDataframe:
    """Tests for validate_dataframe helper."""

    def test_none_raises(self):
        with pytest.raises(EmptyDataError):
            validate_dataframe(None, "test_df")

    def test_non_dataframe_raises(self):
        with pytest.raises(ValidationError, match="Expected DataFrame"):
            validate_dataframe([1, 2, 3], "test_list")

    def test_empty_with_min_rows(self):
        df = pd.DataFrame(columns=["a", "b"])
        with pytest.raises(EmptyDataError):
            validate_dataframe(df, "empty_df", min_rows=1)

    def test_too_few_rows(self):
        df = pd.DataFrame({"a": [1], "b": [2]})
        with pytest.raises(ValidationError, match="at least 5"):
            validate_dataframe(df, "small_df", min_rows=5)

    def test_missing_column(self):
        df = pd.DataFrame({"chr": ["chr1"], "start": [100]})
        with pytest.raises(MissingColumnError):
            validate_dataframe(df, "peaks", required_columns=["chr", "start", "end"])

    def test_valid_passes(self):
        df = pd.DataFrame({"chr": ["chr1"], "start": [100], "end": [200]})
        # Should not raise
        validate_dataframe(df, "peaks", required_columns=["chr", "start", "end"], min_rows=1)

    def test_no_constraints(self):
        df = pd.DataFrame({"x": [1]})
        # Should not raise with no constraints
        validate_dataframe(df, "simple")


class TestValidateNumericParam:
    """Tests for validate_numeric_param helper."""

    def test_below_min(self):
        with pytest.raises(InvalidParameterError):
            validate_numeric_param(-1, "fdr", min_val=0)

    def test_above_max(self):
        with pytest.raises(InvalidParameterError):
            validate_numeric_param(2.0, "fdr", max_val=1.0)

    def test_valid_range(self):
        # Should not raise
        validate_numeric_param(0.05, "fdr", min_val=0, max_val=1)

    def test_boundary_values(self):
        # Exact boundary should be valid
        validate_numeric_param(0, "count", min_val=0)
        validate_numeric_param(1.0, "fdr", max_val=1.0)
