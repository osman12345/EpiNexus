"""
Custom exception classes for EpiNexus.

Provides clear, module-specific error types for better error handling
and debugging throughout the analysis pipeline.
"""


class EpiNexusError(Exception):
    """Base exception for all EpiNexus errors."""
    pass


# ============================================================================
# Input / File errors
# ============================================================================

class FileFormatError(EpiNexusError):
    """Raised when an input file has an unexpected or invalid format."""
    pass


class PeakFileFormatError(FileFormatError):
    """Raised when a BED/narrowPeak/broadPeak file is malformed."""
    pass


class GTFParseError(FileFormatError):
    """Raised when a GTF/GFF annotation file cannot be parsed."""
    pass


class ExpressionFileFormatError(FileFormatError):
    """Raised when an expression data file has an unrecognized format."""
    pass


# ============================================================================
# Data validation errors
# ============================================================================

class ValidationError(EpiNexusError):
    """Raised when input data fails validation checks."""
    pass


class MissingColumnError(ValidationError):
    """Raised when a required column is missing from a DataFrame."""

    def __init__(self, column: str, dataframe_name: str = "DataFrame", available: list = None):
        available_str = f" Available columns: {available}" if available else ""
        super().__init__(
            f"Required column '{column}' not found in {dataframe_name}.{available_str}"
        )
        self.column = column
        self.available = available


class InsufficientSamplesError(ValidationError):
    """Raised when there are too few samples for an analysis."""

    def __init__(self, required: int, actual: int, context: str = "analysis"):
        super().__init__(
            f"Insufficient samples for {context}: need at least {required}, got {actual}"
        )
        self.required = required
        self.actual = actual


class EmptyDataError(ValidationError):
    """Raised when data is empty where it should not be."""

    def __init__(self, data_name: str = "data"):
        super().__init__(f"Empty {data_name} provided where non-empty data is required")
        self.data_name = data_name


class InvalidParameterError(ValidationError):
    """Raised when a parameter value is out of valid range."""

    def __init__(self, param: str, value, valid_range: str = ""):
        msg = f"Invalid value for '{param}': {value}"
        if valid_range:
            msg += f". Expected: {valid_range}"
        super().__init__(msg)
        self.param = param
        self.value = value


# ============================================================================
# Analysis errors
# ============================================================================

class AnalysisError(EpiNexusError):
    """Base class for analysis-specific errors."""
    pass


class DifferentialAnalysisError(AnalysisError):
    """Raised when differential peak analysis fails."""
    pass


class NormalizationError(AnalysisError):
    """Raised when count matrix normalization fails."""
    pass


class PeakAnnotationError(AnalysisError):
    """Raised when peak annotation fails."""
    pass


class MotifAnalysisError(AnalysisError):
    """Raised when TF motif analysis fails."""
    pass


class IntegrationError(AnalysisError):
    """Raised when multi-omics integration fails."""
    pass


class PeakGeneLinkError(AnalysisError):
    """Raised when peak-to-gene linking fails."""
    pass


# ============================================================================
# Pipeline errors
# ============================================================================

class PipelineError(EpiNexusError):
    """Base class for pipeline execution errors."""
    pass


class AlignmentError(PipelineError):
    """Raised when read alignment fails."""
    pass


class PeakCallingError(PipelineError):
    """Raised when peak calling fails."""
    pass


class PipelineConfigError(PipelineError):
    """Raised when pipeline configuration is invalid."""
    pass


# ============================================================================
# Validation helpers
# ============================================================================

def validate_dataframe(
    df,
    name: str = "DataFrame",
    required_columns: list = None,
    min_rows: int = 0,
) -> None:
    """Validate a DataFrame has expected shape and columns.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to validate.
    name : str
        Human-readable name for error messages.
    required_columns : list, optional
        Columns that must be present.
    min_rows : int
        Minimum number of rows required.

    Raises
    ------
    EmptyDataError
        If df is None or empty and min_rows > 0.
    MissingColumnError
        If a required column is missing.
    """
    import pandas as pd

    if df is None:
        raise EmptyDataError(name)

    if not isinstance(df, pd.DataFrame):
        raise ValidationError(f"Expected DataFrame for {name}, got {type(df).__name__}")

    if min_rows > 0 and len(df) < min_rows:
        if len(df) == 0:
            raise EmptyDataError(name)
        raise ValidationError(
            f"{name} has {len(df)} rows but at least {min_rows} are required"
        )

    if required_columns:
        for col in required_columns:
            if col not in df.columns:
                raise MissingColumnError(col, name, available=list(df.columns))


def validate_numeric_param(value, name: str, min_val=None, max_val=None) -> None:
    """Validate a numeric parameter is within acceptable bounds.

    Raises
    ------
    InvalidParameterError
        If the value is out of range.
    """
    if min_val is not None and value < min_val:
        raise InvalidParameterError(name, value, f">= {min_val}")
    if max_val is not None and value > max_val:
        raise InvalidParameterError(name, value, f"<= {max_val}")
