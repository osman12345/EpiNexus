"""
Pydantic schemas for API request/response validation.

Defines schemas for:
- Job management
- Sample metadata
- Analysis configuration
- Results formatting
"""

from datetime import datetime
from typing import Optional, List, Dict, Any, Union
from pydantic import BaseModel, Field
from enum import Enum


# Enums for validation
class JobStatusEnum(str, Enum):
    PENDING = "pending"
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class JobTypeEnum(str, Enum):
    ALIGNMENT = "alignment"
    PEAK_CALLING = "peak_calling"
    DIFFBIND = "diffbind"
    INTEGRATION = "integration"
    QC = "qc"
    BIGWIG = "bigwig"
    FULL_PIPELINE = "full_pipeline"


class HistoneMarkEnum(str, Enum):
    H3K27ac = "H3K27ac"
    H3K27me3 = "H3K27me3"
    H3K4me1 = "H3K4me1"
    H3K4me3 = "H3K4me3"
    H3K9me3 = "H3K9me3"
    H3K36me3 = "H3K36me3"
    OTHER = "other"


class GenomeEnum(str, Enum):
    mm10 = "mm10"
    mm39 = "mm39"
    hg38 = "hg38"
    hg19 = "hg19"
    rn6 = "rn6"
    rn7 = "rn7"
    custom = "custom"


# ============================================================================
# Sample Schemas
# ============================================================================


class SampleBase(BaseModel):
    """Base sample schema."""

    name: str = Field(..., description="Sample name/identifier")
    histone_mark: str = Field(..., description="Histone modification (e.g., H3K27ac)")
    condition: str = Field(..., description="Treatment/condition group")
    replicate: int = Field(default=1, ge=1, description="Replicate number")
    genome: str = Field(default="mm10", description="Reference genome")
    tissue: Optional[str] = None
    cell_type: Optional[str] = None


class SampleCreate(SampleBase):
    """Schema for creating a new sample."""

    fastq_r1: Optional[str] = None
    fastq_r2: Optional[str] = None
    bam_file: Optional[str] = None
    peak_file: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None


class SampleResponse(SampleBase):
    """Schema for sample response."""

    id: str
    bam_file: Optional[str] = None
    peak_file: Optional[str] = None
    bigwig_file: Optional[str] = None
    peak_count: Optional[int] = None
    frip_score: Optional[float] = None
    mapping_rate: Optional[float] = None
    created_at: datetime

    class Config:
        from_attributes = True


class SampleListResponse(BaseModel):
    """Response for listing samples."""

    samples: List[SampleResponse]
    total: int


# ============================================================================
# Comparison Schemas
# ============================================================================


class ComparisonBase(BaseModel):
    """Base comparison schema."""

    name: str = Field(..., description="Comparison name")
    description: Optional[str] = None
    group1: str = Field(..., description="Treatment/case group")
    group2: str = Field(..., description="Control/reference group")
    histone_mark: str = Field(..., description="Histone mark for this comparison")
    genome: str = Field(default="mm10")


class ComparisonCreate(ComparisonBase):
    """Schema for creating a comparison."""

    sample_ids: List[str] = Field(..., description="List of sample IDs to include")
    fdr_threshold: float = Field(default=0.1, ge=0, le=1)
    lfc_threshold: float = Field(default=0.5, ge=0)
    min_overlap: int = Field(default=1, ge=1)
    summit_size: int = Field(default=250, ge=50)


class ComparisonResponse(ComparisonBase):
    """Schema for comparison response."""

    id: str
    total_peaks: Optional[int] = None
    significant_peaks: Optional[int] = None
    gained_peaks: Optional[int] = None
    lost_peaks: Optional[int] = None
    results_dir: Optional[str] = None
    created_at: datetime

    class Config:
        from_attributes = True


# ============================================================================
# Job Schemas
# ============================================================================


class JobBase(BaseModel):
    """Base job schema."""

    name: str = Field(..., description="Job name")
    job_type: JobTypeEnum = Field(..., description="Type of analysis")


class JobCreate(JobBase):
    """Schema for creating a job."""

    config: Dict[str, Any] = Field(default_factory=dict, description="Job configuration")


class JobResponse(JobBase):
    """Schema for job response."""

    id: str
    status: JobStatusEnum
    progress: float = 0.0
    current_step: Optional[str] = None
    created_at: datetime
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    error_message: Optional[str] = None
    results: Optional[Dict[str, Any]] = None

    class Config:
        from_attributes = True


class JobStatus(BaseModel):
    """Schema for job status updates."""

    id: str
    status: JobStatusEnum
    progress: float
    current_step: Optional[str] = None
    error_message: Optional[str] = None


class JobListResponse(BaseModel):
    """Response for listing jobs."""

    jobs: List[JobResponse]
    total: int


# ============================================================================
# Analysis Configuration Schemas
# ============================================================================


class DiffBindConfig(BaseModel):
    """Configuration for DiffBind analysis."""

    comparison_id: str = Field(..., description="Comparison to analyze")

    # Normalization
    normalize_method: str = Field(default="RLE", description="Normalization method (RLE, TMM, lib)")
    use_background: bool = Field(default=True, description="Use background normalization")

    # Peak parameters
    min_overlap: int = Field(default=1, ge=1, description="Minimum samples with peak")
    summit_size: int = Field(default=250, ge=50, description="Summit region size")

    # Significance thresholds
    fdr_threshold: float = Field(default=0.1, ge=0, le=1)
    lfc_threshold: float = Field(default=0.5, ge=0)

    # Gene annotation
    tss_upstream: int = Field(default=3000, ge=0, description="TSS upstream distance")
    tss_downstream: int = Field(default=3000, ge=0, description="TSS downstream distance")

    # Output options
    generate_plots: bool = Field(default=True)
    save_counts: bool = Field(default=True)


class IntegrationConfig(BaseModel):
    """Configuration for multi-mark integration."""

    mark_comparisons: Dict[str, str] = Field(
        ..., description="Map of histone mark to comparison ID (e.g., {'H3K27ac': 'comp_1', 'H3K27me3': 'comp_2'})"
    )

    # Integration type
    integration_type: str = Field(default="two_mark", description="Integration type: two_mark or three_mark")

    # Optional RNA-seq integration
    rnaseq_file: Optional[str] = Field(default=None, description="Path to RNA-seq differential expression results")

    # Thresholds
    peak_fdr: float = Field(default=0.1, ge=0, le=1)
    de_fdr: float = Field(default=0.05, ge=0, le=1)
    de_lfc: float = Field(default=0.3, ge=0)

    # Output
    generate_plots: bool = Field(default=True)


class PipelineConfig(BaseModel):
    """Configuration for full pipeline."""

    # Input type
    input_type: str = Field(..., description="Input type: fastq, bam, or bed")
    sample_sheet: str = Field(..., description="Path to sample sheet CSV")
    genome: str = Field(default="mm10")

    # Alignment options (for FASTQ input)
    align_params: Optional[Dict[str, Any]] = None
    use_spikein: bool = Field(default=True)
    spikein_genome: str = Field(default="ecoli")

    # Peak calling options
    peak_caller: str = Field(default="seacr", description="Peak caller: seacr or macs2")
    peak_params: Optional[Dict[str, Any]] = None

    # Analysis options
    run_diffbind: bool = Field(default=True)
    comparisons: Optional[List[ComparisonCreate]] = None

    # QC options
    generate_bigwig: bool = Field(default=True)
    generate_qc_report: bool = Field(default=True)


class AnalysisConfig(BaseModel):
    """Generic analysis configuration."""

    analysis_type: str
    config: Union[DiffBindConfig, IntegrationConfig, PipelineConfig, Dict[str, Any]]


# ============================================================================
# Results Schemas
# ============================================================================


class PeakResult(BaseModel):
    """Schema for a differential peak."""

    peak_id: str
    chrom: str
    start: int
    end: int
    fold_change: float
    fdr: float
    direction: str  # Gained or Lost
    gene_name: Optional[str] = None
    gene_id: Optional[str] = None
    distance_to_tss: Optional[int] = None
    annotation: Optional[str] = None


class DiffBindResults(BaseModel):
    """Schema for DiffBind analysis results."""

    comparison_name: str
    total_peaks: int
    significant_peaks: int
    gained_peaks: int
    lost_peaks: int
    peaks: List[PeakResult]


class IntegrationResults(BaseModel):
    """Schema for integration analysis results."""

    comparison_name: str
    total_genes: int
    chromatin_states: Dict[str, int]
    high_confidence_targets: int
    genes: List[Dict[str, Any]]


class QCMetrics(BaseModel):
    """Schema for QC metrics."""

    sample_name: str
    total_reads: int
    mapped_reads: int
    mapping_rate: float
    duplicate_rate: float
    peak_count: int
    frip_score: float
    median_fragment_size: Optional[int] = None


# ============================================================================
# File Upload Schemas
# ============================================================================


class FileUpload(BaseModel):
    """Schema for file upload response."""

    filename: str
    file_path: str
    file_size: int
    file_type: str


class SampleSheetUpload(BaseModel):
    """Schema for sample sheet upload."""

    samples: List[SampleCreate]
    validation_errors: Optional[List[str]] = None
    warnings: Optional[List[str]] = None
