"""
SQLAlchemy database models for the Histone Mark Analyzer.

Defines tables for:
- Jobs: Track analysis job status and configuration
- Samples: Store sample metadata and file paths
- Comparisons: Define differential analysis comparisons
- Analysis: Store analysis results and metrics
"""

import uuid
from datetime import datetime
from typing import Optional, Dict, Any, List
from sqlalchemy import (
    Column, String, Integer, Float, Boolean, DateTime,
    ForeignKey, Text, JSON, Enum as SQLEnum, Table
)
from sqlalchemy.orm import relationship, declarative_base
from sqlalchemy.sql import func
import enum

Base = declarative_base()


class JobStatus(enum.Enum):
    """Job status enumeration."""
    PENDING = "pending"
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class JobType(enum.Enum):
    """Type of analysis job."""
    ALIGNMENT = "alignment"
    PEAK_CALLING = "peak_calling"
    DIFFBIND = "diffbind"
    INTEGRATION = "integration"
    QC = "qc"
    BIGWIG = "bigwig"
    FULL_PIPELINE = "full_pipeline"


class DataType(enum.Enum):
    """Input data type."""
    FASTQ = "fastq"
    BAM = "bam"
    BED = "bed"
    BIGWIG = "bigwig"


# Association table for many-to-many relationship between samples and comparisons
sample_comparison = Table(
    'sample_comparison',
    Base.metadata,
    Column('sample_id', String, ForeignKey('samples.id')),
    Column('comparison_id', String, ForeignKey('comparisons.id'))
)


class Job(Base):
    """Job tracking table."""

    __tablename__ = "jobs"

    id = Column(String, primary_key=True, default=lambda: str(uuid.uuid4()))
    name = Column(String, nullable=False)
    job_type = Column(SQLEnum(JobType), nullable=False)
    status = Column(SQLEnum(JobStatus), default=JobStatus.PENDING)

    # Timestamps
    created_at = Column(DateTime, default=func.now())
    started_at = Column(DateTime, nullable=True)
    completed_at = Column(DateTime, nullable=True)

    # Configuration and results
    config = Column(JSON, nullable=True)  # Job configuration parameters
    results = Column(JSON, nullable=True)  # Summary results/metrics
    output_dir = Column(String, nullable=True)  # Results directory path

    # Progress tracking
    progress = Column(Float, default=0.0)  # 0-100
    current_step = Column(String, nullable=True)
    total_steps = Column(Integer, nullable=True)

    # Error handling
    error_message = Column(Text, nullable=True)
    error_traceback = Column(Text, nullable=True)

    # Relationships
    analyses = relationship("Analysis", back_populates="job")

    def __repr__(self):
        return f"<Job(id={self.id}, name={self.name}, status={self.status})>"

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "id": self.id,
            "name": self.name,
            "job_type": self.job_type.value,
            "status": self.status.value,
            "created_at": self.created_at.isoformat() if self.created_at else None,
            "started_at": self.started_at.isoformat() if self.started_at else None,
            "completed_at": self.completed_at.isoformat() if self.completed_at else None,
            "progress": self.progress,
            "current_step": self.current_step,
            "error_message": self.error_message
        }


class Sample(Base):
    """Sample metadata table."""

    __tablename__ = "samples"

    id = Column(String, primary_key=True, default=lambda: str(uuid.uuid4()))
    name = Column(String, nullable=False, index=True)

    # Sample metadata
    histone_mark = Column(String, nullable=False)  # H3K27ac, H3K27me3, H3K4me1, etc.
    condition = Column(String, nullable=False)  # Treatment group
    replicate = Column(Integer, default=1)
    tissue = Column(String, nullable=True)
    cell_type = Column(String, nullable=True)

    # Reference genome
    genome = Column(String, default="mm10")

    # File paths
    fastq_r1 = Column(String, nullable=True)
    fastq_r2 = Column(String, nullable=True)
    bam_file = Column(String, nullable=True)
    bam_index = Column(String, nullable=True)
    peak_file = Column(String, nullable=True)
    bigwig_file = Column(String, nullable=True)

    # QC metrics
    total_reads = Column(Integer, nullable=True)
    mapped_reads = Column(Integer, nullable=True)
    mapping_rate = Column(Float, nullable=True)
    duplicate_rate = Column(Float, nullable=True)
    peak_count = Column(Integer, nullable=True)
    frip_score = Column(Float, nullable=True)  # Fraction of reads in peaks

    # Spike-in calibration
    spikein_reads = Column(Integer, nullable=True)
    scale_factor = Column(Float, nullable=True)

    # Timestamps
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, default=func.now(), onupdate=func.now())

    # Additional metadata
    metadata = Column(JSON, nullable=True)

    # Relationships
    comparisons = relationship(
        "Comparison",
        secondary=sample_comparison,
        back_populates="samples"
    )

    def __repr__(self):
        return f"<Sample(id={self.id}, name={self.name}, mark={self.histone_mark})>"

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "id": self.id,
            "name": self.name,
            "histone_mark": self.histone_mark,
            "condition": self.condition,
            "replicate": self.replicate,
            "genome": self.genome,
            "bam_file": self.bam_file,
            "peak_file": self.peak_file,
            "bigwig_file": self.bigwig_file,
            "peak_count": self.peak_count,
            "frip_score": self.frip_score
        }


class Comparison(Base):
    """Differential analysis comparison definition."""

    __tablename__ = "comparisons"

    id = Column(String, primary_key=True, default=lambda: str(uuid.uuid4()))
    name = Column(String, nullable=False, index=True)
    description = Column(Text, nullable=True)

    # Comparison groups
    group1 = Column(String, nullable=False)  # Treatment/case
    group2 = Column(String, nullable=False)  # Control/reference
    histone_mark = Column(String, nullable=False)

    # Reference genome
    genome = Column(String, default="mm10")

    # Analysis parameters
    fdr_threshold = Column(Float, default=0.1)
    lfc_threshold = Column(Float, default=0.5)
    min_overlap = Column(Integer, default=1)
    summit_size = Column(Integer, default=250)

    # Results summary
    total_peaks = Column(Integer, nullable=True)
    significant_peaks = Column(Integer, nullable=True)
    gained_peaks = Column(Integer, nullable=True)
    lost_peaks = Column(Integer, nullable=True)
    results_dir = Column(String, nullable=True)

    # Timestamps
    created_at = Column(DateTime, default=func.now())

    # Relationships
    samples = relationship(
        "Sample",
        secondary=sample_comparison,
        back_populates="comparisons"
    )
    analyses = relationship("Analysis", back_populates="comparison")

    def __repr__(self):
        return f"<Comparison(id={self.id}, name={self.name}, {self.group1} vs {self.group2})>"

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "id": self.id,
            "name": self.name,
            "description": self.description,
            "group1": self.group1,
            "group2": self.group2,
            "histone_mark": self.histone_mark,
            "genome": self.genome,
            "fdr_threshold": self.fdr_threshold,
            "total_peaks": self.total_peaks,
            "significant_peaks": self.significant_peaks,
            "gained_peaks": self.gained_peaks,
            "lost_peaks": self.lost_peaks
        }


class Analysis(Base):
    """Analysis results storage."""

    __tablename__ = "analyses"

    id = Column(String, primary_key=True, default=lambda: str(uuid.uuid4()))
    analysis_type = Column(String, nullable=False)  # diffbind, integration, qc

    # Foreign keys
    job_id = Column(String, ForeignKey("jobs.id"), nullable=True)
    comparison_id = Column(String, ForeignKey("comparisons.id"), nullable=True)

    # Results
    results_file = Column(String, nullable=True)  # Path to main results CSV
    plots_file = Column(String, nullable=True)  # Path to plots PDF
    summary = Column(JSON, nullable=True)  # Summary statistics

    # Metrics specific to analysis type
    metrics = Column(JSON, nullable=True)

    # Timestamps
    created_at = Column(DateTime, default=func.now())

    # Relationships
    job = relationship("Job", back_populates="analyses")
    comparison = relationship("Comparison", back_populates="analyses")

    def __repr__(self):
        return f"<Analysis(id={self.id}, type={self.analysis_type})>"


# Database initialization functions
def init_db(engine):
    """Initialize database tables."""
    Base.metadata.create_all(bind=engine)


def get_session(engine):
    """Create a new database session."""
    from sqlalchemy.orm import sessionmaker
    Session = sessionmaker(bind=engine)
    return Session()
