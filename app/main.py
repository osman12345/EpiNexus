"""
FastAPI application for the Histone Mark Analyzer.

Provides REST API endpoints for:
- Job management (submit, status, cancel)
- Sample management (upload, list, delete)
- Analysis execution (DiffBind, integration)
- Results retrieval and visualization
"""

import logging
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, HTTPException, Depends, UploadFile, File, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse, JSONResponse
from sqlalchemy import create_engine
from sqlalchemy.orm import Session, sessionmaker

from .config import settings
from .models.database import Base, init_db, Job, Sample, Comparison
from .models.schemas import (
    JobCreate, JobResponse, JobListResponse, JobStatus,
    SampleCreate, SampleResponse, SampleListResponse,
    ComparisonCreate, ComparisonResponse,
    DiffBindConfig, IntegrationConfig
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Database setup
engine = create_engine(
    settings.database_url,
    connect_args={"check_same_thread": False} if "sqlite" in settings.database_url else {}
)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan manager."""
    # Startup
    logger.info("Starting Histone Mark Analyzer API...")
    settings.ensure_directories()
    init_db(engine)
    logger.info("Database initialized")
    yield
    # Shutdown
    logger.info("Shutting down Histone Mark Analyzer API...")


# Create FastAPI application
app = FastAPI(
    title=settings.app_name,
    description="API for analyzing histone modifications from CUT&Tag/ChIP-seq data",
    version=settings.app_version,
    lifespan=lifespan
)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Configure appropriately for production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# Dependency for database session
def get_db():
    """Get database session."""
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


# ============================================================================
# Health Check
# ============================================================================

@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "healthy", "version": settings.app_version}


@app.get("/config")
async def get_config():
    """Get application configuration."""
    from .config import GenomeConfig
    return {
        "supported_genomes": list(GenomeConfig.SUPPORTED_GENOMES.keys()),
        "histone_marks": list(GenomeConfig.HISTONE_MARKS.keys()),
        "default_genome": settings.default_genome,
        "default_fdr": settings.default_fdr_threshold
    }


# ============================================================================
# Sample Endpoints
# ============================================================================

@app.get("/samples", response_model=SampleListResponse)
async def list_samples(
    histone_mark: str = None,
    condition: str = None,
    genome: str = None,
    skip: int = 0,
    limit: int = 100,
    db: Session = Depends(get_db)
):
    """List all samples with optional filtering."""
    query = db.query(Sample)

    if histone_mark:
        query = query.filter(Sample.histone_mark == histone_mark)
    if condition:
        query = query.filter(Sample.condition == condition)
    if genome:
        query = query.filter(Sample.genome == genome)

    total = query.count()
    samples = query.offset(skip).limit(limit).all()

    return SampleListResponse(
        samples=[SampleResponse.model_validate(s) for s in samples],
        total=total
    )


@app.post("/samples", response_model=SampleResponse)
async def create_sample(sample: SampleCreate, db: Session = Depends(get_db)):
    """Create a new sample."""
    db_sample = Sample(
        name=sample.name,
        histone_mark=sample.histone_mark,
        condition=sample.condition,
        replicate=sample.replicate,
        genome=sample.genome,
        tissue=sample.tissue,
        cell_type=sample.cell_type,
        fastq_r1=sample.fastq_r1,
        fastq_r2=sample.fastq_r2,
        bam_file=sample.bam_file,
        peak_file=sample.peak_file,
        metadata=sample.metadata
    )
    db.add(db_sample)
    db.commit()
    db.refresh(db_sample)
    return SampleResponse.model_validate(db_sample)


@app.get("/samples/{sample_id}", response_model=SampleResponse)
async def get_sample(sample_id: str, db: Session = Depends(get_db)):
    """Get a specific sample by ID."""
    sample = db.query(Sample).filter(Sample.id == sample_id).first()
    if not sample:
        raise HTTPException(status_code=404, detail="Sample not found")
    return SampleResponse.model_validate(sample)


@app.delete("/samples/{sample_id}")
async def delete_sample(sample_id: str, db: Session = Depends(get_db)):
    """Delete a sample."""
    sample = db.query(Sample).filter(Sample.id == sample_id).first()
    if not sample:
        raise HTTPException(status_code=404, detail="Sample not found")
    db.delete(sample)
    db.commit()
    return {"message": "Sample deleted"}


@app.post("/samples/upload-sheet")
async def upload_sample_sheet(
    file: UploadFile = File(...),
    db: Session = Depends(get_db)
):
    """Upload a sample sheet CSV and create samples."""
    import pandas as pd
    import io

    content = await file.read()
    df = pd.read_csv(io.BytesIO(content))

    required_cols = ["SampleID", "Condition", "Factor"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise HTTPException(
            status_code=400,
            detail=f"Missing required columns: {missing}"
        )

    created_samples = []
    errors = []

    for idx, row in df.iterrows():
        try:
            sample = Sample(
                name=row["SampleID"],
                condition=row["Condition"],
                histone_mark=row.get("Factor", row.get("histone_mark", "unknown")),
                replicate=int(row.get("Replicate", 1)),
                genome=row.get("Genome", settings.default_genome),
                bam_file=row.get("bamReads", row.get("bam_file")),
                peak_file=row.get("Peaks", row.get("peak_file"))
            )
            db.add(sample)
            created_samples.append(sample.name)
        except Exception as e:
            errors.append(f"Row {idx}: {str(e)}")

    db.commit()

    return {
        "message": f"Created {len(created_samples)} samples",
        "samples": created_samples,
        "errors": errors if errors else None
    }


# ============================================================================
# Comparison Endpoints
# ============================================================================

@app.get("/comparisons")
async def list_comparisons(
    histone_mark: str = None,
    db: Session = Depends(get_db)
):
    """List all comparisons."""
    query = db.query(Comparison)
    if histone_mark:
        query = query.filter(Comparison.histone_mark == histone_mark)
    comparisons = query.all()
    return {"comparisons": [c.to_dict() for c in comparisons]}


@app.post("/comparisons", response_model=ComparisonResponse)
async def create_comparison(
    comparison: ComparisonCreate,
    db: Session = Depends(get_db)
):
    """Create a new comparison."""
    # Verify samples exist
    samples = db.query(Sample).filter(Sample.id.in_(comparison.sample_ids)).all()
    if len(samples) != len(comparison.sample_ids):
        raise HTTPException(status_code=400, detail="Some sample IDs not found")

    db_comparison = Comparison(
        name=comparison.name,
        description=comparison.description,
        group1=comparison.group1,
        group2=comparison.group2,
        histone_mark=comparison.histone_mark,
        genome=comparison.genome,
        fdr_threshold=comparison.fdr_threshold,
        lfc_threshold=comparison.lfc_threshold,
        min_overlap=comparison.min_overlap,
        summit_size=comparison.summit_size
    )
    db_comparison.samples = samples

    db.add(db_comparison)
    db.commit()
    db.refresh(db_comparison)

    return ComparisonResponse.model_validate(db_comparison)


@app.get("/comparisons/{comparison_id}", response_model=ComparisonResponse)
async def get_comparison(comparison_id: str, db: Session = Depends(get_db)):
    """Get a specific comparison."""
    comparison = db.query(Comparison).filter(Comparison.id == comparison_id).first()
    if not comparison:
        raise HTTPException(status_code=404, detail="Comparison not found")
    return ComparisonResponse.model_validate(comparison)


# ============================================================================
# Job Endpoints
# ============================================================================

@app.get("/jobs", response_model=JobListResponse)
async def list_jobs(
    status: str = None,
    job_type: str = None,
    skip: int = 0,
    limit: int = 50,
    db: Session = Depends(get_db)
):
    """List all jobs with optional filtering."""
    from .models.database import JobStatus as DBJobStatus, JobType

    query = db.query(Job).order_by(Job.created_at.desc())

    if status:
        query = query.filter(Job.status == DBJobStatus(status))
    if job_type:
        query = query.filter(Job.job_type == JobType(job_type))

    total = query.count()
    jobs = query.offset(skip).limit(limit).all()

    return JobListResponse(
        jobs=[JobResponse.model_validate(j) for j in jobs],
        total=total
    )


@app.post("/jobs", response_model=JobResponse)
async def create_job(
    job: JobCreate,
    background_tasks: BackgroundTasks,
    db: Session = Depends(get_db)
):
    """Create and submit a new analysis job."""
    from .models.database import JobType, JobStatus as DBJobStatus

    db_job = Job(
        name=job.name,
        job_type=JobType(job.job_type.value),
        status=DBJobStatus.PENDING,
        config=job.config
    )
    db.add(db_job)
    db.commit()
    db.refresh(db_job)

    # Queue the job for background execution
    # background_tasks.add_task(execute_job, db_job.id)

    return JobResponse.model_validate(db_job)


@app.get("/jobs/{job_id}", response_model=JobResponse)
async def get_job(job_id: str, db: Session = Depends(get_db)):
    """Get job status and details."""
    job = db.query(Job).filter(Job.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    return JobResponse.model_validate(job)


@app.post("/jobs/{job_id}/cancel")
async def cancel_job(job_id: str, db: Session = Depends(get_db)):
    """Cancel a running job."""
    from .models.database import JobStatus as DBJobStatus

    job = db.query(Job).filter(Job.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status not in [DBJobStatus.PENDING, DBJobStatus.QUEUED, DBJobStatus.RUNNING]:
        raise HTTPException(status_code=400, detail="Job cannot be cancelled")

    job.status = DBJobStatus.CANCELLED
    db.commit()

    return {"message": "Job cancelled"}


# ============================================================================
# Analysis Endpoints
# ============================================================================

@app.post("/analysis/diffbind")
async def run_diffbind_analysis(
    config: DiffBindConfig,
    background_tasks: BackgroundTasks,
    db: Session = Depends(get_db)
):
    """Submit a DiffBind differential analysis job."""
    from .models.database import JobType, JobStatus as DBJobStatus

    # Verify comparison exists
    comparison = db.query(Comparison).filter(
        Comparison.id == config.comparison_id
    ).first()
    if not comparison:
        raise HTTPException(status_code=404, detail="Comparison not found")

    # Create job
    job = Job(
        name=f"DiffBind: {comparison.name}",
        job_type=JobType.DIFFBIND,
        status=DBJobStatus.PENDING,
        config=config.model_dump()
    )
    db.add(job)
    db.commit()
    db.refresh(job)

    # TODO: Queue for background execution
    # background_tasks.add_task(run_diffbind, job.id, config)

    return {"job_id": job.id, "message": "DiffBind analysis submitted"}


@app.post("/analysis/integration")
async def run_integration_analysis(
    config: IntegrationConfig,
    background_tasks: BackgroundTasks,
    db: Session = Depends(get_db)
):
    """Submit a multi-mark integration analysis job."""
    from .models.database import JobType, JobStatus as DBJobStatus

    # Verify all comparisons exist
    for mark, comp_id in config.mark_comparisons.items():
        comparison = db.query(Comparison).filter(Comparison.id == comp_id).first()
        if not comparison:
            raise HTTPException(
                status_code=404,
                detail=f"Comparison not found for {mark}: {comp_id}"
            )

    # Create job
    job = Job(
        name=f"Integration: {config.integration_type}",
        job_type=JobType.INTEGRATION,
        status=DBJobStatus.PENDING,
        config=config.model_dump()
    )
    db.add(job)
    db.commit()
    db.refresh(job)

    return {"job_id": job.id, "message": "Integration analysis submitted"}


# ============================================================================
# Results Endpoints
# ============================================================================

@app.get("/results/{job_id}")
async def get_results(job_id: str, db: Session = Depends(get_db)):
    """Get analysis results for a completed job."""
    from .models.database import JobStatus as DBJobStatus

    job = db.query(Job).filter(Job.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != DBJobStatus.COMPLETED:
        raise HTTPException(
            status_code=400,
            detail=f"Job is not completed (status: {job.status.value})"
        )

    return {
        "job_id": job.id,
        "job_type": job.job_type.value,
        "results": job.results,
        "output_dir": job.output_dir
    }


@app.get("/results/{job_id}/download/{filename}")
async def download_result_file(
    job_id: str,
    filename: str,
    db: Session = Depends(get_db)
):
    """Download a specific result file."""
    job = db.query(Job).filter(Job.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if not job.output_dir:
        raise HTTPException(status_code=404, detail="No output directory")

    file_path = Path(job.output_dir) / filename
    if not file_path.exists():
        raise HTTPException(status_code=404, detail="File not found")

    return FileResponse(
        path=str(file_path),
        filename=filename,
        media_type="application/octet-stream"
    )


# ============================================================================
# Run with uvicorn
# ============================================================================

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(
        "app.main:app",
        host=settings.api_host,
        port=settings.api_port,
        reload=settings.debug
    )
