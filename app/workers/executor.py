"""
Background job executor for EpiNexus.

Executes analysis jobs (DiffBind, integration, QC) in the background
using FastAPI BackgroundTasks. Updates job status and results in the
database as the analysis progresses.
"""

import logging
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Any


logger = logging.getLogger(__name__)


def execute_job(job_id: str, db_factory):
    """Execute an analysis job in the background.

    Parameters
    ----------
    job_id : str
        The ID of the job to execute.
    db_factory : callable
        A callable that returns a new SQLAlchemy session.
    """
    db = db_factory()
    try:
        from ..models.database import Job, JobStatus, JobType

        job = db.query(Job).filter(Job.id == job_id).first()
        if not job:
            logger.error(f"Job {job_id} not found")
            return

        # Mark as running
        job.status = JobStatus.RUNNING
        job.started_at = datetime.now(timezone.utc)
        job.current_step = "Initializing"
        db.commit()

        logger.info(f"Starting job {job_id}: {job.name} (type={job.job_type.value})")

        # Route to the appropriate handler
        handlers = {
            JobType.DIFFBIND: _run_diffbind,
            JobType.INTEGRATION: _run_integration,
            JobType.QC: _run_qc,
        }

        handler = handlers.get(job.job_type)
        if not handler:
            raise NotImplementedError(f"No handler for job type: {job.job_type.value}")

        results = handler(job, db)

        # Mark as completed
        job.status = JobStatus.COMPLETED
        job.completed_at = datetime.now(timezone.utc)
        job.progress = 100.0
        job.current_step = "Done"
        job.results = results
        db.commit()

        logger.info(f"Job {job_id} completed successfully")

    except Exception as e:
        logger.error(f"Job {job_id} failed: {e}")
        try:
            job = db.query(Job).filter(Job.id == job_id).first()
            if job:
                job.status = JobStatus.FAILED
                job.completed_at = datetime.now(timezone.utc)
                job.error_message = str(e)
                job.error_traceback = traceback.format_exc()
                db.commit()
        except Exception:
            logger.error(f"Failed to update error status for job {job_id}")
    finally:
        db.close()


def _update_progress(job, db, progress: float, step: str):
    """Update job progress in the database."""
    job.progress = progress
    job.current_step = step
    db.commit()


def _run_diffbind(job, db) -> Dict[str, Any]:
    """Execute a DiffBind differential analysis job."""
    from sqlalchemy.orm import joinedload

    from ..core.differential import DifferentialAnalyzer, DifferentialConfig, Sample
    from ..models.database import Comparison

    config = job.config or {}
    comparison_id = config.get("comparison_id")

    # Eager-load samples to avoid N+1 query pattern
    comparison = (
        db.query(Comparison).options(joinedload(Comparison.samples)).filter(Comparison.id == comparison_id).first()
    )
    if not comparison:
        raise ValueError(f"Comparison {comparison_id} not found")

    _update_progress(job, db, 10, "Loading samples")

    # Build sample list from comparison
    samples = []
    for s in comparison.samples:
        samples.append(
            Sample(
                sample_id=s.name,
                condition=s.condition,
                histone_mark=s.histone_mark,
                replicate=s.replicate,
                bam_file=s.bam_file,
                peak_file=s.peak_file,
            )
        )

    if len(samples) < 2:
        raise ValueError("Need at least 2 samples for differential analysis")

    _update_progress(job, db, 20, "Creating consensus peaks")

    analyzer = DifferentialAnalyzer()

    # Create diff config
    diff_config = DifferentialConfig(
        comparison_name=comparison.name,
        group1=comparison.group1,
        group2=comparison.group2,
        min_overlap=config.get("min_overlap", 1),
        summit_extend=config.get("summit_size", 250),
        normalize_method=config.get("normalize_method", "RLE"),
        fdr_threshold=config.get("fdr_threshold", 0.1),
        lfc_threshold=config.get("lfc_threshold", 0.5),
    )

    _update_progress(job, db, 40, "Counting reads in peaks")

    # Create consensus peaks
    consensus = analyzer.create_consensus_peaks(samples, min_overlap=diff_config.min_overlap)

    _update_progress(job, db, 70, "Running differential analysis")

    # Save output
    output_dir = Path(f"results/{job.id}")
    output_dir.mkdir(parents=True, exist_ok=True)
    consensus.to_csv(output_dir / "consensus_peaks.csv", index=False)

    job.output_dir = str(output_dir)

    _update_progress(job, db, 90, "Saving results")

    return {
        "total_peaks": len(consensus),
        "comparison": comparison.name,
        "output_dir": str(output_dir),
    }


def _run_integration(job, db) -> Dict[str, Any]:
    """Execute a multi-mark integration job."""
    from ..core.integration import MarkIntegration, IntegrationConfig

    config = job.config or {}

    _update_progress(job, db, 10, "Loading integration configuration")

    int_config = IntegrationConfig(
        integration_type=config.get("integration_type", "two_mark"),
        h3k27ac_file=config.get("h3k27ac_file"),
        h3k27me3_file=config.get("h3k27me3_file"),
        h3k4me1_file=config.get("h3k4me1_file"),
        rnaseq_file=config.get("rnaseq_file"),
        peak_fdr=config.get("peak_fdr", 0.1),
        de_fdr=config.get("de_fdr", 0.05),
        de_lfc=config.get("de_lfc", 0.3),
    )

    # Set output directory
    output_dir = Path(f"results/{job.id}")
    output_dir.mkdir(parents=True, exist_ok=True)
    int_config.output_dir = str(output_dir)

    _update_progress(job, db, 30, "Running multi-mark integration")

    integrator = MarkIntegration()

    if int_config.integration_type == "three_mark":
        results = integrator.run_three_mark_integration(int_config)
    else:
        results = integrator.run_two_mark_integration(int_config)

    _update_progress(job, db, 80, "Saving results")

    job.output_dir = str(output_dir)

    return {
        "integration_type": results.integration_type,
        "total_genes": results.total_genes,
        "chromatin_states": results.chromatin_states,
        "high_confidence_count": results.high_confidence_count,
        "integrated_file": results.integrated_file,
        "summary": results.summary,
        "output_dir": str(output_dir),
    }


def _run_qc(job, db) -> Dict[str, Any]:
    """Execute a QC analysis job."""
    from ..core.qc import QCAnalyzer

    config = job.config or {}

    _update_progress(job, db, 10, "Initializing QC analysis")

    analyzer = QCAnalyzer()

    sample_name = config.get("sample_name", job.name)
    bam_path = config.get("bam_file")
    peak_path = config.get("peak_file")
    spikein_bam = config.get("spikein_bam")

    _update_progress(job, db, 30, "Running QC metrics")

    metrics = analyzer.run_full_qc(
        sample_name=sample_name,
        bam_path=bam_path,
        peak_path=peak_path,
        spikein_bam=spikein_bam,
    )

    _update_progress(job, db, 70, "Generating QC report")

    # Save report
    output_dir = Path(f"results/{job.id}")
    output_dir.mkdir(parents=True, exist_ok=True)

    report_df = analyzer.generate_qc_report([metrics], str(output_dir / "qc_report.csv"))

    job.output_dir = str(output_dir)

    _update_progress(job, db, 90, "Saving results")

    return {
        "status": "completed",
        "metrics": metrics.to_dict(),
        "report_file": str(output_dir / "qc_report.csv"),
        "output_dir": str(output_dir),
    }
