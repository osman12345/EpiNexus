"""
Batch Processing System

Queue and process multiple samples/analyses:
- Job queue management
- Progress tracking
- Parallel processing
- Result aggregation
"""

import os
import json
import time
import uuid
import threading
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, List, Optional, Callable
from dataclasses import dataclass, asdict, field
from enum import Enum
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd


class JobStatus(Enum):
    """Job execution status."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass
class Job:
    """Represents a batch processing job."""
    id: str
    name: str
    job_type: str
    params: Dict[str, Any]
    status: JobStatus = JobStatus.PENDING
    progress: float = 0.0
    message: str = ""
    created_at: str = ""
    started_at: str = ""
    completed_at: str = ""
    result: Dict[str, Any] = field(default_factory=dict)
    error: str = ""

    def __post_init__(self):
        if not self.created_at:
            self.created_at = datetime.now().isoformat()

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d['status'] = self.status.value
        return d

    @classmethod
    def from_dict(cls, data: Dict) -> 'Job':
        data['status'] = JobStatus(data['status'])
        return cls(**data)


class BatchProcessor:
    """
    Manage batch processing of multiple analyses.

    Supports:
    - Differential analysis across multiple comparisons
    - Multi-sample QC processing
    - Peak calling for multiple samples
    - Annotation of multiple peak sets
    """

    def __init__(
        self,
        max_workers: int = 4,
        jobs_dir: str = "batch_jobs"
    ):
        self.max_workers = max_workers
        self.jobs_dir = Path(jobs_dir)
        self.jobs_dir.mkdir(parents=True, exist_ok=True)

        self.jobs: Dict[str, Job] = {}
        self.executor: Optional[ThreadPoolExecutor] = None
        self.job_handlers: Dict[str, Callable] = {}

        # Register default handlers
        self._register_default_handlers()

        # Load existing jobs
        self._load_jobs()

    def _register_default_handlers(self):
        """Register default job type handlers."""
        self.register_handler("differential_analysis", self._handle_differential)
        self.register_handler("qc_analysis", self._handle_qc)
        self.register_handler("peak_annotation", self._handle_annotation)
        self.register_handler("super_enhancer", self._handle_super_enhancer)

    def register_handler(self, job_type: str, handler: Callable):
        """Register a handler function for a job type."""
        self.job_handlers[job_type] = handler

    def create_job(
        self,
        name: str,
        job_type: str,
        params: Dict[str, Any]
    ) -> Job:
        """Create a new batch job."""

        job = Job(
            id=str(uuid.uuid4())[:8],
            name=name,
            job_type=job_type,
            params=params
        )

        self.jobs[job.id] = job
        self._save_job(job)

        return job

    def submit_batch(
        self,
        jobs: List[Dict[str, Any]],
        batch_name: str = None
    ) -> List[Job]:
        """Submit multiple jobs as a batch."""

        created_jobs = []
        batch_id = str(uuid.uuid4())[:8]

        for i, job_spec in enumerate(jobs):
            job = self.create_job(
                name=job_spec.get('name', f'{batch_name or "Batch"}_{i+1}'),
                job_type=job_spec['job_type'],
                params=job_spec.get('params', {})
            )
            job.params['batch_id'] = batch_id
            job.params['batch_index'] = i
            created_jobs.append(job)

        return created_jobs

    def start_processing(self, job_ids: List[str] = None):
        """Start processing jobs."""

        if self.executor is None:
            self.executor = ThreadPoolExecutor(max_workers=self.max_workers)

        # Get jobs to process
        if job_ids:
            jobs_to_run = [self.jobs[jid] for jid in job_ids if jid in self.jobs]
        else:
            jobs_to_run = [j for j in self.jobs.values() if j.status == JobStatus.PENDING]

        # Submit jobs
        futures = {}
        for job in jobs_to_run:
            future = self.executor.submit(self._run_job, job)
            futures[future] = job

        return len(futures)

    def _run_job(self, job: Job):
        """Execute a single job."""

        job.status = JobStatus.RUNNING
        job.started_at = datetime.now().isoformat()
        job.progress = 0.0
        self._save_job(job)

        try:
            handler = self.job_handlers.get(job.job_type)
            if handler is None:
                raise ValueError(f"Unknown job type: {job.job_type}")

            # Run handler with progress callback
            def update_progress(progress: float, message: str = ""):
                job.progress = progress
                job.message = message
                self._save_job(job)

            result = handler(job.params, update_progress)

            job.status = JobStatus.COMPLETED
            job.progress = 100.0
            job.result = result
            job.completed_at = datetime.now().isoformat()

        except Exception as e:
            job.status = JobStatus.FAILED
            job.error = str(e)
            job.completed_at = datetime.now().isoformat()

        self._save_job(job)
        return job

    def cancel_job(self, job_id: str) -> bool:
        """Cancel a pending job."""
        if job_id in self.jobs:
            job = self.jobs[job_id]
            if job.status == JobStatus.PENDING:
                job.status = JobStatus.CANCELLED
                self._save_job(job)
                return True
        return False

    def get_job(self, job_id: str) -> Optional[Job]:
        """Get job by ID."""
        return self.jobs.get(job_id)

    def get_all_jobs(self) -> List[Job]:
        """Get all jobs sorted by creation time."""
        return sorted(self.jobs.values(), key=lambda j: j.created_at, reverse=True)

    def get_queue_status(self) -> Dict[str, int]:
        """Get counts by status."""
        status_counts = {s.value: 0 for s in JobStatus}
        for job in self.jobs.values():
            status_counts[job.status.value] += 1
        return status_counts

    def clear_completed(self):
        """Remove completed and failed jobs."""
        to_remove = [
            jid for jid, job in self.jobs.items()
            if job.status in (JobStatus.COMPLETED, JobStatus.FAILED, JobStatus.CANCELLED)
        ]
        for jid in to_remove:
            job_file = self.jobs_dir / f"{jid}.json"
            if job_file.exists():
                job_file.unlink()
            del self.jobs[jid]

    def _save_job(self, job: Job):
        """Save job to disk."""
        job_file = self.jobs_dir / f"{job.id}.json"
        with open(job_file, 'w') as f:
            json.dump(job.to_dict(), f, indent=2)

    def _load_jobs(self):
        """Load existing jobs from disk."""
        import logging
        logger = logging.getLogger(__name__)

        for job_file in self.jobs_dir.glob("*.json"):
            try:
                with open(job_file, 'r') as f:
                    data = json.load(f)
                    job = Job.from_dict(data)
                    self.jobs[job.id] = job
            except (json.JSONDecodeError, KeyError, TypeError) as e:
                # Skip corrupted or invalid job files
                logger.warning(f"Failed to load job file {job_file}: {e}")
                continue

    # Default job handlers
    def _handle_differential(self, params: Dict, progress_cb: Callable) -> Dict:
        """Handle differential analysis job."""
        progress_cb(10, "Loading samples...")
        time.sleep(1)  # Simulate processing

        progress_cb(30, "Building count matrix...")
        time.sleep(1)

        progress_cb(50, "Running statistical test...")
        time.sleep(2)

        progress_cb(80, "Generating results...")
        time.sleep(1)

        progress_cb(100, "Complete")

        return {
            "total_peaks": 45000,
            "differential_peaks": 2500,
            "up_regulated": 1300,
            "down_regulated": 1200
        }

    def _handle_qc(self, params: Dict, progress_cb: Callable) -> Dict:
        """Handle QC analysis job."""
        progress_cb(20, "Calculating mapping statistics...")
        time.sleep(1)

        progress_cb(50, "Computing FRiP scores...")
        time.sleep(1)

        progress_cb(80, "Generating QC report...")
        time.sleep(1)

        progress_cb(100, "Complete")

        return {
            "samples_processed": params.get("n_samples", 8),
            "samples_passed": params.get("n_samples", 8) - 1,
            "avg_frip": 0.32
        }

    def _handle_annotation(self, params: Dict, progress_cb: Callable) -> Dict:
        """Handle peak annotation job."""
        progress_cb(25, "Loading gene annotations...")
        time.sleep(1)

        progress_cb(50, "Annotating peaks...")
        time.sleep(2)

        progress_cb(75, "Running enrichment...")
        time.sleep(1)

        progress_cb(100, "Complete")

        return {
            "peaks_annotated": params.get("n_peaks", 30000),
            "promoter_peaks": 12000,
            "enhancer_peaks": 8000
        }

    def _handle_super_enhancer(self, params: Dict, progress_cb: Callable) -> Dict:
        """Handle super-enhancer detection job."""
        progress_cb(20, "Stitching enhancers...")
        time.sleep(1)

        progress_cb(50, "Ranking by signal...")
        time.sleep(1)

        progress_cb(80, "Identifying super-enhancers...")
        time.sleep(1)

        progress_cb(100, "Complete")

        return {
            "total_enhancers": 15000,
            "super_enhancers": 500,
            "typical_enhancers": 14500
        }

    def shutdown(self):
        """Shutdown the executor."""
        if self.executor:
            self.executor.shutdown(wait=True)
            self.executor = None
