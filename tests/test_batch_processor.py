"""Tests for Batch Processing System."""

import pytest
import tempfile
import time
import json
import threading
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import Mock, patch, MagicMock

from app.core.batch_processor import Job, JobStatus, BatchProcessor


class TestJob:
    """Test Job class."""

    def test_job_creation_with_auto_created_at(self):
        """Test creating a job with automatic created_at timestamp."""
        job = Job(
            id="job123",
            name="test_job",
            job_type="qc_analysis",
            params={"sample": "sample1"},
        )

        assert job.id == "job123"
        assert job.name == "test_job"
        assert job.job_type == "qc_analysis"
        assert job.params == {"sample": "sample1"}
        assert job.status == JobStatus.PENDING
        assert job.progress == 0.0
        assert job.created_at != ""
        assert job.started_at == ""
        assert job.completed_at == ""
        assert job.message == ""
        assert job.error == ""
        assert job.result == {}

    def test_job_creation_with_custom_created_at(self):
        """Test creating a job with custom created_at."""
        custom_time = "2024-01-01T12:00:00"
        job = Job(
            id="job123",
            name="test_job",
            job_type="qc_analysis",
            params={},
            created_at=custom_time,
        )

        assert job.created_at == custom_time

    def test_job_status_transitions(self):
        """Test job status transitions."""
        job = Job(
            id="job123",
            name="test_job",
            job_type="qc_analysis",
            params={},
        )

        assert job.status == JobStatus.PENDING

        job.status = JobStatus.RUNNING
        assert job.status == JobStatus.RUNNING

        job.status = JobStatus.COMPLETED
        assert job.status == JobStatus.COMPLETED

    def test_job_to_dict(self):
        """Test converting job to dictionary."""
        job = Job(
            id="job123",
            name="test_job",
            job_type="qc_analysis",
            params={"sample": "sample1"},
            status=JobStatus.RUNNING,
            progress=50.0,
            message="Processing...",
            result={"result_key": "result_value"},
        )

        result = job.to_dict()

        assert isinstance(result, dict)
        assert result["id"] == "job123"
        assert result["name"] == "test_job"
        assert result["job_type"] == "qc_analysis"
        assert result["status"] == "running"  # Enum converted to string
        assert result["progress"] == 50.0
        assert result["message"] == "Processing..."
        assert result["result"] == {"result_key": "result_value"}

    def test_job_from_dict(self):
        """Test creating job from dictionary."""
        data = {
            "id": "job123",
            "name": "test_job",
            "job_type": "qc_analysis",
            "params": {"sample": "sample1"},
            "status": "running",
            "progress": 50.0,
            "message": "Processing...",
            "created_at": "2024-01-01T12:00:00",
            "started_at": "2024-01-01T12:01:00",
            "completed_at": "",
            "result": {},
            "error": "",
        }

        job = Job.from_dict(data)

        assert job.id == "job123"
        assert job.name == "test_job"
        assert job.job_type == "qc_analysis"
        assert job.params == {"sample": "sample1"}
        assert job.status == JobStatus.RUNNING
        assert job.progress == 50.0
        assert job.message == "Processing..."

    def test_job_serialization_roundtrip(self):
        """Test that job can be serialized and deserialized."""
        original = Job(
            id="job123",
            name="test_job",
            job_type="qc_analysis",
            params={"sample": "sample1"},
            status=JobStatus.RUNNING,
            progress=75.0,
            result={"data": "value"},
        )

        # Serialize
        data = original.to_dict()

        # Deserialize
        restored = Job.from_dict(data)

        assert restored.id == original.id
        assert restored.name == original.name
        assert restored.job_type == original.job_type
        assert restored.status == original.status
        assert restored.progress == original.progress
        assert restored.result == original.result


class TestBatchProcessorJobCreation:
    """Test BatchProcessor job creation."""

    def test_create_job(self):
        """Test creating a job."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(jobs_dir=tmpdir)

            job = processor.create_job(
                name="test_job",
                job_type="qc_analysis",
                params={"sample": "sample1"},
            )

            assert job.id is not None
            assert job.name == "test_job"
            assert job.job_type == "qc_analysis"
            assert job.params == {"sample": "sample1"}
            assert job.status == JobStatus.PENDING
            assert job.id in processor.jobs

    def test_create_job_persisted(self):
        """Test that created job is persisted to disk."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(jobs_dir=tmpdir)

            job = processor.create_job(
                name="test_job",
                job_type="qc_analysis",
                params={"sample": "sample1"},
            )

            # Check that job file was created
            job_file = Path(tmpdir) / f"{job.id}.json"
            assert job_file.exists()

            # Load and verify job data
            with open(job_file) as f:
                data = json.load(f)
                assert data["id"] == job.id
                assert data["name"] == "test_job"

    def test_submit_batch(self):
        """Test submitting multiple jobs as a batch."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(jobs_dir=tmpdir)

            jobs_spec = [
                {"name": "job1", "job_type": "qc_analysis", "params": {"sample": "s1"}},
                {"name": "job2", "job_type": "qc_analysis", "params": {"sample": "s2"}},
                {"name": "job3", "job_type": "qc_analysis", "params": {"sample": "s3"}},
            ]

            jobs = processor.submit_batch(jobs_spec, batch_name="test_batch")

            assert len(jobs) == 3
            for job in jobs:
                assert "batch_id" in job.params
                assert "batch_index" in job.params


class TestBatchProcessorHandlers:
    """Test BatchProcessor handler registration."""

    def test_register_handler(self):
        """Test registering a custom handler."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(jobs_dir=tmpdir)

            def custom_handler(params, progress_cb):
                progress_cb(50, "Processing...")
                return {"status": "success"}

            processor.register_handler("custom_job", custom_handler)

            assert "custom_job" in processor.job_handlers
            assert processor.job_handlers["custom_job"] == custom_handler

    def test_default_handlers_registered(self):
        """Test that default handlers are registered on init."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(jobs_dir=tmpdir)

            # Check that default handlers are present
            assert "differential_analysis" in processor.job_handlers
            assert "qc_analysis" in processor.job_handlers
            assert "peak_annotation" in processor.job_handlers
            assert "super_enhancer" in processor.job_handlers


class TestBatchProcessorProcessing:
    """Test BatchProcessor job processing."""

    def test_start_processing_with_simple_handler(self):
        """Test processing a job with a simple handler."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(max_workers=2, jobs_dir=tmpdir)

            # Create a simple handler
            def simple_handler(params, progress_cb):
                progress_cb(50, "Processing...")
                return {"output": "success", "value": params.get("input_value", 42)}

            processor.register_handler("simple", simple_handler)

            # Create and process a job
            job = processor.create_job(
                name="simple_job",
                job_type="simple",
                params={"input_value": 100},
            )

            job_count = processor.start_processing([job.id])
            assert job_count == 1

            # Wait for completion
            time.sleep(1)

            # Check job status
            completed_job = processor.get_job(job.id)
            assert completed_job.status == JobStatus.COMPLETED
            assert completed_job.result["output"] == "success"
            assert completed_job.result["value"] == 100
            assert completed_job.progress == 100.0

    def test_start_processing_multiple_jobs(self):
        """Test processing multiple jobs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(max_workers=2, jobs_dir=tmpdir)

            def simple_handler(params, progress_cb):
                progress_cb(100, "Done")
                return {"job_id": params.get("id")}

            processor.register_handler("simple", simple_handler)

            # Create multiple jobs
            jobs = []
            for i in range(3):
                job = processor.create_job(
                    name=f"job_{i}",
                    job_type="simple",
                    params={"id": f"id_{i}"},
                )
                jobs.append(job.id)

            job_count = processor.start_processing(jobs)
            assert job_count == 3

            # Wait for completion
            time.sleep(2)

            # Check all jobs are processed
            for job_id in jobs:
                job = processor.get_job(job_id)
                assert job.status == JobStatus.COMPLETED

    def test_start_processing_all_pending_jobs(self):
        """Test processing all pending jobs when no job_ids specified."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(max_workers=2, jobs_dir=tmpdir)

            def simple_handler(params, progress_cb):
                progress_cb(100, "Done")
                return {"status": "completed"}

            processor.register_handler("simple", simple_handler)

            # Create multiple jobs
            for i in range(2):
                processor.create_job(
                    name=f"job_{i}",
                    job_type="simple",
                    params={},
                )

            # Start processing without specifying job_ids
            job_count = processor.start_processing()
            assert job_count == 2

            time.sleep(1)

            all_jobs = processor.get_all_jobs()
            assert all(j.status == JobStatus.COMPLETED for j in all_jobs)

    def test_job_progress_callback(self):
        """Test that progress callback is called."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(max_workers=1, jobs_dir=tmpdir)

            progress_calls = []

            def tracking_handler(params, progress_cb):
                for i in range(3):
                    progress_cb(int((i + 1) * 33), f"Step {i + 1}")
                    time.sleep(0.1)
                return {"status": "done"}

            processor.register_handler("tracking", tracking_handler)

            job = processor.create_job(
                name="tracking_job",
                job_type="tracking",
                params={},
            )

            processor.start_processing([job.id])
            time.sleep(1)

            # Reload job to check final state
            completed_job = processor.get_job(job.id)
            assert completed_job.status == JobStatus.COMPLETED
            assert completed_job.progress == 100.0

    def test_handler_exception_sets_failed_status(self):
        """Test that handler exceptions set job to FAILED."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(max_workers=1, jobs_dir=tmpdir)

            def failing_handler(params, progress_cb):
                raise ValueError("Handler failed")

            processor.register_handler("failing", failing_handler)

            job = processor.create_job(
                name="failing_job",
                job_type="failing",
                params={},
            )

            processor.start_processing([job.id])
            time.sleep(0.5)

            failed_job = processor.get_job(job.id)
            assert failed_job.status == JobStatus.FAILED
            assert "Handler failed" in failed_job.error


class TestBatchProcessorJobCancellation:
    """Test job cancellation."""

    def test_cancel_pending_job(self):
        """Test cancelling a pending job."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(jobs_dir=tmpdir)

            job = processor.create_job(
                name="test_job",
                job_type="qc_analysis",
                params={},
            )

            assert job.status == JobStatus.PENDING

            success = processor.cancel_job(job.id)
            assert success is True

            cancelled_job = processor.get_job(job.id)
            assert cancelled_job.status == JobStatus.CANCELLED

    def test_cancel_running_job_fails(self):
        """Test that cancelling a running job returns False."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(jobs_dir=tmpdir)

            job = processor.create_job(
                name="test_job",
                job_type="qc_analysis",
                params={},
            )

            # Manually set to running
            job.status = JobStatus.RUNNING

            success = processor.cancel_job(job.id)
            assert success is False
            assert job.status == JobStatus.RUNNING

    def test_cancel_nonexistent_job(self):
        """Test cancelling a non-existent job."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(jobs_dir=tmpdir)

            success = processor.cancel_job("nonexistent_id")
            assert success is False


class TestBatchProcessorJobRetrieval:
    """Test job retrieval methods."""

    def test_get_job(self):
        """Test retrieving a single job."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(jobs_dir=tmpdir)

            job = processor.create_job(
                name="test_job",
                job_type="qc_analysis",
                params={"sample": "s1"},
            )

            retrieved = processor.get_job(job.id)
            assert retrieved is not None
            assert retrieved.id == job.id
            assert retrieved.name == "test_job"

    def test_get_nonexistent_job(self):
        """Test retrieving a non-existent job."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(jobs_dir=tmpdir)

            result = processor.get_job("nonexistent_id")
            assert result is None

    def test_get_all_jobs(self):
        """Test retrieving all jobs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(jobs_dir=tmpdir)

            # Create multiple jobs
            for i in range(3):
                processor.create_job(
                    name=f"job_{i}",
                    job_type="qc_analysis",
                    params={},
                )

            all_jobs = processor.get_all_jobs()
            assert len(all_jobs) == 3

    def test_get_all_jobs_sorted_by_creation_time(self):
        """Test that all jobs are sorted by creation time (reverse)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(jobs_dir=tmpdir)

            # Create jobs with small delays
            for i in range(3):
                processor.create_job(
                    name=f"job_{i}",
                    job_type="qc_analysis",
                    params={},
                )
                time.sleep(0.01)

            all_jobs = processor.get_all_jobs()
            # Should be in reverse creation order
            assert all_jobs[0].name == "job_2"
            assert all_jobs[1].name == "job_1"
            assert all_jobs[2].name == "job_0"


class TestBatchProcessorQueueStatus:
    """Test queue status reporting."""

    def test_get_queue_status(self):
        """Test getting queue status counts."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(jobs_dir=tmpdir)

            def simple_handler(params, progress_cb):
                progress_cb(100, "Done")
                return {"status": "success"}

            processor.register_handler("simple", simple_handler)

            # Create jobs in different states
            job1 = processor.create_job("job1", "simple", {})
            job2 = processor.create_job("job2", "simple", {})
            job3 = processor.create_job("job3", "simple", {})

            # Process one job
            processor.start_processing([job1.id])
            time.sleep(0.5)

            # Cancel one
            processor.cancel_job(job3.id)

            status = processor.get_queue_status()

            assert status["pending"] >= 1  # job2 should be pending
            assert status["completed"] >= 1  # job1 should be completed
            assert status["cancelled"] == 1


class TestBatchProcessorPersistence:
    """Test job persistence and loading."""

    def test_load_jobs_on_startup(self):
        """Test that jobs are loaded from disk on startup."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create processor and add jobs
            processor1 = BatchProcessor(jobs_dir=tmpdir)
            job1 = processor1.create_job("job1", "qc_analysis", {})
            job2 = processor1.create_job("job2", "qc_analysis", {})

            # Create new processor - should load existing jobs
            processor2 = BatchProcessor(jobs_dir=tmpdir)
            all_jobs = processor2.get_all_jobs()

            assert len(all_jobs) == 2
            assert any(j.id == job1.id for j in all_jobs)
            assert any(j.id == job2.id for j in all_jobs)

    def test_clear_completed_jobs(self):
        """Test clearing completed and failed jobs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(jobs_dir=tmpdir)

            def simple_handler(params, progress_cb):
                progress_cb(100, "Done")
                return {"status": "success"}

            processor.register_handler("simple", simple_handler)

            # Create and complete a job
            job1 = processor.create_job("job1", "simple", {})
            processor.start_processing([job1.id])
            time.sleep(0.5)

            # Create a pending job
            job2 = processor.create_job("job2", "simple", {})

            # Clear completed
            processor.clear_completed()

            remaining = processor.get_all_jobs()
            assert len(remaining) == 1
            assert remaining[0].id == job2.id

            # Job file should be deleted
            job1_file = Path(tmpdir) / f"{job1.id}.json"
            assert not job1_file.exists()


class TestBatchProcessorTimeout:
    """Test timeout behavior."""

    def test_job_timeout_configuration(self):
        """Test that timeout configuration is applied correctly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Set a very short timeout
            processor = BatchProcessor(
                max_workers=1,
                jobs_dir=tmpdir,
                job_timeouts={"slow": 2},  # 2 second timeout
            )

            # Verify that the timeout is configured
            assert processor.job_timeouts["slow"] == 2

            # Test that custom timeouts override defaults
            assert processor.job_timeouts["differential_analysis"] == 3600

            # Test that GLOBAL_TIMEOUT is used as fallback
            assert processor.GLOBAL_TIMEOUT == 7200

    def test_custom_job_timeouts(self):
        """Test custom timeout configuration."""
        with tempfile.TemporaryDirectory() as tmpdir:
            custom_timeouts = {
                "fast_job": 10,
                "slow_job": 60,
            }

            processor = BatchProcessor(
                jobs_dir=tmpdir,
                job_timeouts=custom_timeouts,
            )

            assert processor.job_timeouts["fast_job"] == 10
            assert processor.job_timeouts["slow_job"] == 60
            # Default timeout should still be available
            assert processor.job_timeouts["differential_analysis"] == 3600


class TestBatchProcessorShutdown:
    """Test processor shutdown."""

    def test_shutdown(self):
        """Test graceful shutdown."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(jobs_dir=tmpdir)

            def simple_handler(params, progress_cb):
                progress_cb(100, "Done")
                return {}

            processor.register_handler("simple", simple_handler)

            job = processor.create_job("job", "simple", {})
            processor.start_processing([job.id])

            # Shutdown should not raise
            processor.shutdown()

            # Executor should be None
            assert processor.executor is None


class TestBatchProcessorEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_empty_batch(self):
        """Test submitting an empty batch."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(jobs_dir=tmpdir)

            jobs = processor.submit_batch([], batch_name="empty")

            assert len(jobs) == 0

    def test_job_with_empty_params(self):
        """Test creating a job with empty parameters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(jobs_dir=tmpdir)

            job = processor.create_job(
                name="empty_params_job",
                job_type="qc_analysis",
                params={},
            )

            assert job.params == {}

    def test_handler_with_no_progress_updates(self):
        """Test handler that doesn't call progress callback."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(max_workers=1, jobs_dir=tmpdir)

            def silent_handler(params, progress_cb):
                # Never calls progress_cb
                return {"result": "done"}

            processor.register_handler("silent", silent_handler)

            job = processor.create_job("silent_job", "silent", {})
            processor.start_processing([job.id])
            time.sleep(0.5)

            completed = processor.get_job(job.id)
            assert completed.status == JobStatus.COMPLETED
            assert completed.result["result"] == "done"

    def test_concurrent_job_creation_and_processing(self):
        """Test concurrent job creation and processing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BatchProcessor(max_workers=2, jobs_dir=tmpdir)

            def simple_handler(params, progress_cb):
                time.sleep(0.1)
                progress_cb(100, "Done")
                return {"id": params.get("id")}

            processor.register_handler("simple", simple_handler)

            # Create jobs in one thread
            created_jobs = []
            for i in range(3):
                job = processor.create_job(
                    f"job_{i}",
                    "simple",
                    {"id": i},
                )
                created_jobs.append(job.id)

            # Process them
            processor.start_processing(created_jobs)
            time.sleep(1)

            # All should be completed
            for job_id in created_jobs:
                job = processor.get_job(job_id)
                assert job.status == JobStatus.COMPLETED
