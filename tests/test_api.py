"""
Unit tests for FastAPI REST API endpoints.

Uses the api_client fixture from conftest.py (TestClient with in-memory SQLite).
"""

import pytest


# ============================================================================
# Health & Config
# ============================================================================


class TestHealthEndpoints:
    """Tests for health check and config endpoints."""

    def test_health(self, api_client):
        resp = api_client.get("/health")
        assert resp.status_code == 200
        data = resp.json()
        assert data["status"] == "healthy"
        assert "version" in data

    def test_config(self, api_client):
        resp = api_client.get("/config")
        assert resp.status_code == 200
        data = resp.json()
        assert "supported_genomes" in data
        assert "histone_marks" in data
        assert isinstance(data["supported_genomes"], list)


# ============================================================================
# Sample CRUD
# ============================================================================


class TestSampleEndpoints:
    """Tests for sample create / read / list / delete."""

    @pytest.fixture
    def _sample_payload(self):
        return {
            "name": "sample_A",
            "histone_mark": "H3K27ac",
            "condition": "Treatment",
            "replicate": 1,
            "genome": "mm10",
        }

    def test_create_sample(self, api_client, _sample_payload):
        resp = api_client.post("/samples", json=_sample_payload)
        assert resp.status_code == 200
        data = resp.json()
        assert data["name"] == "sample_A"
        assert "id" in data

    def test_list_samples(self, api_client, _sample_payload):
        api_client.post("/samples", json=_sample_payload)
        resp = api_client.get("/samples")
        assert resp.status_code == 200
        data = resp.json()
        assert data["total"] >= 1
        assert len(data["samples"]) >= 1

    def test_get_sample(self, api_client, _sample_payload):
        create = api_client.post("/samples", json=_sample_payload).json()
        resp = api_client.get(f"/samples/{create['id']}")
        assert resp.status_code == 200
        assert resp.json()["name"] == "sample_A"

    def test_get_sample_not_found(self, api_client):
        resp = api_client.get("/samples/nonexistent-id")
        assert resp.status_code == 404

    def test_delete_sample(self, api_client, _sample_payload):
        create = api_client.post("/samples", json=_sample_payload).json()
        resp = api_client.delete(f"/samples/{create['id']}")
        assert resp.status_code == 200
        # Verify gone
        resp2 = api_client.get(f"/samples/{create['id']}")
        assert resp2.status_code == 404

    def test_delete_sample_not_found(self, api_client):
        resp = api_client.delete("/samples/nonexistent-id")
        assert resp.status_code == 404

    def test_filter_by_histone_mark(self, api_client):
        api_client.post(
            "/samples",
            json={
                "name": "s1",
                "histone_mark": "H3K27ac",
                "condition": "Ctrl",
                "genome": "mm10",
            },
        )
        api_client.post(
            "/samples",
            json={
                "name": "s2",
                "histone_mark": "H3K4me3",
                "condition": "Ctrl",
                "genome": "mm10",
            },
        )
        resp = api_client.get("/samples", params={"histone_mark": "H3K4me3"})
        data = resp.json()
        assert all(s["histone_mark"] == "H3K4me3" for s in data["samples"])


# ============================================================================
# Sample Sheet Upload
# ============================================================================


class TestSampleSheetUpload:
    """Tests for CSV sample sheet upload."""

    def test_upload_valid_sheet(self, api_client):
        csv = "SampleID,Condition,Factor,Replicate\ns1,Treat,H3K27ac,1\ns2,Ctrl,H3K27ac,1\n"
        resp = api_client.post(
            "/samples/upload-sheet",
            files={"file": ("sheet.csv", csv, "text/csv")},
        )
        assert resp.status_code == 200
        data = resp.json()
        assert data["message"].startswith("Created 2")

    def test_upload_missing_columns(self, api_client):
        csv = "Name,Group\ns1,Treat\n"
        resp = api_client.post(
            "/samples/upload-sheet",
            files={"file": ("bad.csv", csv, "text/csv")},
        )
        assert resp.status_code == 400
        assert "Missing required columns" in resp.json()["detail"]


# ============================================================================
# Job CRUD
# ============================================================================


class TestJobEndpoints:
    """Tests for job create / read / list / cancel."""

    @pytest.fixture
    def _job_payload(self):
        return {
            "name": "Test QC Job",
            "job_type": "qc",
            "config": {},
        }

    def test_create_job(self, api_client, _job_payload):
        resp = api_client.post("/jobs", json=_job_payload)
        assert resp.status_code == 200
        data = resp.json()
        assert data["name"] == "Test QC Job"
        assert data["status"] == "pending"

    def test_list_jobs(self, api_client, _job_payload):
        api_client.post("/jobs", json=_job_payload)
        resp = api_client.get("/jobs")
        assert resp.status_code == 200
        assert resp.json()["total"] >= 1

    def test_get_job(self, api_client, _job_payload):
        create = api_client.post("/jobs", json=_job_payload).json()
        resp = api_client.get(f"/jobs/{create['id']}")
        assert resp.status_code == 200

    def test_get_job_not_found(self, api_client):
        resp = api_client.get("/jobs/nonexistent-id")
        assert resp.status_code == 404

    def test_cancel_pending_job(self, api_client, _job_payload):
        create = api_client.post("/jobs", json=_job_payload).json()
        resp = api_client.post(f"/jobs/{create['id']}/cancel")
        assert resp.status_code == 200
        # Verify cancelled
        status = api_client.get(f"/jobs/{create['id']}").json()
        assert status["status"] == "cancelled"

    def test_cancel_nonexistent_job(self, api_client):
        resp = api_client.post("/jobs/nonexistent/cancel")
        assert resp.status_code == 404


# ============================================================================
# Results Endpoints
# ============================================================================


class TestResultsEndpoints:
    """Tests for results retrieval."""

    def test_results_not_found(self, api_client):
        resp = api_client.get("/results/nonexistent")
        assert resp.status_code == 404

    def test_results_not_completed(self, api_client):
        """Requesting results for a pending job returns 400."""
        job = api_client.post(
            "/jobs",
            json={
                "name": "Pending Job",
                "job_type": "qc",
                "config": {},
            },
        ).json()
        resp = api_client.get(f"/results/{job['id']}")
        assert resp.status_code == 400
        assert "not completed" in resp.json()["detail"]

    def test_download_not_found(self, api_client):
        resp = api_client.get("/results/nonexistent/download/file.csv")
        assert resp.status_code == 404


# ============================================================================
# Comparison Endpoints
# ============================================================================


class TestComparisonEndpoints:
    """Tests for comparison creation and retrieval."""

    def _create_samples(self, api_client):
        """Helper to create two samples and return their IDs."""
        ids = []
        for name, cond in [("s1", "Treat"), ("s2", "Ctrl")]:
            resp = api_client.post(
                "/samples",
                json={
                    "name": name,
                    "histone_mark": "H3K27ac",
                    "condition": cond,
                    "genome": "mm10",
                },
            )
            ids.append(resp.json()["id"])
        return ids

    def test_create_comparison(self, api_client):
        sample_ids = self._create_samples(api_client)
        resp = api_client.post(
            "/comparisons",
            json={
                "name": "Test Comparison",
                "group1": "Treat",
                "group2": "Ctrl",
                "histone_mark": "H3K27ac",
                "genome": "mm10",
                "sample_ids": sample_ids,
            },
        )
        assert resp.status_code == 200
        data = resp.json()
        assert data["name"] == "Test Comparison"
        assert "id" in data

    def test_create_comparison_bad_samples(self, api_client):
        resp = api_client.post(
            "/comparisons",
            json={
                "name": "Bad",
                "group1": "A",
                "group2": "B",
                "histone_mark": "H3K27ac",
                "sample_ids": ["fake-id-1", "fake-id-2"],
            },
        )
        assert resp.status_code == 400

    def test_list_comparisons(self, api_client):
        resp = api_client.get("/comparisons")
        assert resp.status_code == 200
        assert "comparisons" in resp.json()

    def test_get_comparison_not_found(self, api_client):
        resp = api_client.get("/comparisons/nonexistent")
        assert resp.status_code == 404
