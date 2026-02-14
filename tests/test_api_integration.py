"""
Integration tests for FastAPI REST API endpoints.

These tests exercise the full HTTP request lifecycle through the FastAPI
TestClient, including middleware, dependency injection, database round-trips,
and response serialization.  They complement the unit-level tests in
test_api.py by testing multi-step workflows and edge cases.
"""

import pytest


# ============================================================================
# Fixtures
# ============================================================================


@pytest.fixture
def _sample(api_client):
    """Create a sample and return its JSON representation."""
    resp = api_client.post(
        "/samples",
        json={
            "name": "integ_sample",
            "histone_mark": "H3K27ac",
            "condition": "Treatment",
            "replicate": 1,
            "genome": "hg38",
        },
    )
    assert resp.status_code == 200
    return resp.json()


@pytest.fixture
def _job(api_client):
    """Create a pending job and return its JSON."""
    resp = api_client.post(
        "/jobs",
        json={
            "name": "Integration Job",
            "job_type": "qc",
            "config": {"foo": "bar"},
        },
    )
    assert resp.status_code == 200
    return resp.json()


# ============================================================================
# Sample lifecycle
# ============================================================================


class TestSampleLifecycle:
    """End-to-end: create → read → update → list → delete → verify gone."""

    def test_full_sample_lifecycle(self, api_client):
        # Create
        resp = api_client.post(
            "/samples",
            json={
                "name": "lifecycle_s1",
                "histone_mark": "H3K4me3",
                "condition": "Control",
                "genome": "mm10",
            },
        )
        assert resp.status_code == 200
        sample = resp.json()
        sample_id = sample["id"]

        # Read back
        resp = api_client.get(f"/samples/{sample_id}")
        assert resp.status_code == 200
        assert resp.json()["histone_mark"] == "H3K4me3"

        # List (should include our sample)
        resp = api_client.get("/samples")
        ids = [s["id"] for s in resp.json()["samples"]]
        assert sample_id in ids

        # Delete
        resp = api_client.delete(f"/samples/{sample_id}")
        assert resp.status_code == 200

        # Verify gone
        resp = api_client.get(f"/samples/{sample_id}")
        assert resp.status_code == 404


class TestSampleValidation:
    """Input validation for sample endpoints."""

    def test_create_sample_missing_name(self, api_client):
        resp = api_client.post(
            "/samples",
            json={
                "histone_mark": "H3K27ac",
                "condition": "Treatment",
                "genome": "hg38",
            },
        )
        assert resp.status_code == 422  # Pydantic validation error

    def test_create_duplicate_samples(self, api_client):
        """Creating two samples with the same name should succeed (names are not unique)."""
        payload = {
            "name": "dup_sample",
            "histone_mark": "H3K27ac",
            "condition": "Ctrl",
            "genome": "hg38",
        }
        r1 = api_client.post("/samples", json=payload)
        r2 = api_client.post("/samples", json=payload)
        assert r1.status_code == 200
        assert r2.status_code == 200
        assert r1.json()["id"] != r2.json()["id"]

    def test_sample_sheet_upload_semicolon_sep(self, api_client):
        """Non-comma CSV should fail gracefully if columns are missing."""
        tsv = "SampleID;Condition;Factor;Replicate\ns1;Treat;H3K27ac;1\n"
        resp = api_client.post(
            "/samples/upload-sheet",
            files={"file": ("sheet.csv", tsv, "text/csv")},
        )
        # Should either parse (if pandas handles ; delimiter) or return 400
        assert resp.status_code in (200, 400)


# ============================================================================
# Job lifecycle
# ============================================================================


class TestJobLifecycle:
    """End-to-end job management."""

    def test_create_and_cancel_job(self, api_client):
        resp = api_client.post(
            "/jobs",
            json={
                "name": "Cancel Me",
                "job_type": "qc",
                "config": {},
            },
        )
        assert resp.status_code == 200
        job_id = resp.json()["id"]

        # Cancel
        resp = api_client.post(f"/jobs/{job_id}/cancel")
        assert resp.status_code == 200

        # Verify cancelled
        resp = api_client.get(f"/jobs/{job_id}")
        assert resp.json()["status"] == "cancelled"

    def test_job_list_pagination(self, api_client):
        """Creating multiple jobs should all appear in the list."""
        ids = []
        for i in range(5):
            r = api_client.post(
                "/jobs",
                json={
                    "name": f"batch_{i}",
                    "job_type": "qc",
                    "config": {},
                },
            )
            ids.append(r.json()["id"])

        resp = api_client.get("/jobs")
        listed_ids = [j["id"] for j in resp.json()["jobs"]]
        for jid in ids:
            assert jid in listed_ids


# ============================================================================
# Comparison workflow
# ============================================================================


class TestComparisonWorkflow:
    """Create samples → create comparison → verify."""

    def test_comparison_with_valid_samples(self, api_client):
        # Create treatment and control samples
        sample_ids = []
        for name, cond in [("t1", "Treated"), ("t2", "Treated"), ("c1", "Control"), ("c2", "Control")]:
            r = api_client.post(
                "/samples",
                json={
                    "name": name,
                    "histone_mark": "H3K27ac",
                    "condition": cond,
                    "genome": "hg38",
                },
            )
            sample_ids.append(r.json()["id"])

        # Create comparison
        resp = api_client.post(
            "/comparisons",
            json={
                "name": "Treated_vs_Control",
                "group1": "Treated",
                "group2": "Control",
                "histone_mark": "H3K27ac",
                "genome": "hg38",
                "sample_ids": sample_ids,
            },
        )
        assert resp.status_code == 200
        comp = resp.json()
        assert comp["name"] == "Treated_vs_Control"

        # Verify it appears in the list
        resp = api_client.get("/comparisons")
        names = [c["name"] for c in resp.json()["comparisons"]]
        assert "Treated_vs_Control" in names


# ============================================================================
# Download / path traversal security
# ============================================================================


class TestDownloadSecurity:
    """Verify path traversal protection on the download endpoint."""

    def test_download_rejects_dot_dot(self, api_client, _job):
        """Filenames containing '..' must be rejected."""
        resp = api_client.get(f"/results/{_job['id']}/download/../../../etc/passwd")
        assert resp.status_code in (400, 404, 422)

    def test_download_rejects_encoded_traversal(self, api_client, _job):
        """Filenames with slashes must be rejected."""
        resp = api_client.get(f"/results/{_job['id']}/download/foo/bar.csv")
        # FastAPI path parameter may not match, or our guard should reject
        assert resp.status_code in (400, 404, 422)

    def test_download_nonexistent_job(self, api_client):
        resp = api_client.get("/results/no-such-id/download/file.csv")
        assert resp.status_code == 404

    def test_download_pending_job(self, api_client, _job):
        """Pending jobs have no output_dir, should return 404."""
        resp = api_client.get(f"/results/{_job['id']}/download/result.csv")
        # The job exists but has no output_dir → 400 ("not completed") or 404
        assert resp.status_code in (400, 404)


# ============================================================================
# Health / config sanity
# ============================================================================


class TestHealthDeep:
    """Deeper health/config checks."""

    def test_health_response_structure(self, api_client):
        resp = api_client.get("/health")
        data = resp.json()
        assert "status" in data
        assert "version" in data
        assert data["status"] == "healthy"

    def test_config_lists_are_nonempty(self, api_client):
        resp = api_client.get("/config")
        data = resp.json()
        assert len(data.get("supported_genomes", [])) > 0
        assert len(data.get("histone_marks", [])) > 0


# ============================================================================
# Error handling
# ============================================================================


class TestErrorHandling:
    """Verify proper error responses for invalid inputs."""

    def test_invalid_json_body(self, api_client):
        resp = api_client.post(
            "/samples",
            content="not json",
            headers={"Content-Type": "application/json"},
        )
        assert resp.status_code == 422

    def test_empty_sample_sheet(self, api_client):
        resp = api_client.post(
            "/samples/upload-sheet",
            files={"file": ("empty.csv", "", "text/csv")},
        )
        assert resp.status_code == 400

    def test_nonexistent_endpoint(self, api_client):
        resp = api_client.get("/api/v99/does-not-exist")
        assert resp.status_code == 404
