"""
Shared test fixtures for EpiNexus test suite.
"""

import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# ============================================================================
# Peak DataFrames
# ============================================================================


@pytest.fixture
def sample_peaks():
    """A small set of peaks for basic testing."""
    return pd.DataFrame({
        "chr": ["chr1", "chr1", "chr2", "chr2", "chr3"],
        "start": [1000, 5000, 2000, 8000, 3000],
        "end": [2000, 6000, 3000, 9000, 4000],
        "name": [f"peak_{i}" for i in range(5)],
        "score": [100, 200, 150, 300, 250],
    })


@pytest.fixture
def sample_peaks_with_signal():
    """Peaks with signal and fold-change columns."""
    np.random.seed(42)
    n = 20
    chroms = np.random.choice(["chr1", "chr2", "chr3"], n)
    starts = np.random.randint(1000, 100000, n)
    ends = starts + np.random.randint(200, 2000, n)
    return pd.DataFrame({
        "chr": chroms,
        "start": starts,
        "end": ends,
        "signal": np.random.uniform(1, 100, n),
        "log2FC": np.random.normal(0, 1.5, n),
        "pvalue": 10 ** -np.random.uniform(0.5, 8, n),
        "FDR": np.random.uniform(0, 1, n),
    })


@pytest.fixture
def overlapping_peaks():
    """Two sets of peaks with known overlaps for testing overlap detection."""
    query = pd.DataFrame({
        "chr": ["chr1", "chr1", "chr2"],
        "start": [100, 500, 200],
        "end": [300, 700, 400],
    })
    subject = pd.DataFrame({
        "chr": ["chr1", "chr1", "chr3"],
        "start": [250, 800, 100],
        "end": [350, 900, 300],
    })
    return query, subject


# ============================================================================
# Temporary files
# ============================================================================


@pytest.fixture
def temp_dir():
    """Provide a temporary directory that gets cleaned up after the test."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_bed_file(temp_dir):
    """Create a temporary BED file for testing."""
    bed_path = temp_dir / "test_peaks.bed"
    bed_content = """chr1\t1000\t2000\tpeak_0\t100\t.
chr1\t5000\t6000\tpeak_1\t200\t.
chr2\t2000\t3000\tpeak_2\t150\t.
chr2\t8000\t9000\tpeak_3\t300\t.
chr3\t3000\t4000\tpeak_4\t250\t."""
    bed_path.write_text(bed_content)
    return bed_path


@pytest.fixture
def sample_narrowpeak_file(temp_dir):
    """Create a temporary narrowPeak file for testing."""
    np_path = temp_dir / "test.narrowPeak"
    np_content = """chr1\t1000\t2000\tpeak_0\t100\t.\t5.5\t3.2\t2.1\t500
chr1\t5000\t6000\tpeak_1\t200\t.\t8.1\t5.4\t4.3\t400
chr2\t2000\t3000\tpeak_2\t150\t.\t6.2\t4.1\t3.0\t450"""
    np_path.write_text(np_content)
    return np_path


@pytest.fixture
def empty_bed_file(temp_dir):
    """Create an empty BED file."""
    bed_path = temp_dir / "empty.bed"
    bed_path.write_text("")
    return bed_path


@pytest.fixture
def malformed_bed_file(temp_dir):
    """Create a BED file with invalid data."""
    bed_path = temp_dir / "malformed.bed"
    bed_content = """chr1\tnot_a_number\t2000
chr2\t3000\talso_bad"""
    bed_path.write_text(bed_content)
    return bed_path


# ============================================================================
# Database fixtures
# ============================================================================


@pytest.fixture
def test_db_url():
    """In-memory SQLite database URL for testing."""
    return "sqlite:///:memory:"


@pytest.fixture
def db_engine(test_db_url):
    """Create a test database engine.

    Uses StaticPool so that all connections share the same in-memory
    SQLite database (otherwise each connection gets its own empty DB).
    """
    from sqlalchemy import create_engine, event
    from sqlalchemy.pool import StaticPool
    from app.models.database import Base

    engine = create_engine(
        test_db_url,
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    Base.metadata.create_all(bind=engine)
    yield engine
    Base.metadata.drop_all(bind=engine)


@pytest.fixture
def db_session(db_engine):
    """Create a test database session."""
    from sqlalchemy.orm import sessionmaker

    Session = sessionmaker(bind=db_engine)
    session = Session()
    yield session
    session.close()


# ============================================================================
# FastAPI test client
# ============================================================================


@pytest.fixture
def api_client(db_engine):
    """Create a FastAPI test client with a test database."""
    from fastapi.testclient import TestClient
    from sqlalchemy.orm import sessionmaker

    from app.main import app, get_db

    TestSession = sessionmaker(bind=db_engine)

    def override_get_db():
        db = TestSession()
        try:
            yield db
        finally:
            db.close()

    app.dependency_overrides[get_db] = override_get_db
    client = TestClient(app)
    yield client
    app.dependency_overrides.clear()
