"""
Unit tests for DiffBind-compatible interface.
"""

import pytest
import pandas as pd
import tempfile
from pathlib import Path
import sys

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.core.diffbind import (
    DiffBindConfig,
    DiffBindResults,
    DiffBindRunner,
    run_diffbind
)


class TestDiffBindConfig:
    """Tests for DiffBindConfig dataclass."""

    def test_config_defaults(self):
        """Test default configuration values."""
        config = DiffBindConfig(
            comparison_name="Test",
            group1="Treatment",
            group2="Control",
            histone_mark="H3K27ac"
        )

        assert config.genome == "hg38"
        assert config.fdr_threshold == 0.05
        assert config.lfc_threshold == 1.0
        assert config.min_overlap == 2
        assert config.summit_size == 250
        assert config.normalize_method == "RLE"

    def test_config_custom_values(self):
        """Test custom configuration values."""
        config = DiffBindConfig(
            comparison_name="Custom",
            group1="A",
            group2="B",
            histone_mark="H3K4me3",
            genome="mm10",
            fdr_threshold=0.01,
            lfc_threshold=2.0,
            normalize_method="TMM"
        )

        assert config.genome == "mm10"
        assert config.fdr_threshold == 0.01
        assert config.lfc_threshold == 2.0
        assert config.normalize_method == "TMM"


class TestDiffBindResults:
    """Tests for DiffBindResults dataclass."""

    def test_results_creation(self):
        """Test basic results creation."""
        results = DiffBindResults(
            comparison_name="Test",
            total_peaks=10000,
            significant_peaks=500,
            gained_peaks=300,
            lost_peaks=200
        )

        assert results.comparison_name == "Test"
        assert results.total_peaks == 10000
        assert results.significant_peaks == 500

    def test_results_summary(self):
        """Test that summary is auto-generated."""
        results = DiffBindResults(
            comparison_name="Test",
            total_peaks=10000,
            significant_peaks=500,
            gained_peaks=300,
            lost_peaks=200
        )

        assert "comparison" in results.summary
        assert results.summary["total_peaks"] == 10000
        assert results.summary["gained"] == 300

    def test_results_with_dataframes(self):
        """Test results with DataFrames."""
        all_peaks = pd.DataFrame({
            "peak_id": ["p1", "p2", "p3"],
            "log2FC": [1.5, -2.0, 0.5],
            "FDR": [0.01, 0.001, 0.2]
        })

        sig_peaks = all_peaks[all_peaks["FDR"] < 0.05]

        results = DiffBindResults(
            comparison_name="Test",
            total_peaks=3,
            significant_peaks=2,
            gained_peaks=1,
            lost_peaks=1,
            all_peaks=all_peaks,
            significant_df=sig_peaks
        )

        assert len(results.all_peaks) == 3
        assert len(results.significant_df) == 2


class TestDiffBindRunner:
    """Tests for DiffBindRunner class."""

    @pytest.fixture
    def runner(self):
        """Create runner instance."""
        return DiffBindRunner()

    @pytest.fixture
    def sample_sheet_data(self):
        """Create sample sheet data."""
        return [
            {
                "SampleID": "treat_1",
                "Condition": "Treatment",
                "Factor": "H3K27ac",
                "Replicate": 1
            },
            {
                "SampleID": "treat_2",
                "Condition": "Treatment",
                "Factor": "H3K27ac",
                "Replicate": 2
            },
            {
                "SampleID": "ctrl_1",
                "Condition": "Control",
                "Factor": "H3K27ac",
                "Replicate": 1
            },
            {
                "SampleID": "ctrl_2",
                "Condition": "Control",
                "Factor": "H3K27ac",
                "Replicate": 2
            }
        ]

    def test_create_sample_sheet(self, runner, sample_sheet_data, tmp_path):
        """Test sample sheet creation."""
        output_path = str(tmp_path / "samples.csv")

        result_path = runner.create_sample_sheet(sample_sheet_data, output_path)

        assert Path(result_path).exists()

        # Read and verify
        df = pd.read_csv(result_path)
        assert len(df) == 4
        assert "SampleID" in df.columns
        assert "Condition" in df.columns

    def test_sample_sheet_validation(self, runner, tmp_path):
        """Test sample sheet validation."""
        incomplete_data = [
            {"SampleID": "s1", "Condition": "A"}
            # Missing Factor and Replicate
        ]

        output_path = str(tmp_path / "bad_samples.csv")

        with pytest.raises(ValueError, match="missing columns"):
            runner.create_sample_sheet(incomplete_data, output_path)


class TestConvenienceFunction:
    """Test the run_diffbind convenience function."""

    def test_function_exists(self):
        """Test that the convenience function is importable."""
        assert callable(run_diffbind)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
