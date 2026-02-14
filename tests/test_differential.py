"""
Unit tests for differential analysis module.
"""

import pytest
import pandas as pd
from pathlib import Path
import sys

# Check for scipy availability
try:
    import scipy

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.core.differential import Sample, DifferentialConfig, DifferentialAnalyzer, DifferentialResults


class TestSample:
    """Tests for Sample dataclass."""

    def test_sample_creation(self):
        """Test basic sample creation."""
        sample = Sample(sample_id="test_sample", condition="Treatment", histone_mark="H3K27ac", replicate=1)
        assert sample.sample_id == "test_sample"
        assert sample.condition == "Treatment"
        assert sample.histone_mark == "H3K27ac"
        assert sample.replicate == 1

    def test_sample_id_string_conversion(self):
        """Test that sample_id is converted to string."""
        sample = Sample(sample_id=123, condition="Control", histone_mark="H3K4me3")
        assert sample.sample_id == "123"
        assert isinstance(sample.sample_id, str)


class TestDifferentialConfig:
    """Tests for DifferentialConfig dataclass."""

    def test_config_defaults(self):
        """Test default configuration values."""
        config = DifferentialConfig(comparison_name="Test", group1="Treatment", group2="Control")
        assert config.fdr_threshold == 0.1
        assert config.lfc_threshold == 0.5
        assert config.min_overlap == 2
        assert config.normalize_method == "RLE"

    def test_config_custom_values(self):
        """Test custom configuration values."""
        config = DifferentialConfig(
            comparison_name="Custom",
            group1="A",
            group2="B",
            fdr_threshold=0.05,
            lfc_threshold=1.0,
            normalize_method="TMM",
        )
        assert config.fdr_threshold == 0.05
        assert config.lfc_threshold == 1.0
        assert config.normalize_method == "TMM"


class TestDifferentialResults:
    """Tests for DifferentialResults dataclass."""

    def test_results_to_dict(self):
        """Test conversion to dictionary."""
        results = DifferentialResults(
            comparison_name="Test", total_peaks=1000, significant_peaks=100, gained_peaks=60, lost_peaks=40
        )
        result_dict = results.to_dict()

        assert result_dict["comparison_name"] == "Test"
        assert result_dict["total_peaks"] == 1000
        assert result_dict["significant_peaks"] == 100
        assert result_dict["gained_peaks"] == 60
        assert result_dict["lost_peaks"] == 40


class TestDifferentialAnalyzer:
    """Tests for DifferentialAnalyzer class."""

    @pytest.fixture
    def analyzer(self):
        """Create analyzer instance."""
        return DifferentialAnalyzer()

    @pytest.fixture
    def temp_peak_file(self, tmp_path):
        """Create a temporary peak file."""
        peak_file = tmp_path / "test_peaks.bed"
        # Create sample narrowPeak data
        data = "chr1\t1000\t2000\tpeak1\t100\t.\t5.0\t3.0\t2.0\t500\n"
        data += "chr1\t5000\t6000\tpeak2\t150\t.\t7.0\t4.0\t3.0\t500\n"
        data += "chr2\t10000\t11000\tpeak3\t200\t.\t10.0\t5.0\t4.0\t500\n"
        peak_file.write_text(data)
        return str(peak_file)

    def test_load_peaks_narrowpeak(self, analyzer, temp_peak_file):
        """Test loading narrowPeak format."""
        peaks = analyzer.load_peaks(temp_peak_file)

        assert len(peaks) == 3
        assert "chrom" in peaks.columns
        assert "start" in peaks.columns
        assert "end" in peaks.columns
        assert "summit" in peaks.columns

    def test_load_peaks_bed(self, analyzer, tmp_path):
        """Test loading simple BED format."""
        bed_file = tmp_path / "test.bed"
        data = "chr1\t1000\t2000\tpeak1\t100\n"
        data += "chr1\t5000\t6000\tpeak2\t150\n"
        bed_file.write_text(data)

        peaks = analyzer.load_peaks(str(bed_file))

        assert len(peaks) == 2
        assert "chrom" in peaks.columns
        assert "summit" in peaks.columns  # Should be computed

    def test_rle_normalization(self, analyzer):
        """Test RLE normalization."""
        # Create test count matrix
        counts = pd.DataFrame(
            {"sample1": [100, 200, 300, 400], "sample2": [110, 220, 280, 420], "sample3": [90, 180, 320, 380]},
            index=["peak1", "peak2", "peak3", "peak4"],
        )

        size_factors = analyzer._rle_normalization(counts)

        assert len(size_factors) == 3
        # Use .values (property) not .values() (method)
        assert all(sf > 0 for sf in size_factors.values)

    def test_tmm_normalization(self, analyzer):
        """Test TMM normalization."""
        counts = pd.DataFrame(
            {"sample1": [100, 200, 300, 400], "sample2": [120, 240, 280, 480], "sample3": [80, 160, 320, 360]},
            index=["peak1", "peak2", "peak3", "peak4"],
        )

        size_factors = analyzer._tmm_normalization(counts)

        assert len(size_factors) == 3
        # Use .values (property) not .values() (method)
        assert all(sf > 0 for sf in size_factors.values)

    def test_merge_peaks_pandas(self, analyzer):
        """Test peak merging with pandas fallback."""
        peaks = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr1", "chr2"],
                "start": [1000, 1100, 5000, 10000],
                "end": [1500, 1600, 5500, 10500],
                "summit": [250, 250, 250, 250],
                "sample": ["s1", "s2", "s1", "s1"],
            }
        )

        consensus = analyzer._merge_peaks_pandas(peaks, min_overlap=1, summit_extend=250)

        assert len(consensus) > 0
        assert "peak_id" in consensus.columns
        assert "n_samples" in consensus.columns


class TestDifferentialFallback:
    """Test fallback differential analysis (when PyDESeq2 not available)."""

    @pytest.fixture
    def analyzer(self):
        return DifferentialAnalyzer()

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy required for fallback differential analysis")
    def test_fallback_differential(self, analyzer):
        """Test scipy-based fallback differential analysis."""
        # Create test data
        counts = pd.DataFrame(
            {
                "treat1": [100, 200, 50, 300, 150],
                "treat2": [120, 180, 60, 280, 160],
                "ctrl1": [50, 210, 100, 150, 140],
                "ctrl2": [60, 190, 90, 160, 150],
            },
            index=["p1", "p2", "p3", "p4", "p5"],
        )

        samples = [
            Sample("treat1", "Treatment", "H3K27ac"),
            Sample("treat2", "Treatment", "H3K27ac"),
            Sample("ctrl1", "Control", "H3K27ac"),
            Sample("ctrl2", "Control", "H3K27ac"),
        ]

        results = analyzer._run_differential_fallback(counts, samples, "Treatment", "Control")

        assert len(results) == 5
        assert "log2FC" in results.columns
        assert "pvalue" in results.columns
        assert "FDR" in results.columns
        assert "direction" in results.columns

        # Check that FDR is properly bounded
        assert all(results["FDR"] >= 0)
        assert all(results["FDR"] <= 1)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
