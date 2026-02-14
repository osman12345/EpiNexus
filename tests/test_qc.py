"""Tests for Quality Control (QC) module."""

import pytest
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import pandas as pd

from app.core.qc import QCMetrics, QCAnalyzer


class TestQCMetrics:
    """Test QCMetrics dataclass."""

    def test_creation_with_defaults(self):
        """Test creating a QCMetrics object with default values."""
        metrics = QCMetrics(sample_name="sample1")

        assert metrics.sample_name == "sample1"
        assert metrics.total_reads == 0
        assert metrics.mapped_reads == 0
        assert metrics.mapping_rate == 0.0
        assert metrics.pass_qc is True
        assert metrics.qc_warnings == []

    def test_creation_with_values(self):
        """Test creating a QCMetrics object with custom values."""
        metrics = QCMetrics(
            sample_name="sample1",
            total_reads=1000000,
            mapped_reads=850000,
            mapping_rate=0.85,
            peak_count=5000,
            frip_score=0.25,
            median_fragment_size=200,
        )

        assert metrics.sample_name == "sample1"
        assert metrics.total_reads == 1000000
        assert metrics.mapped_reads == 850000
        assert metrics.mapping_rate == 0.85
        assert metrics.peak_count == 5000
        assert metrics.frip_score == 0.25
        assert metrics.median_fragment_size == 200

    def test_to_dict(self):
        """Test converting QCMetrics to dictionary."""
        metrics = QCMetrics(
            sample_name="sample1",
            total_reads=1000000,
            mapped_reads=850000,
            mapping_rate=0.85,
            peak_count=5000,
            frip_score=0.25,
            median_fragment_size=200,
            spikein_reads=50000,
            scale_factor=0.2,
            pass_qc=True,
        )

        result = metrics.to_dict()

        assert isinstance(result, dict)
        assert result["sample_name"] == "sample1"
        assert result["total_reads"] == 1000000
        assert result["mapped_reads"] == 850000
        assert result["mapping_rate"] == 0.85
        assert result["peak_count"] == 5000
        assert result["frip_score"] == 0.25
        assert result["median_fragment_size"] == 200
        assert result["spikein_reads"] == 50000
        assert result["scale_factor"] == 0.2
        assert result["pass_qc"] is True

    def test_qc_warnings_list_initialization(self):
        """Test that qc_warnings is properly initialized as an empty list."""
        metrics1 = QCMetrics(sample_name="sample1")
        metrics2 = QCMetrics(sample_name="sample2")

        # Ensure each instance has its own list (not shared)
        metrics1.qc_warnings.append("warning1")
        assert len(metrics1.qc_warnings) == 1
        assert len(metrics2.qc_warnings) == 0


class TestCheckQCThresholds:
    """Test _check_qc_thresholds method."""

    def test_low_mapping_rate_warning(self):
        """Test warning when mapping rate is below threshold."""
        analyzer = QCAnalyzer()
        metrics = QCMetrics(
            sample_name="sample1",
            total_reads=1000000,
            mapped_reads=600000,
            mapping_rate=0.60,  # Below 0.70 threshold
        )

        result = analyzer._check_qc_thresholds(metrics)

        assert result.pass_qc is False
        assert len(result.qc_warnings) == 1
        assert "Low mapping rate" in result.qc_warnings[0]
        assert "60.0%" in result.qc_warnings[0]

    def test_high_duplication_rate_warning(self):
        """Test warning when duplication rate is above threshold."""
        analyzer = QCAnalyzer()
        metrics = QCMetrics(
            sample_name="sample1",
            total_reads=1000000,
            mapped_reads=900000,
            mapping_rate=0.90,
            duplicate_rate=0.40,  # Above 0.30 threshold
        )

        result = analyzer._check_qc_thresholds(metrics)

        assert result.pass_qc is False
        assert len(result.qc_warnings) == 1
        assert "High duplication" in result.qc_warnings[0]
        assert "40.0%" in result.qc_warnings[0]

    def test_zero_frip_with_peaks_warning(self):
        """Test warning when FRiP is zero but peaks exist."""
        analyzer = QCAnalyzer()
        metrics = QCMetrics(
            sample_name="sample1",
            mapping_rate=0.85,
            duplicate_rate=0.10,
            frip_score=0.0,
            peak_count=5000,  # Peaks exist but no reads in them
        )

        result = analyzer._check_qc_thresholds(metrics)

        assert result.pass_qc is False
        assert any("FRiP score is 0.0" in w for w in result.qc_warnings)
        assert any("no reads overlap peaks" in w for w in result.qc_warnings)

    def test_low_frip_warning(self):
        """Test warning when FRiP is low but not zero."""
        analyzer = QCAnalyzer()
        metrics = QCMetrics(
            sample_name="sample1",
            mapping_rate=0.85,
            duplicate_rate=0.10,
            frip_score=0.05,  # Below 0.10 threshold
            peak_count=5000,
        )

        result = analyzer._check_qc_thresholds(metrics)

        assert result.pass_qc is False
        assert any("Low FRiP" in w for w in result.qc_warnings)
        assert any("5.00%" in w for w in result.qc_warnings)

    def test_low_peak_count_warning(self):
        """Test warning when peak count is below threshold."""
        analyzer = QCAnalyzer()
        metrics = QCMetrics(
            sample_name="sample1",
            mapping_rate=0.85,
            duplicate_rate=0.10,
            frip_score=0.15,
            peak_count=500,  # Below 1000 threshold
        )

        result = analyzer._check_qc_thresholds(metrics)

        assert result.pass_qc is False
        assert any("Low peak count" in w for w in result.qc_warnings)
        assert any("500" in w for w in result.qc_warnings)

    def test_small_fragment_size_warning(self):
        """Test warning when fragment size is too small."""
        analyzer = QCAnalyzer()
        metrics = QCMetrics(
            sample_name="sample1",
            mapping_rate=0.85,
            duplicate_rate=0.10,
            frip_score=0.15,
            peak_count=5000,
            median_fragment_size=100,  # Below 150 threshold
        )

        result = analyzer._check_qc_thresholds(metrics)

        assert result.pass_qc is False
        assert any("Small fragments" in w for w in result.qc_warnings)
        assert any("100bp" in w for w in result.qc_warnings)

    def test_large_fragment_size_warning(self):
        """Test warning when fragment size is too large."""
        analyzer = QCAnalyzer()
        metrics = QCMetrics(
            sample_name="sample1",
            mapping_rate=0.85,
            duplicate_rate=0.10,
            frip_score=0.15,
            peak_count=5000,
            median_fragment_size=400,  # Above 300 threshold
        )

        result = analyzer._check_qc_thresholds(metrics)

        assert result.pass_qc is False
        assert any("Large fragments" in w for w in result.qc_warnings)
        assert any("400bp" in w for w in result.qc_warnings)

    def test_optimal_metrics_pass_qc(self):
        """Test that optimal metrics pass QC."""
        analyzer = QCAnalyzer()
        metrics = QCMetrics(
            sample_name="sample1",
            mapping_rate=0.90,
            duplicate_rate=0.10,
            frip_score=0.25,
            peak_count=5000,
            median_fragment_size=200,
        )

        result = analyzer._check_qc_thresholds(metrics)

        assert result.pass_qc is True
        assert len(result.qc_warnings) == 0

    def test_multiple_warnings(self):
        """Test that multiple warnings are accumulated."""
        analyzer = QCAnalyzer()
        metrics = QCMetrics(
            sample_name="sample1",
            mapping_rate=0.60,  # Low
            duplicate_rate=0.40,  # High
            frip_score=0.05,  # Low
            peak_count=500,  # Low
            median_fragment_size=100,  # Small
        )

        result = analyzer._check_qc_thresholds(metrics)

        assert result.pass_qc is False
        assert len(result.qc_warnings) >= 5


class TestCountPeaks:
    """Test count_peaks method."""

    def test_count_peaks_valid_bed(self):
        """Test counting peaks in a valid BED file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            bed_file = Path(tmpdir) / "peaks.bed"

            # Create a simple BED file
            bed_content = """chr1	1000	2000	peak1	100
chr1	3000	4000	peak2	150
chr2	5000	6000	peak3	120
"""
            bed_file.write_text(bed_content)

            analyzer = QCAnalyzer()
            peak_count = analyzer.count_peaks(str(bed_file))

            assert peak_count == 3

    def test_count_peaks_with_comments(self):
        """Test that comment lines are ignored."""
        with tempfile.TemporaryDirectory() as tmpdir:
            bed_file = Path(tmpdir) / "peaks.bed"

            bed_content = """# This is a comment
chr1	1000	2000	peak1	100
# Another comment
chr1	3000	4000	peak2	150
chr2	5000	6000	peak3	120
"""
            bed_file.write_text(bed_content)

            analyzer = QCAnalyzer()
            peak_count = analyzer.count_peaks(str(bed_file))

            assert peak_count == 3

    def test_count_peaks_with_empty_lines(self):
        """Test that empty lines are ignored."""
        with tempfile.TemporaryDirectory() as tmpdir:
            bed_file = Path(tmpdir) / "peaks.bed"

            bed_content = """chr1	1000	2000	peak1	100

chr1	3000	4000	peak2	150

chr2	5000	6000	peak3	120
"""
            bed_file.write_text(bed_content)

            analyzer = QCAnalyzer()
            peak_count = analyzer.count_peaks(str(bed_file))

            assert peak_count == 3

    def test_count_peaks_empty_file(self):
        """Test counting peaks in an empty file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            bed_file = Path(tmpdir) / "peaks.bed"
            bed_file.write_text("")

            analyzer = QCAnalyzer()
            peak_count = analyzer.count_peaks(str(bed_file))

            assert peak_count == 0

    def test_count_peaks_nonexistent_file(self):
        """Test handling of nonexistent file."""
        analyzer = QCAnalyzer()
        peak_count = analyzer.count_peaks("/nonexistent/path/peaks.bed")

        assert peak_count == 0

    def test_count_peaks_large_file(self):
        """Test counting peaks in a large file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            bed_file = Path(tmpdir) / "peaks.bed"

            # Create a file with many peaks
            with open(bed_file, "w") as f:
                for i in range(10000):
                    f.write(f"chr1\t{i * 100}\t{i * 100 + 1000}\tpeak{i}\t{100 + i % 50}\n")

            analyzer = QCAnalyzer()
            peak_count = analyzer.count_peaks(str(bed_file))

            assert peak_count == 10000


class TestGenerateQCReport:
    """Test generate_qc_report method."""

    def test_generate_qc_report_single_sample(self):
        """Test generating a report for a single sample."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = Path(tmpdir) / "qc_report.csv"

            metrics = QCMetrics(
                sample_name="sample1",
                total_reads=1000000,
                mapped_reads=850000,
                mapping_rate=0.85,
                peak_count=5000,
                frip_score=0.25,
                median_fragment_size=200,
                pass_qc=True,
            )

            analyzer = QCAnalyzer()
            df = analyzer.generate_qc_report([metrics], str(output_file))

            # Check output file was created
            assert output_file.exists()

            # Check DataFrame content
            assert len(df) == 1
            assert df.iloc[0]["sample_name"] == "sample1"
            assert df.iloc[0]["qc_status"] == "PASS"

            # Check CSV content
            csv_content = output_file.read_text()
            assert "sample_name" in csv_content
            assert "sample1" in csv_content
            assert "PASS" in csv_content

    def test_generate_qc_report_multiple_samples(self):
        """Test generating a report for multiple samples."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = Path(tmpdir) / "qc_report.csv"

            samples = [
                QCMetrics(
                    sample_name="sample1",
                    total_reads=1000000,
                    mapped_reads=850000,
                    mapping_rate=0.85,
                    peak_count=5000,
                    frip_score=0.25,
                    pass_qc=True,
                ),
                QCMetrics(
                    sample_name="sample2",
                    total_reads=1000000,
                    mapped_reads=600000,
                    mapping_rate=0.60,  # Below threshold
                    peak_count=3000,
                    frip_score=0.15,
                    pass_qc=False,
                    qc_warnings=["Low mapping rate"],
                ),
            ]

            analyzer = QCAnalyzer()
            df = analyzer.generate_qc_report(samples, str(output_file))

            assert len(df) == 2
            assert df.iloc[0]["qc_status"] == "PASS"
            assert df.iloc[1]["qc_status"] == "WARN"

            # Check CSV file exists and is readable
            assert output_file.exists()
            read_df = pd.read_csv(output_file)
            assert len(read_df) == 2

    def test_generate_qc_report_with_spike_in(self):
        """Test generating a report with spike-in data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = Path(tmpdir) / "qc_report.csv"

            metrics = QCMetrics(
                sample_name="sample1",
                total_reads=1000000,
                mapped_reads=850000,
                mapping_rate=0.85,
                peak_count=5000,
                frip_score=0.25,
                spikein_reads=50000,
                scale_factor=0.2,
                pass_qc=True,
            )

            analyzer = QCAnalyzer()
            df = analyzer.generate_qc_report([metrics], str(output_file))

            assert "spikein_reads" in df.columns
            assert "scale_factor" in df.columns
            assert df.iloc[0]["spikein_reads"] == 50000
            assert df.iloc[0]["scale_factor"] == 0.2

    def test_qc_status_mapping(self):
        """Test that qc_status column is correctly mapped."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = Path(tmpdir) / "qc_report.csv"

            metrics_pass = QCMetrics(sample_name="pass_sample", pass_qc=True)
            metrics_warn = QCMetrics(sample_name="warn_sample", pass_qc=False)

            analyzer = QCAnalyzer()
            df = analyzer.generate_qc_report([metrics_pass, metrics_warn], str(output_file))

            assert df.iloc[0]["qc_status"] == "PASS"
            assert df.iloc[1]["qc_status"] == "WARN"


class TestQCAnalyzerInitialization:
    """Test QCAnalyzer initialization."""

    def test_default_initialization(self):
        """Test initializing QCAnalyzer with default paths."""
        analyzer = QCAnalyzer()

        assert analyzer.samtools == "samtools"
        assert analyzer.bedtools == "bedtools"

    def test_custom_initialization(self):
        """Test initializing QCAnalyzer with custom paths."""
        analyzer = QCAnalyzer(samtools_path="/custom/samtools", bedtools_path="/custom/bedtools")

        assert analyzer.samtools == "/custom/samtools"
        assert analyzer.bedtools == "/custom/bedtools"

    def test_thresholds_constants(self):
        """Test that threshold constants are properly defined."""
        analyzer = QCAnalyzer()

        assert analyzer.MIN_MAPPING_RATE == 0.70
        assert analyzer.MAX_DUPLICATE_RATE == 0.30
        assert analyzer.MIN_FRIP == 0.10
        assert analyzer.MIN_PEAK_COUNT == 1000
        assert analyzer.OPTIMAL_FRAGMENT_SIZE == (150, 300)


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_metrics_with_zero_peak_count(self):
        """Test thresholds with zero peak count."""
        analyzer = QCAnalyzer()
        metrics = QCMetrics(
            sample_name="sample1",
            mapping_rate=0.85,
            peak_count=0,  # No peaks
            frip_score=0.0,
        )

        result = analyzer._check_qc_thresholds(metrics)

        # Should not warn about zero FRiP if peak_count is 0
        zero_frip_warnings = [w for w in result.qc_warnings if "FRiP score is 0.0" in w]
        assert len(zero_frip_warnings) == 0

    def test_metrics_with_zero_fragment_size(self):
        """Test thresholds with zero fragment size."""
        analyzer = QCAnalyzer()
        metrics = QCMetrics(
            sample_name="sample1",
            mapping_rate=0.85,
            duplicate_rate=0.10,
            frip_score=0.15,
            peak_count=5000,
            median_fragment_size=0,  # No fragment size calculated
        )

        result = analyzer._check_qc_thresholds(metrics)

        # Should not warn about fragment size if it's 0
        frag_warnings = [w for w in result.qc_warnings if "fragment" in w.lower()]
        assert len(frag_warnings) == 0

    def test_all_metrics_at_threshold_boundaries(self):
        """Test metrics at exact threshold boundaries."""
        analyzer = QCAnalyzer()
        metrics = QCMetrics(
            sample_name="sample1",
            mapping_rate=0.70,  # Exactly at MIN_MAPPING_RATE
            duplicate_rate=0.30,  # Exactly at MAX_DUPLICATE_RATE
            frip_score=0.10,  # Exactly at MIN_FRIP
            peak_count=1000,  # Exactly at MIN_PEAK_COUNT
            median_fragment_size=150,  # At lower optimal bound
        )

        result = analyzer._check_qc_thresholds(metrics)

        # At threshold should pass
        assert result.pass_qc is True
        assert len(result.qc_warnings) == 0

    def test_metrics_just_below_threshold_boundaries(self):
        """Test metrics just below threshold boundaries."""
        analyzer = QCAnalyzer()
        metrics = QCMetrics(
            sample_name="sample1",
            mapping_rate=0.699,  # Just below MIN_MAPPING_RATE
            duplicate_rate=0.300,  # At limit
            frip_score=0.099,  # Just below MIN_FRIP
            peak_count=999,  # Just below MIN_PEAK_COUNT
            median_fragment_size=149,  # Just below lower optimal bound
        )

        result = analyzer._check_qc_thresholds(metrics)

        # Should have warnings
        assert result.pass_qc is False
        assert len(result.qc_warnings) >= 4
