"""
Unit tests for ATAC-seq analysis module.
"""

import pytest
import numpy as np
import pandas as pd
from pathlib import Path
import sys

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Check for scipy availability
try:
    import scipy

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

from app.core.atacseq import ATACQCMetrics, FragmentSizeDistribution, ATACSeqAnalyzer, generate_demo_atac_data


class TestATACQCMetrics:
    """Tests for ATACQCMetrics dataclass."""

    def test_metrics_creation(self):
        """Test basic metrics creation."""
        metrics = ATACQCMetrics(
            total_reads=80000000,
            mapped_reads=76000000,
            mapping_rate=0.95,
            mitochondrial_rate=0.05,
            duplicate_rate=0.10,
            tss_enrichment=8.5,
            frip=0.35,
            nfr_ratio=0.45,
            peak_count=50000,
            median_fragment_size=150,
        )

        assert metrics.total_reads == 80000000
        assert metrics.mapping_rate == 0.95
        assert metrics.tss_enrichment == 8.5


class TestFragmentSizeDistribution:
    """Tests for FragmentSizeDistribution dataclass."""

    def test_distribution_creation(self):
        """Test distribution creation."""
        dist = FragmentSizeDistribution(
            sizes=np.arange(0, 1000, 10),
            counts=np.random.randint(0, 1000, 100),
            nfr_fraction=0.35,
            mono_fraction=0.30,
            di_fraction=0.15,
            tri_fraction=0.05,
        )

        assert len(dist.sizes) == 100
        assert dist.nfr_fraction == 0.35


class TestATACSeqAnalyzer:
    """Tests for ATACSeqAnalyzer class."""

    @pytest.fixture
    def analyzer(self):
        """Create analyzer instance."""
        return ATACSeqAnalyzer()

    def test_init(self, analyzer):
        """Test analyzer initialization."""
        assert analyzer.nfr_range == (0, 100)
        assert analyzer.mono_range == (180, 247)
        assert analyzer.di_range == (315, 473)

    def test_calculate_qc_metrics(self, analyzer):
        """Test QC metrics calculation."""
        bam_stats = {
            "total_reads": 80000000,
            "mapped_reads": 76000000,
            "mitochondrial_reads": 3800000,
            "duplicates": 7600000,
            "frip": 0.35,
            "nfr_ratio": 0.45,
            "median_fragment": 150,
        }

        peaks = pd.DataFrame({"chr": ["chr1"] * 100, "start": range(100), "end": range(100, 200)})

        tss_scores = np.random.uniform(5, 15, 1000)

        metrics = analyzer.calculate_qc_metrics(bam_stats, peaks, tss_scores)

        assert metrics.total_reads == 80000000
        assert metrics.mapped_reads == 76000000
        assert metrics.mapping_rate == 76000000 / 80000000
        assert metrics.peak_count == 100
        assert metrics.frip == 0.35

    def test_analyze_fragment_sizes(self, analyzer):
        """Test fragment size analysis."""
        # Generate realistic fragment sizes
        np.random.seed(42)
        nfr_fragments = np.random.normal(50, 20, 3500)
        mono_fragments = np.random.normal(200, 30, 3000)
        di_fragments = np.random.normal(400, 40, 2000)
        tri_fragments = np.random.normal(580, 30, 1000)

        fragment_sizes = np.concatenate([nfr_fragments, mono_fragments, di_fragments, tri_fragments])
        fragment_sizes = fragment_sizes[fragment_sizes > 0]

        dist = analyzer.analyze_fragment_sizes(fragment_sizes)

        assert len(dist.sizes) > 0
        assert len(dist.counts) > 0
        assert 0 <= dist.nfr_fraction <= 1
        assert 0 <= dist.mono_fraction <= 1

    def test_calculate_tss_enrichment(self, analyzer):
        """Test TSS enrichment calculation."""
        np.random.seed(42)

        # Simulate coverage at TSS with enrichment at center
        n_tss = 100
        window_size = 2000

        # Create profile with peak at center
        x = np.arange(window_size)
        center = window_size // 2
        base_signal = np.ones(window_size)
        peak_signal = 10 * np.exp(-((x - center) ** 2) / (2 * 50**2))

        coverage = np.zeros((n_tss, window_size))
        for i in range(n_tss):
            noise = np.random.uniform(0.5, 1.5, window_size)
            coverage[i] = (base_signal + peak_signal) * noise

        normalized, score = analyzer.calculate_tss_enrichment(coverage, window_size)

        assert len(normalized) == window_size
        assert score > 1  # Should be enriched

    def test_benjamini_hochberg(self, analyzer):
        """Test FDR correction."""
        pvalues = np.array([0.01, 0.04, 0.05, 0.10, 0.20, 0.50])

        fdr = analyzer._benjamini_hochberg(pvalues)

        assert len(fdr) == len(pvalues)
        assert all(f <= 1.0 for f in fdr)
        assert all(f >= 0.0 for f in fdr)

    def test_integrate_with_histones(self, analyzer):
        """Test integration with histone data."""
        atac_peaks = pd.DataFrame(
            {"chr": ["chr1", "chr1", "chr2"], "start": [1000, 5000, 10000], "end": [2000, 6000, 11000]}
        )

        h3k4me3_peaks = pd.DataFrame({"chr": ["chr1", "chr1"], "start": [1100, 3000], "end": [1900, 4000]})

        h3k27ac_peaks = pd.DataFrame({"chr": ["chr1", "chr2"], "start": [1200, 10500], "end": [1800, 10800]})

        histone_peaks = {"H3K4me3": h3k4me3_peaks, "H3K27ac": h3k27ac_peaks}

        results = analyzer.integrate_with_histones(atac_peaks, histone_peaks)

        assert "H3K4me3_overlap" in results.columns
        assert "H3K27ac_overlap" in results.columns
        assert "chromatin_state" in results.columns

    def test_footprint_analysis(self, analyzer):
        """Test TF footprinting analysis."""
        np.random.seed(42)

        # Create mock coverage with footprints
        coverage = np.ones(10000) * 100
        motif_positions = np.array([1000, 3000, 5000, 7000])

        # Add dips at motif positions
        for pos in motif_positions:
            coverage[pos - 10 : pos + 10] = 50

        result = analyzer.footprint_analysis(coverage, motif_positions, window=100)

        assert "positions" in result
        assert "mean_footprint" in result
        assert "footprint_depth" in result
        assert len(result["positions"]) == 200


class TestDemoDataGeneration:
    """Test demo data generation."""

    def test_generate_demo_atac_data(self):
        """Test that demo data generates without errors."""
        peaks, bam_stats = generate_demo_atac_data()

        assert isinstance(peaks, pd.DataFrame)
        assert isinstance(bam_stats, dict)
        assert len(peaks) == 50000
        assert "chr" in peaks.columns
        assert "start" in peaks.columns
        assert "end" in peaks.columns
        assert "summit" in peaks.columns
        assert "total_reads" in bam_stats


class TestDifferentialAccessibility:
    """Tests for differential accessibility analysis."""

    @pytest.fixture
    def analyzer(self):
        return ATACSeqAnalyzer()

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy required for differential accessibility")
    def test_differential_accessibility(self, analyzer):
        """Test differential accessibility analysis."""
        # Create test count matrix
        np.random.seed(42)
        n_peaks = 100

        # Create conditions with some differential peaks
        ctrl_base = np.random.poisson(100, n_peaks)
        treat_base = ctrl_base.copy()
        treat_base[:20] *= 3  # Up in treatment
        treat_base[20:40] //= 3  # Down in treatment

        count_matrix = pd.DataFrame(
            {
                "ctrl_1": ctrl_base + np.random.poisson(10, n_peaks),
                "ctrl_2": ctrl_base + np.random.poisson(10, n_peaks),
                "treat_1": treat_base + np.random.poisson(10, n_peaks),
                "treat_2": treat_base + np.random.poisson(10, n_peaks),
            },
            index=[f"peak_{i}" for i in range(n_peaks)],
        )

        results = analyzer.differential_accessibility(
            count_matrix, condition1=["treat_1", "treat_2"], condition2=["ctrl_1", "ctrl_2"]
        )

        assert len(results) == n_peaks
        assert "log2FoldChange" in results.columns
        assert "pvalue" in results.columns
        assert "padj" in results.columns


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
