"""
Comprehensive tests for the peak gene linking module.

Tests cover:
- PeakGeneLinkingEngine initialization
- Distance-based linking with power-law decay
- ABC (Activity-by-Contact) linking
- Nearest gene linking
- Method dispatching via link_peaks_to_genes
- Edge cases: no nearby genes, identical positions, multiple peaks same gene
"""

import pytest
import pandas as pd
import numpy as np
from app.core.peak_gene_linking import (
    PeakGeneLinkingEngine,
    EnhancerGeneLink,
    load_gene_annotations,
    generate_demo_genes,
)


class TestPeakGeneLinkingEngineInitialization:
    """Test PeakGeneLinkingEngine initialization with default and custom params."""

    def test_init_with_defaults(self):
        """Test initialization with default parameters."""
        engine = PeakGeneLinkingEngine()
        assert engine.max_distance == 1000000  # 1 Mb
        assert engine.activity_column == "signal"
        assert engine.use_hic is False

    def test_init_with_custom_params(self):
        """Test initialization with custom parameters."""
        engine = PeakGeneLinkingEngine(
            max_distance=500000,
            activity_column="h3k27ac",
            use_hic=True,
        )
        assert engine.max_distance == 500000
        assert engine.activity_column == "h3k27ac"
        assert engine.use_hic is True

    def test_init_with_partial_custom_params(self):
        """Test initialization with some custom parameters."""
        engine = PeakGeneLinkingEngine(max_distance=750000)
        assert engine.max_distance == 750000
        assert engine.activity_column == "signal"
        assert engine.use_hic is False


class TestDistanceLinking:
    """Test _distance_linking method with synthetic data and power-law decay."""

    @pytest.fixture
    def simple_peaks_genes(self):
        """Create simple synthetic peaks and genes."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1", "chr1", "chr2"],
                "start": [1000, 50000, 100000],
                "end": [2000, 51000, 101000],
                "peak_id": ["peak_1", "peak_2", "peak_3"],
                "signal": [10.0, 5.0, 8.0],
            }
        )

        genes = pd.DataFrame(
            {
                "chr": ["chr1", "chr1", "chr2"],
                "tss": [5000, 60000, 110000],
                "gene_id": ["ENSG001", "ENSG002", "ENSG003"],
                "gene_symbol": ["GENE1", "GENE2", "GENE3"],
            }
        )

        return peaks, genes

    def test_distance_linking_basic(self, simple_peaks_genes):
        """Test basic distance-based linking."""
        peaks, genes = simple_peaks_genes
        engine = PeakGeneLinkingEngine(activity_column="signal")
        links = engine._distance_linking(peaks, genes, distance_threshold=60000)

        assert len(links) > 0
        assert "distance" in links.columns
        assert "activity" in links.columns
        assert "contact" in links.columns
        assert "method" in links.columns
        assert all(links["method"] == "Distance")

    def test_distance_linking_power_law_contact_score(self, simple_peaks_genes):
        """Verify contact score uses power-law decay with 0.87 exponent."""
        peaks, genes = simple_peaks_genes
        engine = PeakGeneLinkingEngine(activity_column="signal")
        links = engine._distance_linking(peaks, genes, distance_threshold=60000)

        # Check that contact scores are calculated correctly
        # contact = 1.0 / (1.0 + (distance / 5000) ** 0.87)
        for _, link in links.iterrows():
            distance = link["distance"]
            expected_contact = 1.0 / (1.0 + (distance / 5000) ** 0.87)
            assert np.isclose(link["contact"], expected_contact, rtol=1e-5)

    def test_distance_linking_respects_threshold(self, simple_peaks_genes):
        """Test that distance threshold is respected."""
        peaks, genes = simple_peaks_genes
        engine = PeakGeneLinkingEngine(activity_column="signal")

        # Link with small threshold - should get fewer results
        links_small = engine._distance_linking(peaks, genes, distance_threshold=10000)

        # Link with large threshold
        links_large = engine._distance_linking(peaks, genes, distance_threshold=100000)

        assert len(links_small) <= len(links_large)

    def test_distance_linking_same_chromosome_only(self, simple_peaks_genes):
        """Test that linking only occurs on same chromosome."""
        peaks, genes = simple_peaks_genes
        engine = PeakGeneLinkingEngine(activity_column="signal")
        links = engine._distance_linking(peaks, genes, distance_threshold=1000000)

        # Check that peaks and genes are on same chromosome by merging with gene data
        for _, link in links.iterrows():
            gene = genes[genes["gene_id"] == link["gene_id"]].iloc[0]
            assert link["enhancer_chr"] == gene["chr"]

    def test_distance_linking_activity_column(self):
        """Test that custom activity columns are used correctly."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1", "chr1"],
                "start": [1000, 50000],
                "end": [2000, 51000],
                "peak_id": ["peak_1", "peak_2"],
                "h3k27ac": [20.0, 15.0],
            }
        )

        genes = pd.DataFrame(
            {
                "chr": ["chr1"],
                "tss": [5000],
                "gene_id": ["ENSG001"],
                "gene_symbol": ["GENE1"],
            }
        )

        engine = PeakGeneLinkingEngine(activity_column="h3k27ac")
        links = engine._distance_linking(peaks, genes, distance_threshold=60000)

        # Check that h3k27ac values are used as activity
        assert links.iloc[0]["activity"] == 20.0
        assert links.iloc[1]["activity"] == 15.0

    def test_distance_linking_missing_activity_defaults_to_one(self):
        """Test that missing activity column defaults to 1.0."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1"],
                "start": [1000],
                "end": [2000],
                "peak_id": ["peak_1"],
            }
        )

        genes = pd.DataFrame(
            {
                "chr": ["chr1"],
                "tss": [5000],
                "gene_id": ["ENSG001"],
                "gene_symbol": ["GENE1"],
            }
        )

        engine = PeakGeneLinkingEngine(activity_column="missing_column")
        links = engine._distance_linking(peaks, genes, distance_threshold=60000)

        assert len(links) > 0
        assert links.iloc[0]["activity"] == 1.0


class TestABCLinking:
    """Test _abc_linking method with synthetic data and ABC score verification."""

    @pytest.fixture
    def abc_peaks_genes(self):
        """Create synthetic peaks and genes for ABC testing."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1", "chr1", "chr1"],
                "start": [1000, 10000, 20000],
                "end": [2000, 11000, 21000],
                "peak_id": ["peak_1", "peak_2", "peak_3"],
                "signal": [10.0, 5.0, 3.0],
            }
        )

        genes = pd.DataFrame(
            {
                "chr": ["chr1", "chr1"],
                "tss": [5000, 25000],
                "gene_id": ["ENSG001", "ENSG002"],
                "gene_symbol": ["GENE1", "GENE2"],
            }
        )

        return peaks, genes

    def test_abc_linking_basic(self, abc_peaks_genes):
        """Test basic ABC linking."""
        peaks, genes = abc_peaks_genes
        engine = PeakGeneLinkingEngine(activity_column="signal", max_distance=100000)
        links = engine._abc_linking(peaks, genes)

        assert len(links) > 0
        assert "abc_score" in links.columns
        assert "contact" in links.columns
        assert "activity" in links.columns
        assert all(links["method"] == "ABC")

    def test_abc_scores_sum_to_one_per_gene(self, abc_peaks_genes):
        """Verify ABC scores sum to 1.0 for each gene."""
        peaks, genes = abc_peaks_genes
        engine = PeakGeneLinkingEngine(activity_column="signal", max_distance=100000)
        links = engine._abc_linking(peaks, genes)

        # Group by gene and check that ABC scores sum to 1.0
        for gene_id in links["gene_id"].unique():
            gene_links = links[links["gene_id"] == gene_id]
            abc_sum = gene_links["abc_score"].sum()
            # Allow for small floating point errors
            assert np.isclose(abc_sum, 1.0, atol=1e-6), f"ABC scores for {gene_id} sum to {abc_sum}"

    def test_abc_linking_respects_max_distance(self):
        """Test that ABC linking respects max_distance parameter."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1", "chr1"],
                "start": [1000, 1000000],
                "end": [2000, 1001000],
                "peak_id": ["peak_1", "peak_2"],
                "signal": [10.0, 5.0],
            }
        )

        genes = pd.DataFrame(
            {
                "chr": ["chr1"],
                "tss": [5000],
                "gene_id": ["ENSG001"],
                "gene_symbol": ["GENE1"],
            }
        )

        # Small max_distance - should exclude distant peak
        engine_small = PeakGeneLinkingEngine(activity_column="signal", max_distance=100000)
        links_small = engine_small._abc_linking(peaks, genes)

        # Large max_distance - should include distant peak
        engine_large = PeakGeneLinkingEngine(activity_column="signal", max_distance=2000000)
        links_large = engine_large._abc_linking(peaks, genes)

        assert len(links_small) < len(links_large)

    def test_abc_linking_activity_and_contact_product(self, abc_peaks_genes):
        """Verify ABC scores are calculated as (activity × contact) / sum."""
        peaks, genes = abc_peaks_genes
        engine = PeakGeneLinkingEngine(activity_column="signal", max_distance=100000)
        links = engine._abc_linking(peaks, genes)

        # For the first gene, manually verify ABC calculation
        gene1_links = links[links["gene_id"] == "ENSG001"].reset_index(drop=True)

        if len(gene1_links) > 0:
            # Calculate expected ABC scores
            activities = gene1_links["activity"].values
            contacts = gene1_links["contact"].values
            products = activities * contacts
            expected_abc = products / products.sum()

            # Compare with actual ABC scores
            actual_abc = gene1_links["abc_score"].values

            for i in range(len(gene1_links)):
                assert np.isclose(actual_abc[i], expected_abc[i], rtol=1e-5)

    def test_abc_linking_no_nearby_peaks(self):
        """Test ABC linking when gene has no nearby peaks."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1"],
                "start": [1000],
                "end": [2000],
                "peak_id": ["peak_1"],
                "signal": [10.0],
            }
        )

        genes = pd.DataFrame(
            {
                "chr": ["chr1", "chr1"],
                "tss": [5000, 2000000],
                "gene_id": ["ENSG001", "ENSG002"],
                "gene_symbol": ["GENE1", "GENE2"],
            }
        )

        engine = PeakGeneLinkingEngine(activity_column="signal", max_distance=100000)
        links = engine._abc_linking(peaks, genes)

        # Should only have links for GENE1, not GENE2 (too far away)
        assert len(links) > 0
        assert "ENSG002" not in links["gene_id"].values or len(links[links["gene_id"] == "ENSG002"]) == 0


class TestNearestGeneLinking:
    """Test _nearest_gene_linking method."""

    @pytest.fixture
    def nearest_peaks_genes(self):
        """Create synthetic peaks and genes for nearest gene testing."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1", "chr1", "chr2"],
                "start": [10000, 50000, 100000],
                "end": [11000, 51000, 101000],
                "peak_id": ["peak_1", "peak_2", "peak_3"],
                "signal": [10.0, 5.0, 8.0],
            }
        )

        genes = pd.DataFrame(
            {
                "chr": ["chr1", "chr1", "chr1", "chr2", "chr2"],
                "tss": [5000, 15000, 60000, 90000, 110000],
                "gene_id": ["ENSG001", "ENSG002", "ENSG003", "ENSG004", "ENSG005"],
                "gene_symbol": ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"],
            }
        )

        return peaks, genes

    def test_nearest_gene_linking_basic(self, nearest_peaks_genes):
        """Test basic nearest gene linking."""
        peaks, genes = nearest_peaks_genes
        engine = PeakGeneLinkingEngine()
        links = engine._nearest_gene_linking(peaks, genes)

        assert len(links) > 0
        assert len(links) == len(peaks)  # One link per peak
        assert all(links["method"] == "Nearest")

    def test_nearest_gene_each_peak_one_gene(self, nearest_peaks_genes):
        """Verify each peak links to exactly one nearest gene."""
        peaks, genes = nearest_peaks_genes
        engine = PeakGeneLinkingEngine()
        links = engine._nearest_gene_linking(peaks, genes)

        # Group by peak and verify each has exactly one gene
        for peak_id in links["enhancer_id"].unique():
            peak_links = links[links["enhancer_id"] == peak_id]
            assert len(peak_links) == 1

    def test_nearest_gene_is_closest(self, nearest_peaks_genes):
        """Verify that the linked gene is indeed the nearest one."""
        peaks, genes = nearest_peaks_genes
        engine = PeakGeneLinkingEngine()
        links = engine._nearest_gene_linking(peaks, genes)

        for _, link in links.iterrows():
            peak_center = (link["enhancer_start"] + link["enhancer_end"]) // 2
            gene_tss = link["gene_tss"]
            distance = abs(peak_center - gene_tss)

            # Verify this is the closest gene
            peak_chr = link["enhancer_chr"]
            chr_genes = genes[genes["chr"] == peak_chr]

            # Find the actual closest distance
            min_distance = np.abs(chr_genes["tss"].values - peak_center).min()

            assert distance == min_distance

    def test_nearest_gene_no_genes_on_chromosome(self):
        """Test nearest gene linking when no genes exist on peak chromosome."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1", "chr2"],
                "start": [10000, 100000],
                "end": [11000, 101000],
                "peak_id": ["peak_1", "peak_2"],
                "signal": [10.0, 8.0],
            }
        )

        genes = pd.DataFrame(
            {
                "chr": ["chr3"],
                "tss": [5000],
                "gene_id": ["ENSG001"],
                "gene_symbol": ["GENE1"],
            }
        )

        engine = PeakGeneLinkingEngine()
        links = engine._nearest_gene_linking(peaks, genes)

        # No links should be created
        assert len(links) == 0

    def test_nearest_gene_contact_score_one(self, nearest_peaks_genes):
        """Verify contact score is always 1.0 for nearest gene linking."""
        peaks, genes = nearest_peaks_genes
        engine = PeakGeneLinkingEngine()
        links = engine._nearest_gene_linking(peaks, genes)

        assert all(links["contact"] == 1.0)


class TestMethodDispatching:
    """Test link_peaks_to_genes method dispatching."""

    @pytest.fixture
    def test_data(self):
        """Create test peaks and genes."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1", "chr1"],
                "start": [10000, 50000],
                "end": [11000, 51000],
                "peak_id": ["peak_1", "peak_2"],
                "signal": [10.0, 5.0],
            }
        )

        genes = pd.DataFrame(
            {
                "chr": ["chr1", "chr1"],
                "tss": [15000, 60000],
                "gene_id": ["ENSG001", "ENSG002"],
                "gene_symbol": ["GENE1", "GENE2"],
            }
        )

        return peaks, genes

    def test_dispatch_to_abc_method(self, test_data):
        """Test dispatching to ABC method."""
        peaks, genes = test_data
        engine = PeakGeneLinkingEngine(activity_column="signal")
        links = engine.link_peaks_to_genes(peaks, genes, method="abc")

        assert len(links) > 0
        assert all(links["method"] == "ABC")

    def test_dispatch_to_distance_method(self, test_data):
        """Test dispatching to distance method."""
        peaks, genes = test_data
        engine = PeakGeneLinkingEngine(activity_column="signal")
        links = engine.link_peaks_to_genes(peaks, genes, method="distance")

        assert len(links) > 0
        assert all(links["method"] == "Distance")

    def test_dispatch_to_nearest_method(self, test_data):
        """Test dispatching to nearest method."""
        peaks, genes = test_data
        engine = PeakGeneLinkingEngine()
        links = engine.link_peaks_to_genes(peaks, genes, method="nearest")

        assert len(links) > 0
        assert all(links["method"] == "Nearest")

    def test_dispatch_invalid_method_raises_error(self, test_data):
        """Test that invalid method raises ValueError."""
        peaks, genes = test_data
        engine = PeakGeneLinkingEngine()

        with pytest.raises(ValueError, match="Unknown method"):
            engine.link_peaks_to_genes(peaks, genes, method="invalid_method")

    def test_dispatch_returns_dataframe(self, test_data):
        """Test that all methods return DataFrame."""
        peaks, genes = test_data
        engine = PeakGeneLinkingEngine(activity_column="signal")

        for method in ["abc", "distance", "nearest"]:
            links = engine.link_peaks_to_genes(peaks, genes, method=method)
            assert isinstance(links, pd.DataFrame)


class TestEdgeCases:
    """Test edge cases: no nearby genes, identical positions, multiple peaks same gene."""

    def test_no_nearby_genes(self):
        """Test linking when no genes are nearby."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1"],
                "start": [1000],
                "end": [2000],
                "peak_id": ["peak_1"],
                "signal": [10.0],
            }
        )

        genes = pd.DataFrame(
            {
                "chr": ["chr1"],
                "tss": [2000000],
                "gene_id": ["ENSG001"],
                "gene_symbol": ["GENE1"],
            }
        )

        engine = PeakGeneLinkingEngine(activity_column="signal", max_distance=100000)
        links_abc = engine._abc_linking(peaks, genes)
        links_distance = engine._distance_linking(peaks, genes, distance_threshold=100000)

        # ABC linking should skip this gene since peak is too far
        assert len(links_abc) == 0
        # Distance linking also should skip it
        assert len(links_distance) == 0

    def test_identical_peak_and_gene_position(self):
        """Test linking when peak and gene are at identical positions."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1"],
                "start": [10000],
                "end": [11000],
                "peak_id": ["peak_1"],
                "signal": [10.0],
            }
        )

        genes = pd.DataFrame(
            {
                "chr": ["chr1"],
                "tss": [10500],  # Center of peak
                "gene_id": ["ENSG001"],
                "gene_symbol": ["GENE1"],
            }
        )

        engine = PeakGeneLinkingEngine(activity_column="signal")
        links = engine._distance_linking(peaks, genes, distance_threshold=50000)

        assert len(links) > 0
        assert links.iloc[0]["distance"] == 0

    def test_multiple_peaks_same_gene(self):
        """Test linking multiple peaks to the same gene."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1", "chr1", "chr1"],
                "start": [10000, 15000, 20000],
                "end": [11000, 16000, 21000],
                "peak_id": ["peak_1", "peak_2", "peak_3"],
                "signal": [10.0, 5.0, 3.0],
            }
        )

        genes = pd.DataFrame(
            {
                "chr": ["chr1"],
                "tss": [18000],
                "gene_id": ["ENSG001"],
                "gene_symbol": ["GENE1"],
            }
        )

        engine = PeakGeneLinkingEngine(activity_column="signal", max_distance=100000)
        links = engine._abc_linking(peaks, genes)

        # All three peaks should link to the same gene
        assert len(links) == 3
        assert all(links["gene_id"] == "ENSG001")

        # ABC scores should sum to 1.0
        assert np.isclose(links["abc_score"].sum(), 1.0)

    def test_multiple_peaks_same_gene_different_methods(self):
        """Test multiple peaks linking to same gene with different methods."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1", "chr1", "chr1"],
                "start": [10000, 15000, 20000],
                "end": [11000, 16000, 21000],
                "peak_id": ["peak_1", "peak_2", "peak_3"],
                "signal": [10.0, 5.0, 3.0],
            }
        )

        genes = pd.DataFrame(
            {
                "chr": ["chr1"],
                "tss": [18000],
                "gene_id": ["ENSG001"],
                "gene_symbol": ["GENE1"],
            }
        )

        engine = PeakGeneLinkingEngine(activity_column="signal", max_distance=100000)

        # Distance method
        links_distance = engine._distance_linking(peaks, genes, distance_threshold=100000)
        assert len(links_distance) == 3
        assert all(links_distance["gene_id"] == "ENSG001")

        # ABC method
        links_abc = engine._abc_linking(peaks, genes)
        assert len(links_abc) == 3
        assert all(links_abc["gene_id"] == "ENSG001")

        # Nearest method - should only link to GENE1 once per peak
        links_nearest = engine._nearest_gene_linking(peaks, genes)
        assert len(links_nearest) == 3
        assert all(links_nearest["gene_id"] == "ENSG001")

    def test_peak_equidistant_from_two_genes(self):
        """Test peak that is equidistant from two genes."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1"],
                "start": [15000],
                "end": [16000],
                "peak_id": ["peak_1"],
                "signal": [10.0],
            }
        )

        genes = pd.DataFrame(
            {
                "chr": ["chr1", "chr1"],
                "tss": [10000, 20000],  # Equidistant from peak at 15500
                "gene_id": ["ENSG001", "ENSG002"],
                "gene_symbol": ["GENE1", "GENE2"],
            }
        )

        engine = PeakGeneLinkingEngine()
        links = engine._nearest_gene_linking(peaks, genes)

        # Should link to exactly one gene (first encountered with minimum distance)
        assert len(links) == 1

    def test_single_peak_single_gene(self):
        """Test with minimal data: single peak and single gene."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1"],
                "start": [10000],
                "end": [11000],
                "peak_id": ["peak_1"],
                "signal": [10.0],
            }
        )

        genes = pd.DataFrame(
            {
                "chr": ["chr1"],
                "tss": [15000],
                "gene_id": ["ENSG001"],
                "gene_symbol": ["GENE1"],
            }
        )

        engine = PeakGeneLinkingEngine(activity_column="signal")

        links_abc = engine._abc_linking(peaks, genes)
        links_distance = engine._distance_linking(peaks, genes, distance_threshold=50000)
        links_nearest = engine._nearest_gene_linking(peaks, genes)

        # All methods should produce exactly one link
        assert len(links_abc) == 1
        assert len(links_distance) == 1
        assert len(links_nearest) == 1

        # ABC score should be 1.0 (only link)
        assert links_abc.iloc[0]["abc_score"] == 1.0


class TestEnhancerGeneLinkDataclass:
    """Test EnhancerGeneLink dataclass."""

    def test_enhancer_gene_link_creation(self):
        """Test creating EnhancerGeneLink instance."""
        link = EnhancerGeneLink(
            enhancer_id="peak_1",
            enhancer_chr="chr1",
            enhancer_start=10000,
            enhancer_end=11000,
            gene_id="ENSG001",
            gene_symbol="GENE1",
            gene_tss=15000,
            distance=4500,
            activity=10.0,
            contact=0.9,
            abc_score=0.5,
            method="ABC",
        )

        assert link.enhancer_id == "peak_1"
        assert link.gene_symbol == "GENE1"
        assert link.distance == 4500


class TestFilterAndSummarize:
    """Test filter_links and summarize_links methods."""

    @pytest.fixture
    def sample_links(self):
        """Create sample links DataFrame."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1", "chr1", "chr1"],
                "start": [10000, 15000, 20000],
                "end": [11000, 16000, 21000],
                "peak_id": ["peak_1", "peak_2", "peak_3"],
                "signal": [10.0, 5.0, 3.0],
            }
        )

        genes = pd.DataFrame(
            {
                "chr": ["chr1", "chr1"],
                "tss": [18000, 50000],
                "gene_id": ["ENSG001", "ENSG002"],
                "gene_symbol": ["GENE1", "GENE2"],
            }
        )

        engine = PeakGeneLinkingEngine(activity_column="signal", max_distance=100000)
        return engine._abc_linking(peaks, genes)

    def test_filter_links_by_abc_score(self, sample_links):
        """Test filtering links by ABC score."""
        engine = PeakGeneLinkingEngine()
        filtered = engine.filter_links(sample_links, min_abc_score=0.1)

        assert all(filtered["abc_score"] >= 0.1)

    def test_filter_links_by_distance(self, sample_links):
        """Test filtering links by distance."""
        engine = PeakGeneLinkingEngine()
        filtered = engine.filter_links(sample_links, max_distance=10000)

        assert all(filtered["distance"] <= 10000)

    def test_summarize_links_statistics(self, sample_links):
        """Test link summarization statistics."""
        engine = PeakGeneLinkingEngine()
        summary = engine.summarize_links(sample_links)

        assert "total_links" in summary
        assert "unique_enhancers" in summary
        assert "unique_genes" in summary
        assert "mean_abc_score" in summary
        assert "median_distance" in summary

        assert summary["total_links"] == len(sample_links)
        assert isinstance(summary["mean_abc_score"], float)


class TestHelperFunctions:
    """Test helper functions for loading annotations and generating demo data."""

    def test_generate_demo_genes(self):
        """Test demo gene generation."""
        genes = generate_demo_genes()

        assert len(genes) > 0
        assert "gene_id" in genes.columns
        assert "gene_symbol" in genes.columns
        assert "chr" in genes.columns
        assert "tss" in genes.columns
        assert "strand" in genes.columns

    def test_load_gene_annotations_fallback(self):
        """Test that load_gene_annotations falls back to demo genes."""
        genes = load_gene_annotations()

        assert len(genes) > 0
        assert "gene_symbol" in genes.columns
        assert "chr" in genes.columns


class TestPowerLawConsistency:
    """Test power-law consistency across methods."""

    def test_power_law_exponent_consistency(self):
        """Verify 0.87 exponent is used consistently."""
        peaks = pd.DataFrame(
            {
                "chr": ["chr1", "chr1"],
                "start": [10000, 20000],
                "end": [11000, 21000],
                "peak_id": ["peak_1", "peak_2"],
                "signal": [10.0, 5.0],
            }
        )

        genes = pd.DataFrame(
            {
                "chr": ["chr1"],
                "tss": [30000],
                "gene_id": ["ENSG001"],
                "gene_symbol": ["GENE1"],
            }
        )

        engine = PeakGeneLinkingEngine(activity_column="signal", max_distance=100000)

        # Get links from both methods
        links_abc = engine._abc_linking(peaks, genes)
        links_distance = engine._distance_linking(peaks, genes, distance_threshold=100000)

        # Contact scores should match for same peak-gene pairs
        for method_links in [links_abc, links_distance]:
            for _, link in method_links.iterrows():
                distance = link["distance"]
                expected_contact = 1.0 / (1.0 + (distance / 5000) ** 0.87)
                assert np.isclose(link["contact"], expected_contact, rtol=1e-5)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
