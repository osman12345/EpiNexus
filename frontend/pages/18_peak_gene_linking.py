"""
Peak-to-Gene Linking Page

Link enhancers to target genes using ABC model and other methods.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from app.core.peak_gene_linking import PeakGeneLinkingEngine, generate_demo_genes, load_gene_annotations

st.set_page_config(page_title="Peak-Gene Linking - EpiNexus", page_icon="🔗", layout="wide")

# Try to import data manager
try:
    from frontend.components.data_manager import DataManager

    HAS_DATA_MANAGER = True
except ImportError:
    HAS_DATA_MANAGER = False

# Try to import workflow manager
try:
    from frontend.components.workflow_manager import WorkflowManager

    HAS_WORKFLOW_MANAGER = True
except ImportError:
    try:
        from components.workflow_manager import WorkflowManager

        HAS_WORKFLOW_MANAGER = True
    except ImportError:
        HAS_WORKFLOW_MANAGER = False


def has_data():
    """Check if user has loaded peak data."""
    if HAS_DATA_MANAGER:
        peaks = DataManager.get_data("peaks")
        return peaks is not None and len(peaks) > 0
    return len(st.session_state.get("samples", [])) > 0


def render_empty_state():
    """Show empty state when no data is loaded."""
    st.markdown("---")
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        st.markdown(
            """
        <div style="text-align: center; padding: 3rem; background: #f8f9fa;
                    border-radius: 12px; border: 2px dashed #dee2e6;">
            <div style="font-size: 3rem; margin-bottom: 1rem;">🔗</div>
            <h2 style="color: #6c757d; margin-bottom: 0.5rem;">No Data Loaded</h2>
            <p style="color: #6c757d; font-size: 1.1rem;">
                Upload your enhancer peaks to link them to target genes.
            </p>
        </div>
        """,
            unsafe_allow_html=True,
        )
        st.markdown("")
        if st.button("📁 Go to Data & Project", type="primary", use_container_width=True):
            st.switch_page("pages/01_data_project.py")
        st.markdown("")
        st.markdown("**Linking methods:**")
        st.markdown("- ABC Model (Activity-by-Contact)")
        st.markdown("- Distance-based scoring")
        st.markdown("- Nearest gene assignment")


def main():
    st.title("🔗 Peak-to-Gene Linking")
    st.markdown("""
    Predict which genes are regulated by each enhancer using the Activity-by-Contact (ABC)
    model or other linking methods.
    """)

    # Check if data is loaded
    if not has_data():
        render_empty_state()
        return

    # Sidebar settings
    with st.sidebar:
        st.header("Settings")

        method = st.selectbox(
            "Linking Method",
            ["abc", "distance", "nearest"],
            format_func=lambda x: {
                "abc": "ABC Model (Recommended)",
                "distance": "Distance-based",
                "nearest": "Nearest Gene",
            }.get(x),
        )

        max_distance = st.slider("Max Distance (kb)", 10, 2000, 500) * 1000

        if method == "abc":
            min_abc = st.slider("Min ABC Score", 0.0, 0.2, 0.02)
        else:
            min_abc = 0

        st.markdown("---")

        st.markdown("""
        **ABC Model:**
        Score = (Activity × Contact) / Σ

        Activity = Enhancer signal
        Contact = 3D proximity (Hi-C or distance)
        """)

    tab1, tab2, tab3, tab4 = st.tabs(["📤 Input Data", "🔗 Linking Results", "📊 Visualization", "📋 Gene Summary"])

    with tab1:
        peaks, genes = render_input_data()

    with tab2:
        if peaks is not None and genes is not None:
            links = render_linking_results(peaks, genes, method, max_distance, min_abc)
        else:
            st.info("Upload or generate data in the Input Data tab first.")
            links = None

    with tab3:
        if links is not None and len(links) > 0:
            render_visualizations(links)
        else:
            st.info("Run linking first to see visualizations.")

    with tab4:
        if links is not None and len(links) > 0:
            render_gene_summary(links)
        else:
            st.info("Run linking first to see gene summary.")


def render_input_data():
    """Handle input data upload."""
    st.header("Input Data")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Enhancer Peaks")

        # Try to get peaks from session state
        peaks = None
        if HAS_DATA_MANAGER:
            peaks = DataManager.get_data("peaks")
        if peaks is None and "uploaded_peaks" in st.session_state:
            peaks = st.session_state.uploaded_peaks
        if peaks is None and "samples" in st.session_state:
            for sample in st.session_state.samples:
                if "peaks" in sample:
                    peaks = sample["peaks"].copy()
                    break

        if peaks is not None:
            st.success(f"Using {len(peaks)} peaks from session")
        else:
            peak_file = st.file_uploader("Upload peak file", type=["bed", "csv", "tsv"], key="peak_upload")
            if peak_file:
                peaks = pd.read_csv(peak_file, sep="\t", header=None)
                peaks.columns = ["chr", "start", "end", "name", "signal"][: len(peaks.columns)]
                st.success(f"Loaded {len(peaks)} peaks")

        if peaks is not None:
            st.dataframe(peaks.head(10), use_container_width=True)

    with col2:
        st.subheader("Gene Annotations")

        gene_source = st.radio(
            "Gene source", ["Built-in (GENCODE)", "Demo genes", "Upload custom"], horizontal=True, key="gene_source"
        )

        if gene_source == "Built-in (GENCODE)":
            genome = st.selectbox("Genome", ["hg38", "hg19", "mm10", "mm39"], key="gene_genome")
            with st.spinner("Loading gene annotations..."):
                genes = load_gene_annotations(genome)
            st.success(f"Loaded {len(genes)} curated genes ({genome})")

        elif gene_source == "Demo genes":
            genes = generate_demo_genes()
            st.success(f"Loaded {len(genes)} demo genes")

        else:
            gene_file = st.file_uploader("Upload gene annotation", type=["bed", "csv", "tsv", "gtf"], key="gene_upload")
            if gene_file:
                genes = pd.read_csv(gene_file, sep="\t")
                # Standardize column names if needed
                col_map = {
                    "chrom": "chr",
                    "chromosome": "chr",
                    "txStart": "tss",
                    "name2": "gene_symbol",
                    "geneName": "gene_symbol",
                    "name": "gene_id",
                }
                genes = genes.rename(columns={k: v for k, v in col_map.items() if k in genes.columns})
                st.success(f"Loaded {len(genes)} genes")
            else:
                genes = None

        if genes is not None:
            st.dataframe(genes.head(10), use_container_width=True)

    return peaks, genes


def render_linking_results(peaks, genes, method, max_distance, min_abc):
    """Run and display linking results."""
    st.header("Linking Results")

    if st.button("🔗 Run Peak-Gene Linking", type="primary"):
        with st.spinner("Linking peaks to genes..."):
            engine = PeakGeneLinkingEngine(max_distance=max_distance, activity_column="signal")

            links = engine.link_peaks_to_genes(peaks, genes, method=method)

            # Filter links
            if min_abc > 0 and method == "abc":
                links = engine.filter_links(links, min_abc_score=min_abc)

            st.session_state.peak_gene_links = links

            # Summary stats
            summary = engine.summarize_links(links)

            col1, col2, col3, col4 = st.columns(4)
            col1.metric("Total Links", f"{summary['total_links']:,}")
            col2.metric("Unique Enhancers", f"{summary['unique_enhancers']:,}")
            col3.metric("Unique Genes", f"{summary['unique_genes']:,}")
            col4.metric("Mean ABC Score", f"{summary['mean_abc_score']:.3f}")

            # Record workflow step
            if HAS_WORKFLOW_MANAGER:
                WorkflowManager.record_step(
                    step_name="Peak-Gene Linking",
                    tool="EpiNexus ABC Model",
                    parameters={
                        "method": method,
                        "max_distance": max_distance,
                        "min_abc_score": min_abc,
                        "total_links": summary["total_links"],
                        "unique_genes": summary["unique_genes"],
                    },
                    inputs=["peaks", "genes"],
                    outputs=["peak_gene_links"],
                )

    # Display links from session state
    links = st.session_state.get("peak_gene_links", None)

    if links is not None and len(links) > 0:
        st.subheader("Top Links")

        # Filter options
        col1, col2, col3 = st.columns(3)

        with col1:
            gene_filter = st.text_input("Filter by gene", "")

        with col2:
            max_dist_filter = st.number_input("Max distance (kb)", 0, 2000, 500) * 1000

        with col3:
            sort_by = st.selectbox("Sort by", ["abc_score", "distance", "activity"])

        # Apply filters
        filtered = links.copy()
        if gene_filter:
            filtered = filtered[filtered["gene_symbol"].str.contains(gene_filter, case=False)]
        filtered = filtered[filtered["distance"] <= max_dist_filter]
        filtered = filtered.sort_values(sort_by, ascending=(sort_by == "distance"))

        # Display
        display_cols = ["enhancer_id", "gene_symbol", "distance", "activity", "abc_score", "method"]
        st.dataframe(filtered[display_cols].head(100).round(4), use_container_width=True, hide_index=True)

        # Download
        csv = filtered.to_csv(index=False)
        st.download_button("📥 Download Links", csv, "peak_gene_links.csv", "text/csv")

    return links


def render_visualizations(links):
    """Visualize linking results."""
    st.header("Visualizations")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("ABC Score Distribution")

        fig = px.histogram(links, x="abc_score", nbins=50, title="Distribution of ABC Scores")
        fig.add_vline(x=0.02, line_dash="dash", line_color="red", annotation_text="Typical threshold")
        fig.update_layout(height=400)
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("Distance Distribution")

        fig = px.histogram(links, x=links["distance"] / 1000, nbins=50, title="Distance to Gene TSS")
        fig.update_layout(xaxis_title="Distance (kb)", height=400)
        st.plotly_chart(fig, use_container_width=True)

    # ABC vs Distance scatter
    st.subheader("ABC Score vs Distance")

    fig = px.scatter(
        links,
        x=links["distance"] / 1000,
        y="abc_score",
        color="activity",
        hover_data=["gene_symbol", "enhancer_id"],
        opacity=0.6,
        title="ABC Score vs Distance (colored by activity)",
    )
    fig.update_layout(xaxis_title="Distance (kb)", yaxis_title="ABC Score", height=500)
    st.plotly_chart(fig, use_container_width=True)

    # Links per gene
    st.subheader("Enhancer Links per Gene")

    links_per_gene = links.groupby("gene_symbol").size().sort_values(ascending=False).head(20)

    fig = px.bar(
        x=links_per_gene.values,
        y=links_per_gene.index,
        orientation="h",
        title="Top Genes by Number of Linked Enhancers",
    )
    fig.update_layout(xaxis_title="Number of Enhancers", yaxis_title="Gene", height=500)
    st.plotly_chart(fig, use_container_width=True)


def render_gene_summary(links):
    """Gene-centric summary of links."""
    st.header("Gene Summary")

    # Aggregate by gene
    gene_summary = (
        links.groupby("gene_symbol")
        .agg({"enhancer_id": "count", "abc_score": ["sum", "max", "mean"], "distance": "min", "activity": "sum"})
        .round(4)
    )

    gene_summary.columns = ["n_enhancers", "total_abc", "max_abc", "mean_abc", "min_distance", "total_activity"]
    gene_summary = gene_summary.sort_values("total_abc", ascending=False).reset_index()

    st.dataframe(gene_summary.head(50), use_container_width=True, hide_index=True)

    # Gene detail view
    st.subheader("Gene Detail View")

    selected_gene = st.selectbox("Select gene", gene_summary["gene_symbol"].tolist())

    if selected_gene:
        gene_links = links[links["gene_symbol"] == selected_gene].sort_values("abc_score", ascending=False)

        st.markdown(f"**{selected_gene}** has {len(gene_links)} linked enhancers")

        # Show enhancers for this gene
        fig = go.Figure()

        # Gene TSS (center point)
        gene_tss = gene_links["gene_tss"].iloc[0]

        # Plot enhancers relative to TSS
        for _, link in gene_links.iterrows():
            enhancer_center = (link["enhancer_start"] + link["enhancer_end"]) / 2
            relative_pos = (enhancer_center - gene_tss) / 1000  # kb

            fig.add_trace(
                go.Scatter(
                    x=[relative_pos],
                    y=[link["abc_score"]],
                    mode="markers",
                    marker=dict(
                        size=link["activity"] * 2 + 5, color=link["abc_score"], colorscale="Reds", showscale=True
                    ),
                    name=link["enhancer_id"],
                    hovertemplate=f"Distance: {relative_pos:.1f} kb<br>ABC: {link['abc_score']:.3f}",
                )
            )

        fig.add_vline(x=0, line_dash="dash", annotation_text="TSS")

        fig.update_layout(
            title=f"Enhancers linked to {selected_gene}",
            xaxis_title="Distance from TSS (kb)",
            yaxis_title="ABC Score",
            showlegend=False,
            height=400,
        )

        st.plotly_chart(fig, use_container_width=True)


def generate_demo_peaks(n: int = 500, genome: str = "hg38"):
    """
    Generate demo enhancer peak data based on realistic characteristics.

    Creates peaks near known gene loci with realistic signal distributions
    typical of H3K27ac ChIP-seq data at enhancers.
    """
    np.random.seed(42)

    # Known enhancer regions near important genes (hg38 coordinates)
    known_enhancers = [
        # MYC super-enhancer region
        {"chr": "chr8", "region_start": 128700000, "region_end": 128900000, "n_peaks": 20},
        # TP53 regulatory region
        {"chr": "chr17", "region_start": 7500000, "region_end": 7800000, "n_peaks": 15},
        # EGFR regulatory region
        {"chr": "chr7", "region_start": 55000000, "region_end": 55300000, "n_peaks": 12},
        # VEGFA regulatory region
        {"chr": "chr6", "region_start": 43600000, "region_end": 43900000, "n_peaks": 10},
        # CCND1 regulatory region
        {"chr": "chr11", "region_start": 69500000, "region_end": 69800000, "n_peaks": 8},
        # BCL2 regulatory region
        {"chr": "chr18", "region_start": 63000000, "region_end": 63300000, "n_peaks": 10},
    ]

    peaks_list = []
    peak_id = 0

    # Generate peaks in known enhancer regions
    for region in known_enhancers:
        region_size = region["region_end"] - region["region_start"]
        for _ in range(region["n_peaks"]):
            start = region["region_start"] + np.random.randint(0, region_size - 2000)
            width = np.random.randint(500, 3000)
            # Higher signals for super-enhancer regions
            signal = np.random.lognormal(3, 0.8) if "MYC" in str(region) else np.random.lognormal(2, 1)

            peaks_list.append(
                {
                    "chr": region["chr"],
                    "start": start,
                    "end": start + width,
                    "peak_id": f"enhancer_{peak_id}",
                    "signal": signal,
                }
            )
            peak_id += 1

    # Generate random peaks across the genome for the remaining count
    chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    chr_weights = np.array([8, 7, 6, 5, 5, 5, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 3])
    chr_weights = chr_weights / chr_weights.sum()

    remaining = n - len(peaks_list)
    for i in range(remaining):
        chrom = np.random.choice(chromosomes, p=chr_weights)
        start = np.random.randint(1000000, 200000000)
        width = np.random.randint(500, 3000)

        peaks_list.append(
            {
                "chr": chrom,
                "start": start,
                "end": start + width,
                "peak_id": f"enhancer_{peak_id}",
                "signal": np.random.lognormal(2, 1),
            }
        )
        peak_id += 1

    peaks = pd.DataFrame(peaks_list)

    # Ensure end > start
    peaks["end"] = peaks["start"] + np.abs(peaks["end"] - peaks["start"])

    return peaks


if __name__ == "__main__":
    main()
