"""
Transcription Factor ChIP-seq Analysis Page

Dedicated analysis workflow for TF ChIP-seq data:
- Peak characterization
- Motif enrichment
- Target gene identification
- Co-binding analysis
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from app.core.tf_analysis import (
    TFAnalysisEngine,
    TFCoBindingAnalyzer,
    MotifEnrichmentAnalyzer,
    generate_demo_tf_peaks,
    generate_demo_genes,
)

st.set_page_config(page_title="TF ChIP-seq - EpiNexus", page_icon="🔬", layout="wide")

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
    """Check if user has loaded TF ChIP-seq data."""
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
            <div style="font-size: 3rem; margin-bottom: 1rem;">🔬</div>
            <h2 style="color: #6c757d; margin-bottom: 0.5rem;">No Data Loaded</h2>
            <p style="color: #6c757d; font-size: 1.1rem;">
                Upload your TF ChIP-seq peaks to analyze binding patterns.
            </p>
        </div>
        """,
            unsafe_allow_html=True,
        )
        st.markdown("")
        if st.button("📁 Go to Data & Project", type="primary", use_container_width=True):
            st.switch_page("pages/01_data_project.py")
        st.markdown("")
        st.markdown("**TF ChIP-seq analysis:**")
        st.markdown("- Peak characterization")
        st.markdown("- Motif enrichment")
        st.markdown("- Target gene identification")
        st.markdown("- Co-binding analysis")


def main():
    st.title("🔬 Transcription Factor ChIP-seq Analysis")
    st.markdown("""
    Specialized analysis for transcription factor binding data. Upload your TF ChIP-seq
    peaks or use demo data to explore binding patterns, motifs, and target genes.
    """)

    # Check if data is loaded
    if not has_data():
        render_empty_state()
        return

    # Demo data indicator
    render_data_source_indicator()

    tab1, tab2, tab3, tab4, tab5 = st.tabs(
        ["📤 Data Input", "📊 Peak Analysis", "🧬 Motif Enrichment", "🎯 Target Genes", "🤝 Co-binding"]
    )

    with tab1:
        peaks, genes = render_data_input()

    with tab2:
        if peaks is not None:
            render_peak_analysis(peaks)
        else:
            st.info("Load TF peak data first.")

    with tab3:
        if peaks is not None:
            render_motif_analysis(peaks)
        else:
            st.info("Load TF peak data first.")

    with tab4:
        if peaks is not None and genes is not None:
            render_target_analysis(peaks, genes)
        else:
            st.info("Load both peak and gene data first.")

    with tab5:
        render_cobinding_analysis()


def render_data_source_indicator():
    """Show clear indicator of data source."""
    if st.session_state.get("tf_using_demo", True):
        st.info("📌 **Demo Mode**: Viewing simulated TF ChIP-seq data. Upload your own data for real analysis.")
    else:
        st.success("✅ **Your Data**: Analyzing your uploaded TF binding data.")


def render_data_input():
    """Handle TF data input."""
    st.header("Data Input")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("TF Binding Peaks")

        data_source = st.radio(
            "Data source", ["Demo data", "Upload peaks", "From preprocessing"], horizontal=True, key="tf_data_source"
        )

        if data_source == "Demo data":
            st.session_state.tf_using_demo = True

            tf_name = st.selectbox("Select Transcription Factor", ["MYC", "P53", "CTCF", "STAT3", "NFkB", "Custom"])

            if tf_name == "Custom":
                tf_name = st.text_input("Enter TF name", "MyTF")

            st.session_state.tf_name = tf_name
            peaks = generate_demo_tf_peaks(tf_name)
            st.success(f"Loaded {len(peaks):,} demo {tf_name} peaks")

        elif data_source == "Upload peaks":
            st.session_state.tf_using_demo = False

            tf_name = st.text_input("Transcription Factor name", "MyTF")
            st.session_state.tf_name = tf_name

            uploaded = st.file_uploader(
                "Upload peak file (BED/narrowPeak/CSV)", type=["bed", "narrowpeak", "csv", "tsv"]
            )

            if uploaded:
                if uploaded.name.endswith(".csv"):
                    peaks = pd.read_csv(uploaded)
                else:
                    peaks = pd.read_csv(uploaded, sep="\t", header=None)
                    # Assign narrowPeak column names if applicable
                    if len(peaks.columns) >= 10:
                        peaks.columns = [
                            "chr",
                            "start",
                            "end",
                            "peak_id",
                            "score",
                            "strand",
                            "signal",
                            "pvalue",
                            "qvalue",
                            "summit",
                        ][: len(peaks.columns)]
                    else:
                        peaks.columns = ["chr", "start", "end", "peak_id", "signal"][: len(peaks.columns)]

                st.success(f"Loaded {len(peaks):,} peaks")
            else:
                peaks = None

        else:  # From preprocessing
            st.session_state.tf_using_demo = False
            peaks = st.session_state.get("preprocessed_peaks", None)
            if peaks is not None:
                st.success(f"Using {len(peaks):,} peaks from preprocessing")
            else:
                st.warning("No peaks from preprocessing. Run the preprocessing pipeline first.")

        if peaks is not None:
            with st.expander("Preview peak data"):
                st.dataframe(peaks.head(20), use_container_width=True, hide_index=True)

    with col2:
        st.subheader("Gene Annotations")

        gene_source = st.radio("Gene source", ["Demo genes", "Upload genes"], horizontal=True, key="gene_source")

        if gene_source == "Demo genes":
            genes = generate_demo_genes()
            st.success(f"Loaded {len(genes):,} demo genes")
        else:
            gene_file = st.file_uploader("Upload gene annotation (BED/GTF/CSV)", type=["bed", "gtf", "csv", "tsv"])

            if gene_file:
                genes = pd.read_csv(gene_file, sep="\t" if not gene_file.name.endswith(".csv") else ",")
                st.success(f"Loaded {len(genes):,} genes")
            else:
                genes = generate_demo_genes()
                st.info("Using demo genes")

        if genes is not None:
            with st.expander("Preview gene data"):
                st.dataframe(genes.head(20), use_container_width=True, hide_index=True)

    return peaks, genes


def render_peak_analysis(peaks: pd.DataFrame):
    """Comprehensive TF peak analysis."""
    st.header("Peak Analysis")

    tf_name = st.session_state.get("tf_name", "TF")
    engine = TFAnalysisEngine(tf_name)

    genes = generate_demo_genes()
    results = engine.analyze_peaks(peaks, genes)

    # Summary metrics
    summary = results["peak_summary"]
    col1, col2, col3, col4, col5 = st.columns(5)

    col1.metric("Total Peaks", f"{summary['total_peaks']:,}")
    col2.metric("Mean Width", f"{summary['mean_width']:.0f} bp")
    col3.metric("Median Width", f"{summary['median_width']:.0f} bp")
    col4.metric("Chromosomes", summary["chromosomes"])
    col5.metric("Total Coverage", f"{summary['total_coverage_bp'] / 1e6:.1f} Mb")

    st.markdown("---")

    col1, col2 = st.columns(2)

    with col1:
        # Peak width distribution
        st.subheader("Peak Width Distribution")

        widths = peaks["end"] - peaks["start"]
        fig = px.histogram(x=widths, nbins=50, title=f"{tf_name} Peak Width Distribution")
        fig.update_layout(xaxis_title="Peak Width (bp)", yaxis_title="Count", height=400)
        fig.add_vline(
            x=widths.median(), line_dash="dash", line_color="red", annotation_text=f"Median: {widths.median():.0f} bp"
        )
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        # Signal distribution
        st.subheader("Peak Signal Distribution")

        if "signal" in peaks.columns:
            fig = px.histogram(peaks, x="signal", nbins=50, title=f"{tf_name} Peak Signal Distribution")
            fig.update_layout(xaxis_title="Signal (fold enrichment)", yaxis_title="Count", height=400)
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No signal column in peak data")

    # Genomic distribution
    st.subheader("Genomic Distribution")

    col1, col2 = st.columns([2, 1])

    with col1:
        dist = results["genomic_distribution"]
        categories = list(dist["counts"].keys())
        values = list(dist["counts"].values())

        fig = px.pie(
            values=values,
            names=categories,
            title=f"{tf_name} Binding Site Distribution",
            color_discrete_sequence=px.colors.qualitative.Set3,
        )
        fig.update_traces(textposition="inside", textinfo="percent+label")
        fig.update_layout(height=400)
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.markdown("**Distribution Summary:**")
        for cat, count in dist["counts"].items():
            pct = dist["percentages"][cat]
            st.markdown(f"- **{cat}**: {count:,} ({pct}%)")

    # Chromosome distribution
    st.subheader("Peaks per Chromosome")

    chr_counts = peaks.groupby("chr").size().reset_index(name="count")
    # Sort chromosomes naturally
    chr_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chr_counts["chr"] = pd.Categorical(chr_counts["chr"], categories=chr_order, ordered=True)
    chr_counts = chr_counts.sort_values("chr")

    fig = px.bar(chr_counts, x="chr", y="count", title="Peak Distribution Across Chromosomes")
    fig.update_layout(xaxis_title="Chromosome", yaxis_title="Number of Peaks", height=400)
    st.plotly_chart(fig, use_container_width=True)


def render_motif_analysis(peaks: pd.DataFrame):
    """Motif enrichment analysis for TF peaks."""
    st.header("Motif Enrichment Analysis")

    tf_name = st.session_state.get("tf_name", "TF")

    if st.session_state.get("tf_using_demo", True):
        st.info(
            "🔬 **Demo Mode**: Showing simulated motif enrichment results. "
            "With real data, this would run HOMER or MEME-ChIP for de novo motif discovery."
        )

    analyzer = MotifEnrichmentAnalyzer()
    results = analyzer.analyze_motif_enrichment(peaks, tf_name)

    # Primary motif
    st.subheader("Primary Motif")

    primary = results["primary_motif"]

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Consensus", primary["consensus"])
    col2.metric("Peaks with Motif", f"{primary['peaks_with_motif']:,}")
    col3.metric("Enrichment P-value", f"{primary['enrichment_pvalue']:.2e}")
    col4.metric("Centrality Score", f"{primary['centrality_score']:.2f}")

    st.markdown(f"""
    **Interpretation**: The {primary["name"]} ({primary["consensus"]}) was found in
    {primary["peaks_with_motif"] / len(peaks) * 100:.1f}% of {tf_name} binding sites,
    with strong central enrichment (centrality score: {primary["centrality_score"]:.2f}).
    """)

    st.markdown("---")

    # Secondary motifs (potential co-binders)
    st.subheader("Secondary Motifs (Potential Co-binders)")

    secondary_df = pd.DataFrame(results["secondary_motifs"])

    fig = px.bar(
        secondary_df,
        x="tf_name",
        y="peaks_with_motif",
        color=-np.log10(secondary_df["enrichment_pvalue"]),
        color_continuous_scale="Reds",
        title="Secondary Motif Enrichment",
    )
    fig.update_layout(
        xaxis_title="Transcription Factor",
        yaxis_title="Peaks with Motif",
        coloraxis_colorbar_title="-log10(P)",
        height=400,
    )
    st.plotly_chart(fig, use_container_width=True)

    with st.expander("View all secondary motifs"):
        display_df = secondary_df.copy()
        display_df["enrichment_pvalue"] = display_df["enrichment_pvalue"].apply(lambda x: f"{x:.2e}")
        st.dataframe(display_df, use_container_width=True, hide_index=True)

    # De novo motifs
    st.subheader("De Novo Motifs")

    st.markdown("Motifs discovered without prior knowledge:")

    denovo_df = pd.DataFrame(results["de_novo_motifs"])

    col1, col2, col3 = st.columns(3)

    for i, (_, motif) in enumerate(denovo_df.iterrows()):
        with [col1, col2, col3][i]:
            st.markdown(f"**{motif['motif_id']}**")
            st.code(motif["consensus"])
            st.markdown(f"Peaks: {motif['peaks_with_motif']:,}")
            st.markdown(f"P-value: {motif['enrichment_pvalue']:.2e}")
            st.markdown(f"Best match: {motif['best_match']}")


def render_target_analysis(peaks: pd.DataFrame, genes: pd.DataFrame):
    """Target gene identification."""
    st.header("Target Gene Identification")

    tf_name = st.session_state.get("tf_name", "TF")
    engine = TFAnalysisEngine(tf_name)

    # Settings
    col1, col2 = st.columns(2)

    with col1:
        max_distance = st.slider("Maximum distance from TSS (kb)", 10, 500, 100) * 1000

    with col2:
        promoter_def = st.slider("Promoter definition (bp from TSS)", 500, 5000, 2000)

    engine.promoter_distance = promoter_def

    # Identify targets
    targets = engine._identify_target_genes(peaks, genes, max_distance)

    if len(targets) == 0:
        st.warning("No target genes found. Try increasing the maximum distance.")
        return

    # Summary
    st.subheader("Target Gene Summary")

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Total Peak-Gene Links", f"{len(targets):,}")
    col2.metric("Unique Target Genes", f"{targets['gene_symbol'].nunique():,}")
    col3.metric("Promoter-Bound Genes", f"{targets[targets['is_promoter']]['gene_symbol'].nunique():,}")
    col4.metric("Distal Binding", f"{(~targets['is_promoter']).sum():,}")

    st.markdown("---")

    # Distance distribution
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Distance to TSS Distribution")

        fig = px.histogram(
            targets,
            x=targets["distance_to_tss"] / 1000,
            nbins=50,
            color="is_promoter",
            title="Distribution of Binding Sites Relative to TSS",
        )
        fig.update_layout(xaxis_title="Distance to TSS (kb)", yaxis_title="Count", height=400)
        fig.add_vline(x=0, line_dash="dash", line_color="gray")
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("Binding Location")

        location_counts = targets["binding_location"].value_counts()
        fig = px.pie(values=location_counts.values, names=location_counts.index, title="Upstream vs Downstream Binding")
        st.plotly_chart(fig, use_container_width=True)

    # Top target genes
    st.subheader("Top Target Genes")

    gene_summary = targets.groupby("gene_symbol").agg(
        {"peak_signal": ["count", "max", "sum"], "abs_distance": "min", "is_promoter": "any"}
    )
    gene_summary.columns = ["n_peaks", "max_signal", "total_signal", "min_distance", "has_promoter_binding"]
    gene_summary = gene_summary.sort_values("total_signal", ascending=False).reset_index()

    st.dataframe(gene_summary.head(50).round(2), use_container_width=True, hide_index=True)

    # Download
    csv = targets.to_csv(index=False)
    st.download_button("📥 Download Target Genes", csv, f"{tf_name}_target_genes.csv", "text/csv")


def render_cobinding_analysis():
    """Co-binding analysis between TFs."""
    st.header("Co-binding Analysis")

    st.markdown("""
    Analyze co-localization between two transcription factors. Upload data for both TFs
    or use demo data to explore co-binding patterns.
    """)

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("TF 1")

        tf1_source = st.radio("TF1 source", ["Demo", "Upload"], key="tf1_source", horizontal=True)

        if tf1_source == "Demo":
            tf1_name = st.selectbox("Select TF1", ["MYC", "P53", "CTCF"], key="tf1_name")
            tf1_peaks = generate_demo_tf_peaks(tf1_name, n_peaks=3000)
            st.success(f"Loaded {len(tf1_peaks):,} {tf1_name} peaks")
        else:
            tf1_name = st.text_input("TF1 name", "TF1")
            tf1_file = st.file_uploader("Upload TF1 peaks", type=["bed", "csv"], key="tf1_upload")
            if tf1_file:
                tf1_peaks = pd.read_csv(tf1_file, sep="\t" if tf1_file.name.endswith(".bed") else ",")
            else:
                tf1_peaks = None

    with col2:
        st.subheader("TF 2")

        tf2_source = st.radio("TF2 source", ["Demo", "Upload"], key="tf2_source", horizontal=True)

        if tf2_source == "Demo":
            tf2_name = st.selectbox("Select TF2", ["CTCF", "P53", "MYC"], key="tf2_name")
            tf2_peaks = generate_demo_tf_peaks(tf2_name, n_peaks=3000)
            st.success(f"Loaded {len(tf2_peaks):,} {tf2_name} peaks")
        else:
            tf2_name = st.text_input("TF2 name", "TF2")
            tf2_file = st.file_uploader("Upload TF2 peaks", type=["bed", "csv"], key="tf2_upload")
            if tf2_file:
                tf2_peaks = pd.read_csv(tf2_file, sep="\t" if tf2_file.name.endswith(".bed") else ",")
            else:
                tf2_peaks = None

    if tf1_peaks is None or tf2_peaks is None:
        st.info("Load peaks for both TFs to analyze co-binding.")
        return

    st.markdown("---")

    # Run analysis
    if st.button("🔍 Analyze Co-binding", type="primary"):
        with st.spinner("Analyzing co-localization..."):
            analyzer = TFCoBindingAnalyzer(overlap_threshold=100)
            results = analyzer.analyze_cobinding(tf1_peaks, tf2_peaks, tf1_name, tf2_name)

            st.session_state.cobinding_results = results
            st.session_state.tf1_name_cobind = tf1_name
            st.session_state.tf2_name_cobind = tf2_name

            # Record workflow step
            if HAS_WORKFLOW_MANAGER:
                WorkflowManager.record_step(
                    step_name="TF Co-binding Analysis",
                    tool="EpiNexus TF ChIP-seq Module",
                    parameters={
                        "tf1": tf1_name,
                        "tf2": tf2_name,
                        "overlap_threshold": 100,
                        "overlapping_peaks": results["overlapping_peaks"],
                        "fold_enrichment": results["fold_enrichment"],
                    },
                    inputs=["tf1_peaks", "tf2_peaks"],
                    outputs=["cobinding_results"],
                )

    # Display results
    results = st.session_state.get("cobinding_results", None)

    if results:
        tf1_name = st.session_state.get("tf1_name_cobind", "TF1")
        tf2_name = st.session_state.get("tf2_name_cobind", "TF2")

        st.subheader("Co-binding Results")

        col1, col2, col3, col4 = st.columns(4)

        col1.metric(f"{tf1_name} Peaks", f"{results[f'{tf1_name}_peaks']:,}")
        col2.metric(f"{tf2_name} Peaks", f"{results[f'{tf2_name}_peaks']:,}")
        col3.metric("Overlapping", f"{results['overlapping_peaks']:,}")
        col4.metric("Fold Enrichment", f"{results['fold_enrichment']:.1f}x")

        st.markdown("---")

        # Venn-style visualization
        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Overlap Statistics")

            # Create a simple bar chart showing overlap
            categories = [f"{tf1_name} only", "Overlapping", f"{tf2_name} only"]
            values = [results[f"{tf1_name}_only"], results["overlapping_peaks"], results[f"{tf2_name}_only"]]

            fig = px.bar(x=categories, y=values, color=categories, title="Peak Categories")
            fig.update_layout(xaxis_title="", yaxis_title="Number of Peaks", showlegend=False, height=400)
            st.plotly_chart(fig, use_container_width=True)

        with col2:
            st.subheader("Overlap Metrics")

            metrics_df = pd.DataFrame(
                {
                    "Metric": [
                        f"{tf1_name} overlap fraction",
                        f"{tf2_name} overlap fraction",
                        "Jaccard index",
                        "Fold enrichment",
                    ],
                    "Value": [
                        f"{results[f'{tf1_name}_overlap_fraction']:.1%}",
                        f"{results[f'{tf2_name}_overlap_fraction']:.1%}",
                        f"{results['jaccard_index']:.3f}",
                        f"{results['fold_enrichment']:.2f}",
                    ],
                }
            )

            st.dataframe(metrics_df, use_container_width=True, hide_index=True)

            # Interpretation
            jaccard = results["jaccard_index"]
            if jaccard > 0.5:
                st.success(f"✅ **Strong co-binding**: {tf1_name} and {tf2_name} frequently co-localize.")
            elif jaccard > 0.2:
                st.warning(f"⚠️ **Moderate co-binding**: Some overlap between {tf1_name} and {tf2_name}.")
            else:
                st.info(f"ℹ️ **Limited co-binding**: {tf1_name} and {tf2_name} bind largely distinct sites.")


if __name__ == "__main__":
    main()
