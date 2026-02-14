"""
Super-Enhancer Analysis Page

Identify and analyze super-enhancers using the ROSE algorithm.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import sys

from frontend.components.theme import COLORS

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from app.core.super_enhancers import SuperEnhancerDetector

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

st.set_page_config(page_title="Super-Enhancers - EpiNexus", page_icon="⭐", layout="wide")


def main():
    st.title("⭐ Super-Enhancer Analysis")
    st.markdown("""
    Identify super-enhancers using the ROSE algorithm. Super-enhancers are large clusters
    of enhancers that drive cell identity and are often associated with disease.
    """)

    # Sidebar settings
    with st.sidebar:
        st.header("ROSE Parameters")

        stitch_distance = st.slider(
            "Stitch distance (bp)", 5000, 25000, 12500, 500, help="Maximum distance to merge nearby enhancers"
        )

        tss_exclusion = st.slider(
            "TSS exclusion (bp)", 0, 5000, 2500, 250, help="Exclude peaks within this distance of TSS"
        )

        min_size = st.slider("Min enhancer size (bp)", 100, 2000, 500, 100, help="Minimum size for constituent peaks")

        st.markdown("---")
        st.markdown("""
        **Algorithm:**
        1. Stitch nearby H3K27ac peaks
        2. Rank by total signal
        3. Find inflection point
        4. Classify as SE or TE
        """)

    # Main content tabs
    tab1, tab2, tab3, tab4 = st.tabs(["📤 Input Data", "📊 Hockey Stick Plot", "📋 Results Table", "🧬 SE Annotation"])

    with tab1:
        result = render_input_and_run(stitch_distance, tss_exclusion, min_size)

    with tab2:
        if result:
            render_hockey_stick(result)
        else:
            st.info("Run analysis first to see the hockey stick plot.")

    with tab3:
        if result:
            render_results_table(result)
        else:
            st.info("Run analysis first to see results.")

    with tab4:
        if result:
            render_se_annotation(result)
        else:
            st.info("Run analysis first to see annotations.")


def render_input_and_run(stitch_distance, tss_exclusion, min_size):
    """Handle input data and run analysis."""
    st.header("Input Data")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("H3K27ac Peaks")
        st.markdown("Upload H3K27ac peak file (BED format with signal column)")

        uploaded_file = st.file_uploader("Upload peak file", type=["bed", "narrowPeak", "csv", "tsv"], key="peak_file")

    with col2:
        st.subheader("TSS Regions (Optional)")
        st.markdown("Exclude promoter regions from analysis")

        tss_file = st.file_uploader("Upload TSS file", type=["bed", "csv", "tsv"], key="tss_file")

    # Load data
    if uploaded_file:
        peaks = load_peak_file(uploaded_file)
        st.success(f"Loaded {len(peaks)} peaks")
    else:
        peaks = None
        # Show empty state hint
        st.markdown(
            """
        <div style="text-align: center; padding: 2rem; background: #f8f9fa;
                    border-radius: 8px; border: 1px dashed #dee2e6; margin-top: 1rem;">
            <p style="color: #6c757d; margin: 0;">
                Upload your H3K27ac peaks to run super-enhancer detection.<br/>
                <small>Accepted formats: BED, narrowPeak, CSV, TSV</small>
            </p>
        </div>
        """,
            unsafe_allow_html=True,
        )

    # Preview data
    if peaks is not None:
        st.subheader("Peak Preview")
        st.dataframe(peaks.head(10), use_container_width=True)

        # Run analysis
        st.markdown("---")

        if st.button("🚀 Run Super-Enhancer Detection", type="primary"):
            with st.spinner("Running ROSE algorithm..."):
                detector = SuperEnhancerDetector(
                    stitch_distance=stitch_distance, tss_exclusion=tss_exclusion, min_enhancer_size=min_size
                )

                tss_df = None
                if tss_file:
                    tss_df = load_peak_file(tss_file)

                result = detector.detect(peaks, signal_col="signal", tss_regions=tss_df)

                # Store in session state
                st.session_state.se_result = result

                # Show summary
                st.success("Analysis complete!")

                # Record workflow step
                if HAS_WORKFLOW_MANAGER:
                    WorkflowManager.record_step(
                        step_name="Super-Enhancer Detection",
                        tool="EpiNexus ROSE Algorithm",
                        parameters={
                            "stitch_distance": stitch_distance,
                            "tss_exclusion": tss_exclusion,
                            "min_enhancer_size": min_size,
                        },
                        inputs=["peaks"],
                        outputs=["super_enhancers", "typical_enhancers"],
                    )

                col1, col2, col3, col4 = st.columns(4)
                col1.metric("Total Enhancers", result.statistics["total_enhancers"])
                col2.metric("Super-Enhancers", result.statistics["super_enhancer_count"])
                col3.metric("Typical Enhancers", result.statistics["typical_enhancer_count"])
                col4.metric("SE Signal Fraction", f"{result.statistics['se_signal_fraction']:.1%}")

                return result

    return st.session_state.get("se_result", None)


def render_hockey_stick(result):
    """Render the hockey stick plot."""
    st.header("Hockey Stick Plot")

    st.markdown("""
    The hockey stick plot shows all enhancers ranked by signal intensity.
    Super-enhancers are identified at the inflection point where signal
    increases dramatically.
    """)

    # Prepare data
    all_enh = result.all_enhancers.copy()

    # Create hockey stick plot
    fig = go.Figure()

    # Typical enhancers
    typical = all_enh[~all_enh["is_super_enhancer"]]
    fig.add_trace(
        go.Scatter(
            x=typical["rank"],
            y=typical["signal"],
            mode="markers",
            marker=dict(color="#3498db", size=6, opacity=0.6),
            name="Typical Enhancers",
            hovertemplate="Rank: %{x}<br>Signal: %{y:.1f}<extra></extra>",
        )
    )

    # Super-enhancers
    super_e = all_enh[all_enh["is_super_enhancer"]]
    fig.add_trace(
        go.Scatter(
            x=super_e["rank"],
            y=super_e["signal"],
            mode="markers",
            marker=dict(color="#e74c3c", size=8, opacity=0.8),
            name="Super-Enhancers",
            hovertemplate="Rank: %{x}<br>Signal: %{y:.1f}<extra></extra>",
        )
    )

    # Add cutoff line
    fig.add_hline(
        y=result.cutoff_value, line_dash="dash", line_color="gray", annotation_text=f"Cutoff: {result.cutoff_value:.1f}"
    )

    fig.update_layout(
        xaxis_title="Enhancer Rank",
        yaxis_title="H3K27ac Signal",
        height=500,
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
    )

    st.plotly_chart(fig, use_container_width=True)

    # Statistics
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Size Distribution")

        size_data = pd.DataFrame(
            {
                "Class": ["Super-Enhancer"] * len(result.super_enhancers)
                + ["Typical Enhancer"] * len(result.typical_enhancers),
                "Size (kb)": list(result.super_enhancers["size"] / 1000)
                + list(result.typical_enhancers["size"] / 1000),
            }
        )

        fig2 = px.box(size_data, x="Class", y="Size (kb)", color="Class", color_discrete_map=COLORS.SE_MAP)
        fig2.update_layout(height=350, showlegend=False)
        st.plotly_chart(fig2, use_container_width=True)

    with col2:
        st.subheader("Constituent Peaks")

        const_data = pd.DataFrame(
            {
                "Class": ["Super-Enhancer"] * len(result.super_enhancers)
                + ["Typical Enhancer"] * len(result.typical_enhancers),
                "Constituent Peaks": list(result.super_enhancers["constituent_count"])
                + list(result.typical_enhancers["constituent_count"]),
            }
        )

        fig3 = px.box(const_data, x="Class", y="Constituent Peaks", color="Class", color_discrete_map=COLORS.SE_MAP)
        fig3.update_layout(height=350, showlegend=False)
        st.plotly_chart(fig3, use_container_width=True)


def render_results_table(result):
    """Render results table."""
    st.header("Super-Enhancer Results")

    # Filter options
    col1, col2, col3 = st.columns(3)

    with col1:
        show_type = st.selectbox("Show", ["All", "Super-Enhancers only", "Typical Enhancers only"])

    with col2:
        min_signal = st.number_input("Min signal", value=0.0)

    with col3:
        sort_by = st.selectbox("Sort by", ["signal", "size", "rank", "constituent_count"])

    # Filter data
    if show_type == "Super-Enhancers only":
        data = result.super_enhancers.copy()
    elif show_type == "Typical Enhancers only":
        data = result.typical_enhancers.copy()
    else:
        data = result.all_enhancers.copy()

    data = data[data["signal"] >= min_signal]
    data = data.sort_values(sort_by, ascending=False)

    # Display columns
    display_cols = ["chr", "start", "end", "size", "signal", "rank", "constituent_count", "enhancer_class"]
    available_cols = [c for c in display_cols if c in data.columns]

    st.dataframe(data[available_cols], use_container_width=True, hide_index=True)

    # Download
    csv = data.to_csv(index=False)
    st.download_button("📥 Download Results", csv, "super_enhancers.csv", "text/csv")


def render_se_annotation(result):
    """Annotate super-enhancers with nearby genes."""
    st.header("Super-Enhancer Annotation")

    st.markdown("""
    Super-enhancers often regulate cell identity genes. Below are the predicted
    target genes based on proximity.
    """)

    # Demo annotation (in real implementation, would use ChIPseeker or similar)
    se_data = result.super_enhancers.copy()

    if len(se_data) == 0:
        st.warning("No super-enhancers detected.")
        return

    # Add demo gene annotations
    np.random.seed(42)
    demo_genes = [
        "MYC",
        "SOX2",
        "NANOG",
        "OCT4",
        "KLF4",
        "GATA1",
        "RUNX1",
        "ERG",
        "FLI1",
        "TAL1",
        "LMO2",
        "MECOM",
        "HOXA9",
        "MEIS1",
    ]

    se_data["nearest_gene"] = np.random.choice(demo_genes, len(se_data))
    se_data["distance_to_tss"] = np.random.randint(-50000, 50000, len(se_data))

    # Display
    st.dataframe(
        se_data[["chr", "start", "end", "signal", "nearest_gene", "distance_to_tss"]].head(20),
        use_container_width=True,
        hide_index=True,
    )

    # Gene list
    st.subheader("SE-Associated Genes")

    genes = se_data["nearest_gene"].value_counts()
    st.markdown(f"**{len(genes)} unique genes** associated with super-enhancers")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Top SE-associated genes:**")
        for gene, count in genes.head(10).items():
            st.write(f"• {gene}: {count} SEs")

    with col2:
        # Pathway enrichment placeholder
        st.markdown("**Enriched pathways (GO):**")
        pathways = [
            "Cell fate commitment",
            "Stem cell differentiation",
            "Transcription factor activity",
            "Chromatin organization",
        ]
        for p in pathways:
            st.write(f"• {p}")


def generate_demo_peaks():
    """Generate demo peak data."""
    np.random.seed(42)
    n_peaks = 2000

    # Generate peaks with realistic distribution
    # Most peaks have low signal, few have very high (super-enhancers)
    signals = np.concatenate(
        [
            np.random.exponential(10, int(n_peaks * 0.9)),  # Typical enhancers
            np.random.exponential(50, int(n_peaks * 0.1)),  # Strong enhancers
        ]
    )

    chroms = np.random.choice([f"chr{i}" for i in range(1, 23)], n_peaks)
    starts = np.random.randint(1000000, 200000000, n_peaks)

    peaks = pd.DataFrame(
        {"chr": chroms, "start": starts, "end": starts + np.random.randint(500, 3000, n_peaks), "signal": signals}
    )

    return peaks.sort_values(["chr", "start"]).reset_index(drop=True)


def load_peak_file(uploaded_file):
    """Load peak file from upload."""
    name = uploaded_file.name.lower()

    if name.endswith(".csv"):
        df = pd.read_csv(uploaded_file)
    elif name.endswith(".tsv") or name.endswith(".bed") or name.endswith(".narrowpeak"):
        df = pd.read_csv(uploaded_file, sep="\t", header=None)

        # Standard BED columns
        if len(df.columns) >= 6:
            df.columns = ["chr", "start", "end", "name", "score", "strand"] + [
                f"col{i}" for i in range(6, len(df.columns))
            ]
            df["signal"] = df["score"]
        else:
            df.columns = ["chr", "start", "end"] + [f"col{i}" for i in range(3, len(df.columns))]
            df["signal"] = 1.0

    return df


if __name__ == "__main__":
    main()
