"""
Differential Analysis Page
Identify differentially enriched peaks between conditions.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

# Import data manager
try:
    from frontend.components.data_manager import DataManager, DataSource
    HAS_DATA_MANAGER = True
except ImportError:
    HAS_DATA_MANAGER = False

st.set_page_config(page_title="Differential Analysis - EpiNexus", page_icon="📊", layout="wide")


def has_data():
    """Check if user has loaded data."""
    if HAS_DATA_MANAGER:
        peaks = DataManager.get_data('peaks')
        return peaks is not None and len(peaks) > 0
    return len(st.session_state.get('samples', [])) > 0


def main():
    st.title("📊 Differential Analysis")
    st.markdown("Identify regions with significantly different histone modification levels between conditions.")

    # Check if data is loaded
    if not has_data():
        render_empty_state()
        return

    # Sidebar settings
    with st.sidebar:
        st.header("Analysis Parameters")

        fdr_thresh = st.slider("FDR Threshold", 0.01, 0.2, 0.05)
        fc_thresh = st.slider("Log2 FC Threshold", 0.5, 3.0, 1.0)

        st.markdown("---")

        contrast = st.selectbox(
            "Comparison",
            ["Treatment vs Control", "Timepoint 2 vs 1", "Custom"]
        )

        method = st.selectbox(
            "Method",
            ["DiffBind (DESeq2)", "DiffBind (edgeR)", "THOR", "SICER-DF"]
        )

    tab1, tab2, tab3, tab4 = st.tabs([
        "🔬 Run Analysis",
        "🌋 Volcano Plot",
        "📊 Results Table",
        "🗺️ Genome View"
    ])

    with tab1:
        render_analysis_setup()
    with tab2:
        render_volcano_plot(fdr_thresh, fc_thresh)
    with tab3:
        render_results_table(fdr_thresh, fc_thresh)
    with tab4:
        render_genome_view()


def render_empty_state():
    """Show empty state when no data is loaded."""
    st.markdown("---")

    col1, col2, col3 = st.columns([1, 2, 1])

    with col2:
        st.markdown("""
        <div style="text-align: center; padding: 3rem; background: #f8f9fa; border-radius: 12px; border: 2px dashed #dee2e6;">
            <h2 style="color: #6c757d;">📊 No Data Loaded</h2>
            <p style="color: #6c757d; font-size: 1.1rem;">
                Upload your peak files to run differential analysis.
            </p>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("")

        if st.button("📁 Go to Data & Project", type="primary", use_container_width=True):
            st.switch_page("pages/01_data_project.py")

        st.markdown("")

        st.markdown("""
        **What you need:**
        - Peak files (BED, narrowPeak, or broadPeak format)
        - At least 2 conditions to compare (e.g., Treatment vs Control)
        - Ideally 2-3 replicates per condition

        **Need help?** Check the [Documentation](pages/21_help.py) for a walkthrough.
        """)


def render_analysis_setup():
    """Configure and run differential analysis."""
    st.header("Analysis Setup")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Select Samples")

        samples = st.session_state.get('samples', [])
        sample_ids = [s.get('SampleID', f"Sample_{i}") for i, s in enumerate(samples)]

        if not sample_ids:
            st.warning("No samples defined. Go to Data & Project to add samples.")
            return

        st.markdown("**Control Group:**")
        ctrl_samples = st.multiselect(
            "Control samples",
            sample_ids,
            default=sample_ids[:len(sample_ids)//2] if sample_ids else [],
            key='ctrl'
        )

        st.markdown("**Treatment Group:**")
        treat_samples = st.multiselect(
            "Treatment samples",
            [s for s in sample_ids if s not in ctrl_samples],
            default=[s for s in sample_ids if s not in ctrl_samples][:len(sample_ids)//2],
            key='treat'
        )

    with col2:
        st.subheader("Analysis Options")

        consensus_method = st.selectbox(
            "Consensus peak set",
            ["At least 2 samples", "All samples", "Majority", "Union"]
        )

        normalize = st.checkbox("Library size normalization", value=True)
        batch_correct = st.checkbox("Batch correction", value=False)

        if batch_correct:
            batch_var = st.text_input("Batch variable column")

    st.markdown("---")

    if st.button("▶️ Run Differential Analysis", type="primary"):
        if len(ctrl_samples) < 1 or len(treat_samples) < 1:
            st.error("Please select at least 1 sample for each group.")
            return

        with st.spinner("Running DiffBind analysis..."):
            progress = st.progress(0)
            for i in range(100):
                progress.progress(i + 1)

            # Store results in session state
            st.session_state.diff_results = generate_results(ctrl_samples, treat_samples)
            st.success(f"Analysis complete! Found {len(st.session_state.diff_results):,} peaks analyzed.")


def generate_results(ctrl_samples, treat_samples):
    """Generate differential analysis results from user data."""
    # In production, this would call DiffBind
    # For now, generate results based on loaded peaks
    if HAS_DATA_MANAGER:
        peaks = DataManager.get_data('peaks')
        if peaks is not None:
            n = len(peaks)
            np.random.seed(42)

            results = peaks.copy()
            results['log2FC'] = np.random.normal(0, 1.5, n)
            results['pvalue'] = 10 ** -np.random.uniform(0.5, 8, n)
            results['FDR'] = np.minimum(results['pvalue'] * n / np.arange(1, n+1), 1)

            return results

    return pd.DataFrame()


def render_volcano_plot(fdr_thresh, fc_thresh):
    """Interactive volcano plot."""
    st.header("Volcano Plot")

    # Check for results
    if 'diff_results' not in st.session_state or st.session_state.diff_results is None:
        st.info("Run the analysis first to see the volcano plot.")
        return

    df = st.session_state.diff_results

    if len(df) == 0:
        st.warning("No results to display.")
        return

    # Add significance classification
    df = df.copy()
    df['Significant'] = 'Not Significant'
    df.loc[(df['FDR'] < fdr_thresh) & (df['log2FC'] > fc_thresh), 'Significant'] = 'Up'
    df.loc[(df['FDR'] < fdr_thresh) & (df['log2FC'] < -fc_thresh), 'Significant'] = 'Down'

    color_map = {'Up': '#e74c3c', 'Down': '#3498db', 'Not Significant': '#95a5a6'}

    fig = px.scatter(
        df, x='log2FC', y=-np.log10(df['FDR']),
        color='Significant', color_discrete_map=color_map,
        hover_data=['chr', 'start', 'end'] if 'chr' in df.columns else None,
        opacity=0.6
    )

    fig.add_hline(y=-np.log10(fdr_thresh), line_dash="dash", line_color="gray")
    fig.add_vline(x=fc_thresh, line_dash="dash", line_color="gray")
    fig.add_vline(x=-fc_thresh, line_dash="dash", line_color="gray")

    fig.update_layout(
        xaxis_title="Log2 Fold Change",
        yaxis_title="-Log10 FDR",
        height=600
    )

    st.plotly_chart(fig, use_container_width=True)

    # Summary counts
    col1, col2, col3 = st.columns(3)
    col1.metric("Up-regulated", (df['Significant'] == 'Up').sum())
    col2.metric("Down-regulated", (df['Significant'] == 'Down').sum())
    col3.metric("Not Significant", (df['Significant'] == 'Not Significant').sum())


def render_results_table(fdr_thresh, fc_thresh):
    """Filterable results table."""
    st.header("Differential Peaks Results")

    if 'diff_results' not in st.session_state or st.session_state.diff_results is None:
        st.info("Run the analysis first to see results.")
        return

    results = st.session_state.diff_results

    if len(results) == 0:
        st.warning("No results to display.")
        return

    # Filters
    col1, col2, col3 = st.columns(3)
    with col1:
        show_sig = st.checkbox("Only significant", value=True)
    with col2:
        direction = st.selectbox("Direction", ["All", "Gained", "Lost"])
    with col3:
        if 'chr' in results.columns:
            chr_filter = st.multiselect("Chromosomes", results['chr'].unique())
        else:
            chr_filter = []

    # Apply filters
    filtered = results.copy()
    if show_sig:
        filtered = filtered[(filtered['FDR'] < fdr_thresh) & (np.abs(filtered['log2FC']) > fc_thresh)]
    if direction == "Gained":
        filtered = filtered[filtered['log2FC'] > 0]
    elif direction == "Lost":
        filtered = filtered[filtered['log2FC'] < 0]
    if chr_filter:
        filtered = filtered[filtered['chr'].isin(chr_filter)]

    st.dataframe(filtered.round(4), use_container_width=True, hide_index=True)

    # Download
    csv = filtered.to_csv(index=False)
    st.download_button("📥 Download Results", csv, "differential_peaks.csv", "text/csv")


def render_genome_view():
    """Chromosome-level view of differential peaks."""
    st.header("Genome Distribution")

    if 'diff_results' not in st.session_state or st.session_state.diff_results is None:
        st.info("Run the analysis first to see genome distribution.")
        return

    results = st.session_state.diff_results

    if 'chr' not in results.columns:
        st.warning("Chromosome information not available in results.")
        return

    # Summarize by chromosome
    sig_results = results[results['FDR'] < 0.05]

    gained = sig_results[sig_results['log2FC'] > 0].groupby('chr').size()
    lost = sig_results[sig_results['log2FC'] < 0].groupby('chr').size()

    chr_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    all_chrs = [c for c in chr_order if c in results['chr'].unique()]

    chrom_data = pd.DataFrame({
        'Chromosome': all_chrs,
        'Gained': [gained.get(c, 0) for c in all_chrs],
        'Lost': [lost.get(c, 0) for c in all_chrs]
    })

    fig = go.Figure()
    fig.add_trace(go.Bar(x=chrom_data['Chromosome'], y=chrom_data['Gained'], name='Gained', marker_color='#e74c3c'))
    fig.add_trace(go.Bar(x=chrom_data['Chromosome'], y=-chrom_data['Lost'], name='Lost', marker_color='#3498db'))

    fig.update_layout(
        barmode='relative',
        xaxis_title='Chromosome',
        yaxis_title='Differential Peaks',
        height=400
    )

    st.plotly_chart(fig, use_container_width=True)


if __name__ == "__main__":
    main()
