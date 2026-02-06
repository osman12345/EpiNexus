"""
ATAC-seq Analysis Page

Chromatin accessibility analysis and integration with histone modifications.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from app.core.atacseq import ATACSeqAnalyzer, generate_demo_atac_data

st.set_page_config(
    page_title="ATAC-seq Analysis - EpiNexus",
    page_icon="🔓",
    layout="wide"
)

# Try to import data manager
try:
    from frontend.components.data_manager import DataManager
    HAS_DATA_MANAGER = True
except ImportError:
    HAS_DATA_MANAGER = False


def has_data():
    """Check if user has loaded ATAC-seq data."""
    if HAS_DATA_MANAGER:
        peaks = DataManager.get_data('peaks')
        return peaks is not None and len(peaks) > 0
    return len(st.session_state.get('samples', [])) > 0


def render_empty_state():
    """Show empty state when no data is loaded."""
    st.markdown("---")
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        st.markdown("""
        <div style="text-align: center; padding: 3rem; background: #f8f9fa;
                    border-radius: 12px; border: 2px dashed #dee2e6;">
            <div style="font-size: 3rem; margin-bottom: 1rem;">🔓</div>
            <h2 style="color: #6c757d; margin-bottom: 0.5rem;">No Data Loaded</h2>
            <p style="color: #6c757d; font-size: 1.1rem;">
                Upload your ATAC-seq data to analyze chromatin accessibility.
            </p>
        </div>
        """, unsafe_allow_html=True)
        st.markdown("")
        if st.button("📁 Go to Data & Project", type="primary", use_container_width=True):
            st.switch_page("pages/01_data_project.py")
        st.markdown("")
        st.markdown("**ATAC-seq analysis features:**")
        st.markdown("- QC metrics & fragment analysis")
        st.markdown("- Differential accessibility")
        st.markdown("- TF footprinting")
        st.markdown("- Histone mark integration")


def main():
    st.title("🔓 ATAC-seq Analysis")
    st.markdown("""
    Analyze chromatin accessibility data and integrate with histone modifications
    to characterize regulatory elements.
    """)

    # Check if data is loaded
    if not has_data():
        render_empty_state()
        return

    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "📊 QC Metrics",
        "📈 Fragment Analysis",
        "🔬 Differential Accessibility",
        "👣 TF Footprinting",
        "🔗 Histone Integration"
    ])

    with tab1:
        render_qc_metrics()

    with tab2:
        render_fragment_analysis()

    with tab3:
        render_differential()

    with tab4:
        render_footprinting()

    with tab5:
        render_histone_integration()


def render_qc_metrics():
    """Display ATAC-seq QC metrics."""
    st.header("Quality Control Metrics")

    # Load demo data
    peaks, bam_stats = generate_demo_atac_data()
    analyzer = ATACSeqAnalyzer()

    # Calculate TSS enrichment (simulated)
    np.random.seed(42)
    tss_scores = np.random.lognormal(1.5, 0.5, 1000)

    qc = analyzer.calculate_qc_metrics(bam_stats, peaks, tss_scores)

    # Display metrics in cards
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric("Total Reads", f"{qc.total_reads/1e6:.1f}M")
        st.metric("Mapping Rate", f"{qc.mapping_rate:.1%}")

    with col2:
        st.metric("Mitochondrial", f"{qc.mitochondrial_rate:.1%}",
                 delta="Good" if qc.mitochondrial_rate < 0.1 else "High",
                 delta_color="normal" if qc.mitochondrial_rate < 0.1 else "inverse")
        st.metric("Duplicate Rate", f"{qc.duplicate_rate:.1%}")

    with col3:
        st.metric("TSS Enrichment", f"{qc.tss_enrichment:.1f}",
                 delta="Good" if qc.tss_enrichment > 5 else "Low")
        st.metric("FRiP", f"{qc.frip:.1%}")

    with col4:
        st.metric("Peak Count", f"{qc.peak_count:,}")
        st.metric("NFR Ratio", f"{qc.nfr_ratio:.1%}")

    st.markdown("---")

    # QC thresholds guide
    st.subheader("QC Thresholds Guide")

    thresholds = pd.DataFrame({
        'Metric': ['TSS Enrichment', 'FRiP', 'Mitochondrial %', 'NFR Ratio'],
        'Ideal': ['>7', '>0.3', '<5%', '>0.4'],
        'Acceptable': ['5-7', '0.2-0.3', '5-10%', '0.3-0.4'],
        'Poor': ['<5', '<0.2', '>10%', '<0.3']
    })

    st.dataframe(thresholds, use_container_width=True, hide_index=True)


def render_fragment_analysis():
    """Analyze fragment size distribution."""
    st.header("Fragment Size Distribution")

    st.markdown("""
    ATAC-seq fragments show characteristic sizes corresponding to:
    - **Nucleosome-free regions (NFR)**: <100 bp
    - **Mononucleosome**: 180-247 bp
    - **Dinucleosome**: 315-473 bp
    """)

    # Generate demo fragment sizes
    np.random.seed(42)
    fragment_sizes = np.concatenate([
        np.random.normal(50, 20, 5000),    # NFR
        np.random.normal(200, 25, 3000),   # Mono
        np.random.normal(400, 40, 1500),   # Di
        np.random.normal(600, 50, 500)     # Tri
    ])
    fragment_sizes = fragment_sizes[fragment_sizes > 0]

    analyzer = ATACSeqAnalyzer()
    frag_dist = analyzer.analyze_fragment_sizes(fragment_sizes.astype(int))

    col1, col2 = st.columns(2)

    with col1:
        # Fragment size histogram
        fig = go.Figure()

        fig.add_trace(go.Histogram(
            x=fragment_sizes,
            nbinsx=100,
            name='Fragments',
            marker_color='#3498db'
        ))

        # Add region annotations
        fig.add_vrect(x0=0, x1=100, fillcolor="#27ae60", opacity=0.2,
                     annotation_text="NFR", annotation_position="top")
        fig.add_vrect(x0=180, x1=247, fillcolor="#f39c12", opacity=0.2,
                     annotation_text="Mono")
        fig.add_vrect(x0=315, x1=473, fillcolor="#9b59b6", opacity=0.2,
                     annotation_text="Di")

        fig.update_layout(
            title="Fragment Size Distribution",
            xaxis_title="Fragment Size (bp)",
            yaxis_title="Count",
            height=400
        )

        st.plotly_chart(fig, use_container_width=True)

    with col2:
        # Pie chart of fractions
        fractions = pd.DataFrame({
            'Region': ['NFR (<100bp)', 'Mono (180-247bp)', 'Di (315-473bp)', 'Other'],
            'Fraction': [
                frag_dist.nfr_fraction,
                frag_dist.mono_fraction,
                frag_dist.di_fraction,
                1 - frag_dist.nfr_fraction - frag_dist.mono_fraction - frag_dist.di_fraction
            ]
        })

        fig = px.pie(fractions, values='Fraction', names='Region',
                    color_discrete_sequence=['#27ae60', '#f39c12', '#9b59b6', '#95a5a6'])
        fig.update_layout(title="Fragment Fraction by Region", height=400)

        st.plotly_chart(fig, use_container_width=True)


def render_differential():
    """Differential accessibility analysis."""
    st.header("Differential Accessibility")

    col1, col2 = st.columns([1, 3])

    with col1:
        st.subheader("Settings")

        fdr_thresh = st.slider("FDR Threshold", 0.01, 0.2, 0.05)
        fc_thresh = st.slider("Log2 FC Threshold", 0.5, 3.0, 1.0)

        method = st.selectbox(
            "Method",
            ["DESeq2", "edgeR", "Wilcoxon"]
        )

        if st.button("Run Analysis", type="primary"):
            st.session_state.atac_diff_run = True

    with col2:
        # Generate demo differential results
        np.random.seed(42)
        n = 2000

        diff_results = pd.DataFrame({
            'peak_id': [f'peak_{i}' for i in range(n)],
            'log2FoldChange': np.random.normal(0, 1.5, n),
            'pvalue': 10 ** -np.random.uniform(0, 6, n),
            'padj': 10 ** -np.random.uniform(0, 4, n),
            'baseMean': np.random.lognormal(5, 1.5, n)
        })

        # Volcano plot
        diff_results['Significant'] = (
            (diff_results['padj'] < fdr_thresh) &
            (np.abs(diff_results['log2FoldChange']) > fc_thresh)
        )

        diff_results['Direction'] = 'Not Significant'
        diff_results.loc[diff_results['Significant'] & (diff_results['log2FoldChange'] > 0), 'Direction'] = 'More Accessible'
        diff_results.loc[diff_results['Significant'] & (diff_results['log2FoldChange'] < 0), 'Direction'] = 'Less Accessible'

        color_map = {
            'More Accessible': '#27ae60',
            'Less Accessible': '#e74c3c',
            'Not Significant': '#95a5a6'
        }

        fig = px.scatter(
            diff_results,
            x='log2FoldChange',
            y=-np.log10(diff_results['padj']),
            color='Direction',
            color_discrete_map=color_map,
            opacity=0.6,
            hover_data=['peak_id', 'baseMean']
        )

        fig.add_hline(y=-np.log10(fdr_thresh), line_dash='dash', line_color='gray')
        fig.add_vline(x=fc_thresh, line_dash='dash', line_color='gray')
        fig.add_vline(x=-fc_thresh, line_dash='dash', line_color='gray')

        fig.update_layout(
            title="Differential Accessibility Volcano Plot",
            xaxis_title="Log2 Fold Change",
            yaxis_title="-Log10 Adjusted P-value",
            height=500
        )

        st.plotly_chart(fig, use_container_width=True)

        # Summary
        col1, col2, col3 = st.columns(3)
        col1.metric("More Accessible", (diff_results['Direction'] == 'More Accessible').sum())
        col2.metric("Less Accessible", (diff_results['Direction'] == 'Less Accessible').sum())
        col3.metric("Not Significant", (diff_results['Direction'] == 'Not Significant').sum())


def render_footprinting():
    """TF footprinting analysis."""
    st.header("TF Footprinting Analysis")

    st.markdown("""
    Transcription factor footprints appear as local dips in ATAC-seq signal
    where bound TFs protect DNA from transposase cutting.
    """)

    col1, col2 = st.columns([1, 2])

    with col1:
        tf_select = st.selectbox(
            "Select Transcription Factor",
            ["CTCF", "SP1", "NRF1", "REST", "GATA1", "MAX"]
        )

        window_size = st.slider("Window Size (bp)", 50, 200, 100)

    with col2:
        # Generate demo footprint
        np.random.seed(42)
        x = np.arange(-window_size, window_size)

        # Background signal with dip at center (footprint)
        footprint_depths = {'CTCF': 0.4, 'SP1': 0.25, 'NRF1': 0.3, 'REST': 0.35, 'GATA1': 0.2, 'MAX': 0.15}
        depth = footprint_depths.get(tf_select, 0.3)

        # Flank signal
        signal = np.ones_like(x, dtype=float)
        # Add footprint dip
        footprint_width = 15
        footprint = depth * np.exp(-0.5 * (x / footprint_width) ** 2)
        signal = signal - footprint
        # Add noise
        signal += np.random.normal(0, 0.05, len(x))

        # Add flanking peaks (characteristic of ATAC)
        flank_peaks = 0.3 * np.exp(-0.5 * ((x - 30) / 10) ** 2) + 0.3 * np.exp(-0.5 * ((x + 30) / 10) ** 2)
        signal += flank_peaks

        fig = go.Figure()

        fig.add_trace(go.Scatter(
            x=x,
            y=signal,
            mode='lines',
            fill='tozeroy',
            name=f'{tf_select} Footprint',
            line=dict(color='#3498db', width=2)
        ))

        # Add motif annotation
        fig.add_vrect(x0=-10, x1=10, fillcolor="rgba(231, 76, 60, 0.3)",
                     annotation_text=f"{tf_select} motif")

        fig.update_layout(
            title=f"{tf_select} Aggregate Footprint",
            xaxis_title="Distance from Motif Center (bp)",
            yaxis_title="Normalized ATAC Signal",
            height=400
        )

        st.plotly_chart(fig, use_container_width=True)

        # Footprint statistics
        st.markdown(f"""
        **Footprint Statistics for {tf_select}:**
        - Footprint depth: {depth:.1%}
        - Motifs analyzed: 5,234
        - Footprint score: {depth * 100:.1f}
        """)


def render_histone_integration():
    """Integrate ATAC-seq with histone modifications."""
    st.header("Integration with Histone Modifications")

    st.markdown("""
    Classify accessible chromatin regions by their histone modification profiles
    to identify promoters, enhancers, and other regulatory elements.
    """)

    # Generate demo integration data
    np.random.seed(42)
    n_peaks = 1000

    integration_data = pd.DataFrame({
        'peak_id': [f'peak_{i}' for i in range(n_peaks)],
        'chr': np.random.choice([f'chr{i}' for i in range(1, 23)], n_peaks),
        'H3K4me3': np.random.choice([True, False], n_peaks, p=[0.25, 0.75]),
        'H3K27ac': np.random.choice([True, False], n_peaks, p=[0.4, 0.6]),
        'H3K4me1': np.random.choice([True, False], n_peaks, p=[0.35, 0.65]),
        'H3K27me3': np.random.choice([True, False], n_peaks, p=[0.15, 0.85])
    })

    # Classify regions
    def classify(row):
        if row['H3K4me3'] and row['H3K27ac']:
            return 'Active Promoter'
        elif row['H3K4me1'] and row['H3K27ac']:
            return 'Active Enhancer'
        elif row['H3K4me1'] and not row['H3K27ac']:
            return 'Poised Enhancer'
        elif row['H3K4me3'] and row['H3K27me3']:
            return 'Bivalent'
        elif row['H3K27me3']:
            return 'Repressed'
        else:
            return 'Other Accessible'

    integration_data['Chromatin_State'] = integration_data.apply(classify, axis=1)

    col1, col2 = st.columns(2)

    with col1:
        # Pie chart of chromatin states
        state_counts = integration_data['Chromatin_State'].value_counts()

        fig = px.pie(
            values=state_counts.values,
            names=state_counts.index,
            title="Chromatin State Classification",
            color_discrete_sequence=px.colors.qualitative.Set2
        )
        fig.update_layout(height=400)
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        # Heatmap of mark combinations
        mark_cols = ['H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K27me3']

        # Group by chromatin state
        state_mark_means = integration_data.groupby('Chromatin_State')[mark_cols].mean()

        fig = px.imshow(
            state_mark_means,
            color_continuous_scale='RdYlBu_r',
            title="Histone Mark Enrichment by State"
        )
        fig.update_layout(height=400)
        st.plotly_chart(fig, use_container_width=True)

    # Summary table
    st.subheader("Classification Summary")
    summary = integration_data['Chromatin_State'].value_counts().reset_index()
    summary.columns = ['Chromatin State', 'Count']
    summary['Percentage'] = (summary['Count'] / summary['Count'].sum() * 100).round(1)

    st.dataframe(summary, use_container_width=True, hide_index=True)


if __name__ == "__main__":
    main()
