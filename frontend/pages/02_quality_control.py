"""
Quality Control Page
QC metrics and visualizations for ChIP-seq data.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

st.set_page_config(page_title="Quality Control - EpiNexus", page_icon="✅", layout="wide")

# Try to import data manager
try:
    from frontend.components.data_manager import DataManager
    HAS_DATA_MANAGER = True
except ImportError:
    HAS_DATA_MANAGER = False


def has_data():
    """Check if user has loaded data."""
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
            <div style="font-size: 3rem; margin-bottom: 1rem;">✅</div>
            <h2 style="color: #6c757d; margin-bottom: 0.5rem;">No Data Loaded</h2>
            <p style="color: #6c757d; font-size: 1.1rem;">
                Upload your samples to view quality control metrics.
            </p>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("")

        if st.button("📁 Go to Data & Project", type="primary", use_container_width=True):
            st.switch_page("pages/01_data_project.py")

        st.markdown("")
        st.markdown("**What you need:**")
        st.markdown("- BAM files from alignment")
        st.markdown("- Or peak files (BED format)")
        st.markdown("- Sample metadata")


def main():
    st.title("✅ Quality Control")
    st.markdown("Assess the quality of your ChIP-seq samples and identify potential issues.")

    # Check if data is loaded
    if not has_data():
        render_empty_state()
        return

    tab1, tab2, tab3, tab4 = st.tabs([
        "📊 Sample Summary",
        "📈 Peak Metrics",
        "🔗 Correlation",
        "📉 Fragment Analysis"
    ])

    with tab1:
        render_sample_summary()
    with tab2:
        render_peak_metrics()
    with tab3:
        render_correlation_analysis()
    with tab4:
        render_fragment_analysis()


def render_sample_summary():
    """Overall sample quality summary."""
    st.header("Sample Quality Summary")

    # Demo data
    qc_data = pd.DataFrame({
        'Sample': ['H3K4me3_Ctrl_R1', 'H3K4me3_Ctrl_R2', 'H3K4me3_Treat_R1', 'H3K4me3_Treat_R2',
                   'H3K27ac_Ctrl_R1', 'H3K27ac_Ctrl_R2', 'H3K27ac_Treat_R1', 'H3K27ac_Treat_R2'],
        'Total Reads': [45e6, 42e6, 48e6, 44e6, 52e6, 49e6, 51e6, 47e6],
        'Mapped %': [95.2, 94.8, 95.5, 94.9, 96.1, 95.7, 95.9, 95.3],
        'Peaks': [28500, 27800, 31200, 29900, 42100, 40800, 45600, 43200],
        'FRiP': [0.32, 0.30, 0.35, 0.33, 0.28, 0.27, 0.31, 0.29],
        'NSC': [1.8, 1.7, 1.9, 1.8, 1.6, 1.5, 1.7, 1.6],
        'RSC': [1.2, 1.1, 1.3, 1.2, 1.0, 0.9, 1.1, 1.0],
        'Status': ['✅', '✅', '✅', '✅', '✅', '⚠️', '✅', '✅']
    })

    # Summary metrics
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Total Samples", len(qc_data))
    col2.metric("Passing QC", (qc_data['Status'] == '✅').sum())
    col3.metric("Avg FRiP", f"{qc_data['FRiP'].mean():.2%}")
    col4.metric("Avg Peaks", f"{qc_data['Peaks'].mean():,.0f}")

    st.dataframe(qc_data, use_container_width=True, hide_index=True)

    # QC thresholds
    with st.expander("QC Thresholds"):
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("""
            **Recommended Thresholds:**
            - Mapped reads: >90%
            - FRiP (broad marks): >0.1
            - FRiP (narrow marks): >0.2
            - NSC: >1.1
            - RSC: >0.8
            """)
        with col2:
            st.markdown("""
            **Status Legend:**
            - ✅ Pass: Meets all thresholds
            - ⚠️ Warning: Below one threshold
            - ❌ Fail: Multiple issues
            """)


def render_peak_metrics():
    """Peak-level quality metrics."""
    st.header("Peak Quality Metrics")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Peak Count Distribution")
        samples = ['Ctrl_R1', 'Ctrl_R2', 'Treat_R1', 'Treat_R2']
        marks = ['H3K4me3', 'H3K27ac', 'H3K27me3']

        peak_counts = pd.DataFrame({
            'Sample': samples * 3,
            'Mark': [m for m in marks for _ in samples],
            'Peaks': [28500, 27800, 31200, 29900,
                     42100, 40800, 45600, 43200,
                     15200, 14800, 18100, 17200]
        })

        fig = px.bar(peak_counts, x='Sample', y='Peaks', color='Mark', barmode='group')
        fig.update_layout(height=400)
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("Peak Width Distribution")
        np.random.seed(42)
        widths = {
            'H3K4me3': np.random.lognormal(6, 0.5, 1000),
            'H3K27ac': np.random.lognormal(6.5, 0.6, 1000),
            'H3K27me3': np.random.lognormal(8, 0.8, 1000)
        }

        fig = go.Figure()
        for mark, w in widths.items():
            fig.add_trace(go.Histogram(x=w, name=mark, opacity=0.7, nbinsx=50))
        fig.update_layout(barmode='overlay', xaxis_title='Peak Width (bp)',
                         yaxis_title='Count', height=400)
        st.plotly_chart(fig, use_container_width=True)


def render_correlation_analysis():
    """Sample correlation analysis."""
    st.header("Sample Correlation")

    st.subheader("Peak Signal Correlation Matrix")

    samples = ['H3K4me3_C1', 'H3K4me3_C2', 'H3K4me3_T1', 'H3K4me3_T2',
               'H3K27ac_C1', 'H3K27ac_C2', 'H3K27ac_T1', 'H3K27ac_T2']

    np.random.seed(42)
    corr = np.eye(8)
    # High correlation within mark/condition
    for i in range(0, 8, 2):
        corr[i, i+1] = corr[i+1, i] = 0.95 + np.random.uniform(-0.03, 0.03)
    # Moderate correlation within mark
    for i in range(4):
        corr[i, (i+2)%4] = corr[(i+2)%4, i] = 0.7 + np.random.uniform(-0.1, 0.1)
    for i in range(4, 8):
        corr[i, 4+(i+2)%4] = corr[4+(i+2)%4, i] = 0.7 + np.random.uniform(-0.1, 0.1)

    fig = px.imshow(corr, x=samples, y=samples, color_continuous_scale='RdBu_r',
                    zmin=0, zmax=1, aspect='auto')
    fig.update_layout(height=500)
    st.plotly_chart(fig, use_container_width=True)


def render_fragment_analysis():
    """Fragment size and library complexity."""
    st.header("Fragment Analysis")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Fragment Size Distribution")
        np.random.seed(42)
        sizes = np.concatenate([
            np.random.normal(150, 20, 5000),  # Mononucleosome
            np.random.normal(310, 30, 1500),  # Dinucleosome
        ])

        fig = go.Figure(data=go.Histogram(x=sizes, nbinsx=100))
        fig.update_layout(xaxis_title='Fragment Size (bp)', yaxis_title='Count',
                         height=400)
        fig.add_vline(x=150, line_dash="dash", annotation_text="Mono")
        fig.add_vline(x=310, line_dash="dash", annotation_text="Di")
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("Library Complexity")
        complexity = pd.DataFrame({
            'Sample': ['Ctrl_R1', 'Ctrl_R2', 'Treat_R1', 'Treat_R2'],
            'NRF': [0.82, 0.79, 0.85, 0.81],
            'PBC1': [0.91, 0.88, 0.93, 0.90],
            'PBC2': [15.2, 12.1, 18.5, 14.3]
        })

        fig = px.bar(complexity.melt(id_vars='Sample'), x='Sample', y='value',
                    color='variable', barmode='group', facet_col='variable')
        fig.update_layout(height=400)
        st.plotly_chart(fig, use_container_width=True)


if __name__ == "__main__":
    main()
