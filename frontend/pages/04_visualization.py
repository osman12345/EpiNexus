"""
Visualization Page
Interactive visualizations for histone modification data.
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

st.set_page_config(page_title="Visualization - EpiNexus", page_icon="📈", layout="wide")

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
            <div style="font-size: 3rem; margin-bottom: 1rem;">📈</div>
            <h2 style="color: #6c757d; margin-bottom: 0.5rem;">No Data Loaded</h2>
            <p style="color: #6c757d; font-size: 1.1rem;">
                Upload your data to generate visualizations.
            </p>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("")

        if st.button("📁 Go to Data & Project", type="primary", use_container_width=True):
            st.switch_page("pages/01_data_project.py")

        st.markdown("")
        st.markdown("**Available visualizations:**")
        st.markdown("- Signal heatmaps at genomic features")
        st.markdown("- Average profile plots")
        st.markdown("- Genomic distribution analysis")
        st.markdown("- Multi-track signal browser")


def main():
    st.title("📈 Visualization")
    st.markdown("Explore your histone modification data with interactive visualizations.")

    # Check if data is loaded
    if not has_data():
        render_empty_state()
        return

    tab1, tab2, tab3, tab4 = st.tabs([
        "🔥 Heatmaps",
        "📊 Profile Plots",
        "🥧 Genomic Distribution",
        "📉 Signal Tracks"
    ])

    with tab1:
        render_heatmaps()
    with tab2:
        render_profile_plots()
    with tab3:
        render_genomic_distribution()
    with tab4:
        render_signal_tracks()


def render_heatmaps():
    """Peak signal heatmaps."""
    st.header("Signal Heatmaps")

    col1, col2 = st.columns([1, 3])

    with col1:
        st.subheader("Settings")

        region_type = st.selectbox(
            "Center on",
            ["TSS", "Peak Summit", "Peak Center", "TES"]
        )

        window = st.slider("Window (kb)", 1, 20, 5)

        sort_by = st.selectbox(
            "Sort by",
            ["Signal intensity", "Cluster", "Gene expression", "Peak width"]
        )

        n_clusters = st.slider("Clusters", 1, 10, 4)

        colorscale = st.selectbox(
            "Color scale",
            ["Reds", "Blues", "Viridis", "RdBu", "YlOrRd"]
        )

    with col2:
        st.subheader("H3K27ac Signal at TSS ±5kb")

        # Generate demo heatmap data
        np.random.seed(42)
        n_genes = 200
        n_bins = 100

        # Create signal patterns
        signal = np.zeros((n_genes, n_bins))
        for i in range(n_genes):
            center = 50 + np.random.randint(-10, 10)
            width = np.random.randint(10, 30)
            height = np.random.uniform(0.5, 3)
            x = np.arange(n_bins)
            signal[i] = height * np.exp(-0.5 * ((x - center) / width) ** 2)
            signal[i] += np.random.normal(0, 0.1, n_bins)

        # Sort by total signal
        signal = signal[signal.sum(axis=1).argsort()[::-1]]

        fig = go.Figure(data=go.Heatmap(
            z=signal,
            colorscale=colorscale,
            showscale=True
        ))

        fig.update_layout(
            xaxis_title=f"Distance from {region_type} (kb)",
            yaxis_title="Genes/Peaks",
            height=500
        )

        # Update x-axis to show kb
        fig.update_xaxes(
            tickvals=[0, 25, 50, 75, 100],
            ticktext=[f"-{window}", f"-{window//2}", "0", f"+{window//2}", f"+{window}"]
        )

        st.plotly_chart(fig, use_container_width=True)


def render_profile_plots():
    """Average signal profile plots."""
    st.header("Average Signal Profiles")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("TSS Profile by Condition")

        x = np.linspace(-5, 5, 100)

        fig = go.Figure()

        # Control - typical H3K4me3 pattern
        ctrl_signal = 2 * np.exp(-0.5 * (x / 0.8) ** 2) + 0.3
        fig.add_trace(go.Scatter(x=x, y=ctrl_signal, name='Control', line=dict(color='#3498db', width=2)))

        # Treatment - enhanced signal
        treat_signal = 2.8 * np.exp(-0.5 * (x / 0.7) ** 2) + 0.4
        fig.add_trace(go.Scatter(x=x, y=treat_signal, name='Treatment', line=dict(color='#e74c3c', width=2)))

        fig.add_vline(x=0, line_dash="dash", line_color="gray")
        fig.update_layout(
            xaxis_title="Distance from TSS (kb)",
            yaxis_title="Average Signal",
            height=400
        )
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("Multi-mark Profile")

        fig = go.Figure()

        marks = {
            'H3K4me3': (2.5, 0.6, '#e74c3c'),
            'H3K27ac': (1.8, 1.2, '#f39c12'),
            'H3K4me1': (0.8, 2.0, '#27ae60'),
            'H3K27me3': (0.3, 3.0, '#9b59b6')
        }

        for mark, (height, width, color) in marks.items():
            signal = height * np.exp(-0.5 * (x / width) ** 2)
            fig.add_trace(go.Scatter(x=x, y=signal, name=mark, line=dict(color=color, width=2)))

        fig.add_vline(x=0, line_dash="dash", line_color="gray")
        fig.update_layout(
            xaxis_title="Distance from TSS (kb)",
            yaxis_title="Average Signal",
            height=400
        )
        st.plotly_chart(fig, use_container_width=True)


def render_genomic_distribution():
    """Genomic feature distribution of peaks."""
    st.header("Genomic Distribution")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Peak Annotation")

        # Pie chart of genomic features
        features = ['Promoter', 'Exon', 'Intron', 'Intergenic', "5' UTR", "3' UTR"]
        values = [35, 8, 28, 22, 3, 4]

        fig = px.pie(values=values, names=features, hole=0.4,
                    color_discrete_sequence=px.colors.qualitative.Set2)
        fig.update_layout(height=400)
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("Distance to TSS")

        np.random.seed(42)
        distances = np.concatenate([
            np.random.normal(0, 500, 2000),
            np.random.normal(-2000, 1000, 800),
            np.random.uniform(-50000, 50000, 1200)
        ])

        fig = go.Figure(data=go.Histogram(x=distances, nbinsx=100))
        fig.update_layout(
            xaxis_title="Distance to nearest TSS (bp)",
            yaxis_title="Peak count",
            height=400
        )
        fig.add_vline(x=0, line_dash="dash", line_color="red")
        st.plotly_chart(fig, use_container_width=True)

    st.subheader("Feature Overlap by Mark")

    marks = ['H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K27me3', 'H3K36me3']
    features = ['Promoter', 'Enhancer', 'Gene Body', 'Intergenic']

    # Feature percentages for each mark
    data = {
        'H3K4me3': [75, 10, 5, 10],
        'H3K27ac': [40, 45, 8, 7],
        'H3K4me1': [15, 60, 15, 10],
        'H3K27me3': [25, 5, 20, 50],
        'H3K36me3': [5, 5, 80, 10]
    }

    df = pd.DataFrame(data, index=features).T

    fig = px.bar(df, barmode='stack', color_discrete_sequence=px.colors.qualitative.Set2)
    fig.update_layout(
        xaxis_title="Histone Mark",
        yaxis_title="Percentage of Peaks",
        height=400,
        legend_title="Genomic Feature"
    )
    st.plotly_chart(fig, use_container_width=True)


def render_signal_tracks():
    """Signal track visualization."""
    st.header("Signal Tracks")

    # Region selector
    col1, col2, col3 = st.columns([2, 1, 1])
    with col1:
        gene = st.selectbox("Select gene", ["MYC", "BRCA1", "TP53", "EGFR", "GAPDH"])
    with col2:
        upstream = st.number_input("Upstream (kb)", 1, 50, 10)
    with col3:
        downstream = st.number_input("Downstream (kb)", 1, 50, 10)

    st.markdown("---")

    # Create multi-track figure
    fig = make_subplots(rows=5, cols=1, shared_xaxes=True, vertical_spacing=0.02,
                       row_heights=[0.2, 0.2, 0.2, 0.2, 0.2],
                       subplot_titles=['H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K27me3', 'RNA-seq'])

    np.random.seed(42)
    x = np.linspace(0, 100, 500)

    # Generate track data
    tracks = {
        'H3K4me3': (np.exp(-0.5 * ((x - 20) / 5) ** 2) * 3, '#e74c3c'),
        'H3K27ac': (np.exp(-0.5 * ((x - 20) / 8) ** 2) * 2.5 + np.exp(-0.5 * ((x - 60) / 10) ** 2) * 2, '#f39c12'),
        'H3K4me1': (np.exp(-0.5 * ((x - 60) / 12) ** 2) * 2, '#27ae60'),
        'H3K27me3': (np.ones_like(x) * 0.2 + np.random.normal(0, 0.05, len(x)), '#9b59b6'),
        'RNA-seq': (np.where((x > 20) & (x < 80), np.random.uniform(1, 3, len(x)), 0.1), '#3498db')
    }

    for i, (name, (y, color)) in enumerate(tracks.items(), 1):
        fig.add_trace(go.Scatter(x=x, y=y, fill='tozeroy', name=name,
                                line=dict(color=color, width=1)), row=i, col=1)

    fig.update_layout(height=700, showlegend=False)
    fig.update_xaxes(title_text="Position (kb)", row=5, col=1)

    st.plotly_chart(fig, use_container_width=True)


if __name__ == "__main__":
    main()
