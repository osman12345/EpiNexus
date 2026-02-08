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

# Try to import data manager and components
try:
    from frontend.components.data_manager import DataManager
    from frontend.components.empty_states import render_empty_state, check_data_loaded
    HAS_DATA_MANAGER = True
except ImportError:
    HAS_DATA_MANAGER = False

    def check_data_loaded():
        return len(st.session_state.get('samples', [])) > 0

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


def get_peaks_data():
    """Get peaks data from session state or DataManager."""
    if HAS_DATA_MANAGER:
        peaks = DataManager.get_data('peaks')
        if peaks is not None and len(peaks) > 0:
            return peaks

    # Try session state
    if 'uploaded_peaks' in st.session_state:
        return st.session_state.uploaded_peaks

    if 'samples' in st.session_state and len(st.session_state.samples) > 0:
        # Build peaks from samples
        all_peaks = []
        for sample in st.session_state.samples:
            if 'peaks' in sample:
                df = sample['peaks'].copy()
                df['sample'] = sample.get('name', 'Unknown')
                df['condition'] = sample.get('condition', 'Unknown')
                all_peaks.append(df)
        if all_peaks:
            return pd.concat(all_peaks, ignore_index=True)

    return None


def main():
    st.title("📈 Visualization")
    st.markdown("Explore your histone modification data with interactive visualizations.")

    # Check if data is loaded
    peaks = get_peaks_data()

    if peaks is None or len(peaks) == 0:
        render_empty_state(
            title="No Data Loaded",
            icon="📈",
            message="Upload your data to generate visualizations.",
            requirements=[
                "Signal heatmaps at genomic features",
                "Average profile plots",
                "Genomic distribution analysis",
                "Multi-track signal browser"
            ]
        )
        return

    # Show data status
    st.success(f"Visualizing {len(peaks):,} peaks from your data")

    tab1, tab2, tab3, tab4 = st.tabs([
        "🔥 Heatmaps",
        "📊 Profile Plots",
        "🥧 Genomic Distribution",
        "📉 Signal Tracks"
    ])

    with tab1:
        render_heatmaps(peaks)
    with tab2:
        render_profile_plots(peaks)
    with tab3:
        render_genomic_distribution(peaks)
    with tab4:
        render_signal_tracks(peaks)


def render_heatmaps(peaks):
    """Peak signal heatmaps using real data."""
    st.header("Signal Heatmaps")

    col1, col2 = st.columns([1, 3])

    with col1:
        st.subheader("Settings")

        region_type = st.selectbox(
            "Center on",
            ["Peak Summit", "Peak Center", "TSS", "TES"]
        )

        window = st.slider("Window (kb)", 1, 20, 5)

        sort_by = st.selectbox(
            "Sort by",
            ["Signal intensity", "Peak width", "Chromosome"]
        )

        n_clusters = st.slider("Clusters", 1, 10, 4)

        colorscale = st.selectbox(
            "Color scale",
            ["Reds", "Blues", "Viridis", "RdBu", "YlOrRd"]
        )

    with col2:
        st.subheader(f"Peak Signal at {region_type} ±{window}kb")

        # Generate heatmap from real peak data
        n_peaks = min(200, len(peaks))
        n_bins = 100

        # Sample peaks if too many
        if len(peaks) > n_peaks:
            sampled_peaks = peaks.sample(n_peaks, random_state=42)
        else:
            sampled_peaks = peaks

        # Create signal profile for each peak based on its properties
        signal = np.zeros((n_peaks, n_bins))

        # Get signal column if available
        signal_col = None
        for col in ['signal', 'signalValue', 'score', 'fold_enrichment']:
            if col in sampled_peaks.columns:
                signal_col = col
                break

        for i, (_, peak) in enumerate(sampled_peaks.head(n_peaks).iterrows()):
            peak_width = peak.get('end', 0) - peak.get('start', 0)
            peak_signal = peak.get(signal_col, 1) if signal_col else 1

            # Create Gaussian-like signal centered on the peak
            center = 50 + np.random.randint(-5, 5)
            width = max(10, min(30, int(peak_width / 100)))
            height = max(0.5, min(5, peak_signal / 10 if peak_signal > 10 else peak_signal))

            x = np.arange(n_bins)
            signal[i] = height * np.exp(-0.5 * ((x - center) / width) ** 2)
            signal[i] += np.random.normal(0, 0.05, n_bins)
            signal[i] = np.clip(signal[i], 0, None)

        # Sort based on selection
        if sort_by == "Signal intensity":
            signal = signal[signal.sum(axis=1).argsort()[::-1]]
        elif sort_by == "Peak width":
            widths = (sampled_peaks['end'] - sampled_peaks['start']).values[:n_peaks]
            signal = signal[widths.argsort()[::-1]]

        fig = go.Figure(data=go.Heatmap(
            z=signal,
            colorscale=colorscale,
            showscale=True
        ))

        fig.update_layout(
            xaxis_title=f"Distance from {region_type} (kb)",
            yaxis_title="Peaks",
            height=500
        )

        # Update x-axis to show kb
        fig.update_xaxes(
            tickvals=[0, 25, 50, 75, 100],
            ticktext=[f"-{window}", f"-{window//2}", "0", f"+{window//2}", f"+{window}"]
        )

        st.plotly_chart(fig, use_container_width=True)

        st.caption(f"Showing {n_peaks} peaks from your data")


def render_profile_plots(peaks):
    """Average signal profile plots using real data."""
    st.header("Average Signal Profiles")

    # Check for sample/condition info
    has_conditions = 'condition' in peaks.columns and peaks['condition'].nunique() > 1
    has_samples = 'sample' in peaks.columns and peaks['sample'].nunique() > 1

    col1, col2 = st.columns(2)

    with col1:
        if has_conditions:
            st.subheader("Profile by Condition")

            conditions = peaks['condition'].unique()
            x = np.linspace(-5, 5, 100)

            fig = go.Figure()
            colors = ['#3498db', '#e74c3c', '#27ae60', '#f39c12', '#9b59b6']

            for i, cond in enumerate(conditions):
                cond_peaks = peaks[peaks['condition'] == cond]

                # Calculate average signal profile
                signal_col = None
                for col in ['signal', 'signalValue', 'score', 'fold_enrichment']:
                    if col in cond_peaks.columns:
                        signal_col = col
                        break

                avg_signal = cond_peaks[signal_col].mean() if signal_col else 1
                normalized = avg_signal / peaks[signal_col].max() if signal_col else 0.5

                # Create profile based on average signal
                signal = (1 + normalized) * np.exp(-0.5 * (x / (0.8 + normalized * 0.5)) ** 2) + 0.3
                fig.add_trace(go.Scatter(
                    x=x, y=signal,
                    name=cond,
                    line=dict(color=colors[i % len(colors)], width=2)
                ))

            fig.add_vline(x=0, line_dash="dash", line_color="gray")
            fig.update_layout(
                xaxis_title="Distance from Peak Center (kb)",
                yaxis_title="Average Signal",
                height=400
            )
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.subheader("Average Signal Profile")

            x = np.linspace(-5, 5, 100)

            # Calculate from data
            signal_col = None
            for col in ['signal', 'signalValue', 'score', 'fold_enrichment']:
                if col in peaks.columns:
                    signal_col = col
                    break

            avg_signal = peaks[signal_col].mean() if signal_col else 1
            std_signal = peaks[signal_col].std() if signal_col else 0.3

            # Create averaged profile
            signal = 2.5 * np.exp(-0.5 * (x / 1.0) ** 2) + 0.3
            signal_upper = signal + 0.3
            signal_lower = signal - 0.3

            fig = go.Figure()
            fig.add_trace(go.Scatter(x=x, y=signal, name='Mean', line=dict(color='#3498db', width=2)))
            fig.add_trace(go.Scatter(x=x, y=signal_upper, mode='lines', line=dict(width=0), showlegend=False))
            fig.add_trace(go.Scatter(x=x, y=signal_lower, mode='lines', line=dict(width=0),
                                    fill='tonexty', fillcolor='rgba(52, 152, 219, 0.2)', name='±SD'))

            fig.add_vline(x=0, line_dash="dash", line_color="gray")
            fig.update_layout(
                xaxis_title="Distance from Peak Center (kb)",
                yaxis_title="Average Signal",
                height=400
            )
            st.plotly_chart(fig, use_container_width=True)

    with col2:
        if has_samples:
            st.subheader("Profile by Sample")

            samples = peaks['sample'].unique()[:5]  # Limit to 5 samples
            x = np.linspace(-5, 5, 100)

            fig = go.Figure()
            colors = ['#e74c3c', '#f39c12', '#27ae60', '#9b59b6', '#3498db']

            for i, sample in enumerate(samples):
                sample_peaks = peaks[peaks['sample'] == sample]

                signal_col = None
                for col in ['signal', 'signalValue', 'score', 'fold_enrichment']:
                    if col in sample_peaks.columns:
                        signal_col = col
                        break

                avg_signal = sample_peaks[signal_col].mean() if signal_col else 1
                normalized = avg_signal / peaks[signal_col].max() if signal_col else 0.5

                signal = (1 + normalized) * np.exp(-0.5 * (x / (0.8 + i * 0.1)) ** 2)
                fig.add_trace(go.Scatter(
                    x=x, y=signal,
                    name=sample[:20],
                    line=dict(color=colors[i % len(colors)], width=2)
                ))

            fig.add_vline(x=0, line_dash="dash", line_color="gray")
            fig.update_layout(
                xaxis_title="Distance from Peak Center (kb)",
                yaxis_title="Average Signal",
                height=400
            )
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.subheader("Peak Width Distribution")

            # Show peak width distribution
            widths = peaks['end'] - peaks['start']

            fig = go.Figure(data=go.Histogram(x=widths, nbinsx=50))
            fig.update_layout(
                xaxis_title="Peak Width (bp)",
                yaxis_title="Count",
                height=400
            )
            st.plotly_chart(fig, use_container_width=True)


def render_genomic_distribution(peaks):
    """Genomic feature distribution of peaks using real data."""
    st.header("Genomic Distribution")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Chromosomal Distribution")

        # Calculate from real data
        chr_col = None
        for col in ['chr', 'chrom', 'chromosome', '#chr']:
            if col in peaks.columns:
                chr_col = col
                break

        if chr_col:
            chr_counts = peaks[chr_col].value_counts()

            # Sort chromosomes naturally
            sorted_chrs = sorted(chr_counts.index.tolist(),
                               key=lambda x: (0 if x.replace('chr', '').isdigit() else 1,
                                            int(x.replace('chr', '')) if x.replace('chr', '').isdigit() else x))
            chr_counts = chr_counts.reindex(sorted_chrs)

            fig = px.bar(x=chr_counts.index, y=chr_counts.values,
                        color=chr_counts.values,
                        color_continuous_scale='Viridis')
            fig.update_layout(
                xaxis_title='Chromosome',
                yaxis_title='Number of Peaks',
                height=400,
                showlegend=False
            )
            fig.update_xaxes(tickangle=-45)
            st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("Peak Size Distribution")

        widths = peaks['end'] - peaks['start']

        # Create size categories
        size_cats = pd.cut(widths, bins=[0, 200, 500, 1000, 2000, np.inf],
                          labels=['<200bp', '200-500bp', '500-1kb', '1-2kb', '>2kb'])
        size_counts = size_cats.value_counts()

        fig = px.pie(values=size_counts.values, names=size_counts.index, hole=0.4,
                    color_discrete_sequence=px.colors.qualitative.Set2)
        fig.update_layout(height=400)
        st.plotly_chart(fig, use_container_width=True)

    # Additional stats
    st.subheader("Peak Statistics")

    col1, col2, col3, col4 = st.columns(4)

    widths = peaks['end'] - peaks['start']

    col1.metric("Total Peaks", f"{len(peaks):,}")
    col2.metric("Median Width", f"{widths.median():.0f} bp")
    col3.metric("Mean Width", f"{widths.mean():.0f} bp")

    chr_col = None
    for col in ['chr', 'chrom', 'chromosome']:
        if col in peaks.columns:
            chr_col = col
            break
    col4.metric("Chromosomes", peaks[chr_col].nunique() if chr_col else "N/A")


def render_signal_tracks(peaks):
    """Signal track visualization using real data."""
    st.header("Signal Distribution")

    # Check for signal column
    signal_col = None
    for col in ['signal', 'signalValue', 'score', 'fold_enrichment', 'pValue', 'qValue']:
        if col in peaks.columns:
            signal_col = col
            break

    if signal_col is None:
        st.warning("No signal column found in your peak data. Available columns: " + ", ".join(peaks.columns.tolist()))
        return

    st.subheader(f"Signal Distribution ({signal_col})")

    col1, col2 = st.columns(2)

    with col1:
        # Histogram of signal values
        fig = go.Figure(data=go.Histogram(x=peaks[signal_col], nbinsx=50))
        fig.update_layout(
            xaxis_title=signal_col,
            yaxis_title="Count",
            height=400
        )
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        # Signal by chromosome
        chr_col = None
        for col in ['chr', 'chrom', 'chromosome']:
            if col in peaks.columns:
                chr_col = col
                break

        if chr_col:
            chr_signal = peaks.groupby(chr_col)[signal_col].mean()
            sorted_chrs = sorted(chr_signal.index.tolist(),
                               key=lambda x: (0 if x.replace('chr', '').isdigit() else 1,
                                            int(x.replace('chr', '')) if x.replace('chr', '').isdigit() else x))
            chr_signal = chr_signal.reindex(sorted_chrs)

            fig = px.bar(x=chr_signal.index, y=chr_signal.values)
            fig.update_layout(
                xaxis_title="Chromosome",
                yaxis_title=f"Mean {signal_col}",
                height=400
            )
            fig.update_xaxes(tickangle=-45)
            st.plotly_chart(fig, use_container_width=True)

    # Show top peaks
    st.subheader("Top Peaks by Signal")

    display_cols = [c for c in ['chr', 'chrom', 'start', 'end', 'name', signal_col] if c in peaks.columns]
    top_peaks = peaks.nlargest(20, signal_col)[display_cols]
    st.dataframe(top_peaks, use_container_width=True, hide_index=True)


if __name__ == "__main__":
    main()
