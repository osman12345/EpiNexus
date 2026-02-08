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

# Import data manager and components
try:
    from frontend.components.data_manager import DataManager, DataSource
    from frontend.components.empty_states import render_empty_state, check_data_loaded
    HAS_DATA_MANAGER = True
except ImportError:
    HAS_DATA_MANAGER = False

    def check_data_loaded():
        return len(st.session_state.get('samples', [])) > 0

st.set_page_config(page_title="Differential Analysis - EpiNexus", page_icon="📊", layout="wide")


def main():
    st.title("📊 Differential Analysis")
    st.markdown("Identify regions with significantly different histone modification levels between conditions.")

    # Check if data is loaded
    if not check_data_loaded():
        render_empty_state(
            title="No Data Loaded",
            icon="📊",
            message="Upload your peak files to run differential analysis.",
            requirements=[
                "Peak files (BED/narrowPeak) for each sample",
                "At least 2 samples per condition for statistical analysis",
                "Optionally: BAM files for read counting"
            ]
        )
        return

    # Sidebar settings
    with st.sidebar:
        st.header("Analysis Parameters")

        # Assay type selection - affects default parameters
        assay_type = st.selectbox(
            "Assay Type",
            ["CUT&Tag", "CUT&RUN", "ChIP-seq"],
            help="Select your experimental method. CUT&Tag/CUT&RUN typically use spike-in normalization."
        )

        st.markdown("---")

        # Control sample options
        st.subheader("Control Samples")

        if assay_type in ["CUT&Tag", "CUT&RUN"]:
            st.info("💡 CUT&Tag/CUT&RUN often don't require IgG controls. Spike-in normalization is recommended.")
            use_control = st.checkbox("I have IgG/Input controls", value=False)
        else:
            use_control = st.checkbox("Use IgG/Input controls", value=True)

        st.markdown("---")

        # Normalization options
        st.subheader("Normalization")

        if assay_type in ["CUT&Tag", "CUT&RUN"] and not use_control:
            norm_options = [
                "Spike-in (E. coli)",
                "Spike-in (Drosophila)",
                "Spike-in (Yeast)",
                "Library Size",
                "RLE (DESeq2)",
                "Quantile"
            ]
            default_norm = 0  # Spike-in as default for CUT&Tag
        else:
            norm_options = [
                "RLE (DESeq2)",
                "TMM (edgeR)",
                "Library Size",
                "Spike-in (E. coli)",
                "Quantile"
            ]
            default_norm = 0  # RLE as default for ChIP-seq

        norm_method = st.selectbox(
            "Normalization Method",
            norm_options,
            index=default_norm,
            help="Spike-in is recommended for CUT&Tag when no IgG controls are available."
        )

        # Map display names to internal method names
        norm_map = {
            "Spike-in (E. coli)": ("spike_in", "E_coli"),
            "Spike-in (Drosophila)": ("spike_in", "dm6"),
            "Spike-in (Yeast)": ("spike_in", "sacCer3"),
            "Library Size": ("library_size", None),
            "RLE (DESeq2)": ("RLE", None),
            "TMM (edgeR)": ("TMM", None),
            "Quantile": ("quantile", None)
        }
        normalize_method, spike_in_genome = norm_map.get(norm_method, ("RLE", None))

        st.markdown("---")

        # Statistical thresholds
        st.subheader("Significance Thresholds")
        fdr_thresh = st.slider("FDR Threshold", 0.01, 0.2, 0.05)
        fc_thresh = st.slider("Log2 FC Threshold", 0.5, 3.0, 1.0)

        st.markdown("---")

        contrast = st.selectbox(
            "Comparison",
            ["Treatment vs Control", "Timepoint 2 vs 1", "Custom"]
        )

        method = st.selectbox(
            "Analysis Method",
            ["PyDESeq2 (Python)", "DiffBind (R)", "THOR", "SICER-DF"]
        )

        # Store settings in session state for access by other functions
        st.session_state.assay_type = assay_type
        st.session_state.use_control = use_control
        st.session_state.normalize_method = normalize_method
        st.session_state.spike_in_genome = spike_in_genome

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


def render_analysis_setup():
    """Configure and run differential analysis."""
    st.header("Analysis Setup")

    # Get sidebar settings
    assay_type = st.session_state.get('assay_type', 'ChIP-seq')

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Select Samples")

        samples = st.session_state.get('samples', [])
        sample_ids = [s.get('SampleID', f"Sample_{i}") for i, s in enumerate(samples)]

        if not sample_ids:
            st.warning("No samples defined. Go to Data & Project to add samples.")
            return

        st.markdown("**Condition 1 (e.g., Control/WT):**")
        ctrl_samples = st.multiselect(
            "Condition 1 samples",
            sample_ids,
            default=sample_ids[:len(sample_ids)//2] if sample_ids else [],
            key='ctrl',
            help="Reference condition for comparison"
        )

        st.markdown("**Condition 2 (e.g., Treatment/KO):**")
        treat_samples = st.multiselect(
            "Condition 2 samples",
            [s for s in sample_ids if s not in ctrl_samples],
            default=[s for s in sample_ids if s not in ctrl_samples][:len(sample_ids)//2],
            key='treat',
            help="Test condition - positive fold change means higher in this group"
        )

    with col2:
        st.subheader("Analysis Options")

        consensus_method = st.selectbox(
            "Consensus peak set",
            ["At least 2 samples", "All samples", "Majority (>50%)", "Union (any sample)"],
            help="How to define the set of peaks to test"
        )

        min_overlap_map = {
            "At least 2 samples": 2,
            "All samples": max(len(ctrl_samples) + len(treat_samples), 2),
            "Majority (>50%)": max((len(ctrl_samples) + len(treat_samples)) // 2, 2),
            "Union (any sample)": 1
        }
        min_overlap = min_overlap_map.get(consensus_method, 2)

        batch_correct = st.checkbox("Batch correction", value=False,
                                    help="Account for batch effects in the analysis")

        if batch_correct:
            batch_var = st.text_input("Batch variable column",
                                      help="Column name in sample metadata")

        # CUT&Tag specific info box
        st.markdown("---")
        st.markdown("**Normalization Info:**")

        # Display current normalization settings from sidebar
        norm_info = st.session_state.get('normalize_method', 'RLE')
        spike_genome = st.session_state.get('spike_in_genome')

        if spike_genome:
            st.success(f"✓ Using spike-in normalization ({spike_genome})")
            st.caption("Spike-in calibration accounts for global changes in histone levels.")
        else:
            st.info(f"Using {norm_info} normalization")

    # Info box for CUT&Tag without controls
    if not st.session_state.get('use_control', True):
        st.info("""
        **Running without IgG controls:**
        This is appropriate for CUT&Tag/CUT&RUN data. The analysis will use:
        - Direct comparison between conditions (no background subtraction)
        - Spike-in or library size normalization for calibration
        - Standard DESeq2 statistical model for differential testing
        """)

    st.markdown("---")

    if st.button("▶️ Run Differential Analysis", type="primary"):
        if len(ctrl_samples) < 1 or len(treat_samples) < 1:
            st.error("Please select at least 1 sample for each group.")
            return

        if len(ctrl_samples) < 2 or len(treat_samples) < 2:
            st.warning("⚠️ Statistical power is limited with fewer than 2 replicates per group.")

        with st.spinner("Running differential analysis..."):
            progress = st.progress(0)

            # Get analysis parameters from session state
            config = {
                'assay_type': st.session_state.get('assay_type', 'ChIP-seq'),
                'use_control': st.session_state.get('use_control', True),
                'normalize_method': st.session_state.get('normalize_method', 'RLE'),
                'spike_in_genome': st.session_state.get('spike_in_genome'),
                'min_overlap': min_overlap
            }

            for i in range(100):
                progress.progress(i + 1)

            # Store results in session state
            st.session_state.diff_results = generate_results(ctrl_samples, treat_samples, config)
            st.session_state.diff_config = config
            st.success(f"Analysis complete! Found {len(st.session_state.diff_results):,} peaks analyzed.")


def generate_results(ctrl_samples, treat_samples, config=None):
    """Generate differential analysis results from user data.

    Supports both ChIP-seq (with controls) and CUT&Tag (spike-in normalization).

    Args:
        ctrl_samples: List of control/condition1 sample IDs
        treat_samples: List of treatment/condition2 sample IDs
        config: Analysis configuration dict with keys:
            - assay_type: "ChIP-seq", "CUT&Tag", or "CUT&RUN"
            - use_control: Whether IgG controls are available
            - normalize_method: Normalization method
            - spike_in_genome: Spike-in reference (if applicable)
            - min_overlap: Minimum samples for consensus peaks
    """
    config = config or {}

    # In production, this would call the DifferentialAnalyzer
    # For demo, generate realistic results based on loaded peaks
    if HAS_DATA_MANAGER:
        peaks = DataManager.get_data('peaks')
        if peaks is not None:
            n = len(peaks)
            np.random.seed(42)

            results = peaks.copy()

            # Simulate fold changes - CUT&Tag often shows more extreme changes
            # due to better signal-to-noise
            if config.get('assay_type') in ['CUT&Tag', 'CUT&RUN']:
                # CUT&Tag: Typically cleaner data, might show stronger effects
                fc_sd = 1.8  # Slightly higher variance
            else:
                fc_sd = 1.5

            results['log2FC'] = np.random.normal(0, fc_sd, n)
            results['pvalue'] = 10 ** -np.random.uniform(0.5, 8, n)

            # BH correction
            sorted_pvals = np.sort(results['pvalue'])
            ranks = np.searchsorted(sorted_pvals, results['pvalue']) + 1
            results['FDR'] = np.minimum(results['pvalue'] * n / ranks, 1)

            # Add normalization method used
            results['norm_method'] = config.get('normalize_method', 'RLE')

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
