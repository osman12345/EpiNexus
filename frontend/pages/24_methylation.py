"""
DNA Methylation Analysis Module

Comprehensive analysis of DNA methylation data from bisulfite sequencing
(WGBS, RRBS, EM-seq) and methylation arrays (450K, EPIC).

Copyright (c) 2026 EpiNexus Contributors
SPDX-License-Identifier: AGPL-3.0-or-later OR Commercial
"""

import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import subprocess
import tempfile
import sys
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

# Import workflow manager for step recording
try:
    from frontend.components.workflow_manager import WorkflowManager
    HAS_WORKFLOW_MANAGER = True
except ImportError:
    HAS_WORKFLOW_MANAGER = False

# =============================================================================
# CONFIGURATION
# =============================================================================

SUPPORTED_FORMATS = {
    "Methylation Calls": [".bedGraph", ".cov", ".bismark.cov", ".methylKit"],
    "Alignments": [".bam"],
    "Arrays": [".idat", ".csv"],
}

ANALYSIS_TYPES = [
    "Bisulfite Sequencing (WGBS/RRBS)",
    "Enzymatic Methylation (EM-seq)",
    "Methylation Array (450K/EPIC)",
]

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def check_tool(tool_name: str) -> bool:
    """Check if a command-line tool is available."""
    try:
        subprocess.run([tool_name, "--version"], capture_output=True, check=False)
        return True
    except FileNotFoundError:
        return False

def parse_bismark_coverage(file_path: str) -> pd.DataFrame:
    """Parse Bismark coverage file format."""
    # Format: chr, start, end, methylation%, count_methylated, count_unmethylated
    df = pd.read_csv(file_path, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'meth_pct', 'count_meth', 'count_unmeth'])
    df['coverage'] = df['count_meth'] + df['count_unmeth']
    return df

def parse_bedgraph(file_path: str) -> pd.DataFrame:
    """Parse bedGraph methylation file."""
    df = pd.read_csv(file_path, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'meth_pct'])
    return df

def calculate_methylation_stats(df: pd.DataFrame) -> dict:
    """Calculate summary statistics for methylation data."""
    stats = {
        'total_cpgs': len(df),
        'mean_methylation': df['meth_pct'].mean(),
        'median_methylation': df['meth_pct'].median(),
        'std_methylation': df['meth_pct'].std(),
    }

    if 'coverage' in df.columns:
        stats['mean_coverage'] = df['coverage'].mean()
        stats['median_coverage'] = df['coverage'].median()
        stats['cpgs_10x'] = (df['coverage'] >= 10).sum()
        stats['cpgs_5x'] = (df['coverage'] >= 5).sum()

    # Methylation level distribution
    stats['hypomethylated'] = (df['meth_pct'] < 20).sum()
    stats['intermediate'] = ((df['meth_pct'] >= 20) & (df['meth_pct'] <= 80)).sum()
    stats['hypermethylated'] = (df['meth_pct'] > 80).sum()

    return stats

def generate_demo_methylation_data(n_cpgs: int = 50000) -> pd.DataFrame:
    """Generate demo methylation data for testing."""
    np.random.seed(42)

    chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    chr_sizes = {f'chr{i}': 250_000_000 - i * 5_000_000 for i in range(1, 23)}
    chr_sizes['chrX'] = 155_000_000
    chr_sizes['chrY'] = 60_000_000

    data = []
    for _ in range(n_cpgs):
        chrom = np.random.choice(chromosomes, p=[0.08]*22 + [0.04, 0.02])
        pos = np.random.randint(1, chr_sizes.get(chrom, 100_000_000))

        # Bimodal methylation distribution (typical for CpGs)
        if np.random.random() < 0.7:
            meth = np.clip(np.random.beta(8, 2) * 100, 0, 100)  # High methylation
        else:
            meth = np.clip(np.random.beta(2, 8) * 100, 0, 100)  # Low methylation

        coverage = max(1, int(np.random.exponential(15)))
        count_meth = int(round(coverage * meth / 100))
        count_unmeth = coverage - count_meth

        data.append({
            'chr': chrom,
            'start': pos,
            'end': pos + 1,
            'meth_pct': meth,
            'count_meth': count_meth,
            'count_unmeth': count_unmeth,
            'coverage': coverage,
        })

    return pd.DataFrame(data).sort_values(['chr', 'start']).reset_index(drop=True)

# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

def plot_methylation_distribution(df: pd.DataFrame, title: str = "Methylation Distribution"):
    """Plot histogram of methylation levels."""
    fig = px.histogram(
        df, x='meth_pct',
        nbins=50,
        title=title,
        labels={'meth_pct': 'Methylation (%)'},
        color_discrete_sequence=['#1f77b4']
    )
    fig.update_layout(
        xaxis_title="Methylation Level (%)",
        yaxis_title="Number of CpGs",
        showlegend=False
    )
    return fig

def plot_coverage_distribution(df: pd.DataFrame, title: str = "Coverage Distribution"):
    """Plot histogram of coverage depth."""
    if 'coverage' not in df.columns:
        return None

    fig = px.histogram(
        df[df['coverage'] <= 100],  # Cap at 100x for visualization
        x='coverage',
        nbins=50,
        title=title,
        color_discrete_sequence=['#2ca02c']
    )
    fig.update_layout(
        xaxis_title="Coverage Depth",
        yaxis_title="Number of CpGs",
        showlegend=False
    )
    return fig

def plot_chromosome_methylation(df: pd.DataFrame, title: str = "Methylation by Chromosome"):
    """Plot average methylation per chromosome."""
    chr_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    chr_meth = df.groupby('chr')['meth_pct'].mean().reindex(chr_order).dropna()

    fig = px.bar(
        x=chr_meth.index,
        y=chr_meth.values,
        title=title,
        labels={'x': 'Chromosome', 'y': 'Mean Methylation (%)'},
        color=chr_meth.values,
        color_continuous_scale='RdYlBu_r'
    )
    fig.update_layout(showlegend=False, coloraxis_showscale=False)
    return fig

def plot_methylation_by_level(stats: dict, title: str = "CpG Methylation Categories"):
    """Plot pie chart of methylation categories."""
    labels = ['Hypomethylated (<20%)', 'Intermediate (20-80%)', 'Hypermethylated (>80%)']
    values = [stats['hypomethylated'], stats['intermediate'], stats['hypermethylated']]
    colors = ['#3498db', '#f39c12', '#e74c3c']

    fig = go.Figure(data=[go.Pie(
        labels=labels,
        values=values,
        marker_colors=colors,
        hole=0.4
    )])
    fig.update_layout(title=title)
    return fig

# =============================================================================
# DMR ANALYSIS
# =============================================================================

def find_dmrs_simple(df1: pd.DataFrame, df2: pd.DataFrame,
                     min_diff: float = 20, min_cpgs: int = 3,
                     window_size: int = 1000) -> pd.DataFrame:
    """
    Simple DMR detection by comparing two samples.

    For production use, consider methylKit, DSS, or dmrseq.
    """
    # Merge on position
    merged = pd.merge(
        df1[['chr', 'start', 'meth_pct']],
        df2[['chr', 'start', 'meth_pct']],
        on=['chr', 'start'],
        suffixes=('_1', '_2')
    )

    merged['diff'] = merged['meth_pct_1'] - merged['meth_pct_2']

    # Simple windowed approach
    dmrs = []
    for chrom in merged['chr'].unique():
        chr_data = merged[merged['chr'] == chrom].sort_values('start')

        i = 0
        while i < len(chr_data):
            window_start = chr_data.iloc[i]['start']
            window_end = window_start + window_size

            window_data = chr_data[
                (chr_data['start'] >= window_start) &
                (chr_data['start'] < window_end)
            ]

            if len(window_data) >= min_cpgs:
                mean_diff = window_data['diff'].mean()
                if abs(mean_diff) >= min_diff:
                    dmrs.append({
                        'chr': chrom,
                        'start': int(window_data['start'].min()),
                        'end': int(window_data['start'].max()) + 1,
                        'n_cpgs': len(window_data),
                        'mean_diff': mean_diff,
                        'mean_meth_1': window_data['meth_pct_1'].mean(),
                        'mean_meth_2': window_data['meth_pct_2'].mean(),
                        'direction': 'hyper' if mean_diff > 0 else 'hypo'
                    })

            i += max(1, len(window_data))

    return pd.DataFrame(dmrs)

# =============================================================================
# MAIN PAGE
# =============================================================================

def main():
    st.title("🧬 DNA Methylation Analysis")
    st.markdown("Analyze bisulfite sequencing and methylation array data")

    # Initialize session state
    if 'methylation_samples' not in st.session_state:
        st.session_state.methylation_samples = {}
    if 'methylation_dmrs' not in st.session_state:
        st.session_state.methylation_dmrs = None

    # Tool availability
    tools = {
        'bismark': check_tool('bismark'),
        'methyldackel': check_tool('MethylDackel'),
    }

    # Tabs
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "📥 Data Input", "📊 QC & Stats", "🔍 DMR Analysis",
        "📈 Visualization", "🔗 Integration"
    ])

    # =========================================================================
    # TAB 1: DATA INPUT
    # =========================================================================
    with tab1:
        st.subheader("Load Methylation Data")

        col1, col2 = st.columns(2)

        with col1:
            data_type = st.selectbox(
                "Data Type",
                ANALYSIS_TYPES
            )

            input_method = st.radio(
                "Input Method",
                ["Upload Files", "Demo Data", "Specify Path"]
            )

        with col2:
            st.markdown("**Supported Formats:**")
            for fmt_type, extensions in SUPPORTED_FORMATS.items():
                st.caption(f"• {fmt_type}: {', '.join(extensions)}")

        st.markdown("---")

        if input_method == "Demo Data":
            if st.button("Load Demo Methylation Data", type="primary"):
                with st.spinner("Generating demo data..."):
                    # Generate two samples for comparison
                    demo1 = generate_demo_methylation_data(30000)
                    demo2 = generate_demo_methylation_data(30000)
                    # Add some differential signal
                    demo2.loc[demo2.index[:5000], 'meth_pct'] *= 0.5

                    st.session_state.methylation_samples = {
                        'Control': demo1,
                        'Treatment': demo2
                    }
                    st.success("Loaded 2 demo samples: Control and Treatment")
                    st.rerun()

        elif input_method == "Upload Files":
            uploaded_files = st.file_uploader(
                "Upload methylation files",
                type=['cov', 'bedgraph', 'txt', 'csv'],
                accept_multiple_files=True,
                help="Bismark .cov or bedGraph format"
            )

            if uploaded_files:
                for uploaded_file in uploaded_files:
                    sample_name = st.text_input(
                        f"Sample name for {uploaded_file.name}",
                        value=Path(uploaded_file.name).stem,
                        key=f"name_{uploaded_file.name}"
                    )

                if st.button("Process Uploaded Files"):
                    with st.spinner("Processing files..."):
                        for uploaded_file in uploaded_files:
                            with tempfile.NamedTemporaryFile(delete=False, suffix='.cov') as tmp:
                                tmp.write(uploaded_file.read())
                                tmp_path = tmp.name

                            try:
                                if '.cov' in uploaded_file.name or 'bismark' in uploaded_file.name.lower():
                                    df = parse_bismark_coverage(tmp_path)
                                else:
                                    df = parse_bedgraph(tmp_path)

                                sample_name = Path(uploaded_file.name).stem
                                st.session_state.methylation_samples[sample_name] = df
                                st.success(f"Loaded {sample_name}: {len(df):,} CpGs")
                            except Exception as e:
                                st.error(f"Error parsing {uploaded_file.name}: {e}")
                            finally:
                                Path(tmp_path).unlink(missing_ok=True)

        else:  # Specify Path
            file_path = st.text_input("Path to methylation file")
            sample_name = st.text_input("Sample name")

            if st.button("Load File") and file_path and sample_name:
                if Path(file_path).exists():
                    with st.spinner("Loading..."):
                        try:
                            if '.cov' in file_path:
                                df = parse_bismark_coverage(file_path)
                            else:
                                df = parse_bedgraph(file_path)
                            st.session_state.methylation_samples[sample_name] = df
                            st.success(f"Loaded {len(df):,} CpGs")
                        except Exception as e:
                            st.error(f"Error: {e}")
                else:
                    st.error("File not found")

        # Show loaded samples
        if st.session_state.methylation_samples:
            st.markdown("---")
            st.subheader("Loaded Samples")

            sample_info = []
            for name, df in st.session_state.methylation_samples.items():
                sample_info.append({
                    'Sample': name,
                    'CpGs': f"{len(df):,}",
                    'Mean Meth': f"{df['meth_pct'].mean():.1f}%",
                    'Chromosomes': df['chr'].nunique()
                })

            st.dataframe(pd.DataFrame(sample_info), use_container_width=True, hide_index=True)

            if st.button("Clear All Samples", type="secondary"):
                st.session_state.methylation_samples = {}
                st.rerun()

    # =========================================================================
    # TAB 2: QC & STATS
    # =========================================================================
    with tab2:
        st.subheader("Quality Control & Statistics")

        if not st.session_state.methylation_samples:
            st.info("Load methylation data first")
        else:
            selected_sample = st.selectbox(
                "Select Sample",
                list(st.session_state.methylation_samples.keys()),
                key="qc_sample"
            )

            if selected_sample:
                df = st.session_state.methylation_samples[selected_sample]
                stats = calculate_methylation_stats(df)

                # Summary metrics
                st.markdown("### Summary Statistics")

                col1, col2, col3, col4 = st.columns(4)

                with col1:
                    st.metric("Total CpGs", f"{stats['total_cpgs']:,}")
                with col2:
                    st.metric("Mean Methylation", f"{stats['mean_methylation']:.1f}%")
                with col3:
                    st.metric("Median Methylation", f"{stats['median_methylation']:.1f}%")
                with col4:
                    if 'mean_coverage' in stats:
                        st.metric("Mean Coverage", f"{stats['mean_coverage']:.1f}x")
                    else:
                        st.metric("Std Dev", f"{stats['std_methylation']:.1f}%")

                # Coverage stats if available
                if 'cpgs_10x' in stats:
                    st.markdown("### Coverage Statistics")
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("CpGs ≥5x", f"{stats['cpgs_5x']:,}",
                                 f"{100*stats['cpgs_5x']/stats['total_cpgs']:.1f}%")
                    with col2:
                        st.metric("CpGs ≥10x", f"{stats['cpgs_10x']:,}",
                                 f"{100*stats['cpgs_10x']/stats['total_cpgs']:.1f}%")
                    with col3:
                        st.metric("Median Coverage", f"{stats['median_coverage']:.1f}x")

                # Visualizations
                st.markdown("### Distributions")

                col1, col2 = st.columns(2)

                with col1:
                    fig = plot_methylation_distribution(df, f"{selected_sample} - Methylation")
                    st.plotly_chart(fig, use_container_width=True)

                with col2:
                    fig = plot_methylation_by_level(stats, f"{selected_sample} - Categories")
                    st.plotly_chart(fig, use_container_width=True)

                if 'coverage' in df.columns:
                    col1, col2 = st.columns(2)
                    with col1:
                        fig = plot_coverage_distribution(df, f"{selected_sample} - Coverage")
                        st.plotly_chart(fig, use_container_width=True)
                    with col2:
                        fig = plot_chromosome_methylation(df, f"{selected_sample} - By Chromosome")
                        st.plotly_chart(fig, use_container_width=True)

    # =========================================================================
    # TAB 3: DMR ANALYSIS
    # =========================================================================
    with tab3:
        st.subheader("Differentially Methylated Regions (DMRs)")

        if len(st.session_state.methylation_samples) < 2:
            st.info("Load at least 2 samples for DMR analysis")
        else:
            samples = list(st.session_state.methylation_samples.keys())

            col1, col2 = st.columns(2)
            with col1:
                sample1 = st.selectbox("Sample 1 (Reference)", samples, key="dmr_s1")
            with col2:
                sample2 = st.selectbox("Sample 2 (Comparison)",
                                       [s for s in samples if s != sample1], key="dmr_s2")

            st.markdown("### Parameters")
            col1, col2, col3 = st.columns(3)

            with col1:
                min_diff = st.slider("Min Methylation Difference (%)", 10, 50, 20)
            with col2:
                min_cpgs = st.slider("Min CpGs per DMR", 1, 10, 3)
            with col3:
                window_size = st.slider("Window Size (bp)", 500, 5000, 1000, step=100)

            if st.button("Find DMRs", type="primary"):
                with st.spinner("Detecting DMRs..."):
                    df1 = st.session_state.methylation_samples[sample1]
                    df2 = st.session_state.methylation_samples[sample2]

                    dmrs = find_dmrs_simple(df1, df2, min_diff, min_cpgs, window_size)
                    st.session_state.methylation_dmrs = dmrs

                    if len(dmrs) > 0:
                        st.success(f"Found {len(dmrs)} DMRs")

                        # Record workflow step
                        if HAS_WORKFLOW_MANAGER:
                            hyper = (dmrs['direction'] == 'hyper').sum()
                            hypo = (dmrs['direction'] == 'hypo').sum()
                            WorkflowManager.record_step(
                                step_type="dmr_analysis",
                                parameters={
                                    'min_diff': min_diff,
                                    'min_cpgs': min_cpgs,
                                    'window_size': window_size,
                                    'sample1': sample1,
                                    'sample2': sample2,
                                },
                                output_metadata={
                                    'total_dmrs': len(dmrs),
                                    'hypermethylated': int(hyper),
                                    'hypomethylated': int(hypo),
                                }
                            )
                    else:
                        st.warning("No DMRs found with current parameters")

            # Display DMRs
            if st.session_state.methylation_dmrs is not None and len(st.session_state.methylation_dmrs) > 0:
                dmrs = st.session_state.methylation_dmrs

                st.markdown("### DMR Summary")

                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Total DMRs", len(dmrs))
                with col2:
                    hyper = (dmrs['direction'] == 'hyper').sum()
                    st.metric("Hypermethylated", hyper)
                with col3:
                    hypo = (dmrs['direction'] == 'hypo').sum()
                    st.metric("Hypomethylated", hypo)

                # DMR direction plot
                fig = px.histogram(dmrs, x='mean_diff', nbins=30,
                                   color='direction',
                                   color_discrete_map={'hyper': '#e74c3c', 'hypo': '#3498db'},
                                   title="DMR Methylation Differences")
                fig.update_layout(xaxis_title="Methylation Difference (%)",
                                 yaxis_title="Number of DMRs")
                st.plotly_chart(fig, use_container_width=True)

                # DMR table
                st.markdown("### DMR Table")
                st.dataframe(
                    dmrs.round(2),
                    use_container_width=True,
                    hide_index=True
                )

                # Export
                csv = dmrs.to_csv(index=False)
                st.download_button(
                    "Download DMRs (CSV)",
                    csv,
                    "dmrs.csv",
                    "text/csv"
                )

    # =========================================================================
    # TAB 4: VISUALIZATION
    # =========================================================================
    with tab4:
        st.subheader("Methylation Visualization")

        if not st.session_state.methylation_samples:
            st.info("Load methylation data first")
        else:
            viz_type = st.selectbox(
                "Visualization Type",
                ["Sample Comparison", "Chromosome View", "Methylation Heatmap"]
            )

            if viz_type == "Sample Comparison":
                if len(st.session_state.methylation_samples) >= 2:
                    samples = list(st.session_state.methylation_samples.keys())

                    # Compare distributions
                    fig = go.Figure()
                    for sample_name in samples:
                        df = st.session_state.methylation_samples[sample_name]
                        fig.add_trace(go.Histogram(
                            x=df['meth_pct'],
                            name=sample_name,
                            opacity=0.7,
                            nbinsx=50
                        ))

                    fig.update_layout(
                        title="Methylation Distribution Comparison",
                        xaxis_title="Methylation (%)",
                        yaxis_title="Count",
                        barmode='overlay'
                    )
                    st.plotly_chart(fig, use_container_width=True)

                    # Summary table
                    comparison_data = []
                    for name, df in st.session_state.methylation_samples.items():
                        stats = calculate_methylation_stats(df)
                        comparison_data.append({
                            'Sample': name,
                            'CpGs': stats['total_cpgs'],
                            'Mean': f"{stats['mean_methylation']:.1f}%",
                            'Median': f"{stats['median_methylation']:.1f}%",
                            'Hypo (<20%)': f"{100*stats['hypomethylated']/stats['total_cpgs']:.1f}%",
                            'Hyper (>80%)': f"{100*stats['hypermethylated']/stats['total_cpgs']:.1f}%",
                        })

                    st.dataframe(pd.DataFrame(comparison_data), use_container_width=True, hide_index=True)
                else:
                    st.info("Load at least 2 samples for comparison")

            elif viz_type == "Chromosome View":
                selected = st.selectbox(
                    "Select Sample",
                    list(st.session_state.methylation_samples.keys()),
                    key="chr_view_sample"
                )

                df = st.session_state.methylation_samples[selected]
                fig = plot_chromosome_methylation(df, f"{selected} - Methylation by Chromosome")
                st.plotly_chart(fig, use_container_width=True)

                # Chromosome-level stats
                chr_stats = df.groupby('chr').agg({
                    'meth_pct': ['mean', 'std', 'count']
                }).round(2)
                chr_stats.columns = ['Mean Meth (%)', 'Std Dev', 'CpG Count']
                chr_stats = chr_stats.reset_index()
                st.dataframe(chr_stats, use_container_width=True, hide_index=True)

            elif viz_type == "Methylation Heatmap":
                if len(st.session_state.methylation_samples) >= 2:
                    st.markdown("Correlation heatmap of sample methylation profiles")

                    # Compute correlations (simplified - using chromosome means)
                    samples = list(st.session_state.methylation_samples.keys())
                    chr_order = [f'chr{i}' for i in range(1, 23)]

                    matrix_data = []
                    for name in samples:
                        df = st.session_state.methylation_samples[name]
                        chr_means = df.groupby('chr')['meth_pct'].mean()
                        matrix_data.append(chr_means.reindex(chr_order).fillna(0).values)

                    corr_matrix = np.corrcoef(matrix_data)

                    fig = px.imshow(
                        corr_matrix,
                        x=samples,
                        y=samples,
                        color_continuous_scale='RdBu_r',
                        zmin=-1, zmax=1,
                        title="Sample Correlation (Chromosome-level)"
                    )
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    st.info("Load at least 2 samples for heatmap")

    # =========================================================================
    # TAB 5: INTEGRATION
    # =========================================================================
    with tab5:
        st.subheader("Epigenome Integration")
        st.markdown("Integrate DNA methylation with other epigenomic data")

        integration_type = st.selectbox(
            "Integration Type",
            [
                "Methylation × Histone Marks",
                "Methylation × Accessibility (ATAC-seq)",
                "Methylation × Gene Expression",
                "Export for External Tools"
            ]
        )

        if integration_type == "Methylation × Histone Marks":
            st.markdown("""
            ### Correlate Methylation with Histone Modifications

            Common patterns:
            - **H3K4me3** (promoters): Typically anti-correlated with DNA methylation
            - **H3K27me3** (repression): Can co-occur with DNA methylation at silenced regions
            - **H3K27ac** (enhancers): Generally anti-correlated with methylation

            **To integrate:**
            1. Load your ChIP-seq peaks in the Core Analysis section
            2. Run peak annotation
            3. Return here to correlate with methylation levels
            """)

            if st.session_state.methylation_samples and hasattr(st.session_state, 'peak_data'):
                st.success("Both methylation and peak data available!")
                if st.button("Run Integration Analysis"):
                    st.info("Integration analysis would run here...")
            else:
                missing = []
                if not st.session_state.methylation_samples:
                    missing.append("methylation data")
                if not hasattr(st.session_state, 'peak_data'):
                    missing.append("ChIP-seq peaks")
                st.warning(f"Missing: {', '.join(missing)}")

        elif integration_type == "Methylation × Accessibility (ATAC-seq)":
            st.markdown("""
            ### Correlate Methylation with Chromatin Accessibility

            DNA methylation typically **inversely correlates** with accessibility:
            - Open chromatin regions tend to be hypomethylated
            - Closed/silenced regions are often hypermethylated

            This integration helps identify:
            - Regulatory regions with coordinated epigenetic states
            - Regions with discordant signals (potential regulatory switches)
            """)

        elif integration_type == "Methylation × Gene Expression":
            st.markdown("""
            ### Correlate Methylation with Expression

            Promoter methylation is classically associated with gene silencing:
            - Hypermethylated promoters → Reduced expression
            - Hypomethylated promoters → Active transcription

            Gene body methylation can have complex relationships with expression.

            **Requirements:**
            - RNA-seq gene expression data (load in Expression Integration page)
            - Gene annotations for your genome
            """)

        elif integration_type == "Export for External Tools":
            st.markdown("### Export Data for External Analysis")

            if st.session_state.methylation_samples:
                export_sample = st.selectbox(
                    "Sample to Export",
                    list(st.session_state.methylation_samples.keys()),
                    key="export_sample"
                )

                export_format = st.selectbox(
                    "Export Format",
                    ["bedGraph", "methylKit", "BED (with score)"]
                )

                if st.button("Export"):
                    df = st.session_state.methylation_samples[export_sample]

                    if export_format == "bedGraph":
                        export_df = df[['chr', 'start', 'end', 'meth_pct']]
                        csv_data = export_df.to_csv(sep='\t', index=False, header=False)
                        filename = f"{export_sample}.bedGraph"
                    elif export_format == "methylKit":
                        export_df = df[['chr', 'start', 'meth_pct', 'count_meth', 'count_unmeth']].copy()
                        export_df.columns = ['chr', 'start', 'meth', 'coverage.meth', 'coverage.unmeth']
                        csv_data = export_df.to_csv(sep='\t', index=False)
                        filename = f"{export_sample}.methylKit.txt"
                    else:  # BED
                        export_df = df[['chr', 'start', 'end']].copy()
                        export_df['name'] = '.'
                        export_df['score'] = df['meth_pct'].astype(int)
                        csv_data = export_df.to_csv(sep='\t', index=False, header=False)
                        filename = f"{export_sample}.bed"

                    st.download_button(
                        f"Download {filename}",
                        csv_data,
                        filename,
                        "text/plain"
                    )
            else:
                st.info("Load methylation data first")

if __name__ == "__main__":
    main()
