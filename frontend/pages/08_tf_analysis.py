"""
Transcription Factor Analysis Page

Provides:
- TF motif enrichment analysis
- TF binding site visualization
- TF-target gene prediction
- Integration with histone marks and expression
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
import sys

# Add parent paths
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

st.set_page_config(
    page_title="TF Analysis - EpiNexus",
    page_icon="🧬",
    layout="wide"
)

# Try to import data manager
try:
    from frontend.components.data_manager import DataManager
    HAS_DATA_MANAGER = True
except ImportError:
    HAS_DATA_MANAGER = False

# Initialize session state
if 'tf_results' not in st.session_state:
    st.session_state.tf_results = None
if 'enrichment_results' not in st.session_state:
    st.session_state.enrichment_results = None


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
            <div style="font-size: 3rem; margin-bottom: 1rem;">🔬</div>
            <h2 style="color: #6c757d; margin-bottom: 0.5rem;">No Data Loaded</h2>
            <p style="color: #6c757d; font-size: 1.1rem;">
                Upload your peak files to run TF motif analysis.
            </p>
        </div>
        """, unsafe_allow_html=True)
        st.markdown("")
        if st.button("📁 Go to Data & Project", type="primary", use_container_width=True):
            st.switch_page("pages/01_data_project.py")
        st.markdown("")
        st.markdown("**TF Analysis features:**")
        st.markdown("- Motif enrichment (HOMER, MEME)")
        st.markdown("- TF binding site visualization")
        st.markdown("- TF-target gene prediction")


def main():
    st.title("🔬 Transcription Factor Analysis")
    st.markdown("""
    Analyze transcription factor binding motifs and their relationship with
    histone modifications and gene expression.
    """)

    # Check if data is loaded
    if not has_data():
        render_empty_state()
        return

    # Sidebar configuration
    with st.sidebar:
        st.header("Analysis Settings")

        analysis_type = st.selectbox(
            "Analysis Type",
            [
                "Motif Enrichment",
                "TF Binding Site Scan",
                "TF-Target Prediction",
                "TF + Histone Integration"
            ]
        )

        st.markdown("---")

        st.subheader("Thresholds")
        motif_pval = st.slider("Motif P-value", 1e-6, 0.01, 1e-4, format="%.1e")
        enrichment_fdr = st.slider("Enrichment FDR", 0.01, 0.2, 0.05)

    # Main content based on analysis type
    if analysis_type == "Motif Enrichment":
        render_motif_enrichment()
    elif analysis_type == "TF Binding Site Scan":
        render_binding_site_scan()
    elif analysis_type == "TF-Target Prediction":
        render_target_prediction()
    else:
        render_tf_histone_integration()


def render_motif_enrichment():
    """Render motif enrichment analysis section."""
    st.header("TF Motif Enrichment Analysis")

    st.markdown("""
    Find transcription factor motifs enriched in your peak set compared to background.
    This helps identify TFs that may regulate genes near your differential peaks.
    """)

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Target Peaks")
        target_file = st.file_uploader(
            "Upload target peaks (BED/CSV)",
            type=['bed', 'csv', 'tsv'],
            key='target_peaks',
            help="Peaks where you want to find enriched motifs (e.g., gained H3K27ac peaks)"
        )

        if target_file:
            target_df = load_peak_file(target_file)
            st.success(f"Loaded {len(target_df)} target peaks")
            st.dataframe(target_df.head(), use_container_width=True)

    with col2:
        st.subheader("Background Peaks")
        bg_option = st.radio(
            "Background selection",
            ["All peaks", "Upload custom", "Genome-wide (shuffled)"]
        )

        if bg_option == "Upload custom":
            bg_file = st.file_uploader(
                "Upload background peaks",
                type=['bed', 'csv', 'tsv'],
                key='bg_peaks'
            )

    st.markdown("---")

    # Reference genome for sequence extraction
    st.subheader("Reference Genome")
    col1, col2 = st.columns(2)

    with col1:
        genome = st.selectbox(
            "Select genome",
            ["mm10", "mm39", "hg38", "hg19", "rn6"]
        )

    with col2:
        genome_fasta = st.text_input(
            "Path to genome FASTA (optional)",
            placeholder="/path/to/genome.fa"
        )

    # Motif database selection
    st.subheader("Motif Database")
    motif_source = st.radio(
        "Motif source",
        ["Built-in JASPAR (common TFs)", "Upload custom motifs", "Full JASPAR database"]
    )

    if motif_source == "Upload custom":
        motif_file = st.file_uploader(
            "Upload motifs (JASPAR format)",
            type=['txt', 'jaspar', 'meme']
        )

    # Run analysis
    st.markdown("---")
    if st.button("Run Motif Enrichment", type="primary"):
        with st.spinner("Analyzing motif enrichment..."):
            # Simulated results for demo
            results = simulate_enrichment_results()
            st.session_state.enrichment_results = results
            st.success("Analysis complete!")

    # Display results
    if st.session_state.enrichment_results is not None:
        display_enrichment_results(st.session_state.enrichment_results)


def render_binding_site_scan():
    """Render TF binding site scanning section."""
    st.header("TF Binding Site Scanning")

    st.markdown("""
    Scan your peaks for specific TF binding motifs and visualize their locations.
    """)

    col1, col2 = st.columns([2, 1])

    with col1:
        peaks_file = st.file_uploader(
            "Upload peaks to scan",
            type=['bed', 'csv'],
            key='scan_peaks'
        )

    with col2:
        # TF selection
        st.subheader("Select TFs")
        tf_options = [
            "AP-1 (FOS/JUN)", "NF-κB", "GATA1/2", "ETS1",
            "MEF2A", "SP1", "CREB1", "TP53", "All available"
        ]
        selected_tfs = st.multiselect("TFs to scan", tf_options, default=["AP-1 (FOS/JUN)"])

    if st.button("Scan for Binding Sites", type="primary"):
        with st.spinner("Scanning..."):
            # Demo results
            scan_results = simulate_scan_results()
            display_scan_results(scan_results)


def render_target_prediction():
    """Render TF-target gene prediction section."""
    st.header("TF-Target Gene Prediction")

    st.markdown("""
    Predict target genes regulated by transcription factors based on:
    - TF binding site locations
    - Distance to transcription start sites
    - Histone mark context
    - Gene expression changes
    """)

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("TF Binding Data")
        tf_file = st.file_uploader(
            "Upload TF peaks (ChIP-seq/CUT&Tag)",
            type=['bed', 'csv', 'narrowPeak'],
            key='tf_peaks'
        )

        tf_name = st.text_input("TF Name", placeholder="e.g., NFKB1")

    with col2:
        st.subheader("Expression Data (Optional)")
        expr_file = st.file_uploader(
            "Upload DE results",
            type=['csv', 'tsv', 'xlsx'],
            key='expr_data'
        )

    st.subheader("Prediction Parameters")
    col1, col2, col3 = st.columns(3)

    with col1:
        max_distance = st.number_input(
            "Max distance from TSS (bp)",
            value=100000,
            step=10000
        )

    with col2:
        promoter_dist = st.number_input(
            "Promoter region (bp)",
            value=3000,
            step=500
        )

    with col3:
        min_confidence = st.selectbox(
            "Minimum confidence",
            ["All", "Low", "Medium", "High"]
        )

    if st.button("Predict Targets", type="primary"):
        with st.spinner("Predicting TF targets..."):
            targets = simulate_target_predictions()
            display_target_predictions(targets)


def render_tf_histone_integration():
    """Render TF + histone mark integration section."""
    st.header("TF + Histone Mark Integration")

    st.markdown("""
    Integrate TF binding with histone modification data to understand:
    - TF binding at active vs repressed regions
    - Co-occupancy patterns
    - Regulatory logic
    """)

    # File uploads
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("TF Data")
        tf_file = st.file_uploader(
            "TF binding peaks",
            type=['bed', 'csv'],
            key='tf_int'
        )

    with col2:
        st.subheader("Histone Data")
        h3k27ac_file = st.file_uploader("H3K27ac differential peaks", type=['csv'])
        h3k27me3_file = st.file_uploader("H3K27me3 differential peaks", type=['csv'])
        h3k4me1_file = st.file_uploader("H3K4me1 differential peaks", type=['csv'])

    st.subheader("Expression Data (Optional)")
    expr_file = st.file_uploader(
        "Gene expression DE results",
        type=['csv', 'tsv'],
        key='expr_int'
    )

    if st.button("Run Integration Analysis", type="primary"):
        with st.spinner("Integrating data..."):
            results = simulate_integration_results()
            display_integration_results(results)


# =============================================================================
# Helper Functions
# =============================================================================

def load_peak_file(uploaded_file) -> pd.DataFrame:
    """Load peak file from upload."""
    if uploaded_file.name.endswith('.bed'):
        df = pd.read_csv(uploaded_file, sep='\t', header=None)
        df.columns = ['chrom', 'start', 'end'] + [f'col_{i}' for i in range(3, len(df.columns))]
    else:
        df = pd.read_csv(uploaded_file)
    return df


def simulate_enrichment_results() -> pd.DataFrame:
    """Generate simulated enrichment results for demo."""
    tfs = ['AP-1 (FOS)', 'NF-κB', 'MEF2A', 'GATA2', 'ETS1', 'SP1', 'CREB1', 'TP53', 'KLF4', 'STAT3']
    np.random.seed(42)

    data = {
        'TF': tfs,
        'Motif_ID': [f'MA{i:04d}.1' for i in range(len(tfs))],
        'Target_Count': np.random.randint(50, 500, len(tfs)),
        'Target_Percent': np.random.uniform(5, 40, len(tfs)),
        'Background_Percent': np.random.uniform(1, 15, len(tfs)),
        'Fold_Enrichment': np.random.uniform(1.5, 8, len(tfs)),
        'P_value': 10 ** np.random.uniform(-10, -2, len(tfs)),
        'FDR': 10 ** np.random.uniform(-8, -1, len(tfs))
    }

    return pd.DataFrame(data).sort_values('FDR')


def simulate_scan_results() -> pd.DataFrame:
    """Generate simulated scan results."""
    np.random.seed(42)
    n = 100

    data = {
        'peak_id': [f'peak_{i:05d}' for i in range(n)],
        'chrom': np.random.choice(['chr1', 'chr2', 'chr3', 'chr5', 'chr7'], n),
        'start': np.random.randint(1000000, 50000000, n),
        'tf_motifs': [';'.join(np.random.choice(['AP-1', 'NF-κB', 'GATA', 'ETS'],
                     np.random.randint(0, 4))) for _ in range(n)],
        'n_motifs': np.random.randint(0, 5, n),
        'best_score': np.random.uniform(0.7, 1.0, n)
    }

    return pd.DataFrame(data)


def simulate_target_predictions() -> pd.DataFrame:
    """Generate simulated target predictions."""
    genes = ['Nppa', 'Nppb', 'Myh7', 'Myh6', 'Acta1', 'Tnni3', 'Ryr2', 'Atp2a2',
             'Pln', 'Camk2d', 'Col1a1', 'Col3a1', 'Fn1', 'Postn', 'Tgfb1']
    np.random.seed(42)

    data = {
        'Gene': genes,
        'Distance_to_TSS': np.random.randint(-50000, 50000, len(genes)),
        'Binding_Score': np.random.uniform(5, 20, len(genes)),
        'Expression_log2FC': np.random.uniform(-3, 3, len(genes)),
        'Expression_FDR': 10 ** np.random.uniform(-6, -1, len(genes)),
        'Confidence': np.random.choice(['High', 'Medium', 'Low'], len(genes), p=[0.3, 0.5, 0.2]),
        'Histone_Context': np.random.choice(['Active_enhancer', 'Active_promoter', 'Poised', 'Repressed'], len(genes))
    }

    return pd.DataFrame(data).sort_values('Expression_FDR')


def simulate_integration_results() -> dict:
    """Generate simulated integration results."""
    return {
        'total_genes': 5234,
        'tf_bound_genes': 1523,
        'de_genes': 892,
        'integrated': 412,
        'concordant_activation': 234,
        'concordant_repression': 98,
        'discordant': 80,
        'cooccupancy': pd.DataFrame({
            'Histone_Mark': ['H3K27ac', 'H3K27me3', 'H3K4me1'],
            'TF_Overlap': [456, 123, 389],
            'TF_Overlap_Percent': [29.9, 8.1, 25.5]
        })
    }


# =============================================================================
# Display Functions
# =============================================================================

def display_enrichment_results(results: pd.DataFrame):
    """Display motif enrichment results."""
    st.subheader("Enrichment Results")

    # Summary metrics
    col1, col2, col3, col4 = st.columns(4)
    sig_results = results[results['FDR'] < 0.05]

    with col1:
        st.metric("TFs Tested", len(results))
    with col2:
        st.metric("Significant (FDR<0.05)", len(sig_results))
    with col3:
        st.metric("Top Enrichment", f"{results['Fold_Enrichment'].max():.1f}x")
    with col4:
        st.metric("Best FDR", f"{results['FDR'].min():.2e}")

    # Enrichment plot
    st.subheader("Motif Enrichment Plot")

    fig = go.Figure()

    # Color by significance
    colors = ['#E41A1C' if fdr < 0.05 else '#999999' for fdr in results['FDR']]

    fig.add_trace(go.Bar(
        x=results['TF'],
        y=results['Fold_Enrichment'],
        marker_color=colors,
        text=[f"FDR: {fdr:.2e}" for fdr in results['FDR']],
        hovertemplate="<b>%{x}</b><br>Fold Enrichment: %{y:.2f}<br>%{text}<extra></extra>"
    ))

    fig.update_layout(
        title="TF Motif Enrichment in Target Peaks",
        xaxis_title="Transcription Factor",
        yaxis_title="Fold Enrichment",
        template="plotly_white",
        height=500
    )

    fig.add_hline(y=1, line_dash="dash", line_color="gray")

    st.plotly_chart(fig, use_container_width=True)

    # Results table
    st.subheader("Detailed Results")
    st.dataframe(
        results.style.format({
            'Target_Percent': '{:.1f}%',
            'Background_Percent': '{:.1f}%',
            'Fold_Enrichment': '{:.2f}',
            'P_value': '{:.2e}',
            'FDR': '{:.2e}'
        }),
        use_container_width=True
    )

    # Download
    csv = results.to_csv(index=False)
    st.download_button(
        "📥 Download Results (CSV)",
        csv,
        "motif_enrichment_results.csv",
        "text/csv"
    )


def display_scan_results(results: pd.DataFrame):
    """Display binding site scan results."""
    st.subheader("Binding Site Scan Results")

    # Summary
    peaks_with_motifs = len(results[results['n_motifs'] > 0])

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Peaks Scanned", len(results))
    with col2:
        st.metric("Peaks with Motifs", peaks_with_motifs)
    with col3:
        st.metric("Coverage", f"{peaks_with_motifs/len(results)*100:.1f}%")

    # Distribution plot
    fig = px.histogram(
        results,
        x='n_motifs',
        title="Distribution of Motifs per Peak",
        labels={'n_motifs': 'Number of TF Motifs'},
        template="plotly_white"
    )
    st.plotly_chart(fig, use_container_width=True)

    # Table
    st.dataframe(results.head(50), use_container_width=True)


def display_target_predictions(results: pd.DataFrame):
    """Display TF target predictions."""
    st.subheader("Predicted TF Targets")

    # Summary
    high_conf = len(results[results['Confidence'] == 'High'])

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total Predictions", len(results))
    with col2:
        st.metric("High Confidence", high_conf)
    with col3:
        st.metric("DE Genes Targeted", len(results[results['Expression_FDR'] < 0.05]))

    # Scatter plot: distance vs expression
    fig = px.scatter(
        results,
        x='Distance_to_TSS',
        y='Expression_log2FC',
        color='Confidence',
        size='Binding_Score',
        hover_data=['Gene', 'Histone_Context'],
        title="TF Targets: Distance to TSS vs Expression Change",
        labels={
            'Distance_to_TSS': 'Distance to TSS (bp)',
            'Expression_log2FC': 'Expression Log2 Fold Change'
        },
        color_discrete_map={'High': '#E41A1C', 'Medium': '#FF7F00', 'Low': '#999999'},
        template="plotly_white"
    )

    fig.add_hline(y=0, line_dash="dash", line_color="gray")
    fig.add_vline(x=0, line_dash="dash", line_color="gray")

    st.plotly_chart(fig, use_container_width=True)

    # Table
    st.dataframe(
        results.style.format({
            'Distance_to_TSS': '{:,}',
            'Binding_Score': '{:.1f}',
            'Expression_log2FC': '{:.2f}',
            'Expression_FDR': '{:.2e}'
        }),
        use_container_width=True
    )


def display_integration_results(results: dict):
    """Display TF + histone integration results."""
    st.subheader("Integration Summary")

    # Metrics
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Total Genes", results['total_genes'])
    with col2:
        st.metric("TF-Bound Genes", results['tf_bound_genes'])
    with col3:
        st.metric("DE Genes", results['de_genes'])
    with col4:
        st.metric("Integrated Targets", results['integrated'])

    # Category breakdown
    st.subheader("Regulatory Categories")

    categories = {
        'Concordant Activation': results['concordant_activation'],
        'Concordant Repression': results['concordant_repression'],
        'Discordant': results['discordant']
    }

    fig = px.pie(
        values=list(categories.values()),
        names=list(categories.keys()),
        title="Distribution of Regulatory Categories",
        color_discrete_sequence=['#4DAF4A', '#E41A1C', '#984EA3'],
        template="plotly_white"
    )
    st.plotly_chart(fig, use_container_width=True)

    # Co-occupancy
    st.subheader("TF + Histone Co-occupancy")
    cooc = results['cooccupancy']

    fig = px.bar(
        cooc,
        x='Histone_Mark',
        y='TF_Overlap_Percent',
        text='TF_Overlap',
        title="TF Binding Overlap with Histone Marks",
        labels={'TF_Overlap_Percent': '% TF Peaks Overlapping'},
        template="plotly_white"
    )
    fig.update_traces(textposition='outside')
    st.plotly_chart(fig, use_container_width=True)


if __name__ == "__main__":
    main()
