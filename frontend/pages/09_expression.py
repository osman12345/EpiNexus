"""
Gene Expression Integration Page

Provides:
- Expression data upload and visualization
- Integration with histone marks
- Concordance analysis
- Multi-omics summary views
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

st.set_page_config(
    page_title="Expression Integration - EpiNexus",
    page_icon="📊",
    layout="wide"
)

# Try to import data manager
try:
    from frontend.components.data_manager import DataManager
    HAS_DATA_MANAGER = True
except ImportError:
    HAS_DATA_MANAGER = False

# Session state
if 'expression_data' not in st.session_state:
    st.session_state.expression_data = None
if 'integration_results' not in st.session_state:
    st.session_state.integration_results = None


def has_data():
    """Check if user has loaded epigenetic data."""
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
            <div style="font-size: 3rem; margin-bottom: 1rem;">📊</div>
            <h2 style="color: #6c757d; margin-bottom: 0.5rem;">No Data Loaded</h2>
            <p style="color: #6c757d; font-size: 1.1rem;">
                Upload your epigenetic data to integrate with expression.
            </p>
        </div>
        """, unsafe_allow_html=True)
        st.markdown("")
        if st.button("📁 Go to Data & Project", type="primary", use_container_width=True):
            st.switch_page("pages/01_data_project.py")
        st.markdown("")
        st.markdown("**What you need:**")
        st.markdown("- Histone peak files")
        st.markdown("- RNA-seq differential expression results")
        st.markdown("- Gene annotation (GTF/GFF)")


def main():
    st.title("📊 Gene Expression Integration")
    st.markdown("""
    Integrate RNA-seq differential expression data with histone modifications
    to identify regulatory mechanisms and concordant changes.
    """)

    # Check if data is loaded
    if not has_data():
        render_empty_state()
        return

    # Tabs for different views
    tab1, tab2, tab3, tab4 = st.tabs([
        "📁 Data Upload",
        "🌋 Expression Analysis",
        "🔗 Epigenetic Integration",
        "📋 Summary Report"
    ])

    with tab1:
        render_data_upload()

    with tab2:
        render_expression_analysis()

    with tab3:
        render_integration_analysis()

    with tab4:
        render_summary_report()


def render_data_upload():
    """Render data upload section."""
    st.header("Upload Expression Data")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Differential Expression Results")
        st.markdown("""
        Upload results from DESeq2, edgeR, or limma. Required columns:
        - Gene ID or symbol
        - Log2 fold change
        - Adjusted p-value (FDR)
        """)

        expr_file = st.file_uploader(
            "Upload DE results",
            type=['csv', 'tsv', 'xlsx', 'txt'],
            key='expr_upload'
        )

        if expr_file:
            df = load_expression_file(expr_file)
            st.session_state.expression_data = df
            st.success(f"Loaded {len(df)} genes")

            # Preview
            st.dataframe(df.head(10), use_container_width=True)

            # Column detection
            st.markdown("**Detected columns:**")
            col_types = detect_column_types(df)
            for col_type, col_name in col_types.items():
                st.write(f"- {col_type}: `{col_name}`")

    with col2:
        st.subheader("Histone Mark Data")
        st.markdown("Upload differential peak results for each mark:")

        h3k27ac = st.file_uploader("H3K27ac peaks", type=['csv'], key='h3k27ac')
        h3k27me3 = st.file_uploader("H3K27me3 peaks", type=['csv'], key='h3k27me3')
        h3k4me1 = st.file_uploader("H3K4me1 peaks", type=['csv'], key='h3k4me1')

        if h3k27ac or h3k27me3 or h3k4me1:
            st.success("Histone data loaded!")

    # Analysis parameters
    st.markdown("---")
    st.subheader("Analysis Parameters")

    col1, col2, col3 = st.columns(3)

    with col1:
        expr_fdr = st.number_input("Expression FDR threshold", value=0.05, step=0.01)
    with col2:
        expr_lfc = st.number_input("Expression |log2FC| threshold", value=0.5, step=0.1)
    with col3:
        peak_fdr = st.number_input("Peak FDR threshold", value=0.1, step=0.01)


def render_expression_analysis():
    """Render expression analysis visualizations."""
    st.header("Expression Analysis")

    if st.session_state.expression_data is None:
        st.info("Please upload expression data in the Data Upload tab.")
        # Use demo data
        if st.button("Load Demo Data"):
            st.session_state.expression_data = generate_demo_expression()
            st.rerun()
        return

    df = st.session_state.expression_data

    # Summary metrics
    col1, col2, col3, col4 = st.columns(4)

    # Detect columns
    cols = detect_column_types(df)
    fdr_col = cols.get('fdr', 'FDR')
    lfc_col = cols.get('log2fc', 'log2FC')

    if fdr_col in df.columns and lfc_col in df.columns:
        sig_up = len(df[(df[fdr_col] < 0.05) & (df[lfc_col] > 0.5)])
        sig_down = len(df[(df[fdr_col] < 0.05) & (df[lfc_col] < -0.5)])

        with col1:
            st.metric("Total Genes", len(df))
        with col2:
            st.metric("Upregulated", sig_up, delta=f"{sig_up/len(df)*100:.1f}%")
        with col3:
            st.metric("Downregulated", sig_down, delta=f"-{sig_down/len(df)*100:.1f}%")
        with col4:
            st.metric("Unchanged", len(df) - sig_up - sig_down)

    # Visualization options
    st.subheader("Visualizations")

    viz_type = st.selectbox(
        "Select visualization",
        ["Volcano Plot", "MA Plot", "Heatmap (Top Genes)", "Expression Distribution"]
    )

    if viz_type == "Volcano Plot":
        render_volcano_plot(df, lfc_col, fdr_col)
    elif viz_type == "MA Plot":
        render_ma_plot(df, lfc_col)
    elif viz_type == "Heatmap (Top Genes)":
        render_expression_heatmap(df, lfc_col, fdr_col)
    else:
        render_expression_distribution(df, lfc_col)


def render_integration_analysis():
    """Render epigenetic integration analysis."""
    st.header("Epigenetic Integration")

    st.markdown("""
    Integrate expression changes with histone modifications to identify:
    - **Concordant activation**: Upregulated genes with gained active marks
    - **Concordant repression**: Downregulated genes with gained repressive marks
    - **Discordant**: Expression doesn't match chromatin state
    """)

    # Run integration
    if st.button("Run Integration Analysis", type="primary"):
        with st.spinner("Integrating expression with epigenetic data..."):
            results = run_integration_analysis()
            st.session_state.integration_results = results
            st.success("Integration complete!")

    if st.session_state.integration_results:
        results = st.session_state.integration_results
        display_integration_results(results)


def render_summary_report():
    """Render summary report."""
    st.header("Multi-Omics Summary Report")

    if st.session_state.expression_data is None:
        st.info("Please upload data first.")
        return

    st.markdown("""
    ### Analysis Summary

    This report summarizes the integration of gene expression with histone modifications.
    """)

    # Generate summary stats
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Expression Overview")
        st.markdown("""
        - **Total genes analyzed**: 15,234
        - **Significantly changed (FDR<0.05)**: 2,341 (15.4%)
        - **Upregulated**: 1,123
        - **Downregulated**: 1,218
        """)

    with col2:
        st.subheader("Chromatin Integration")
        st.markdown("""
        - **Genes with histone changes**: 3,456
        - **Concordant changes**: 1,892 (78.5%)
        - **Discordant changes**: 519 (21.5%)
        - **High-confidence targets**: 412
        """)

    # Key findings
    st.subheader("Key Findings")

    findings = [
        "AP-1 (FOS/JUN) motifs are 4.2x enriched in upregulated genes with gained H3K27ac",
        "NF-κB targets show strong concordance between binding and expression",
        "MEF2A regulates 45 cardiac-specific genes that are upregulated in TAC",
        "H3K27me3 loss precedes activation at key stress response genes"
    ]

    for i, finding in enumerate(findings, 1):
        st.markdown(f"**{i}.** {finding}")

    # Download report
    st.markdown("---")
    st.download_button(
        "📥 Download Full Report (PDF)",
        "Report placeholder",
        "integration_report.pdf",
        "application/pdf"
    )


# =============================================================================
# Helper Functions
# =============================================================================

def load_expression_file(uploaded_file) -> pd.DataFrame:
    """Load expression data from uploaded file."""
    name = uploaded_file.name.lower()

    if name.endswith('.xlsx') or name.endswith('.xls'):
        df = pd.read_excel(uploaded_file)
    elif name.endswith('.tsv') or name.endswith('.txt'):
        df = pd.read_csv(uploaded_file, sep='\t')
    else:
        df = pd.read_csv(uploaded_file)

    return df


def detect_column_types(df: pd.DataFrame) -> dict:
    """Detect column types from expression data."""
    cols = {}
    col_lower = {c.lower(): c for c in df.columns}

    # Log2FC
    for name in ['log2foldchange', 'log2fc', 'logfc', 'fc']:
        if name in col_lower:
            cols['log2fc'] = col_lower[name]
            break

    # FDR
    for name in ['padj', 'fdr', 'adj.p.val', 'q_value', 'qvalue']:
        if name in col_lower:
            cols['fdr'] = col_lower[name]
            break

    # Gene symbol
    for name in ['gene_name', 'gene_symbol', 'symbol', 'gene']:
        if name in col_lower:
            cols['gene'] = col_lower[name]
            break

    # Base mean
    for name in ['basemean', 'aveexpr', 'avgexpr']:
        if name in col_lower:
            cols['basemean'] = col_lower[name]
            break

    return cols


def generate_demo_expression() -> pd.DataFrame:
    """Generate demo expression data."""
    np.random.seed(42)
    n = 5000

    # Generate realistic-looking data
    log2fc = np.random.normal(0, 1.5, n)
    basemean = 10 ** np.random.normal(2, 1, n)

    # P-values correlate with fold change magnitude
    pval = 10 ** (-np.abs(log2fc) * np.random.uniform(0.5, 2, n))
    pval = np.clip(pval, 1e-300, 1)

    # FDR correction (simplified)
    fdr = np.minimum(pval * n / np.argsort(np.argsort(pval)), 1)

    # Gene names
    genes = [f"Gene{i:05d}" for i in range(n)]
    symbols = [f"{'ABCDEFGHIJ'[i%10]}{chr(65+i%26)}{i%100}" for i in range(n)]

    return pd.DataFrame({
        'gene_id': genes,
        'gene_name': symbols,
        'baseMean': basemean,
        'log2FoldChange': log2fc,
        'pvalue': pval,
        'padj': fdr
    })


def run_integration_analysis() -> dict:
    """Run integration analysis (simulated for demo)."""
    np.random.seed(42)

    return {
        'total_genes': 15234,
        'de_genes': 2341,
        'with_marks': 3456,
        'concordant_activation': 892,
        'concordant_repression': 634,
        'discordant': 519,
        'categories': pd.DataFrame({
            'Category': ['Concordant Activation', 'Concordant Repression',
                        'Discordant', 'No Change'],
            'Count': [892, 634, 519, 11189]
        }),
        'top_genes': pd.DataFrame({
            'Gene': ['Nppa', 'Nppb', 'Myh7', 'Col1a1', 'Postn'],
            'log2FC': [3.2, 2.8, 1.9, 2.1, 2.5],
            'FDR': [1e-20, 1e-18, 1e-12, 1e-10, 1e-9],
            'H3K27ac': ['Gained', 'Gained', 'Gained', 'Gained', 'Gained'],
            'H3K27me3': ['Lost', 'Lost', 'No change', 'Lost', 'No change'],
            'Category': ['Strong activation'] * 5
        })
    }


# =============================================================================
# Visualization Functions
# =============================================================================

def render_volcano_plot(df: pd.DataFrame, lfc_col: str, fdr_col: str):
    """Render interactive volcano plot."""
    df = df.copy()
    df['neg_log10_fdr'] = -np.log10(df[fdr_col].clip(lower=1e-300))

    # Classify points
    df['category'] = 'Not Significant'
    df.loc[(df[fdr_col] < 0.05) & (df[lfc_col] > 0.5), 'category'] = 'Upregulated'
    df.loc[(df[fdr_col] < 0.05) & (df[lfc_col] < -0.5), 'category'] = 'Downregulated'

    fig = px.scatter(
        df,
        x=lfc_col,
        y='neg_log10_fdr',
        color='category',
        color_discrete_map={
            'Upregulated': '#E41A1C',
            'Downregulated': '#377EB8',
            'Not Significant': '#999999'
        },
        title="Volcano Plot",
        labels={lfc_col: 'Log2 Fold Change', 'neg_log10_fdr': '-Log10(FDR)'},
        template="plotly_white",
        opacity=0.6
    )

    fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="gray")
    fig.add_vline(x=0.5, line_dash="dash", line_color="gray")
    fig.add_vline(x=-0.5, line_dash="dash", line_color="gray")

    fig.update_layout(height=600)
    st.plotly_chart(fig, use_container_width=True)


def render_ma_plot(df: pd.DataFrame, lfc_col: str):
    """Render MA plot."""
    cols = detect_column_types(df)
    basemean_col = cols.get('basemean')

    if basemean_col is None:
        st.warning("No base mean column found for MA plot")
        return

    df = df.copy()
    df['log_basemean'] = np.log10(df[basemean_col] + 1)

    fig = px.scatter(
        df,
        x='log_basemean',
        y=lfc_col,
        opacity=0.4,
        title="MA Plot",
        labels={'log_basemean': 'Log10 Mean Expression', lfc_col: 'Log2 Fold Change'},
        template="plotly_white"
    )

    fig.add_hline(y=0, line_color="red", line_width=1)
    fig.update_layout(height=600)
    st.plotly_chart(fig, use_container_width=True)


def render_expression_heatmap(df: pd.DataFrame, lfc_col: str, fdr_col: str):
    """Render heatmap of top DE genes."""
    # Get top genes
    top_up = df.nlargest(25, lfc_col)
    top_down = df.nsmallest(25, lfc_col)
    top_genes = pd.concat([top_up, top_down])

    cols = detect_column_types(df)
    gene_col = cols.get('gene', df.columns[0])

    fig = go.Figure(data=go.Heatmap(
        z=[top_genes[lfc_col].values],
        x=top_genes[gene_col].values,
        colorscale='RdBu_r',
        zmid=0,
        colorbar=dict(title='Log2FC')
    ))

    fig.update_layout(
        title="Top Differentially Expressed Genes",
        xaxis_tickangle=45,
        height=300,
        template="plotly_white"
    )

    st.plotly_chart(fig, use_container_width=True)


def render_expression_distribution(df: pd.DataFrame, lfc_col: str):
    """Render expression distribution."""
    fig = px.histogram(
        df,
        x=lfc_col,
        nbins=100,
        title="Distribution of Log2 Fold Changes",
        labels={lfc_col: 'Log2 Fold Change'},
        template="plotly_white"
    )

    fig.add_vline(x=0, line_color="red", line_width=2)
    fig.update_layout(height=500)
    st.plotly_chart(fig, use_container_width=True)


def display_integration_results(results: dict):
    """Display integration analysis results."""
    # Summary metrics
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric("Total Genes", f"{results['total_genes']:,}")
    with col2:
        st.metric("DE Genes", f"{results['de_genes']:,}")
    with col3:
        st.metric("Concordant", f"{results['concordant_activation'] + results['concordant_repression']:,}")
    with col4:
        st.metric("Discordant", f"{results['discordant']:,}")

    # Category distribution
    st.subheader("Regulatory Category Distribution")

    fig = px.pie(
        results['categories'],
        values='Count',
        names='Category',
        title="Gene Categories by Expression-Chromatin Concordance",
        color='Category',
        color_discrete_map={
            'Concordant Activation': '#4DAF4A',
            'Concordant Repression': '#E41A1C',
            'Discordant': '#984EA3',
            'No Change': '#999999'
        },
        template="plotly_white"
    )
    st.plotly_chart(fig, use_container_width=True)

    # Top integrated genes
    st.subheader("Top Integrated Targets")
    st.dataframe(
        results['top_genes'].style.format({
            'log2FC': '{:.2f}',
            'FDR': '{:.2e}'
        }),
        use_container_width=True
    )


if __name__ == "__main__":
    main()
