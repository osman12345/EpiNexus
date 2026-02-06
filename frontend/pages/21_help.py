"""
Help & Support Page

Provides:
- Documentation & Tutorials
- Contact/Feedback info
- FAQ
- About information
"""

import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

st.set_page_config(page_title="Help - EpiNexus", page_icon="❓", layout="wide")


def main():
    st.title("❓ Help & Support")

    tab1, tab2, tab3, tab4 = st.tabs([
        "📖 Documentation",
        "💬 Contact & Feedback",
        "❓ FAQ",
        "ℹ️ About"
    ])

    with tab1:
        render_documentation()

    with tab2:
        render_contact()

    with tab3:
        render_faq()

    with tab4:
        render_about()


def render_documentation():
    """Documentation and tutorials."""
    st.header("📖 Documentation & Tutorials")

    # Tutorial selection in sidebar-like expanders
    section = st.selectbox(
        "Select topic",
        [
            "🚀 Getting Started",
            "📤 Data Upload & QC",
            "📊 Differential Analysis",
            "📈 Visualization",
            "⭐ Super-Enhancers",
            "🔬 Multi-omics Integration",
            "💡 Tips & Best Practices"
        ]
    )

    st.markdown("---")

    if section == "🚀 Getting Started":
        render_getting_started_docs()
    elif section == "📤 Data Upload & QC":
        render_data_qc_docs()
    elif section == "📊 Differential Analysis":
        render_differential_docs()
    elif section == "📈 Visualization":
        render_visualization_docs()
    elif section == "⭐ Super-Enhancers":
        render_se_docs()
    elif section == "🔬 Multi-omics Integration":
        render_multiomics_docs()
    elif section == "💡 Tips & Best Practices":
        render_tips_docs()


def render_getting_started_docs():
    """Getting started documentation."""
    st.subheader("🚀 Getting Started with EpiNexus")

    st.markdown("""
    ## What is EpiNexus?

    EpiNexus is an integrated platform for:
    - **Histone modification analysis** (ChIP-seq, CUT&Tag)
    - **Transcription factor binding** (motif analysis, target genes)
    - **Chromatin accessibility** (ATAC-seq)
    - **Multi-omics data integration**

    ## Quick Start Guide

    ### Step 1: Create a Project
    Go to **📁 Data & Project** → **Project** tab to create a new analysis project.

    ### Step 2: Load Your Data
    Choose your entry point:
    - **From FASTQ**: Full preprocessing pipeline (QC → Alignment → Peaks)
    - **From Peaks**: Upload existing BED/narrowPeak files directly

    ### Step 3: Configure Experiment
    In the **Samples & Config** tab:
    1. Select your **Technique** (ChIP-seq, CUT&Tag, ATAC-seq)
    2. Select your **Target** (Histones or Transcription Factors)
    3. Define your samples and conditions

    ### Step 4: Run Analysis
    Navigate to analysis modules:
    - **Differential** - Compare conditions
    - **TF ChIP-seq** - TF-specific analysis
    - **Multi-Mark** - Integrate multiple histone marks

    ## Supported Assays

    | Technique | Targets | Peak Type |
    |-----------|---------|-----------|
    | ChIP-seq | Histones, TFs | Narrow/Broad |
    | CUT&Tag | Histones, TFs | Narrow |
    | CUT&RUN | Histones, TFs | Narrow |
    | ATAC-seq | Accessibility | Narrow |

    ## Data Entry Points

    | Input | Best For |
    |-------|----------|
    | FASTQ | Full pipeline from raw data |
    | BAM | Peak calling + analysis |
    | Peaks | Direct analysis |
    """)


def render_data_qc_docs():
    """Data and QC documentation."""
    st.subheader("📤 Data Upload & Quality Control")

    st.markdown("""
    ## File Formats

    **Peak files:** BED, narrowPeak, broadPeak, MACS2 XLS, CSV/TSV
    **Alignment:** BAM (with index)
    **Raw data:** FASTQ, FASTQ.GZ
    **Signal tracks:** BigWig, bedGraph

    ## Key Quality Metrics

    | Metric | What it measures | Good value |
    |--------|------------------|------------|
    | **FRiP** | Fraction of Reads in Peaks | >0.1 (broad), >0.2 (narrow) |
    | **NSC** | Normalized Strand Cross-correlation | >1.1 |
    | **RSC** | Relative Strand Cross-correlation | >0.8 |
    | **NRF** | Non-Redundant Fraction | >0.8 |

    ## FRiP Score Interpretation
    """)

    # FRiP visualization
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=['Excellent', 'Good', 'Acceptable', 'Poor'],
        y=[0.35, 0.25, 0.15, 0.05],
        marker_color=['#27ae60', '#2ecc71', '#f39c12', '#e74c3c']
    ))
    fig.update_layout(title='FRiP Score Thresholds', yaxis_title='FRiP Score', height=300)
    st.plotly_chart(fig, use_container_width=True)

    st.markdown("""
    ## Troubleshooting Poor QC

    | Issue | Possible Cause | Solution |
    |-------|---------------|----------|
    | Low FRiP | Weak enrichment | Optimize ChIP protocol |
    | Low correlation | Batch effects | Include batch correction |
    | Few peaks | Stringent threshold | Adjust peak calling params |
    | Many peaks | Relaxed threshold | Use more stringent FDR |
    """)


def render_differential_docs():
    """Differential analysis documentation."""
    st.subheader("📊 Differential Analysis")

    st.markdown("""
    ## The DiffBind Workflow

    ```
    1. Load peak files for all samples
           ↓
    2. Create consensus peak set
           ↓
    3. Count reads in consensus peaks
           ↓
    4. Normalize counts
           ↓
    5. Statistical testing (DESeq2/edgeR)
           ↓
    6. Filter by FDR and fold change
    ```

    ## Understanding Volcano Plots
    """)

    # Demo volcano plot
    np.random.seed(42)
    n = 1000
    log2fc = np.random.normal(0, 1.5, n)
    fdr = 10 ** (-np.abs(log2fc) * np.random.uniform(0.5, 2, n))

    demo_df = pd.DataFrame({
        'log2FC': log2fc,
        'FDR': fdr,
        'Significant': ['Yes' if f < 0.05 and abs(l) > 1 else 'No' for f, l in zip(fdr, log2fc)]
    })

    fig = px.scatter(demo_df, x='log2FC', y=-np.log10(demo_df['FDR']),
                    color='Significant', color_discrete_map={'Yes': '#e74c3c', 'No': '#95a5a6'})
    fig.add_hline(y=-np.log10(0.05), line_dash='dash', line_color='gray')
    fig.add_vline(x=1, line_dash='dash', line_color='gray')
    fig.add_vline(x=-1, line_dash='dash', line_color='gray')
    fig.update_layout(height=400)
    st.plotly_chart(fig, use_container_width=True)

    st.markdown("""
    - **X-axis**: Log2 fold change (positive = gained, negative = lost)
    - **Y-axis**: Statistical significance (-log10 FDR)
    - **Thresholds**: Typically |log2FC| > 1 and FDR < 0.05

    ## Choosing Thresholds

    | Scenario | log2FC | FDR | Rationale |
    |----------|--------|-----|-----------|
    | Exploratory | 0.5 | 0.1 | Capture more candidates |
    | Standard | 1.0 | 0.05 | Balanced approach |
    | Stringent | 1.5 | 0.01 | High confidence only |
    """)


def render_visualization_docs():
    """Visualization documentation."""
    st.subheader("📈 Visualization")

    st.markdown("""
    ## Heatmaps

    Show signal intensity across many regions/samples.
    """)

    # Demo heatmap
    np.random.seed(42)
    heatmap_data = np.random.randn(20, 50)
    fig = px.imshow(heatmap_data, color_continuous_scale='RdBu_r',
                   labels={'x': 'Position', 'y': 'Genes'})
    fig.update_layout(height=250)
    st.plotly_chart(fig, use_container_width=True)

    st.markdown("""
    **Best practices:**
    - Sort rows by signal or cluster
    - Use appropriate color scale (diverging for +/-)
    - Center on relevant features (TSS, peak summit)

    ## Profile Plots

    Show average signal across a region type.
    """)

    # Demo profile
    x = np.linspace(-5, 5, 100)
    ctrl = 2 * np.exp(-0.5 * (x / 0.8) ** 2) + 0.3
    treat = 2.8 * np.exp(-0.5 * (x / 0.7) ** 2) + 0.4

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=ctrl, name='Control', line=dict(color='#3498db')))
    fig.add_trace(go.Scatter(x=x, y=treat, name='Treatment', line=dict(color='#e74c3c')))
    fig.add_vline(x=0, line_dash='dash')
    fig.update_layout(xaxis_title='Distance from TSS (kb)', yaxis_title='Signal', height=300)
    st.plotly_chart(fig, use_container_width=True)


def render_se_docs():
    """Super-enhancer documentation."""
    st.subheader("⭐ Super-Enhancers")

    st.markdown("""
    ## What are Super-Enhancers?

    Super-enhancers (SEs) are large clusters of enhancers that:
    - **Drive cell identity** - Control lineage-specific genes
    - **Show exceptional signal** - 10-100x higher than typical enhancers
    - **Are disease-relevant** - Enriched for disease variants

    ## The ROSE Algorithm

    ```
    1. Start with H3K27ac peaks (enhancer mark)
           ↓
    2. Exclude promoters (within 2.5kb of TSS)
           ↓
    3. Stitch nearby peaks (within 12.5kb)
           ↓
    4. Rank all enhancers by signal
           ↓
    5. Find inflection point (hockey stick)
           ↓
    6. Classify: above = SE, below = typical
    ```

    ## Hockey Stick Plot
    """)

    # Demo hockey stick
    np.random.seed(42)
    n = 500
    signal = np.sort(np.concatenate([
        np.random.exponential(10, int(n * 0.9)),
        np.random.exponential(80, int(n * 0.1))
    ]))
    cutoff_idx = int(n * 0.92)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=np.arange(n)[:cutoff_idx], y=signal[:cutoff_idx],
        mode='markers', marker=dict(color='#3498db', size=5), name='Typical Enhancers'))
    fig.add_trace(go.Scatter(x=np.arange(n)[cutoff_idx:], y=signal[cutoff_idx:],
        mode='markers', marker=dict(color='#e74c3c', size=8), name='Super-Enhancers'))
    fig.add_hline(y=signal[cutoff_idx], line_dash='dash')
    fig.update_layout(xaxis_title='Enhancer Rank', yaxis_title='H3K27ac Signal', height=350)
    st.plotly_chart(fig, use_container_width=True)

    st.markdown("""
    ## Key Parameters

    | Parameter | Default | Effect |
    |-----------|---------|--------|
    | Stitch distance | 12.5 kb | Higher = larger SEs |
    | TSS exclusion | 2.5 kb | Removes promoters |
    """)


def render_multiomics_docs():
    """Multi-omics documentation."""
    st.subheader("🔬 Multi-omics Integration")

    st.markdown("""
    ## Data Types

    | Data | What it shows | Example |
    |------|---------------|---------|
    | **Histone ChIP** | Chromatin state | H3K27ac = active enhancer |
    | **TF ChIP/Motif** | Regulator binding | MYC binding sites |
    | **RNA-seq** | Gene expression | Upregulated genes |

    ## Integration Strategies

    **1. Concordance Analysis** - Do histone changes match expression changes?
    """)

    # Concordance demo
    concordance = pd.DataFrame({
        'Histone Change': ['H3K27ac↑', 'H3K27ac↓', 'H3K27me3↑', 'H3K27me3↓'],
        'Expr Up': [425, 38, 25, 180],
        'Expr Down': [52, 312, 85, 42],
        'No Change': [180, 145, 102, 222]
    })

    fig = px.bar(concordance.melt(id_vars='Histone Change'),
                x='Histone Change', y='value', color='variable', barmode='group')
    fig.update_layout(height=300, yaxis_title='Number of Genes')
    st.plotly_chart(fig, use_container_width=True)

    st.markdown("""
    **2. TF-Target Networks** - Link TFs to targets through enhancers

    **3. Regulatory Modules** - Find co-regulated gene groups
    """)


def render_tips_docs():
    """Tips and best practices."""
    st.subheader("💡 Tips & Best Practices")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("""
        ## ✅ Do

        **Data Quality**
        - Always check QC first
        - Use biological replicates (2-3+)
        - Include input controls

        **Analysis**
        - Start with default settings
        - Document your parameters
        - Validate key findings (IGV, qPCR)
        """)

    with col2:
        st.markdown("""
        ## ❌ Avoid

        **Common Pitfalls**
        - P-hacking (adjusting thresholds for desired results)
        - Over-interpretation (correlation ≠ causation)
        - Ignoring batch effects
        - Too few replicates
        """)


def render_contact():
    """Contact and feedback section - simplified."""
    st.header("💬 Contact & Feedback")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("📧 Get in Touch")

        st.markdown("""
        **Report Issues / Bugs**
        - [GitHub Issues](https://github.com/your-repo/epinexus/issues)

        **Feature Requests**
        - [GitHub Discussions](https://github.com/your-repo/epinexus/discussions)

        **General Questions**
        - Email: support@epinexus.org

        **Community**
        - [GitHub Repository](https://github.com/your-repo/epinexus)
        """)

    with col2:
        st.subheader("🔗 Resources")

        st.markdown("""
        **Documentation**
        - [Online Docs](https://epinexus.readthedocs.io)
        - [API Reference](https://epinexus.readthedocs.io/api)

        **Source Code**
        - [GitHub](https://github.com/your-repo/epinexus)

        **Related Tools**
        - [DiffBind](https://bioconductor.org/packages/DiffBind)
        - [MACS2](https://github.com/macs3-project/MACS)
        - [deepTools](https://deeptools.readthedocs.io)
        """)


def render_faq():
    """Frequently asked questions."""
    st.header("❓ Frequently Asked Questions")

    faqs = [
        {
            "q": "What file formats are supported?",
            "a": """
**Peak files:** BED, narrowPeak, broadPeak, MACS2 XLS, CSV/TSV
**Alignment:** BAM (with index)
**Raw data:** FASTQ, FASTQ.GZ
**Signal tracks:** BigWig, bedGraph
**Expression:** CSV with gene_id, log2FC, padj columns
            """
        },
        {
            "q": "What's the difference between ChIP-seq and CUT&Tag?",
            "a": """
Both can profile histones and transcription factors, but:

- **ChIP-seq**: Uses crosslinking and sonication, higher background, needs more cells
- **CUT&Tag/CUT&RUN**: Enzyme-based, lower background, works with fewer cells

EpiNexus supports both with appropriate peak calling settings.
            """
        },
        {
            "q": "How do I analyze transcription factors vs histones?",
            "a": """
1. Go to **Data & Project** → **Samples & Config**
2. Select your **Technique** (ChIP-seq, CUT&Tag, etc.)
3. Select your **Target** (Histone Modifications or Transcription Factors)

For TFs, use the dedicated **TF ChIP-seq** page for motif analysis and target genes.
            """
        },
        {
            "q": "What is the demo data?",
            "a": """
Demo data is **simulated** data for demonstration purposes. It shows realistic patterns but is not from real experiments.

Look for the **📌 Demo Mode** indicator at the top of pages. Upload your own data to switch to real analysis.
            """
        },
        {
            "q": "How do I run the preprocessing pipeline?",
            "a": """
1. Go to **Data & Project** → **From FASTQ** tab
2. Select your technique and target
3. Upload or define your sample sheet
4. Configure pipeline options (aligner, peak caller)
5. Click **Run Pipeline**

**Note:** Requires external tools (Bowtie2, MACS2, etc.) installed on your system.
            """
        },
        {
            "q": "How do I save and load projects?",
            "a": """
**Save:** Go to **Data & Project** → **Project** tab → Click "Download Project File"

**Load:** Same tab → Upload the `.epinexus` file
            """
        },
        {
            "q": "What are good QC thresholds?",
            "a": """
| Metric | Good Value |
|--------|------------|
| FRiP | >0.1 (broad), >0.2 (narrow) |
| NSC | >1.1 |
| RSC | >0.8 |
| Replicate correlation | >0.9 |
            """
        },
        {
            "q": "How do I interpret a volcano plot?",
            "a": """
- **X-axis**: Log2 fold change (positive = gained, negative = lost)
- **Y-axis**: -log10 FDR (higher = more significant)
- **Thresholds**: Typically |log2FC| > 1 and FDR < 0.05
- **Red points**: Significantly changed regions
            """
        },
        {
            "q": "What are super-enhancers?",
            "a": """
Super-enhancers are large clusters of enhancers that:
- Drive cell identity genes
- Show 10-100x higher H3K27ac signal than typical enhancers
- Are often associated with disease variants

Use the **Super-Enhancers** page to identify them with the ROSE algorithm.
            """
        },
        {
            "q": "Can I use public data from ENCODE?",
            "a": """
Yes! Go to **ENCODE Data** in the sidebar to browse and download public ChIP-seq and ATAC-seq datasets.
            """
        },
    ]

    for faq in faqs:
        with st.expander(faq["q"]):
            st.markdown(faq["a"])


def render_about():
    """About EpiNexus."""
    st.header("ℹ️ About EpiNexus")

    col1, col2 = st.columns([2, 1])

    with col1:
        st.markdown("""
        ## EpiNexus

        **Comprehensive Epigenomics Analysis Platform**

        EpiNexus is an integrated platform for analyzing ChIP-seq, CUT&Tag, CUT&RUN,
        and ATAC-seq data. It provides end-to-end workflows from raw FASTQ files to
        publication-ready figures.

        ### Key Features

        - 🧬 **Multi-assay support**: ChIP-seq, CUT&Tag, CUT&RUN, ATAC-seq
        - 🔬 **Histone & TF analysis**: Specialized workflows for both
        - 📊 **Differential binding**: DiffBind integration
        - ⭐ **Super-enhancers**: ROSE algorithm
        - 🧬 **Multi-omics**: RNA-seq integration, GWAS overlap
        - 📈 **Interactive visualization**: Plotly-based charts
        - 🔬 **Genome browser**: IGV.js integration

        ### Citation

        If you use EpiNexus in your research, please cite:
        ```
        EpiNexus: A comprehensive epigenomics analysis platform
        [Your Name et al., 2026]
        ```
        """)

    with col2:
        st.markdown("""
        ### Version

        **Version:** 1.0.0

        ### Built With

        - Python 3.9+
        - Streamlit
        - Plotly
        - Pandas / NumPy

        ### Bioinformatics

        - DiffBind
        - MACS2
        - Bowtie2
        - deepTools

        ### License

        MIT License
        """)

    st.markdown("---")
    st.markdown("""
    **Acknowledgments:** ENCODE Consortium, Bioconductor, nf-core, and the bioinformatics community.
    """)


if __name__ == "__main__":
    main()
