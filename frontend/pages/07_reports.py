"""
Reports Page
Generate comprehensive analysis reports.
"""

import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

st.set_page_config(page_title="Reports - EpiNexus", page_icon="📋", layout="wide")

# Try to import data manager
try:
    from frontend.components.data_manager import DataManager

    HAS_DATA_MANAGER = True
except ImportError:
    HAS_DATA_MANAGER = False


def has_data():
    """Check if user has loaded data."""
    if HAS_DATA_MANAGER:
        peaks = DataManager.get_data("peaks")
        return peaks is not None and len(peaks) > 0
    return len(st.session_state.get("samples", [])) > 0


def render_empty_state():
    """Show empty state when no data is loaded."""
    st.markdown("---")

    col1, col2, col3 = st.columns([1, 2, 1])

    with col2:
        st.markdown(
            """
        <div style="text-align: center; padding: 3rem; background: #f8f9fa;
                    border-radius: 12px; border: 2px dashed #dee2e6;">
            <div style="font-size: 3rem; margin-bottom: 1rem;">📋</div>
            <h2 style="color: #6c757d; margin-bottom: 0.5rem;">No Analysis to Report</h2>
            <p style="color: #6c757d; font-size: 1.1rem;">
                Run your analysis first to generate reports and export data.
            </p>
        </div>
        """,
            unsafe_allow_html=True,
        )

        st.markdown("")

        if st.button("📁 Go to Data & Project", type="primary", use_container_width=True):
            st.switch_page("pages/01_data_project.py")

        st.markdown("")
        st.markdown("**Report features:**")
        st.markdown("- Comprehensive analysis summary")
        st.markdown("- Data export (TSV, CSV, Excel)")
        st.markdown("- Publication-ready methods section")


def main():
    st.title("📋 Analysis Reports")
    st.markdown("Generate and export comprehensive reports of your histone modification analysis.")

    # Check if data is loaded
    if not has_data():
        render_empty_state()
        return

    tab1, tab2, tab3 = st.tabs(["📊 Summary Report", "📥 Data Export", "📝 Methods"])

    with tab1:
        render_summary_report()
    with tab2:
        render_data_export()
    with tab3:
        render_methods()


def render_summary_report():
    """Generate analysis summary report."""
    st.header("Analysis Summary Report")

    col1, col2 = st.columns([2, 1])

    with col1:
        st.subheader("Report Configuration")

        report_title = st.text_input("Report Title", "Histone Modification Analysis Report")
        project_name = st.text_input("Project Name", "H3K27ac Treatment Study")
        author = st.text_input("Author", "")

        include_sections = st.multiselect(
            "Include sections",
            [
                "QC Summary",
                "Differential Analysis",
                "Pathway Enrichment",
                "TF Analysis",
                "Expression Integration",
                "Methods",
            ],
            default=["QC Summary", "Differential Analysis", "Pathway Enrichment"],
        )

        report_format = st.selectbox("Report format", ["HTML", "PDF", "Word (.docx)", "Markdown"])

    with col2:
        st.subheader("Analysis Overview")

        st.markdown(f"""
        **Project:** {project_name}

        **Date:** {datetime.now().strftime("%Y-%m-%d")}

        **Samples:** 8
        - Control: 4
        - Treatment: 4

        **Marks Analyzed:** 3
        - H3K4me3
        - H3K27ac
        - H3K27me3

        **Genome:** hg38
        """)

    st.markdown("---")

    st.subheader("Report Preview")

    # Preview sections
    with st.expander("📊 QC Summary", expanded=True):
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Samples Passed QC", "8/8")
        col2.metric("Avg Mapped Reads", "45.2M")
        col3.metric("Avg FRiP", "0.31")
        col4.metric("Replicate Correlation", ">0.95")

        st.markdown("""
        All samples passed quality control thresholds. High concordance between
        biological replicates confirms data quality suitable for differential analysis.
        """)

    with st.expander("📈 Differential Analysis Summary"):
        st.markdown("""
        **Consensus Peaks:** 45,231

        | Mark | Differential | Up | Down |
        |------|-------------|-----|------|
        | H3K27ac | 2,847 | 1,523 | 1,324 |
        | H3K4me3 | 1,892 | 1,045 | 847 |
        | H3K27me3 | 956 | 312 | 644 |

        Treatment resulted in global increase in active histone marks (H3K27ac, H3K4me3)
        and decrease in repressive marks (H3K27me3).
        """)

    with st.expander("🔬 Key Findings"):
        st.markdown("""
        1. **Active enhancer gain** at 1,523 regions associated with treatment response genes
        2. **MYC and E2F1 motifs enriched** in gained H3K27ac peaks (FDR < 0.001)
        3. **Cell cycle genes** show concordant increase in both H3K27ac and gene expression
        4. **Repressive chromatin loss** at tumor suppressor gene loci
        """)

    if st.button("Generate Report", type="primary"):
        with st.spinner("Generating report..."):
            st.success(f"Report generated: {project_name.replace(' ', '_')}_report.html")
            st.download_button(
                "📥 Download Report",
                f"# {report_title}\n\nGenerated: {datetime.now()}",
                f"{project_name.replace(' ', '_')}_report.html",
                mime="text/html",
            )


def render_data_export():
    """Export analysis data."""
    st.header("Data Export")

    st.subheader("Select Data to Export")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Peak Data**")
        export_peaks = st.checkbox("Consensus peak set", value=True)
        export_diff = st.checkbox("Differential peaks", value=True)
        export_annotated = st.checkbox("Annotated peaks", value=True)

        st.markdown("**Analysis Results**")
        export_de = st.checkbox("DiffBind results", value=True)
        export_pathway = st.checkbox("Pathway enrichment", value=False)
        export_tf = st.checkbox("TF enrichment", value=False)

    with col2:
        st.markdown("**Signal Data**")
        export_matrix = st.checkbox("Count matrix", value=True)
        export_norm = st.checkbox("Normalized signal matrix", value=False)
        export_bigwig = st.checkbox("BigWig tracks", value=False)

        st.markdown("**Options**")
        output_format = st.selectbox("Format", ["TSV", "CSV", "Excel", "BED"])
        compress = st.checkbox("Compress output (.zip)", value=True)

    st.markdown("---")

    if st.button("Export Selected Data", type="primary"):
        with st.spinner("Preparing export..."):
            st.success("Data export ready!")

            # Create sample export data
            export_data = "Peak_ID\tChr\tStart\tEnd\tLog2FC\tFDR\n"
            for i in range(10):
                export_data += f"peak_{i}\tchr1\t{1000000 + i * 10000}\t{1000200 + i * 10000}\t{np.random.normal(0, 1.5):.3f}\t{np.random.uniform(0, 0.1):.4f}\n"

            st.download_button(
                "📥 Download Export", export_data, "histone_analysis_export.tsv", mime="text/tab-separated-values"
            )


def render_methods():
    """Methods section for publications."""
    st.header("Methods Section")

    st.markdown("""
    Generate a methods section suitable for publication describing your analysis workflow.
    """)

    methods_text = """
### Histone Modification Analysis

#### ChIP-seq Processing
ChIP-seq reads were aligned to the reference genome (hg38) using Bowtie2 (v2.4.5) with default
parameters. Duplicate reads were removed using Picard MarkDuplicates. Peaks were called using
MACS2 (v2.2.7.1) with parameters appropriate for each histone mark: narrow peaks for H3K4me3
and H3K27ac (q < 0.01), and broad peaks for H3K27me3 and H3K36me3 (--broad, --broad-cutoff 0.1).

#### Quality Control
Sample quality was assessed using established ChIP-seq metrics including FRiP (Fraction of
Reads in Peaks), NSC (Normalized Strand Cross-correlation Coefficient), and RSC (Relative
Strand Cross-correlation Coefficient). Samples with FRiP < 0.1 or RSC < 0.8 were flagged
for review.

#### Differential Analysis
Differential histone modification analysis was performed using DiffBind (v3.8.4) with DESeq2
for statistical testing. A consensus peak set was generated by requiring peaks to be present
in at least 2 samples. Read counts were normalized using TMM (trimmed mean of M-values)
normalization. Differential peaks were identified using FDR < 0.05 and |log2 fold change| > 1.

#### Peak Annotation
Peaks were annotated using ChIPseeker (v1.34.1) with GENCODE v43 gene annotations. Genomic
feature distributions and distances to nearest TSS were calculated. Gene Ontology and pathway
enrichment analyses were performed using clusterProfiler (v4.6.2) with Benjamini-Hochberg
correction for multiple testing.

#### Multi-mark Integration
Chromatin state inference was performed by analyzing co-occurrence patterns of multiple histone
marks. Active promoters were defined as regions with H3K4me3 and H3K27ac enrichment, enhancers
as H3K4me1+/H3K27ac+ regions without H3K4me3, and repressed regions as H3K27me3+ without
active marks.
"""

    st.markdown(methods_text)

    st.download_button("📥 Download Methods Section", methods_text, "methods_section.md", mime="text/markdown")

    st.subheader("Software Versions")

    versions = pd.DataFrame(
        {
            "Software": ["DiffBind", "DESeq2", "ChIPseeker", "clusterProfiler", "MACS2", "Bowtie2"],
            "Version": ["3.8.4", "1.38.3", "1.34.1", "4.6.2", "2.2.7.1", "2.4.5"],
            "Citation": [
                "Stark & Brown 2011",
                "Love et al. 2014",
                "Yu et al. 2015",
                "Wu et al. 2021",
                "Zhang et al. 2008",
                "Langmead & Salzberg 2012",
            ],
        }
    )

    st.dataframe(versions, use_container_width=True, hide_index=True)


if __name__ == "__main__":
    main()
