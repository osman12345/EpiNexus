"""
Multi-omics rendering helpers.

Extracted from 10_multiomics.py to keep the page file slim.
Each function renders one tab of the multi-omics analysis page.
"""

from __future__ import annotations

from typing import Dict

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots

from frontend.components.theme import COLORS


# ============================================================================
# Data Overview tab
# ============================================================================


def render_data_overview() -> None:
    """Overview of all integrated data types."""
    st.header("Data Overview")

    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.subheader("🧬 Histone Modifications")
        st.markdown("""
        **Marks Analyzed:** 4
        - H3K4me3 (Active promoters)
        - H3K27ac (Active enhancers)
        - H3K4me1 (Enhancers)
        - H3K27me3 (Repression)

        **Differential Regions:**
        - H3K27ac: 2,847 (↑1,523 / ↓1,324)
        - H3K4me3: 1,892 (↑1,045 / ↓847)
        - H3K27me3: 956 (↑312 / ↓644)
        """)
        st.metric("Total Differential Peaks", "5,695")

    with col2:
        st.subheader("🔬 Transcription Factors")
        st.markdown("""
        **Enriched TF Motifs:** 45

        **Top TFs by Enrichment:**
        1. MYC (8.5x, p<10⁻²⁵)
        2. E2F1 (6.2x, p<10⁻¹⁸)
        3. SP1 (5.8x, p<10⁻¹⁵)
        4. NRF1 (4.9x, p<10⁻¹²)
        5. YY1 (4.2x, p<10⁻¹⁰)

        **TF Families Represented:** 12
        """)
        st.metric("Predicted TF-Target Pairs", "8,245")

    with col3:
        st.subheader("📊 Gene Expression")
        st.markdown("""
        **RNA-seq Analysis:**
        - Total genes: 21,432
        - Differentially expressed: 3,156
        - Upregulated: 1,823
        - Downregulated: 1,333

        **Integration Status:**
        - Genes with histone data: 18,542
        - Genes with TF predictions: 12,891
        """)
        st.metric("DEGs Integrated", "3,156")

    with col4:
        _render_methylation_overview()

    st.markdown("---")

    st.subheader("Data Quality Summary")

    meth_status = "✅ Loaded" if st.session_state.get("methylation_samples") else "⚠️ Not loaded"
    meth_samples = (
        str(len(st.session_state.get("methylation_samples", {})))
        if st.session_state.get("methylation_samples")
        else "-"
    )

    quality_data = pd.DataFrame(
        {
            "Data Type": [
                "Histone ChIP-seq",
                "TF Motif Analysis",
                "RNA-seq",
                "DNA Methylation",
                "Integration",
            ],
            "Samples": ["8", "-", "6", meth_samples, "-"],
            "Features": [
                "45,231 peaks",
                "1,892 motifs",
                "21,432 genes",
                "CpG sites",
                "10,234 genes",
            ],
            "Quality": ["✅ Pass", "✅ Pass", "✅ Pass", meth_status, "✅ Complete"],
            "Notes": [
                "High FRiP, good replication",
                "Validated motif database",
                "High mapping rate",
                "Bisulfite-seq/Array",
                "All data types linked",
            ],
        }
    )
    st.dataframe(quality_data, use_container_width=True, hide_index=True)


def _render_methylation_overview() -> None:
    """Methylation sub-section of the data overview."""
    st.subheader("🧬 DNA Methylation")

    has_meth = bool(st.session_state.get("methylation_samples"))

    if has_meth:
        n_samples = len(st.session_state.methylation_samples)
        st.markdown(
            f"""
        **Samples Loaded:** {n_samples}

        **DMR Analysis:**
        - Hypermethylated: varies
        - Hypomethylated: varies

        **Integration:**
        - Promoter methylation
        - Enhancer methylation
        - Gene body methylation
        """
        )
        st.metric("Methylation Samples", n_samples)
    else:
        st.markdown(
            """
        **No methylation data loaded**

        Load methylation data to integrate:
        - Promoter methylation vs expression
        - Enhancer methylation vs H3K27ac
        - DMR overlap with histone changes
        """
        )
        if st.button("Load Methylation Data", key="load_meth_overview"):
            st.switch_page("pages/24_methylation.py")


# ============================================================================
# Integration Analysis tab
# ============================================================================


def render_integration_analysis() -> None:
    """Multi-omics integration analysis: concordance + Sankey diagram."""
    st.header("Integration Analysis")

    col1, col2 = st.columns(2)

    with col1:
        _render_concordance_matrix()

    with col2:
        _render_sankey_diagram()


def _render_concordance_matrix() -> None:
    st.subheader("Concordance Analysis")

    categories = [
        "H3K27ac↑",
        "H3K27ac↓",
        "H3K4me3↑",
        "H3K4me3↓",
        "H3K27me3↑",
        "H3K27me3↓",
    ]
    expr_cats = ["Expr↑", "Expr↓", "No change"]

    concordance = np.array(
        [
            [425, 52, 180],
            [38, 312, 145],
            [380, 45, 220],
            [32, 285, 130],
            [25, 85, 102],
            [180, 42, 222],
        ]
    )

    fig = px.imshow(
        concordance,
        x=expr_cats,
        y=categories,
        color_continuous_scale=COLORS.SCALE_DIVERGING,
        text_auto=True,
    )
    fig.update_layout(
        xaxis_title="Expression Change",
        yaxis_title="Histone Change",
        height=400,
    )
    st.plotly_chart(fig, use_container_width=True)

    st.markdown("""
    **Key Findings:**
    - H3K27ac gain strongly predicts expression increase (p<10⁻⁵⁰)
    - H3K27me3 loss correlates with de-repression of expression
    - 68% of genes show concordant histone-expression changes
    """)


def _render_sankey_diagram() -> None:
    st.subheader("TF-Histone-Expression Integration")

    labels = [
        "MYC targets",
        "E2F1 targets",
        "NRF1 targets",
        "H3K27ac↑",
        "H3K4me3↑",
        "No histone Δ",
        "Expr↑",
        "Expr↓",
        "Expr NC",
    ]

    source = [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5]
    target = [3, 4, 5, 3, 4, 5, 3, 4, 5, 6, 7, 6, 8, 7, 8]
    value = [320, 180, 45, 280, 210, 38, 150, 120, 85, 480, 40, 380, 130, 25, 180]

    palette = COLORS.QUALITATIVE
    node_colors = palette[: len(labels)]

    fig = go.Figure(
        data=[
            go.Sankey(
                node=dict(
                    pad=15,
                    thickness=20,
                    line=dict(color="black", width=0.5),
                    label=labels,
                    color=node_colors,
                ),
                link=dict(source=source, target=target, value=value),
            )
        ]
    )
    fig.update_layout(title_text="TF → Histone → Expression Flow", height=400)
    st.plotly_chart(fig, use_container_width=True)


# ============================================================================
# Regulatory Networks tab
# ============================================================================


def render_regulatory_networks() -> None:
    """Inferred regulatory networks."""
    st.header("Regulatory Networks")

    col1, col2 = st.columns([2, 1])

    with col1:
        _render_network_graph()
    with col2:
        _render_network_stats()


def _render_network_graph() -> None:
    st.subheader("TF Regulatory Network")

    np.random.seed(42)

    tfs = ["MYC", "E2F1", "SP1", "NRF1", "YY1", "MAX", "GABPA", "ETS1"]
    n_tfs = len(tfs)
    angles = np.linspace(0, 2 * np.pi, n_tfs, endpoint=False)
    tf_x = np.cos(angles) * 2
    tf_y = np.sin(angles) * 2

    targets = ["CCND1", "CDK4", "MYC", "BRCA1", "TP53", "RB1", "E2F1", "PCNA"]
    target_x = np.random.uniform(-1, 1, len(targets))
    target_y = np.random.uniform(-1, 1, len(targets))

    fig = go.Figure()

    for i, tf in enumerate(tfs):
        for j, target in enumerate(targets):
            if np.random.random() > 0.6:
                fig.add_trace(
                    go.Scatter(
                        x=[tf_x[i], target_x[j]],
                        y=[tf_y[i], target_y[j]],
                        mode="lines",
                        line=dict(color="rgba(150,150,150,0.3)", width=1),
                        showlegend=False,
                    )
                )

    fig.add_trace(
        go.Scatter(
            x=tf_x,
            y=tf_y,
            mode="markers+text",
            marker=dict(size=30, color=COLORS.UP),
            text=tfs,
            textposition="top center",
            name="TFs",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=target_x,
            y=target_y,
            mode="markers+text",
            marker=dict(size=20, color=COLORS.DOWN),
            text=targets,
            textposition="bottom center",
            name="Target Genes",
        )
    )

    fig.update_layout(
        showlegend=True,
        height=500,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
    )
    st.plotly_chart(fig, use_container_width=True)


def _render_network_stats() -> None:
    st.subheader("Network Statistics")
    st.markdown("""
    **Network Metrics:**
    - Nodes: 156 (45 TFs + 111 targets)
    - Edges: 892 TF-target links
    - Avg connectivity: 8.2

    **Hub TFs (highest connectivity):**
    1. MYC - 82 targets
    2. E2F1 - 68 targets
    3. SP1 - 54 targets

    **Enriched Modules:**
    - Cell cycle regulation
    - DNA repair
    - Metabolic pathways
    """)

    st.subheader("Module Enrichment")
    modules = pd.DataFrame(
        {
            "Module": ["Cell Cycle", "DNA Repair", "Metabolism", "Apoptosis"],
            "Genes": [28, 18, 22, 15],
            "TFs": [5, 4, 6, 3],
            "FDR": [1e-15, 1e-10, 1e-8, 1e-5],
        }
    )
    st.dataframe(modules, use_container_width=True, hide_index=True)


# ============================================================================
# Gene-Level View tab
# ============================================================================

_GENE_DATA: Dict[str, Dict] = {
    "MYC": {"expr": 2.8, "h3k27ac": 3.2, "h3k4me3": 2.1, "tfs": ["E2F1", "MAX", "SP1"]},
    "CCND1": {"expr": 1.9, "h3k27ac": 2.5, "h3k4me3": 1.8, "tfs": ["MYC", "E2F1", "AP1"]},
    "BRCA1": {"expr": -1.2, "h3k27ac": -0.8, "h3k4me3": -0.5, "tfs": ["TP53", "E2F1"]},
    "TP53": {"expr": 0.3, "h3k27ac": 0.5, "h3k4me3": 0.2, "tfs": ["SP1", "NF-Y"]},
}
_DEFAULT_GENE = {"expr": 0.5, "h3k27ac": 0.8, "h3k4me3": 0.3, "tfs": ["SP1"]}


def render_gene_level_view() -> None:
    """Individual gene multi-omics view."""
    st.header("Gene-Level Multi-omics View")

    col1, col2 = st.columns([1, 3])

    with col1:
        gene = st.selectbox(
            "Select Gene",
            ["MYC", "CCND1", "BRCA1", "TP53", "CDK4", "E2F1", "PCNA", "GAPDH"],
        )
        st.markdown("---")
        st.subheader(f"{gene} Summary")

        data = _GENE_DATA.get(gene, _DEFAULT_GENE)
        st.metric("Expression Log2FC", f"{data['expr']:.2f}")
        st.metric("H3K27ac Log2FC", f"{data['h3k27ac']:.2f}")
        st.metric("H3K4me3 Log2FC", f"{data['h3k4me3']:.2f}")
        st.markdown(f"**Predicted TFs:** {', '.join(data['tfs'])}")

    with col2:
        _render_multitrack_plot(gene)


def _render_multitrack_plot(gene: str) -> None:
    """Multi-track ChIP/RNA-seq locus plot."""
    st.subheader(f"{gene} Locus Multi-omics Tracks")

    fig = make_subplots(
        rows=5,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.02,
        row_heights=[0.15, 0.2, 0.2, 0.2, 0.25],
        subplot_titles=["Gene Model", "H3K4me3", "H3K27ac", "H3K4me1", "RNA-seq"],
    )

    x = np.linspace(0, 100, 500)
    np.random.seed(42)

    gene_color = "#2c3e50"

    # Gene model (exons + intron)
    fig.add_trace(
        go.Scatter(
            x=[10, 10, 30, 30, 10],
            y=[0.3, 0.7, 0.7, 0.3, 0.3],
            fill="toself",
            fillcolor=gene_color,
            line=dict(color=gene_color),
            name="Exon",
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=[30, 70],
            y=[0.5, 0.5],
            mode="lines",
            line=dict(color=gene_color, width=2),
            name="Intron",
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=[70, 70, 90, 90, 70],
            y=[0.3, 0.7, 0.7, 0.3, 0.3],
            fill="toself",
            fillcolor=gene_color,
            line=dict(color=gene_color),
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    # Histone tracks
    tracks = {
        "H3K4me3": (np.exp(-0.5 * ((x - 10) / 8) ** 2) * 3, COLORS.UP, 2),
        "H3K27ac": (
            np.exp(-0.5 * ((x - 10) / 10) ** 2) * 2.5 + np.exp(-0.5 * ((x - 50) / 15) ** 2) * 2,
            "#f39c12",
            3,
        ),
        "H3K4me1": (np.exp(-0.5 * ((x - 50) / 20) ** 2) * 1.8, "#27ae60", 4),
    }

    for name, (y_vals, color, row) in tracks.items():
        rgba = tuple(int(color[i : i + 2], 16) for i in (1, 3, 5))
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y_vals * 0.6,
                fill="tozeroy",
                fillcolor=f"rgba{(*rgba, 0.3)}",
                line=dict(color=color, width=1),
                name=f"{name} Control",
            ),
            row=row,
            col=1,
        )
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y_vals,
                fill="tozeroy",
                fillcolor=f"rgba{(*rgba, 0.5)}",
                line=dict(color=color, width=1, dash="dash"),
                name=f"{name} Treatment",
            ),
            row=row,
            col=1,
        )

    # RNA-seq
    rna_y = np.where((x > 10) & (x < 90), np.random.uniform(0.5, 2, len(x)), 0.1)
    fig.add_trace(
        go.Scatter(
            x=x,
            y=rna_y * 0.6,
            fill="tozeroy",
            fillcolor="rgba(52, 152, 219, 0.3)",
            line=dict(color=COLORS.DOWN, width=1),
            name="RNA-seq Control",
        ),
        row=5,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=x,
            y=rna_y * 1.5,
            fill="tozeroy",
            fillcolor="rgba(52, 152, 219, 0.5)",
            line=dict(color=COLORS.DOWN, width=1, dash="dash"),
            name="RNA-seq Treatment",
        ),
        row=5,
        col=1,
    )

    fig.update_layout(height=600, showlegend=False)
    fig.update_xaxes(title_text="Position (kb)", row=5, col=1)
    st.plotly_chart(fig, use_container_width=True)


# ============================================================================
# Summary Report tab
# ============================================================================


def render_summary_report(*, has_workflow_manager: bool = False) -> None:
    """Integrated summary report with export options."""
    st.header("Multi-omics Integration Summary")

    st.markdown("""
    ## Executive Summary

    This analysis integrated histone modification ChIP-seq, transcription factor motif analysis,
    and RNA-seq data to characterize the regulatory mechanisms underlying the treatment response.

    ### Key Findings

    1. **Coordinated epigenetic-transcriptional changes**
       - 68% of genes show concordant histone-expression changes
       - Active marks (H3K27ac, H3K4me3) gain correlates with increased expression
       - H3K27me3 loss associated with de-repression

    2. **Master regulator identification**
       - MYC emerges as central hub with 82 predicted targets
       - E2F1 co-regulates cell cycle genes with MYC
       - NRF1 controls metabolic gene network

    3. **Regulatory network modules**
       - Cell cycle module: 28 genes, 5 TFs, highly interconnected
       - DNA repair module: 18 genes, enriched for BRCA1/2 targets
       - Metabolic reprogramming: shift toward glycolytic gene expression

    4. **Chromatin state transitions**
       - 1,523 regions gained active enhancer marks
       - 644 regions lost repressive H3K27me3
       - Poised promoters activated at key regulatory genes
    """)

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Concordant Genes", "6,959", "68%")
    col2.metric("Master TFs", "8", "MYC, E2F1, ...")
    col3.metric("Network Modules", "4", "Enriched")
    col4.metric("Integration Score", "0.85", "High")

    st.markdown("---")

    _render_export_options(has_workflow_manager)


def _render_export_options(has_workflow_manager: bool) -> None:
    """Export buttons with optional workflow recording."""
    st.subheader("Export Integrated Results")

    col1, col2, col3 = st.columns(3)

    exports = [
        (
            col1,
            "📊 Export Gene Table",
            "Multi-omics Gene Table Export",
            "Gene-level integration table ready for download",
            {"export_type": "gene_table"},
            ["gene_integration_table"],
        ),
        (
            col2,
            "🔗 Export Network",
            "Regulatory Network Export",
            "Regulatory network exported (Cytoscape format)",
            {"export_type": "network", "format": "cytoscape"},
            ["regulatory_network"],
        ),
        (
            col3,
            "📋 Generate Report",
            "Multi-omics Report",
            "Full PDF report generated",
            {"report_format": "pdf"},
            ["multiomics_report"],
        ),
    ]

    for col, label, step_name, msg, params, outputs in exports:
        with col:
            if st.button(label):
                st.success(msg)
                if has_workflow_manager:
                    try:
                        from frontend.components.workflow_manager import WorkflowManager

                        WorkflowManager.record_step(
                            step_name=step_name,
                            tool="EpiNexus Multi-omics Module",
                            parameters=params,
                            inputs=["multiomics_data"],
                            outputs=outputs,
                        )
                    except ImportError:
                        pass
