"""
Annotation Page
Annotate peaks with genomic features and functional information.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

st.set_page_config(page_title="Annotation - EpiNexus", page_icon="🏷️", layout="wide")

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
            <div style="font-size: 3rem; margin-bottom: 1rem;">🏷️</div>
            <h2 style="color: #6c757d; margin-bottom: 0.5rem;">No Data Loaded</h2>
            <p style="color: #6c757d; font-size: 1.1rem;">
                Upload your peak files to annotate them with genomic features.
            </p>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("")

        if st.button("📁 Go to Data & Project", type="primary", use_container_width=True):
            st.switch_page("pages/01_data_project.py")

        st.markdown("")
        st.markdown("**Available annotations:**")
        st.markdown("- Gene annotation (nearest TSS)")
        st.markdown("- Pathway & GO enrichment")
        st.markdown("- Regulatory element overlap")
        st.markdown("- TF motif enrichment")


def main():
    st.title("🏷️ Peak Annotation")
    st.markdown("Annotate peaks with genomic features, nearby genes, and functional information.")

    # Check if data is loaded
    if not has_data():
        render_empty_state()
        return

    tab1, tab2, tab3, tab4 = st.tabs([
        "🧬 Gene Annotation",
        "🔗 Pathway Analysis",
        "📍 Regulatory Elements",
        "🧪 Functional Enrichment"
    ])

    with tab1:
        render_gene_annotation()
    with tab2:
        render_pathway_analysis()
    with tab3:
        render_regulatory_elements()
    with tab4:
        render_functional_enrichment()


def render_gene_annotation():
    """Annotate peaks with nearest genes."""
    st.header("Gene Annotation")

    col1, col2 = st.columns([1, 2])

    with col1:
        st.subheader("Annotation Settings")

        tss_distance = st.slider("TSS distance threshold (kb)", 1, 100, 10)

        annotation_db = st.selectbox(
            "Gene annotation database",
            ["GENCODE v43", "RefSeq", "UCSC Known Genes", "Ensembl"]
        )

        include_types = st.multiselect(
            "Gene types to include",
            ["protein_coding", "lncRNA", "miRNA", "pseudogene", "snRNA"],
            default=["protein_coding", "lncRNA"]
        )

        if st.button("Run Annotation", type="primary"):
            st.success("Annotation complete!")

    with col2:
        st.subheader("Annotated Peaks")

        # Demo data
        np.random.seed(42)
        n = 100

        annotated = pd.DataFrame({
            'Peak_ID': [f'peak_{i}' for i in range(n)],
            'Chr': np.random.choice(['chr1', 'chr2', 'chr3', 'chr4'], n),
            'Start': np.random.randint(1000000, 200000000, n),
            'End': lambda: 0,
            'Nearest_Gene': np.random.choice(['MYC', 'BRCA1', 'TP53', 'EGFR', 'KRAS', 'PTEN', 'CDK4', 'MDM2'], n),
            'Distance_TSS': np.random.randint(-50000, 50000, n),
            'Gene_Type': np.random.choice(['protein_coding', 'lncRNA'], n, p=[0.85, 0.15]),
            'Feature': np.random.choice(['Promoter', 'Intron', 'Exon', 'Intergenic'], n, p=[0.35, 0.3, 0.1, 0.25])
        })
        annotated['End'] = annotated['Start'] + np.random.randint(200, 3000, n)

        st.dataframe(annotated, use_container_width=True, hide_index=True)

        # Summary stats
        st.markdown("**Summary:**")
        col1, col2, col3 = st.columns(3)
        col1.metric("Peaks annotated", len(annotated))
        col2.metric("Unique genes", annotated['Nearest_Gene'].nunique())
        col3.metric("At promoters", f"{(annotated['Feature'] == 'Promoter').mean():.1%}")


def render_pathway_analysis():
    """Pathway and GO enrichment analysis."""
    st.header("Pathway Analysis")

    col1, col2 = st.columns([1, 2])

    with col1:
        st.subheader("Analysis Settings")

        gene_set = st.selectbox(
            "Gene set database",
            ["MSigDB Hallmarks", "KEGG Pathways", "Reactome", "GO Biological Process",
             "GO Molecular Function", "GO Cellular Component"]
        )

        background = st.selectbox(
            "Background genes",
            ["All genes", "Expressed genes", "Custom list"]
        )

        method = st.selectbox(
            "Statistical method",
            ["Hypergeometric", "Fisher's exact", "GSEA"]
        )

    with col2:
        st.subheader("Enriched Pathways")

        # Demo pathway results
        pathways = pd.DataFrame({
            'Pathway': [
                'HALLMARK_E2F_TARGETS',
                'HALLMARK_MYC_TARGETS_V1',
                'HALLMARK_G2M_CHECKPOINT',
                'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                'HALLMARK_DNA_REPAIR',
                'HALLMARK_MTORC1_SIGNALING',
                'HALLMARK_UNFOLDED_PROTEIN_RESPONSE',
                'HALLMARK_GLYCOLYSIS'
            ],
            'Overlap': ['45/200', '38/200', '32/200', '28/200', '25/200', '22/200', '20/200', '18/200'],
            'P-value': [1.2e-15, 3.4e-12, 8.9e-10, 2.1e-8, 5.6e-7, 1.8e-6, 4.2e-5, 9.3e-5],
            'FDR': [4.8e-14, 6.8e-11, 1.2e-8, 2.1e-7, 4.5e-6, 1.2e-5, 2.4e-4, 4.7e-4]
        })

        # Bar plot of enrichment
        fig = px.bar(pathways, x=-np.log10(pathways['FDR']), y='Pathway',
                    orientation='h', color=-np.log10(pathways['FDR']),
                    color_continuous_scale='Reds')
        fig.update_layout(
            xaxis_title='-Log10 FDR',
            yaxis_title='',
            height=400,
            showlegend=False
        )
        st.plotly_chart(fig, use_container_width=True)

        st.dataframe(pathways, use_container_width=True, hide_index=True)


def render_regulatory_elements():
    """Overlap with regulatory element databases."""
    st.header("Regulatory Element Overlap")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("ENCODE cCRE Overlap")

        # Demo cCRE data
        ccre_types = ['PLS (Promoter)', 'pELS (proximal Enhancer)', 'dELS (distal Enhancer)',
                     'DNase-H3K4me3', 'CTCF-only']
        overlaps = [1250, 890, 2100, 320, 180]

        fig = px.bar(x=ccre_types, y=overlaps, color=ccre_types,
                    color_discrete_sequence=px.colors.qualitative.Set2)
        fig.update_layout(
            xaxis_title='cCRE Type',
            yaxis_title='Overlapping Peaks',
            showlegend=False,
            height=400
        )
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("Chromatin State (ChromHMM)")

        states = ['Active TSS', 'Flanking TSS', 'Strong Enhancer', 'Weak Enhancer',
                 'Repressed', 'Heterochromatin', 'Quiescent']
        percentages = [35, 15, 25, 12, 5, 3, 5]

        fig = px.pie(values=percentages, names=states, hole=0.4,
                    color_discrete_sequence=px.colors.qualitative.Set3)
        fig.update_layout(height=400)
        st.plotly_chart(fig, use_container_width=True)

    st.subheader("Overlap with Published Datasets")

    st.markdown("Compare your peaks with published ChIP-seq datasets from ENCODE/Roadmap:")

    datasets = pd.DataFrame({
        'Dataset': ['ENCODE H3K27ac K562', 'ENCODE H3K4me3 K562', 'Roadmap H3K27ac Blood',
                   'ENCODE CTCF K562', 'ENCODE H3K27me3 K562'],
        'Total Peaks': [45000, 32000, 52000, 28000, 18000],
        'Overlap': [12500, 8200, 15800, 3200, 1500],
        'Percent': [27.8, 25.6, 30.4, 11.4, 8.3],
        'P-value': [1.2e-120, 3.4e-85, 8.9e-150, 2.1e-15, 0.045]
    })

    st.dataframe(datasets, use_container_width=True, hide_index=True)


def render_functional_enrichment():
    """Functional enrichment analysis."""
    st.header("Functional Enrichment")

    enrichment_type = st.selectbox(
        "Enrichment type",
        ["TF Binding Sites", "Disease Associations", "Drug Targets", "Protein Domains"]
    )

    if enrichment_type == "TF Binding Sites":
        st.subheader("Enriched TF Motifs")

        tf_results = pd.DataFrame({
            'TF': ['MYC', 'E2F1', 'SP1', 'NRF1', 'YY1', 'ETS1', 'GABPA', 'CREB1'],
            'Motif': ['CACGTG', 'TTTCGCGC', 'GGGCGG', 'GCGCATGCGC', 'CAAGATGGCG', 'GGAA', 'CGGAAGT', 'TGACGTCA'],
            'Enrichment': [8.5, 6.2, 5.8, 4.9, 4.2, 3.8, 3.5, 3.2],
            'P-value': [1e-25, 1e-18, 1e-15, 1e-12, 1e-10, 1e-8, 1e-7, 1e-6],
            'Targets': [245, 189, 312, 156, 134, 98, 87, 76]
        })

        col1, col2 = st.columns(2)

        with col1:
            fig = px.bar(tf_results, x='Enrichment', y='TF', orientation='h',
                        color='Enrichment', color_continuous_scale='Reds')
            fig.update_layout(height=400, xaxis_title='Fold Enrichment', yaxis_title='')
            st.plotly_chart(fig, use_container_width=True)

        with col2:
            st.dataframe(tf_results, use_container_width=True, hide_index=True)

    elif enrichment_type == "Disease Associations":
        st.subheader("Disease-Gene Associations (DisGeNET)")

        disease_results = pd.DataFrame({
            'Disease': ['Breast Cancer', 'Colorectal Cancer', 'Leukemia', 'Lung Cancer', 'Melanoma'],
            'Associated Genes': [45, 38, 32, 28, 22],
            'Enrichment': [4.2, 3.8, 3.5, 3.1, 2.8],
            'FDR': [1e-12, 1e-10, 1e-8, 1e-6, 1e-5]
        })

        fig = px.bar(disease_results, x='Disease', y='Enrichment', color='Associated Genes',
                    color_continuous_scale='Blues')
        fig.update_layout(height=400)
        st.plotly_chart(fig, use_container_width=True)


if __name__ == "__main__":
    main()
