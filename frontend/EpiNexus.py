"""
EpiNexus - Comprehensive Epigenomics Analysis Platform

Main application entry point with programmatic navigation.

Copyright (c) 2026 EpiNexus Contributors
SPDX-License-Identifier: AGPL-3.0-or-later OR Commercial
"""

import streamlit as st
from pathlib import Path
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

# Page configuration - MUST be first Streamlit command
st.set_page_config(
    page_title="EpiNexus - Epigenomics Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# =============================================================================
# INITIALIZE SESSION STATE
# =============================================================================

def init_session_state():
    defaults = {
        "current_project": None,
        "samples": [],
        "comparisons": [],
        "jobs": [],
        "using_demo_data": True,
        "data_loaded": False,
        "selected_genome": "hg38",
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value

init_session_state()

# =============================================================================
# HOME PAGE CONTENT (must be defined before st.Page references it)
# =============================================================================

def render_home_dashboard():
    """Render the main dashboard / home page."""

    st.markdown("""
    <h1 style="background: linear-gradient(90deg, #1f77b4, #7b2cbf);
               -webkit-background-clip: text; -webkit-text-fill-color: transparent;">
        🧬 EpiNexus
    </h1>
    <p style="font-size: 1.2rem; color: #666;">
        Comprehensive Epigenomics Analysis Platform for ChIP-seq, CUT&Tag, and ATAC-seq
    </p>
    """, unsafe_allow_html=True)

    # Quick start cards
    st.markdown("### 🚀 Quick Start")

    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.markdown("#### 📁 New Project")
        st.caption("Create project & load data")
        if st.button("Get Started", key="quick_project", use_container_width=True):
            st.switch_page("pages/01_data_project.py")

    with col2:
        st.markdown("#### 🔬 TF Analysis")
        st.caption("Transcription factor binding")
        if st.button("TF ChIP-seq", key="quick_tf", use_container_width=True):
            st.switch_page("pages/20_tf_chipseq.py")

    with col3:
        st.markdown("#### 📊 Differential")
        st.caption("Compare conditions")
        if st.button("Run Analysis", key="quick_diff", use_container_width=True):
            st.switch_page("pages/03_differential.py")

    with col4:
        st.markdown("#### ❓ Help")
        st.caption("Docs, FAQ, support")
        if st.button("Get Help", key="quick_help", use_container_width=True):
            st.switch_page("pages/21_help.py")

    st.markdown("---")

    # Supported assays
    st.markdown("### 🧬 Supported Assays")

    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("""
        **ChIP-seq**
        - Histone Marks (H3K27ac, H3K4me3...)
        - Transcription Factors (MYC, P53...)
        """)

    with col2:
        st.markdown("""
        **CUT&Tag / CUT&RUN**
        - Low-input profiling
        - Spike-in normalization
        """)

    with col3:
        st.markdown("""
        **ATAC-seq**
        - Chromatin accessibility
        - TF footprinting
        """)

    st.markdown("---")

    # Analysis capabilities
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("""
        ### 📊 Core Analysis
        - ✅ Differential binding (DiffBind)
        - ✅ Peak annotation (ChIPseeker)
        - ✅ Super-enhancer detection (ROSE)
        - ✅ Motif enrichment
        """)

    with col2:
        st.markdown("""
        ### 🔗 Integration
        - ✅ Multi-mark integration
        - ✅ RNA-seq correlation
        - ✅ Peak-to-gene linking (ABC)
        - ✅ GWAS variant overlap
        """)

    # Footer
    st.markdown("---")
    st.caption("EpiNexus v1.0.0 · © 2026 EpiNexus Contributors · AGPL-3.0 / Commercial")


# =============================================================================
# PROGRAMMATIC NAVIGATION - Replaces automatic page discovery
# =============================================================================

pages_dir = Path(__file__).parent / "pages"

# Define pages
pages = {
    "Getting Started": [
        st.Page(render_home_dashboard, title="Home", icon="🏠", default=True),
        st.Page(str(pages_dir / "01_data_project.py"), title="Data & Project", icon="📁"),
        st.Page(str(pages_dir / "00_workflow.py"), title="Workflow Guide", icon="🔄"),
    ],
    "Preprocessing": [
        st.Page(str(pages_dir / "22_alignment.py"), title="Alignment", icon="🧬"),
        st.Page(str(pages_dir / "23_peak_calling.py"), title="Peak Calling", icon="📍"),
    ],
    "Quality Control": [
        st.Page(str(pages_dir / "02_quality_control.py"), title="QC Dashboard", icon="📋"),
        st.Page(str(pages_dir / "16_atacseq.py"), title="ATAC-seq QC", icon="🔓"),
    ],
    "Core Analysis": [
        st.Page(str(pages_dir / "03_differential.py"), title="Differential", icon="📈"),
        st.Page(str(pages_dir / "05_annotation.py"), title="Annotation", icon="🏷️"),
        st.Page(str(pages_dir / "12_super_enhancers.py"), title="Super-Enhancers", icon="⭐"),
    ],
    "TF & Motif": [
        st.Page(str(pages_dir / "20_tf_chipseq.py"), title="TF ChIP-seq", icon="🔬"),
        st.Page(str(pages_dir / "17_motif_viewer.py"), title="Motif Viewer", icon="🔤"),
        st.Page(str(pages_dir / "08_tf_analysis.py"), title="TF Networks", icon="🕸️"),
    ],
    "Integration": [
        st.Page(str(pages_dir / "06_multimark.py"), title="Multi-Mark", icon="🔗"),
        st.Page(str(pages_dir / "09_expression.py"), title="RNA-seq", icon="📊"),
        st.Page(str(pages_dir / "10_multiomics.py"), title="Multi-omics", icon="🧬"),
        st.Page(str(pages_dir / "18_peak_gene_linking.py"), title="Peak-Gene Links", icon="🔗"),
        st.Page(str(pages_dir / "19_gwas_overlap.py"), title="GWAS Overlap", icon="🧬"),
    ],
    "Visualization": [
        st.Page(str(pages_dir / "04_visualization.py"), title="Charts", icon="📊"),
    ],
    "Data & Reports": [
        st.Page(str(pages_dir / "14_encode_integration.py"), title="ENCODE Data", icon="🏛️"),
        st.Page(str(pages_dir / "07_reports.py"), title="Reports", icon="📑"),
        st.Page(str(pages_dir / "13_batch_processing.py"), title="Batch Processing", icon="📦"),
    ],
    "Help": [
        st.Page(str(pages_dir / "21_help.py"), title="Help & Docs", icon="❓"),
    ],
}

# Create navigation
pg = st.navigation(pages)

# =============================================================================
# CUSTOM SIDEBAR CONTENT
# =============================================================================

with st.sidebar:
    # Quick stats
    col1, col2 = st.columns(2)
    with col1:
        st.caption(f"📊 {len(st.session_state.samples)} samples")
    with col2:
        st.caption(f"⚙️ {len(st.session_state.jobs)} jobs")

    st.markdown("---")

    # Genome selector at bottom
    st.selectbox(
        "Reference Genome",
        ["hg38", "hg19", "mm10", "mm39", "dm6"],
        key="selected_genome",
    )

    st.markdown("---")
    st.caption("v1.0.0 · [Docs](https://epinexus.readthedocs.io) · [GitHub](https://github.com/epinexus)")
    st.caption("© 2026 EpiNexus Contributors")

# =============================================================================
# RUN THE APP
# =============================================================================

pg.run()
