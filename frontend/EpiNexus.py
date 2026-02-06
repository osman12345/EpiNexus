"""
EpiNexus - Comprehensive Epigenomics Analysis Platform

Main application entry point with categorized navigation and global state management.
"""

import streamlit as st
import requests
from pathlib import Path
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

# Page configuration
st.set_page_config(
    page_title="EpiNexus - Epigenomics Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for modern styling
st.markdown("""
<style>
    /* Main header styling */
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        background: linear-gradient(90deg, #1f77b4, #7b2cbf);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin-bottom: 0.5rem;
    }
    .sub-header {
        font-size: 1.2rem;
        color: #666;
        margin-bottom: 2rem;
    }

    /* Navigation category styling */
    .nav-category {
        font-size: 0.85rem;
        font-weight: 600;
        color: #1f77b4;
        margin-top: 1rem;
        margin-bottom: 0.5rem;
        padding: 0.3rem 0;
        border-bottom: 1px solid #e0e0e0;
    }

    /* Feature cards */
    .feature-card {
        background: linear-gradient(135deg, #f8f9fa 0%, #ffffff 100%);
        border-radius: 12px;
        padding: 1.5rem;
        text-align: center;
        border: 1px solid #e0e0e0;
        box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        transition: all 0.3s ease;
    }
    .feature-card:hover {
        box-shadow: 0 4px 8px rgba(0,0,0,0.1);
        transform: translateY(-2px);
    }
    .feature-icon {
        font-size: 2.5rem;
        margin-bottom: 0.5rem;
    }
    .feature-title {
        font-size: 1.1rem;
        font-weight: 600;
        margin-bottom: 0.5rem;
        color: #333;
    }
    .feature-desc {
        font-size: 0.9rem;
        color: #666;
    }

    /* Demo mode indicator */
    .demo-indicator {
        background: linear-gradient(90deg, #fff3cd, #ffeeba);
        border: 1px solid #ffc107;
        border-radius: 8px;
        padding: 0.5rem 1rem;
        font-size: 0.85rem;
        margin-bottom: 1rem;
    }

    /* Status badges */
    .badge-connected {
        background: #d4edda;
        color: #155724;
        padding: 0.25rem 0.5rem;
        border-radius: 4px;
        font-size: 0.8rem;
    }
    .badge-demo {
        background: #fff3cd;
        color: #856404;
        padding: 0.25rem 0.5rem;
        border-radius: 4px;
        font-size: 0.8rem;
    }

    /* Hide default Streamlit sidebar navigation */
    section[data-testid="stSidebar"] > div > div > div > div > ul {
        display: none;
    }
</style>
""", unsafe_allow_html=True)

# Initialize session state
def init_session_state():
    defaults = {
        "api_url": "http://localhost:8000",
        "current_project": None,
        "samples": [],
        "comparisons": [],
        "jobs": [],
        "using_demo_data": True,
        "data_loaded": False,
        "selected_genome": "hg38",
        "theme": "light"
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value

init_session_state()


# Navigation structure - Organized into categories
# File names match actual files in pages/ directory
NAVIGATION = {
    "🚀 Getting Started": {
        "pages": [
            {"name": "Home", "file": None, "icon": "🏠", "desc": "Dashboard & Quick Start"},
            {"name": "Data & Project", "file": "pages/01_data_project.py", "icon": "📁", "desc": "Load data, manage projects"},
        ]
    },
    "📊 Quality Control": {
        "pages": [
            {"name": "QC Dashboard", "file": "pages/02_quality_control.py", "icon": "📋", "desc": "Quality metrics"},
            {"name": "ATAC-seq QC", "file": "pages/16_atacseq.py", "icon": "🔓", "desc": "ATAC-specific QC"},
        ]
    },
    "🔬 Core Analysis": {
        "pages": [
            {"name": "Differential", "file": "pages/03_differential.py", "icon": "📈", "desc": "Differential analysis"},
            {"name": "Annotation", "file": "pages/05_annotation.py", "icon": "🏷️", "desc": "Annotate peaks"},
            {"name": "Super-Enhancers", "file": "pages/12_super_enhancers.py", "icon": "⭐", "desc": "ROSE algorithm"},
        ]
    },
    "🔬 TF & Motif Analysis": {
        "pages": [
            {"name": "TF ChIP-seq", "file": "pages/20_tf_chipseq.py", "icon": "🔬", "desc": "TF binding analysis"},
            {"name": "Motif Viewer", "file": "pages/17_motif_viewer.py", "icon": "🔤", "desc": "Sequence logos"},
            {"name": "TF Networks", "file": "pages/08_tf_analysis.py", "icon": "🕸️", "desc": "TF networks"},
        ]
    },
    "🧬 Integration": {
        "pages": [
            {"name": "Multi-Mark", "file": "pages/06_multimark.py", "icon": "🔗", "desc": "Histone integration"},
            {"name": "RNA-seq", "file": "pages/09_expression.py", "icon": "📊", "desc": "Expression data"},
            {"name": "Multi-omics", "file": "pages/10_multiomics.py", "icon": "🧬", "desc": "Combine modalities"},
            {"name": "Peak-Gene Links", "file": "pages/18_peak_gene_linking.py", "icon": "🔗", "desc": "ABC model"},
            {"name": "GWAS Overlap", "file": "pages/19_gwas_overlap.py", "icon": "🧬", "desc": "Disease variants"},
        ]
    },
    "📊 Visualization": {
        "pages": [
            {"name": "Charts", "file": "pages/04_visualization.py", "icon": "📊", "desc": "Interactive plots"},
            {"name": "Genome Browser", "file": "pages/11_genome_browser.py", "icon": "🔬", "desc": "IGV.js viewer"},
        ]
    },
    "📦 Data & Reports": {
        "pages": [
            {"name": "ENCODE Data", "file": "pages/14_encode_integration.py", "icon": "🏛️", "desc": "Public data"},
            {"name": "Reports", "file": "pages/07_reports.py", "icon": "📑", "desc": "Generate reports"},
            {"name": "Batch Processing", "file": "pages/13_batch_processing.py", "icon": "📦", "desc": "Bulk analysis"},
        ]
    },
    "❓ Help": {
        "pages": [
            {"name": "Help & Support", "file": "pages/21_help.py", "icon": "❓", "desc": "Docs, FAQ, contact"},
        ]
    }
}


def check_api_connection():
    """Check if the API is available."""
    try:
        response = requests.get(f"{st.session_state.api_url}/health", timeout=2)
        return response.status_code == 200
    except:
        return False


def render_sidebar():
    """Render the enhanced sidebar with categorized navigation."""
    with st.sidebar:
        # Logo and title
        st.markdown("""
        <div style="text-align: center; margin-bottom: 1rem;">
            <h1 style="margin: 0; font-size: 1.8rem;">🧬 EpiNexus</h1>
            <p style="margin: 0; font-size: 0.85rem; color: #666;">Epigenomics Analysis</p>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("---")

        # Data status indicator
        if st.session_state.using_demo_data:
            st.markdown("""
            <div style="background: #fff3cd; border: 1px solid #ffc107; border-radius: 6px;
                        padding: 0.5rem; text-align: center; margin-bottom: 1rem;">
                <span style="font-size: 0.85rem;">📌 Demo Mode</span>
            </div>
            """, unsafe_allow_html=True)
        else:
            st.markdown("""
            <div style="background: #d4edda; border: 1px solid #28a745; border-radius: 6px;
                        padding: 0.5rem; text-align: center; margin-bottom: 1rem;">
                <span style="font-size: 0.85rem;">✅ Your Data</span>
            </div>
            """, unsafe_allow_html=True)

        # Navigation categories
        for category, content in NAVIGATION.items():
            st.markdown(f'<div class="nav-category">{category}</div>', unsafe_allow_html=True)

            for page in content["pages"]:
                if page["file"]:
                    if st.button(
                        f"{page['icon']} {page['name']}",
                        key=f"nav_{page['name']}",
                        use_container_width=True,
                        help=page['desc']
                    ):
                        st.switch_page(page["file"])
                else:
                    # Home page - no switch needed
                    if st.button(
                        f"{page['icon']} {page['name']}",
                        key=f"nav_{page['name']}",
                        use_container_width=True,
                        help=page['desc'],
                        type="primary" if page["name"] == "Home" else "secondary"
                    ):
                        st.rerun()

        st.markdown("---")

        # Quick stats
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Samples", len(st.session_state.samples), label_visibility="collapsed")
            st.caption("Samples")
        with col2:
            st.metric("Jobs", len(st.session_state.jobs), label_visibility="collapsed")
            st.caption("Jobs")

        # Genome selector
        st.selectbox(
            "Reference Genome",
            ["hg38", "hg19", "mm10", "mm39", "dm6"],
            key="selected_genome",
            help="Select reference genome for analysis"
        )

        st.markdown("---")
        st.caption("EpiNexus v1.0.0")
        st.caption("[Documentation](https://github.com) | [Report Issue](https://github.com)")


def render_home_dashboard():
    """Render the main dashboard / home page."""

    # Header
    st.markdown('<p class="main-header">🧬 EpiNexus</p>', unsafe_allow_html=True)
    st.markdown(
        '<p class="sub-header">Comprehensive Epigenomics Analysis Platform for ChIP-seq, CUT&Tag, and ATAC-seq</p>',
        unsafe_allow_html=True
    )

    # Demo mode notice
    if st.session_state.using_demo_data:
        st.info("📌 **Demo Mode**: You're viewing the platform with simulated data. "
               "Upload your own data or use the preprocessing pipeline to analyze real experiments.")

    # Quick start cards
    st.markdown("### 🚀 Quick Start")

    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.markdown("""
        <div class="feature-card">
            <div class="feature-icon">📁</div>
            <div class="feature-title">New Project</div>
            <div class="feature-desc">Create project & load data</div>
        </div>
        """, unsafe_allow_html=True)
        if st.button("Get Started", key="quick_project", use_container_width=True):
            st.switch_page("pages/01_data_project.py")

    with col2:
        st.markdown("""
        <div class="feature-card">
            <div class="feature-icon">🔬</div>
            <div class="feature-title">TF Analysis</div>
            <div class="feature-desc">Transcription factor binding</div>
        </div>
        """, unsafe_allow_html=True)
        if st.button("TF ChIP-seq", key="quick_tf", use_container_width=True):
            st.switch_page("pages/20_tf_chipseq.py")

    with col3:
        st.markdown("""
        <div class="feature-card">
            <div class="feature-icon">📊</div>
            <div class="feature-title">Differential</div>
            <div class="feature-desc">Compare conditions</div>
        </div>
        """, unsafe_allow_html=True)
        if st.button("Run Analysis", key="quick_diff", use_container_width=True):
            st.switch_page("pages/03_differential.py")

    with col4:
        st.markdown("""
        <div class="feature-card">
            <div class="feature-icon">❓</div>
            <div class="feature-title">Help</div>
            <div class="feature-desc">Docs, FAQ, support</div>
        </div>
        """, unsafe_allow_html=True)
        if st.button("Get Help", key="quick_help", use_container_width=True):
            st.switch_page("pages/21_help.py")

    st.markdown("---")

    # Supported assays
    st.markdown("### 🧬 Supported Assays")

    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("""
        **ChIP-seq**
        - 🧬 Histone Marks (H3K27ac, H3K4me3, H3K27me3...)
        - 🔬 Transcription Factors (MYC, P53, CTCF...)
        - Narrow & broad peak detection
        """)

    with col2:
        st.markdown("""
        **CUT&Tag / CUT&RUN**
        - Low-input profiling
        - Spike-in normalization
        - SEACR peak calling
        """)

    with col3:
        st.markdown("""
        **ATAC-seq**
        - Chromatin accessibility
        - Nucleosome positioning
        - TF footprinting
        """)

    st.markdown("---")

    # Analysis capabilities
    st.markdown("### 📊 Analysis Capabilities")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("""
        **Core Analysis**
        - ✅ Differential binding (DiffBind)
        - ✅ Peak annotation (ChIPseeker)
        - ✅ Super-enhancer detection (ROSE)
        - ✅ Motif enrichment analysis
        - ✅ TF binding analysis
        """)

    with col2:
        st.markdown("""
        **Integration & Visualization**
        - ✅ Multi-mark integration
        - ✅ RNA-seq correlation
        - ✅ Peak-to-gene linking (ABC model)
        - ✅ GWAS variant overlap
        - ✅ Interactive genome browser
        """)

    st.markdown("---")

    # Recent jobs
    st.markdown("### 📋 Recent Jobs")

    if st.session_state.jobs:
        for job in st.session_state.jobs[-5:]:
            status_icon = {
                "completed": "✅",
                "running": "🔄",
                "pending": "⏳",
                "failed": "❌"
            }.get(job.get("status", "pending"), "⏳")

            col1, col2, col3 = st.columns([1, 3, 1])
            with col1:
                st.write(status_icon)
            with col2:
                st.write(f"**{job.get('name', 'Unnamed')}**")
            with col3:
                st.write(job.get("status", "pending"))
    else:
        st.info("No jobs yet. Start an analysis to see your jobs here!")

    # Data entry points comparison
    st.markdown("---")
    st.markdown("### 📥 Data Entry Points")

    st.markdown("""
    | Entry Point | Input | What You Get |
    |-------------|-------|--------------|
    | **FASTQ** | Raw sequencing files | Full pipeline: QC → Alignment → Peak calling → Analysis |
    | **BAM** | Aligned reads | Peak calling → Differential analysis → Visualization |
    | **Peaks** | BED/narrowPeak files | Differential analysis → Integration → Visualization |
    """)


def main():
    """Main application."""
    render_sidebar()
    render_home_dashboard()


if __name__ == "__main__":
    main()
