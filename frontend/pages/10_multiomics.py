"""
Integrated Multi-omics Analysis Page

Comprehensive integration of:
- Histone modifications
- Transcription factor binding
- Gene expression data
- DNA methylation

Copyright (c) 2026 EpiNexus Contributors
SPDX-License-Identifier: AGPL-3.0-or-later OR Commercial
"""

import streamlit as st
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

st.set_page_config(
    page_title="Multi-omics Integration - EpiNexus",
    page_icon="🔬",
    layout="wide"
)

# Try to import data manager
try:
    from frontend.components.data_manager import DataManager
    HAS_DATA_MANAGER = True
except ImportError:
    HAS_DATA_MANAGER = False

# Try to import workflow manager
try:
    from frontend.components.workflow_manager import WorkflowManager
    HAS_WORKFLOW_MANAGER = True
except ImportError:
    try:
        from components.workflow_manager import WorkflowManager
        HAS_WORKFLOW_MANAGER = True
    except ImportError:
        HAS_WORKFLOW_MANAGER = False

# Import rendering helpers from the extracted module
from frontend.components.multiomics import (
    render_data_overview,
    render_integration_analysis,
    render_regulatory_networks,
    render_gene_level_view,
    render_summary_report,
)


def has_data() -> bool:
    """Check if user has loaded data."""
    if HAS_DATA_MANAGER:
        peaks = DataManager.get_data('peaks')
        return peaks is not None and len(peaks) > 0
    return len(st.session_state.get('samples', [])) > 0


def render_empty_state() -> None:
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
                Upload your multi-omics data to run integrative analysis.
            </p>
        </div>
        """, unsafe_allow_html=True)
        st.markdown("")
        if st.button("📁 Go to Data & Project", type="primary", use_container_width=True):
            st.switch_page("pages/01_data_project.py")
        st.markdown("")
        st.markdown("**Multi-omics integration:**")
        st.markdown("- Histone modifications")
        st.markdown("- TF binding data")
        st.markdown("- Gene expression")
        st.markdown("- DNA methylation")
        st.markdown("- Regulatory network inference")


def main() -> None:
    st.title("🔬 Integrated Multi-omics Analysis")
    st.markdown("""
    Comprehensive integration of histone modifications, transcription factor binding,
    gene expression, and DNA methylation data to uncover regulatory mechanisms.
    """)

    if not has_data():
        render_empty_state()
        return

    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "📊 Data Overview",
        "🔗 Integration Analysis",
        "🎯 Regulatory Networks",
        "🧬 Gene-Level View",
        "📋 Summary Report"
    ])

    with tab1:
        render_data_overview()
    with tab2:
        render_integration_analysis()
    with tab3:
        render_regulatory_networks()
    with tab4:
        render_gene_level_view()
    with tab5:
        render_summary_report(has_workflow_manager=HAS_WORKFLOW_MANAGER)


if __name__ == "__main__":
    main()
