"""
Analysis Workflow Page
Clear step-by-step guide through the analysis pipeline.
"""

import streamlit as st
import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

st.set_page_config(page_title="Workflow - EpiNexus", page_icon="🔄", layout="wide")

# Try to import data manager
try:
    from frontend.components.data_manager import DataManager
    HAS_DATA_MANAGER = True
except ImportError:
    HAS_DATA_MANAGER = False


def get_workflow_status():
    """Get current status of each workflow step."""
    status = {
        'project': False,
        'alignment': False,
        'peak_calling': False,
        'data': False,
        'qc': False,
        'analysis': False,
        'results': False
    }

    # Check project
    if st.session_state.get('project_name'):
        status['project'] = True

    # Check alignment
    if st.session_state.get('aligned_bams'):
        status['alignment'] = True

    # Check peak calling
    if st.session_state.get('called_peaks') or st.session_state.get('uploaded_peaks'):
        status['peak_calling'] = True

    # Check data (peaks loaded)
    if HAS_DATA_MANAGER:
        peaks = DataManager.get_data('peaks')
        if peaks is not None and len(peaks) > 0:
            status['data'] = True
    elif st.session_state.get('samples'):
        status['data'] = True
    elif status['peak_calling']:
        status['data'] = True

    # Check QC (assume done if data is loaded)
    if status['data']:
        status['qc'] = st.session_state.get('qc_completed', False)

    # Check analysis
    if 'diff_results' in st.session_state and st.session_state.diff_results is not None:
        status['analysis'] = True

    # Check results
    if status['analysis']:
        status['results'] = True

    return status


def main():
    st.title("🔄 Analysis Workflow")
    st.markdown("Follow these steps to complete your epigenomics analysis.")

    # Get current status
    status = get_workflow_status()

    # Calculate progress (exclude optional preprocessing from percentage)
    core_steps = ['project', 'data', 'qc', 'analysis', 'results']
    completed = sum(1 for k in core_steps if status.get(k, False))
    total = len(core_steps)
    progress = completed / total

    # Progress bar
    st.progress(progress)
    st.markdown(f"**Progress: {completed}/{total} core steps completed**")

    # Show preprocessing status
    if status['alignment'] or status['peak_calling']:
        st.success(f"✅ Preprocessing: Alignment {'✓' if status['alignment'] else '○'} | Peak Calling {'✓' if status['peak_calling'] else '○'}")

    st.markdown("---")

    # Workflow steps
    render_workflow_steps(status)

    # Quick actions sidebar
    with st.sidebar:
        st.header("Quick Actions")

        if not status['project']:
            if st.button("🆕 Start New Project", use_container_width=True):
                st.switch_page("pages/01_data_project.py")
        elif not status['data']:
            if st.button("📤 Upload Data", use_container_width=True):
                st.session_state.workflow_step = 2
                st.switch_page("pages/01_data_project.py")
        elif not status['qc']:
            if st.button("✅ Run QC", use_container_width=True):
                st.switch_page("pages/02_quality_control.py")
        elif not status['analysis']:
            if st.button("📊 Run Analysis", use_container_width=True):
                st.switch_page("pages/03_differential.py")
        else:
            if st.button("📥 Export Results", use_container_width=True):
                st.switch_page("pages/07_reports.py")

        st.markdown("---")

        # Current project info
        st.subheader("Current Project")
        st.write(f"**Name:** {st.session_state.get('project_name', 'Not set')}")
        st.write(f"**Assay:** {st.session_state.get('assay_type', 'Not set')}")
        st.write(f"**Samples:** {len(st.session_state.get('samples', []))}")


def render_workflow_steps(status):
    """Render the workflow steps with status indicators."""

    # Step 1: Project Setup
    render_step(
        number=1,
        title="Project Setup",
        description="Create a new project or load an existing one. Define your assay type and reference genome.",
        status='completed' if status['project'] else 'current' if not any(status.values()) else 'pending',
        page="pages/01_data_project.py",
        details={
            "What you'll do": [
                "Name your project",
                "Select assay type (CUT&Tag, ChIP-seq, etc.)",
                "Choose reference genome (hg38, mm10, etc.)"
            ],
            "Time estimate": "~1 minute"
        }
    )

    # Optional preprocessing section
    st.markdown("### 🧬 Preprocessing (Optional)")
    st.caption("Skip these steps if you already have peak files (BED/narrowPeak)")

    # Step 2a: Alignment
    render_step(
        number="2a",
        title="Alignment (FASTQ → BAM)",
        description="Align raw sequencing reads to the reference genome using Bowtie2 or BWA.",
        status='completed' if status['alignment'] else 'current' if status['project'] and not status['data'] else 'pending',
        page="pages/22_alignment.py",
        details={
            "Supported aligners": [
                "Bowtie2 (ChIP-seq, CUT&Tag)",
                "BWA-MEM (ATAC-seq)"
            ],
            "Requirements": [
                "FASTQ files (single or paired-end)",
                "Reference genome index",
                "Bowtie2 or BWA installed locally"
            ],
            "Time estimate": "~30-60 min per sample"
        }
    )

    # Step 2b: Peak Calling
    render_step(
        number="2b",
        title="Peak Calling (BAM → Peaks)",
        description="Identify regions of enrichment using MACS2 or SEACR.",
        status='completed' if status['peak_calling'] else 'current' if status['alignment'] else 'pending',
        page="pages/23_peak_calling.py",
        details={
            "Supported callers": [
                "MACS2 (ChIP-seq, ATAC-seq)",
                "SEACR (CUT&Tag, CUT&RUN)"
            ],
            "Peak types": [
                "Narrow peaks (TFs, H3K4me3)",
                "Broad peaks (H3K27me3, H3K36me3)"
            ],
            "Time estimate": "~10-30 min per sample"
        }
    )

    st.markdown("### 📊 Core Analysis")

    # Step 2: Data Upload (or use peaks from peak calling)
    render_step(
        number=3,
        title="Load Peak Data",
        description="Upload peak files or use peaks from the peak calling step.",
        status='completed' if status['data'] else 'current' if status['project'] or status['peak_calling'] else 'pending',
        page="pages/01_data_project.py",
        details={
            "Supported formats": [
                "BED, narrowPeak, broadPeak",
                "MACS2 output files",
                "CSV/TSV with chr, start, end"
            ],
            "Time estimate": "~2-5 minutes"
        }
    )

    # Step 4: Quality Control
    render_step(
        number=4,
        title="Quality Control",
        description="Assess data quality with peak metrics, signal distribution, and sample correlation.",
        status='completed' if status['qc'] else 'current' if status['data'] else 'pending',
        page="pages/02_quality_control.py",
        details={
            "QC metrics": [
                "Peak count and width distribution",
                "Signal strength analysis",
                "Chromosomal distribution",
                "Sample correlation matrix"
            ],
            "Time estimate": "~2-3 minutes"
        }
    )

    # Step 5: Differential Analysis
    render_step(
        number=5,
        title="Differential Analysis",
        description="Identify peaks with significant differences between conditions using PyDESeq2.",
        status='completed' if status['analysis'] else 'current' if status['qc'] or status['data'] else 'pending',
        page="pages/03_differential.py",
        details={
            "Analysis options": [
                "Spike-in or RLE normalization",
                "Configurable FDR and fold-change thresholds",
                "Batch correction support",
                "Multiple comparison methods"
            ],
            "Time estimate": "~5-10 minutes"
        }
    )

    # Step 6: Results & Export
    render_step(
        number=6,
        title="Results & Export",
        description="Visualize results, generate reports, and export data for publication.",
        status='completed' if status['results'] else 'current' if status['analysis'] else 'pending',
        page="pages/07_reports.py",
        details={
            "Output formats": [
                "Interactive volcano plots",
                "Publication-ready figures",
                "CSV/Excel data tables",
                "PDF reports"
            ],
            "Time estimate": "~2-5 minutes"
        }
    )

    st.markdown("---")

    # Optional advanced analyses
    st.header("🔬 Optional: Advanced Analyses")
    st.markdown("After completing the core workflow, explore these advanced features:")

    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("### Super-Enhancers")
        st.caption("Identify super-enhancers using the ROSE algorithm")
        if st.button("Go to Super-Enhancers →", key="se"):
            st.switch_page("pages/12_super_enhancers.py")

    with col2:
        st.markdown("### Multi-Mark Integration")
        st.caption("Combine multiple histone marks for chromatin state analysis")
        if st.button("Go to Multi-Mark →", key="mm"):
            st.switch_page("pages/06_multimark.py")

    with col3:
        st.markdown("### Expression Integration")
        st.caption("Correlate epigenetic changes with gene expression")
        if st.button("Go to Expression →", key="expr"):
            st.switch_page("pages/09_expression.py")

    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("### TF Analysis")
        st.caption("Motif enrichment and TF binding analysis")
        if st.button("Go to TF Analysis →", key="tf"):
            st.switch_page("pages/08_tf_analysis.py")

    with col2:
        st.markdown("### GWAS Overlap")
        st.caption("Link peaks to disease-associated variants")
        if st.button("Go to GWAS →", key="gwas"):
            st.switch_page("pages/19_gwas_overlap.py")

    with col3:
        st.markdown("### ENCODE Comparison")
        st.caption("Compare with public reference epigenomes")
        if st.button("Go to ENCODE →", key="encode"):
            st.switch_page("pages/14_encode_integration.py")


def render_step(number, title, description, status, page, details):
    """Render a single workflow step."""

    # Status icons and colors
    if status == 'completed':
        icon = "✅"
        color = "#27ae60"
        status_text = "Completed"
    elif status == 'current':
        icon = "🔵"
        color = "#3498db"
        status_text = "Current Step"
    else:
        icon = "⚪"
        color = "#95a5a6"
        status_text = "Pending"

    # Create the step card
    with st.container():
        col1, col2 = st.columns([3, 1])

        with col1:
            st.markdown(f"### {icon} Step {number}: {title}")
            st.markdown(description)

            # Expandable details
            with st.expander("Details"):
                for key, values in details.items():
                    st.markdown(f"**{key}:**")
                    if isinstance(values, list):
                        for v in values:
                            st.markdown(f"- {v}")
                    else:
                        st.markdown(values)

        with col2:
            st.markdown(f"**Status:** {status_text}")

            if status == 'current':
                if st.button(f"Start Step {number} →", key=f"step_{number}", type="primary"):
                    st.switch_page(page)
            elif status == 'completed':
                if st.button(f"Review →", key=f"step_{number}"):
                    st.switch_page(page)
            else:
                st.button(f"Locked", key=f"step_{number}", disabled=True)

        st.markdown("---")


if __name__ == "__main__":
    main()
