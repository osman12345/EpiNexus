"""
Data & Project Page - Clean, Intuitive Design

Simplified workflow:
1. Create/Load Project
2. Upload Data
3. Configure Samples
"""

import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import json
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

# Import data manager
try:
    from frontend.components.data_manager import (
        DataManager, DataSource, generate_demo_peaks, generate_demo_samples
    )
    HAS_DATA_MANAGER = True
except ImportError:
    HAS_DATA_MANAGER = False

st.set_page_config(page_title="Data & Project - EpiNexus", page_icon="📁", layout="wide")


def main():
    st.title("📁 Data & Project")

    # Initialize session state
    if 'workflow_step' not in st.session_state:
        st.session_state.workflow_step = 1
    if 'samples' not in st.session_state:
        st.session_state.samples = []
    if 'using_demo_data' not in st.session_state:
        st.session_state.using_demo_data = True

    # Workflow progress indicator
    render_workflow_progress()

    st.markdown("---")

    # Main content based on workflow step
    step = st.session_state.workflow_step

    if step == 1:
        render_step1_project()
    elif step == 2:
        render_step2_upload()
    elif step == 3:
        render_step3_configure()
    else:
        render_summary()


def render_workflow_progress():
    """Show workflow progress as a horizontal stepper."""
    step = st.session_state.workflow_step

    steps = [
        ("1. Project", "Create or load"),
        ("2. Data", "Upload files"),
        ("3. Configure", "Set up samples"),
        ("4. Ready", "Start analysis")
    ]

    cols = st.columns(len(steps))

    for i, (title, subtitle) in enumerate(steps, 1):
        with cols[i-1]:
            if i < step:
                # Completed
                st.markdown(f"### ✅ {title}")
                st.caption(subtitle)
            elif i == step:
                # Current
                st.markdown(f"### 🔵 {title}")
                st.caption(f"**{subtitle}**")
            else:
                # Future
                st.markdown(f"### ⚪ {title}")
                st.caption(subtitle)


# =============================================================================
# STEP 1: Project Setup
# =============================================================================

def render_step1_project():
    """Step 1: Create or load a project."""
    st.header("Step 1: Project Setup")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("🆕 Create New Project")

        project_name = st.text_input(
            "Project Name",
            value=st.session_state.get('project_name', 'My_Project'),
            key="new_project_name"
        )

        assay_type = st.selectbox(
            "Assay Type",
            ["CUT&Tag", "CUT&RUN", "ChIP-seq", "ATAC-seq"],
            help="Select your experimental method"
        )

        target_type = st.selectbox(
            "Target",
            ["Histone Modifications", "Transcription Factors", "Chromatin Accessibility"],
            index=0 if assay_type != "ATAC-seq" else 2
        )

        genome = st.selectbox("Reference Genome", ["hg38", "hg19", "mm10", "mm39", "dm6"])

        if st.button("Create Project", type="primary", use_container_width=True):
            st.session_state.project_name = project_name
            st.session_state.assay_type = assay_type
            st.session_state.target_type = target_type
            st.session_state.selected_genome = genome
            st.session_state.project_created = datetime.now().isoformat()
            st.session_state.samples = []
            st.session_state.using_demo_data = False
            st.session_state.workflow_step = 2

            if HAS_DATA_MANAGER:
                DataManager.clear_data()

            st.rerun()

    with col2:
        st.subheader("📂 Load Existing Project")

        uploaded = st.file_uploader(
            "Upload .epinexus file",
            type=['json', 'epinexus'],
            key="load_project"
        )

        if uploaded:
            try:
                data = json.load(uploaded)
                st.session_state.project_name = data.get('name', 'Loaded Project')
                st.session_state.assay_type = data.get('assay_type', 'CUT&Tag')
                st.session_state.target_type = data.get('target_type', 'Histone Modifications')
                st.session_state.selected_genome = data.get('genome', 'hg38')
                st.session_state.samples = data.get('samples', [])
                st.session_state.using_demo_data = False
                st.session_state.workflow_step = 2

                st.success(f"✅ Loaded: **{st.session_state.project_name}**")
                st.rerun()
            except Exception as e:
                st.error(f"Error: {str(e)}")

        st.markdown("---")

        # Quick start with demo
        st.subheader("⚡ Quick Start")
        if st.button("Load Demo Data", use_container_width=True):
            st.session_state.project_name = "Demo_Project"
            st.session_state.assay_type = "CUT&Tag"
            st.session_state.target_type = "Histone Modifications"
            st.session_state.selected_genome = "hg38"
            st.session_state.using_demo_data = True

            if HAS_DATA_MANAGER:
                demo_peaks = generate_demo_peaks()
                DataManager.load_data('peaks', demo_peaks, DataSource.DEMO, "Demo Peaks")
                st.session_state.samples = generate_demo_samples()

            st.session_state.workflow_step = 4
            st.rerun()


# =============================================================================
# STEP 2: Upload Data
# =============================================================================

def render_step2_upload():
    """Step 2: Upload peak files or raw data."""
    st.header("Step 2: Upload Your Data")

    # Show current project info
    col1, col2, col3 = st.columns(3)
    col1.info(f"**Project:** {st.session_state.get('project_name', 'Untitled')}")
    col2.info(f"**Assay:** {st.session_state.get('assay_type', 'CUT&Tag')}")
    col3.info(f"**Genome:** {st.session_state.get('selected_genome', 'hg38')}")

    st.markdown("---")

    # Two options: Peaks or FASTQ
    upload_type = st.radio(
        "What data do you have?",
        ["📊 Peak files (BED/narrowPeak)", "🧬 Raw FASTQ files"],
        horizontal=True
    )

    if "Peak files" in upload_type:
        render_peak_upload()
    else:
        render_fastq_upload()

    # Navigation
    st.markdown("---")
    col1, col2, col3 = st.columns([1, 2, 1])

    with col1:
        if st.button("← Back"):
            st.session_state.workflow_step = 1
            st.rerun()

    with col3:
        # Check if data was uploaded
        has_data = False
        if HAS_DATA_MANAGER:
            peaks = DataManager.get_data('peaks')
            has_data = peaks is not None and len(peaks) > 0

        if has_data or st.session_state.get('fastq_samples'):
            if st.button("Next →", type="primary"):
                st.session_state.workflow_step = 3
                st.rerun()
        else:
            st.button("Next →", disabled=True, help="Upload data first")


def render_peak_upload():
    """Upload peak files."""
    st.subheader("Upload Peak Files")

    col1, col2 = st.columns([2, 1])

    with col1:
        uploaded_files = st.file_uploader(
            "Select peak files",
            type=['bed', 'narrowPeak', 'broadPeak', 'csv', 'tsv'],
            accept_multiple_files=True,
            help="BED, narrowPeak, broadPeak, or CSV with chr/start/end columns"
        )

        if uploaded_files:
            all_peaks = []
            for f in uploaded_files:
                try:
                    # Parse file
                    if f.name.endswith('.csv'):
                        df = pd.read_csv(f)
                    elif f.name.endswith('.tsv'):
                        df = pd.read_csv(f, sep='\t')
                    else:
                        # BED/narrowPeak format
                        df = pd.read_csv(f, sep='\t', header=None)
                        cols = ['chr', 'start', 'end', 'name', 'score', 'strand',
                               'signalValue', 'pValue', 'qValue', 'summit']
                        df.columns = cols[:len(df.columns)]

                    df['source_file'] = f.name
                    all_peaks.append(df)
                    st.success(f"✅ **{f.name}**: {len(df):,} peaks")

                except Exception as e:
                    st.error(f"❌ {f.name}: {str(e)}")

            if all_peaks and HAS_DATA_MANAGER:
                combined = pd.concat(all_peaks, ignore_index=True)
                DataManager.load_data('peaks', combined, DataSource.USER_UPLOAD, "User Peaks")
                st.session_state.using_demo_data = False

                # Auto-detect samples from filenames
                detected_samples = []
                for f in uploaded_files:
                    sample_id = Path(f.name).stem
                    detected_samples.append({
                        'SampleID': sample_id,
                        'Condition': 'Control' if 'ctrl' in sample_id.lower() or 'wt' in sample_id.lower() else 'Treatment',
                        'Factor': 'H3K27ac',
                        'PeakFile': f.name
                    })
                st.session_state.samples = detected_samples

    with col2:
        st.markdown("**Supported formats:**")
        st.markdown("""
        - BED (chr, start, end)
        - narrowPeak (MACS2)
        - broadPeak (MACS2)
        - CSV/TSV with headers
        """)

        st.markdown("---")

        st.markdown("**Peak stats:**")
        if HAS_DATA_MANAGER:
            peaks = DataManager.get_data('peaks')
            if peaks is not None and len(peaks) > 0:
                st.metric("Total Peaks", f"{len(peaks):,}")
                if 'chr' in peaks.columns:
                    st.metric("Chromosomes", peaks['chr'].nunique())


def render_fastq_upload():
    """FASTQ pipeline setup."""
    st.subheader("FASTQ Processing Pipeline")

    st.warning("⚠️ FASTQ processing requires external tools (Bowtie2, MACS2, etc.)")

    # Sample sheet upload
    sample_sheet = st.file_uploader(
        "Upload sample sheet (CSV)",
        type=['csv', 'tsv'],
        help="Columns: SampleID, Condition, R1_path, R2_path (optional)"
    )

    if sample_sheet:
        try:
            df = pd.read_csv(sample_sheet)
            st.session_state.fastq_samples = df.to_dict('records')
            st.dataframe(df, use_container_width=True)
            st.success(f"✅ {len(df)} samples loaded")
        except Exception as e:
            st.error(str(e))

    # Pipeline options
    with st.expander("Pipeline Options"):
        col1, col2 = st.columns(2)
        with col1:
            st.selectbox("Aligner", ["Bowtie2", "BWA-MEM"])
            st.checkbox("Paired-end", value=True)

            if st.session_state.get('assay_type') in ['CUT&Tag', 'CUT&RUN']:
                st.checkbox("Spike-in normalization", value=True)

        with col2:
            st.selectbox("Peak Caller", ["MACS2", "SEACR", "SICER"])
            st.selectbox("Peak Type", ["Narrow", "Broad"])


# =============================================================================
# STEP 3: Configure Samples
# =============================================================================

def render_step3_configure():
    """Step 3: Configure sample metadata."""
    st.header("Step 3: Configure Samples")

    samples = st.session_state.get('samples', [])

    if not samples:
        st.warning("No samples detected. Add samples manually below.")

    # Editable sample table
    st.subheader("Sample Metadata")

    # Create DataFrame for editing
    if samples:
        df = pd.DataFrame(samples)
    else:
        df = pd.DataFrame(columns=['SampleID', 'Condition', 'Factor', 'Replicate'])

    # Add sample form
    with st.expander("➕ Add New Sample", expanded=not samples):
        col1, col2, col3, col4 = st.columns(4)

        with col1:
            new_id = st.text_input("Sample ID", key="add_sample_id")
        with col2:
            new_cond = st.selectbox("Condition", ["Control", "Treatment", "KO", "WT"], key="add_cond")
        with col3:
            new_factor = st.selectbox(
                "Factor/Mark",
                ["H3K27ac", "H3K4me3", "H3K4me1", "H3K27me3", "H3K36me3", "H3K9me3", "Other"],
                key="add_factor"
            )
        with col4:
            new_rep = st.number_input("Replicate", 1, 10, 1, key="add_rep")

        if st.button("Add Sample"):
            if new_id:
                new_sample = {
                    'SampleID': new_id,
                    'Condition': new_cond,
                    'Factor': new_factor,
                    'Replicate': new_rep
                }
                st.session_state.samples.append(new_sample)
                st.rerun()

    # Show and edit existing samples
    if samples:
        edited_df = st.data_editor(
            df,
            num_rows="dynamic",
            use_container_width=True,
            column_config={
                "SampleID": st.column_config.TextColumn("Sample ID", required=True),
                "Condition": st.column_config.SelectboxColumn(
                    "Condition",
                    options=["Control", "Treatment", "KO", "WT"],
                    required=True
                ),
                "Factor": st.column_config.SelectboxColumn(
                    "Factor/Mark",
                    options=["H3K27ac", "H3K4me3", "H3K4me1", "H3K27me3", "H3K36me3", "H3K9me3"],
                    required=True
                ),
                "Replicate": st.column_config.NumberColumn("Replicate", min_value=1, max_value=10)
            }
        )

        # Update session state with edited data
        st.session_state.samples = edited_df.to_dict('records')

        # Group summary
        st.markdown("---")
        st.subheader("Group Summary")

        conditions = edited_df.groupby('Condition').size()
        col1, col2 = st.columns(2)
        for i, (cond, count) in enumerate(conditions.items()):
            col = col1 if i % 2 == 0 else col2
            col.metric(cond, f"{count} samples")

    # Navigation
    st.markdown("---")
    col1, col2, col3 = st.columns([1, 2, 1])

    with col1:
        if st.button("← Back"):
            st.session_state.workflow_step = 2
            st.rerun()

    with col3:
        if samples and len(samples) >= 2:
            if st.button("Finish Setup →", type="primary"):
                st.session_state.workflow_step = 4
                st.session_state.using_demo_data = False
                st.rerun()
        else:
            st.button("Finish Setup →", disabled=True, help="Add at least 2 samples")


# =============================================================================
# STEP 4: Summary / Ready
# =============================================================================

def render_summary():
    """Summary view - ready for analysis."""
    st.header("✅ Project Ready")

    # Project summary
    col1, col2, col3, col4 = st.columns(4)

    col1.metric("Project", st.session_state.get('project_name', 'Untitled'))
    col2.metric("Assay", st.session_state.get('assay_type', 'CUT&Tag'))
    col3.metric("Samples", len(st.session_state.get('samples', [])))

    if HAS_DATA_MANAGER:
        peaks = DataManager.get_data('peaks')
        n_peaks = len(peaks) if peaks is not None else 0
        col4.metric("Peaks", f"{n_peaks:,}")

    # Data status
    is_demo = st.session_state.get('using_demo_data', True)
    if is_demo:
        st.info("📌 **Demo Mode** - Using sample data for demonstration")
    else:
        st.success("✅ **Your Data** - Ready for analysis")

    st.markdown("---")

    # Sample table
    st.subheader("Samples")
    samples = st.session_state.get('samples', [])
    if samples:
        st.dataframe(pd.DataFrame(samples), use_container_width=True, hide_index=True)

    st.markdown("---")

    # Quick actions
    st.subheader("Next Steps")

    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("### 📊 Quality Control")
        st.caption("Check sample quality metrics")
        if st.button("Go to QC →", use_container_width=True):
            st.switch_page("pages/02_quality_control.py")

    with col2:
        st.markdown("### 📈 Differential Analysis")
        st.caption("Find differential peaks")
        if st.button("Go to Differential →", use_container_width=True):
            st.switch_page("pages/03_differential.py")

    with col3:
        st.markdown("### 💾 Save Project")
        st.caption("Download project file")

        project_data = {
            'name': st.session_state.get('project_name', 'Untitled'),
            'assay_type': st.session_state.get('assay_type', 'CUT&Tag'),
            'target_type': st.session_state.get('target_type', 'Histone Modifications'),
            'genome': st.session_state.get('selected_genome', 'hg38'),
            'samples': st.session_state.get('samples', []),
            'saved': datetime.now().isoformat()
        }

        st.download_button(
            "💾 Download .epinexus",
            json.dumps(project_data, indent=2),
            f"{project_data['name']}.epinexus",
            "application/json",
            use_container_width=True
        )

    # Edit options
    st.markdown("---")
    col1, col2 = st.columns([1, 4])
    with col1:
        if st.button("✏️ Edit Project"):
            st.session_state.workflow_step = 1
            st.rerun()


if __name__ == "__main__":
    main()
