"""
Unified Data & Project Page

Combines:
- Project management (new/load/save)
- Data input options (FASTQ pipeline, peak upload)
- Sample and assay configuration
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
        DataManager, DataSource, render_data_status_indicator,
        generate_demo_peaks, generate_demo_samples
    )
    HAS_DATA_MANAGER = True
except ImportError:
    HAS_DATA_MANAGER = False

st.set_page_config(page_title="Data & Project - EpiNexus", page_icon="📁", layout="wide")


def main():
    st.title("📁 Data & Project")

    # Current status bar
    render_status_bar()

    # Main tabs
    tab1, tab2, tab3, tab4 = st.tabs([
        "🆕 Project",
        "⚙️ From FASTQ",
        "📤 From Peaks",
        "📋 Samples & Config"
    ])

    with tab1:
        render_project_tab()

    with tab2:
        render_fastq_tab()

    with tab3:
        render_peaks_tab()

    with tab4:
        render_config_tab()


def render_status_bar():
    """Show current data/project status."""
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        project = st.session_state.get('project_name', 'Untitled')
        st.metric("Project", project)

    with col2:
        n_samples = len(st.session_state.get('samples', []))
        st.metric("Samples", n_samples)

    with col3:
        if HAS_DATA_MANAGER:
            peaks = DataManager.get_data('peaks')
            n_peaks = len(peaks) if peaks is not None else 0
        else:
            n_peaks = 0
        st.metric("Peaks", f"{n_peaks:,}" if n_peaks > 0 else "—")

    with col4:
        status = "📌 Demo" if st.session_state.get('using_demo_data', True) else "✅ Your Data"
        st.metric("Status", status)

    st.markdown("---")


# =============================================================================
# TAB 1: Project Management
# =============================================================================

def render_project_tab():
    """Project creation, loading, and saving."""
    st.header("Project Management")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("🆕 New Project")

        project_name = st.text_input("Project Name", value="My_EpiNexus_Project")
        description = st.text_area("Description (optional)", height=80)

        if st.button("Create New Project", type="primary", use_container_width=True):
            st.session_state.project_name = project_name
            st.session_state.project_description = description
            st.session_state.project_created = datetime.now().isoformat()
            st.session_state.samples = []
            st.session_state.using_demo_data = True

            if HAS_DATA_MANAGER:
                DataManager.clear_data()

            st.success(f"✅ Created new project: **{project_name}**")
            st.rerun()

    with col2:
        st.subheader("📂 Load Existing")

        uploaded_project = st.file_uploader(
            "Load project file (.epinexus or .json)",
            type=['json', 'epinexus']
        )

        if uploaded_project:
            try:
                project_data = json.load(uploaded_project)
                st.session_state.project_name = project_data.get('name', 'Loaded Project')
                st.session_state.project_description = project_data.get('description', '')
                st.session_state.samples = project_data.get('samples', [])
                st.session_state.config = project_data.get('config', {})
                st.session_state.using_demo_data = False

                st.success(f"✅ Loaded project: **{st.session_state.project_name}**")
            except Exception as e:
                st.error(f"Error loading project: {str(e)}")

    st.markdown("---")

    # Save current project
    st.subheader("💾 Save Current Project")

    col1, col2 = st.columns([2, 1])

    with col1:
        if st.button("📥 Download Project File", use_container_width=True):
            project_data = {
                'name': st.session_state.get('project_name', 'Untitled'),
                'description': st.session_state.get('project_description', ''),
                'created': st.session_state.get('project_created', datetime.now().isoformat()),
                'saved': datetime.now().isoformat(),
                'samples': st.session_state.get('samples', []),
                'config': st.session_state.get('config', {}),
                'genome': st.session_state.get('selected_genome', 'hg38')
            }

            st.download_button(
                "💾 Download .epinexus file",
                json.dumps(project_data, indent=2),
                f"{project_data['name']}.epinexus",
                "application/json"
            )

    with col2:
        if st.button("🗑️ Clear All Data", type="secondary"):
            if HAS_DATA_MANAGER:
                DataManager.clear_data()
            st.session_state.samples = []
            st.session_state.using_demo_data = True
            st.rerun()


# =============================================================================
# TAB 2: FASTQ Pipeline
# =============================================================================

def render_fastq_tab():
    """FASTQ preprocessing pipeline."""
    st.header("Start from FASTQ Files")

    st.info("🧬 Process raw sequencing data through the full pipeline: QC → Trimming → Alignment → Peak Calling")

    # Step 1: Technique & Target
    st.subheader("1️⃣ Experiment Setup")

    col1, col2, col3 = st.columns(3)

    with col1:
        technique = st.selectbox(
            "Technique",
            ["ChIP-seq", "CUT&Tag", "CUT&RUN", "ATAC-seq"],
            help="How was the experiment done?"
        )
        st.session_state.technique = technique

    with col2:
        if technique == "ATAC-seq":
            st.info("🔓 Chromatin accessibility")
            target = "Accessibility"
        else:
            target = st.selectbox(
                "Target",
                ["Histone Modifications", "Transcription Factors"],
                help="What are you profiling?"
            )
        st.session_state.target = target

    with col3:
        genome = st.selectbox(
            "Reference Genome",
            ["hg38", "hg19", "mm10", "mm39", "dm6"]
        )
        st.session_state.selected_genome = genome

    # Step 2: Sample Sheet
    st.markdown("---")
    st.subheader("2️⃣ Sample Sheet")

    sample_method = st.radio(
        "How would you like to define samples?",
        ["Upload CSV", "Manual entry", "Auto-detect from folder"],
        horizontal=True
    )

    if sample_method == "Upload CSV":
        st.markdown("""
        **Required columns:** SampleID, Condition, R1_path, R2_path (if paired-end)

        **Optional:** Replicate, Control (for ChIP-seq input samples)
        """)

        sample_csv = st.file_uploader("Upload sample sheet", type=['csv', 'tsv'])
        if sample_csv:
            df = pd.read_csv(sample_csv)
            st.session_state.samples = df.to_dict('records')
            st.session_state.using_demo_data = False
            st.dataframe(df, use_container_width=True, hide_index=True)

    elif sample_method == "Manual entry":
        render_manual_sample_entry()

    else:  # Auto-detect
        folder_path = st.text_input("Path to FASTQ folder")
        if folder_path and st.button("🔍 Scan for FASTQ files"):
            st.info("Would scan folder for .fastq.gz files and auto-generate sample sheet")

    # Step 3: Pipeline Options
    st.markdown("---")
    st.subheader("3️⃣ Pipeline Options")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Alignment**")
        aligner = st.selectbox("Aligner", ["Bowtie2", "BWA-MEM", "STAR"])
        paired_end = st.checkbox("Paired-end", value=True)

        if technique in ["CUT&Tag", "CUT&RUN"]:
            spike_in = st.checkbox("Spike-in normalization", value=True)
            if spike_in:
                spike_genome = st.selectbox("Spike-in", ["E. coli", "S. cerevisiae"])

    with col2:
        st.markdown("**Peak Calling**")
        if target == "Histone Modifications":
            peak_caller = st.selectbox("Peak caller", ["MACS2", "SICER", "SEACR"])
            peak_type = st.selectbox("Peak type", ["Narrow", "Broad", "Auto-detect"])
        elif target == "Transcription Factors":
            peak_caller = st.selectbox("Peak caller", ["MACS2", "HOMER"])
            st.caption("Using narrow peaks for TF binding sites")
        else:  # ATAC-seq
            peak_caller = st.selectbox("Peak caller", ["MACS2", "Genrich"])
            st.caption("Optimized for accessibility data")

    # Run button
    st.markdown("---")

    col1, col2 = st.columns([1, 2])
    with col1:
        if st.button("▶️ Run Pipeline", type="primary", use_container_width=True):
            st.warning("⚠️ Pipeline execution requires external tools (Bowtie2, MACS2, etc.)")
            with st.spinner("This would run the full preprocessing pipeline..."):
                st.session_state.using_demo_data = False
                st.success("Pipeline would process samples and generate peaks")

    with col2:
        if st.button("📋 Show Commands (Dry Run)", use_container_width=True):
            st.code(f"""
# Example commands that would be run:
bowtie2 -x {genome} -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz | samtools view -bS - > sample.bam
samtools sort sample.bam -o sample.sorted.bam
samtools index sample.sorted.bam
macs2 callpeak -t sample.sorted.bam -f BAMPE -g hs -n sample --outdir peaks/
            """, language="bash")


# =============================================================================
# TAB 3: Upload Peaks
# =============================================================================

def render_peaks_tab():
    """Upload existing peak files."""
    st.header("Upload Peak Files")

    st.info("📤 Already have peak files? Upload them directly to skip preprocessing.")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Peak Files")

        st.markdown("""
        **Supported formats:**
        - BED (.bed)
        - narrowPeak (.narrowPeak)
        - broadPeak (.broadPeak)
        - MACS2 output (.xls)
        - CSV/TSV with chr, start, end columns
        """)

        uploaded_peaks = st.file_uploader(
            "Upload peak files",
            type=['bed', 'narrowPeak', 'broadPeak', 'xls', 'csv', 'tsv'],
            accept_multiple_files=True,
            key='peak_upload'
        )

        if uploaded_peaks:
            for f in uploaded_peaks:
                try:
                    if f.name.endswith(('.csv', '.tsv')):
                        sep = '\t' if f.name.endswith('.tsv') else ','
                        df = pd.read_csv(f, sep=sep)
                    else:
                        df = pd.read_csv(f, sep='\t', header=None)
                        if len(df.columns) >= 3:
                            df.columns = ['chr', 'start', 'end'] + [f'col_{i}' for i in range(3, len(df.columns))]

                    st.success(f"✅ **{f.name}**: {len(df):,} peaks")
                    st.session_state.using_demo_data = False

                    if HAS_DATA_MANAGER:
                        DataManager.load_data(
                            data_type='peaks',
                            data=df,
                            source=DataSource.USER_UPLOAD,
                            name=f.name
                        )

                except Exception as e:
                    st.error(f"Error: {f.name} - {str(e)}")

    with col2:
        st.subheader("Signal Tracks (Optional)")

        st.markdown("""
        **For visualization:**
        - BigWig (.bw, .bigwig)
        - bedGraph (.bedgraph)
        """)

        uploaded_bigwig = st.file_uploader(
            "Upload signal tracks",
            type=['bw', 'bigwig', 'bedgraph'],
            accept_multiple_files=True,
            key='bigwig_upload'
        )

        if uploaded_bigwig:
            for f in uploaded_bigwig:
                st.success(f"✅ **{f.name}** ready for visualization")

        st.markdown("---")

        # Demo data option
        st.subheader("Or Use Demo Data")
        if st.button("📌 Load Demo Peaks", use_container_width=True):
            if HAS_DATA_MANAGER:
                demo_peaks = generate_demo_peaks()
                DataManager.load_data(
                    data_type='peaks',
                    data=demo_peaks,
                    source=DataSource.DEMO,
                    name="Demo Peaks"
                )
            st.session_state.using_demo_data = True
            st.success("Demo peaks loaded!")
            st.rerun()


# =============================================================================
# TAB 4: Samples & Configuration
# =============================================================================

def render_config_tab():
    """Sample sheet and analysis configuration."""
    st.header("Samples & Configuration")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("📋 Sample Sheet")

        # Show current samples
        samples = st.session_state.get('samples', [])

        if samples:
            st.dataframe(pd.DataFrame(samples), use_container_width=True, hide_index=True)

            if st.button("🗑️ Clear Samples"):
                st.session_state.samples = []
                st.rerun()
        else:
            st.info("No samples defined yet. Add samples below or upload in the FASTQ/Peaks tabs.")

        # Quick add
        st.markdown("**Quick Add Sample:**")
        with st.form("add_sample"):
            c1, c2 = st.columns(2)
            with c1:
                sample_id = st.text_input("Sample ID")
                condition = st.selectbox("Condition", ["Control", "Treatment", "Other"])
            with c2:
                factor = st.text_input("Factor/Mark", "H3K27ac")
                replicate = st.number_input("Replicate", 1, 10, 1)

            if st.form_submit_button("Add Sample"):
                if sample_id:
                    new_sample = {
                        'SampleID': sample_id,
                        'Condition': condition,
                        'Factor': factor,
                        'Replicate': replicate
                    }
                    if 'samples' not in st.session_state:
                        st.session_state.samples = []
                    st.session_state.samples.append(new_sample)
                    st.session_state.using_demo_data = False
                    st.rerun()

    with col2:
        st.subheader("⚙️ Analysis Settings")

        # Technique & Target (if not already set)
        technique = st.selectbox(
            "Technique",
            ["ChIP-seq", "CUT&Tag", "CUT&RUN", "ATAC-seq"],
            key="config_technique"
        )

        if technique != "ATAC-seq":
            target = st.selectbox(
                "Target Type",
                ["Histone Modifications", "Transcription Factors"],
                key="config_target"
            )

            if target == "Histone Modifications":
                marks = st.multiselect(
                    "Histone Marks",
                    ["H3K4me3", "H3K4me1", "H3K27ac", "H3K27me3", "H3K36me3",
                     "H3K9me3", "H3K9ac", "H3K79me2"],
                    default=["H3K27ac"]
                )
                st.session_state.histone_marks = marks
            else:
                tfs = st.text_input("TF Names (comma-separated)", "MYC, P53")
                st.session_state.transcription_factors = [t.strip() for t in tfs.split(',')]

        st.markdown("---")

        # Reference genome
        genome = st.selectbox(
            "Reference Genome",
            ["hg38", "hg19", "mm10", "mm39", "dm6"],
            key="config_genome"
        )
        st.session_state.selected_genome = genome

        # Save config
        if st.button("💾 Save Configuration", type="primary", use_container_width=True):
            st.session_state.config = {
                'technique': technique,
                'target': target if technique != "ATAC-seq" else "accessibility",
                'genome': genome
            }
            st.success("✅ Configuration saved!")


def render_manual_sample_entry():
    """Form for manual sample entry."""
    st.markdown("**Add samples one by one:**")

    with st.form("manual_sample"):
        col1, col2, col3, col4 = st.columns(4)

        with col1:
            sample_id = st.text_input("Sample ID", key="man_id")
        with col2:
            condition = st.selectbox("Condition", ["Control", "Treatment"], key="man_cond")
        with col3:
            r1_path = st.text_input("R1 FASTQ path", key="man_r1")
        with col4:
            r2_path = st.text_input("R2 FASTQ path", key="man_r2")

        submitted = st.form_submit_button("Add Sample")

        if submitted and sample_id:
            new_sample = {
                'SampleID': sample_id,
                'Condition': condition,
                'R1': r1_path,
                'R2': r2_path
            }
            if 'samples' not in st.session_state:
                st.session_state.samples = []
            st.session_state.samples.append(new_sample)
            st.session_state.using_demo_data = False
            st.success(f"Added: {sample_id}")

    # Show current samples
    if st.session_state.get('samples'):
        st.markdown("**Current samples:**")
        st.dataframe(pd.DataFrame(st.session_state.samples), use_container_width=True, hide_index=True)


if __name__ == "__main__":
    main()
