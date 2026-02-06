"""
Batch Processing Page

Queue and manage multiple analysis jobs.
"""

import streamlit as st
import pandas as pd
import time
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from app.core.batch_processor import BatchProcessor, JobStatus

st.set_page_config(
    page_title="Batch Processing - EpiNexus",
    page_icon="⚙️",
    layout="wide"
)

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
            <div style="font-size: 3rem; margin-bottom: 1rem;">⚙️</div>
            <h2 style="color: #6c757d; margin-bottom: 0.5rem;">No Data Loaded</h2>
            <p style="color: #6c757d; font-size: 1.1rem;">
                Upload your data to run batch analysis jobs.
            </p>
        </div>
        """, unsafe_allow_html=True)
        st.markdown("")
        if st.button("📁 Go to Data & Project", type="primary", use_container_width=True):
            st.switch_page("pages/01_data_project.py")
        st.markdown("")
        st.markdown("**Batch processing features:**")
        st.markdown("- Parallel sample processing")
        st.markdown("- Job queue management")
        st.markdown("- Progress monitoring")


@st.cache_resource
def get_batch_processor():
    return BatchProcessor(max_workers=4, jobs_dir="batch_jobs")


def main():
    st.title("⚙️ Batch Processing")
    st.markdown("Queue and process multiple samples or analyses in parallel.")

    # Check if data is loaded
    if not has_data():
        render_empty_state()
        return

    bp = get_batch_processor()

    # Tabs
    tab1, tab2, tab3, tab4 = st.tabs([
        "➕ Create Jobs",
        "📋 Job Queue",
        "📊 Results",
        "⚙️ Settings"
    ])

    with tab1:
        render_create_jobs(bp)

    with tab2:
        render_job_queue(bp)

    with tab3:
        render_results(bp)

    with tab4:
        render_settings(bp)


def render_create_jobs(bp: BatchProcessor):
    """Create new batch jobs."""
    st.header("Create Batch Jobs")

    job_type = st.selectbox(
        "Analysis Type",
        [
            "differential_analysis",
            "qc_analysis",
            "peak_annotation",
            "super_enhancer"
        ],
        format_func=lambda x: {
            "differential_analysis": "🔬 Differential Analysis",
            "qc_analysis": "✅ Quality Control",
            "peak_annotation": "🏷️ Peak Annotation",
            "super_enhancer": "⭐ Super-Enhancer Detection"
        }.get(x, x)
    )

    st.markdown("---")

    if job_type == "differential_analysis":
        render_differential_batch(bp)
    elif job_type == "qc_analysis":
        render_qc_batch(bp)
    elif job_type == "peak_annotation":
        render_annotation_batch(bp)
    elif job_type == "super_enhancer":
        render_se_batch(bp)


def render_differential_batch(bp: BatchProcessor):
    """Create differential analysis batch."""
    st.subheader("Differential Analysis Batch")

    st.markdown("Run multiple comparisons across different conditions or histone marks.")

    # Upload sample sheet
    st.markdown("**Define Comparisons:**")

    n_comparisons = st.number_input("Number of comparisons", 1, 20, 3)

    comparisons = []
    for i in range(n_comparisons):
        with st.expander(f"Comparison {i+1}", expanded=(i == 0)):
            col1, col2, col3 = st.columns(3)

            with col1:
                name = st.text_input("Name", f"Comparison_{i+1}", key=f"name_{i}")
            with col2:
                control = st.text_input("Control group", "Control", key=f"ctrl_{i}")
            with col3:
                treatment = st.text_input("Treatment group", "Treatment", key=f"treat_{i}")

            mark = st.selectbox(
                "Histone mark",
                ["H3K27ac", "H3K4me3", "H3K27me3", "H3K4me1"],
                key=f"mark_{i}"
            )

            comparisons.append({
                "name": name,
                "control": control,
                "treatment": treatment,
                "mark": mark
            })

    if st.button("🚀 Submit Batch", type="primary"):
        jobs = []
        for comp in comparisons:
            jobs.append({
                "name": f"Diff: {comp['name']} ({comp['mark']})",
                "job_type": "differential_analysis",
                "params": comp
            })

        created = bp.submit_batch(jobs, "Differential Analysis")
        st.success(f"Created {len(created)} jobs!")

        # Start processing
        bp.start_processing([j.id for j in created])
        st.info("Jobs submitted to queue. Check the Job Queue tab for progress.")


def render_qc_batch(bp: BatchProcessor):
    """Create QC analysis batch."""
    st.subheader("Quality Control Batch")

    uploaded_files = st.file_uploader(
        "Upload sample files (BAM or peak files)",
        type=['bam', 'bed', 'narrowPeak'],
        accept_multiple_files=True
    )

    if uploaded_files:
        st.info(f"Selected {len(uploaded_files)} files for QC")

        # QC parameters
        col1, col2 = st.columns(2)
        with col1:
            min_frip = st.slider("Minimum FRiP threshold", 0.05, 0.5, 0.1)
        with col2:
            min_reads = st.number_input("Minimum mapped reads (M)", 1, 100, 10)

        if st.button("🚀 Run QC Batch", type="primary"):
            jobs = []
            for f in uploaded_files:
                jobs.append({
                    "name": f"QC: {f.name}",
                    "job_type": "qc_analysis",
                    "params": {
                        "file": f.name,
                        "min_frip": min_frip,
                        "min_reads": min_reads * 1e6
                    }
                })

            created = bp.submit_batch(jobs, "QC Analysis")
            bp.start_processing([j.id for j in created])
            st.success(f"Started QC for {len(created)} samples!")
    else:
        # Demo mode
        st.markdown("**Or use demo samples:**")
        n_samples = st.slider("Number of demo samples", 2, 20, 8)

        if st.button("🚀 Run Demo QC Batch", type="primary"):
            jobs = []
            for i in range(n_samples):
                jobs.append({
                    "name": f"QC: Sample_{i+1}",
                    "job_type": "qc_analysis",
                    "params": {"n_samples": 1, "sample_id": i+1}
                })

            created = bp.submit_batch(jobs, "Demo QC")
            bp.start_processing([j.id for j in created])
            st.success(f"Started QC for {len(created)} demo samples!")


def render_annotation_batch(bp: BatchProcessor):
    """Create annotation batch."""
    st.subheader("Peak Annotation Batch")

    st.markdown("Annotate multiple peak sets with genomic features and nearby genes.")

    peak_files = st.file_uploader(
        "Upload peak files",
        type=['bed', 'narrowPeak', 'csv'],
        accept_multiple_files=True
    )

    # Annotation options
    col1, col2 = st.columns(2)
    with col1:
        genome = st.selectbox("Genome", ["hg38", "hg19", "mm10"])
        tss_distance = st.slider("TSS distance (kb)", 1, 100, 10)
    with col2:
        annotation_db = st.selectbox("Gene database", ["GENCODE", "RefSeq", "UCSC"])
        include_go = st.checkbox("Include GO enrichment", value=True)

    n_files = len(peak_files) if peak_files else 5

    if st.button("🚀 Run Annotation Batch", type="primary"):
        jobs = []
        for i in range(n_files):
            name = peak_files[i].name if peak_files else f"PeakSet_{i+1}"
            jobs.append({
                "name": f"Annotate: {name}",
                "job_type": "peak_annotation",
                "params": {
                    "genome": genome,
                    "tss_distance": tss_distance * 1000,
                    "database": annotation_db,
                    "go_enrichment": include_go
                }
            })

        created = bp.submit_batch(jobs, "Annotation")
        bp.start_processing([j.id for j in created])
        st.success(f"Started annotation for {len(created)} peak sets!")


def render_se_batch(bp: BatchProcessor):
    """Create super-enhancer batch."""
    st.subheader("Super-Enhancer Detection Batch")

    st.markdown("Run SE detection across multiple samples or conditions.")

    n_samples = st.number_input("Number of samples", 1, 50, 6)

    col1, col2 = st.columns(2)
    with col1:
        stitch = st.slider("Stitch distance (kb)", 5, 25, 12)
    with col2:
        tss_excl = st.slider("TSS exclusion (kb)", 0, 5, 2)

    if st.button("🚀 Run SE Detection Batch", type="primary"):
        jobs = []
        for i in range(n_samples):
            jobs.append({
                "name": f"SE: Sample_{i+1}",
                "job_type": "super_enhancer",
                "params": {
                    "stitch_distance": stitch * 1000,
                    "tss_exclusion": tss_excl * 1000,
                    "sample_id": i + 1
                }
            })

        created = bp.submit_batch(jobs, "Super-Enhancers")
        bp.start_processing([j.id for j in created])
        st.success(f"Started SE detection for {len(created)} samples!")


def render_job_queue(bp: BatchProcessor):
    """Display job queue."""
    st.header("Job Queue")

    # Refresh button
    col1, col2, col3 = st.columns([1, 1, 2])
    with col1:
        if st.button("🔄 Refresh"):
            st.rerun()
    with col2:
        auto_refresh = st.checkbox("Auto-refresh (5s)")

    # Queue status
    status = bp.get_queue_status()
    col1, col2, col3, col4, col5 = st.columns(5)
    col1.metric("Pending", status['pending'])
    col2.metric("Running", status['running'])
    col3.metric("Completed", status['completed'])
    col4.metric("Failed", status['failed'])
    col5.metric("Cancelled", status['cancelled'])

    st.markdown("---")

    # Job list
    jobs = bp.get_all_jobs()

    if not jobs:
        st.info("No jobs in queue. Create jobs in the 'Create Jobs' tab.")
        return

    for job in jobs:
        status_emoji = {
            JobStatus.PENDING: "⏳",
            JobStatus.RUNNING: "🔄",
            JobStatus.COMPLETED: "✅",
            JobStatus.FAILED: "❌",
            JobStatus.CANCELLED: "🚫"
        }

        with st.container():
            col1, col2, col3, col4 = st.columns([3, 1, 1, 1])

            with col1:
                st.markdown(f"**{status_emoji[job.status]} {job.name}**")
                st.caption(f"ID: {job.id} | Type: {job.job_type}")

                if job.status == JobStatus.RUNNING:
                    st.progress(job.progress / 100, text=job.message)

            with col2:
                st.markdown(f"**{job.status.value}**")
                if job.progress > 0:
                    st.caption(f"{job.progress:.0f}%")

            with col3:
                st.caption(f"Created: {job.created_at[:16]}")
                if job.completed_at:
                    st.caption(f"Done: {job.completed_at[:16]}")

            with col4:
                if job.status == JobStatus.PENDING:
                    if st.button("Cancel", key=f"cancel_{job.id}"):
                        bp.cancel_job(job.id)
                        st.rerun()

                if job.status == JobStatus.FAILED:
                    if st.button("Retry", key=f"retry_{job.id}"):
                        job.status = JobStatus.PENDING
                        bp.start_processing([job.id])
                        st.rerun()

            st.markdown("---")

    # Auto-refresh
    if auto_refresh:
        time.sleep(5)
        st.rerun()


def render_results(bp: BatchProcessor):
    """Display job results."""
    st.header("Batch Results")

    completed_jobs = [j for j in bp.get_all_jobs() if j.status == JobStatus.COMPLETED]

    if not completed_jobs:
        st.info("No completed jobs yet.")
        return

    # Results summary
    st.subheader("Completed Jobs")

    for job in completed_jobs:
        with st.expander(f"✅ {job.name}", expanded=False):
            col1, col2 = st.columns(2)

            with col1:
                st.markdown(f"**Job ID:** {job.id}")
                st.markdown(f"**Type:** {job.job_type}")
                st.markdown(f"**Completed:** {job.completed_at[:16]}")

            with col2:
                if job.result:
                    st.markdown("**Results:**")
                    for key, value in job.result.items():
                        st.write(f"• {key}: {value}")

    # Aggregate results
    st.subheader("Aggregate Summary")

    if any(j.job_type == "differential_analysis" for j in completed_jobs):
        diff_jobs = [j for j in completed_jobs if j.job_type == "differential_analysis"]

        agg_data = []
        for j in diff_jobs:
            if j.result:
                agg_data.append({
                    "Job": j.name,
                    "Total Peaks": j.result.get("total_peaks", 0),
                    "Differential": j.result.get("differential_peaks", 0),
                    "Up": j.result.get("up_regulated", 0),
                    "Down": j.result.get("down_regulated", 0)
                })

        if agg_data:
            st.dataframe(pd.DataFrame(agg_data), use_container_width=True)


def render_settings(bp: BatchProcessor):
    """Batch processing settings."""
    st.header("Settings")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Processing")
        max_workers = st.slider("Max parallel jobs", 1, 16, 4)
        st.caption("Number of jobs to run simultaneously")

    with col2:
        st.subheader("Storage")
        st.markdown(f"**Jobs directory:** `batch_jobs/`")

        n_jobs = len(bp.jobs)
        st.markdown(f"**Total jobs:** {n_jobs}")

        if st.button("🗑️ Clear Completed Jobs"):
            bp.clear_completed()
            st.success("Cleared completed, failed, and cancelled jobs.")
            st.rerun()


if __name__ == "__main__":
    main()
