"""
Alignment/Mapping Page

Align FASTQ reads to reference genome using bowtie2 or BWA.
FASTQ → BAM
"""

import streamlit as st
import pandas as pd
from pathlib import Path
import subprocess
import shutil
import tempfile
import os
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

# Import workflow manager for step recording
try:
    from frontend.components.workflow_manager import WorkflowManager

    HAS_WORKFLOW_MANAGER = True
except ImportError:
    HAS_WORKFLOW_MANAGER = False

st.set_page_config(page_title="Alignment - EpiNexus", page_icon="🧬", layout="wide")


def check_tool_installed(tool: str) -> bool:
    """Check if a command-line tool is installed."""
    return shutil.which(tool) is not None


def get_available_tools():
    """Get list of available alignment tools."""
    tools = {}
    tools["bowtie2"] = check_tool_installed("bowtie2")
    tools["bwa"] = check_tool_installed("bwa")
    tools["samtools"] = check_tool_installed("samtools")
    return tools


def run_alignment(
    fastq_r1: str,
    fastq_r2: str,
    output_bam: str,
    aligner: str,
    index_path: str,
    threads: int = 4,
    extra_params: str = "",
) -> dict:
    """
    Run alignment using specified aligner.

    Returns dict with status, stdout, stderr, and output file.
    """
    result = {"success": False, "stdout": "", "stderr": "", "output_bam": output_bam}

    try:
        if aligner == "bowtie2":
            # Bowtie2 alignment
            if fastq_r2:  # Paired-end
                cmd = f"bowtie2 -x {index_path} -1 {fastq_r1} -2 {fastq_r2} -p {threads} {extra_params}"
            else:  # Single-end
                cmd = f"bowtie2 -x {index_path} -U {fastq_r1} -p {threads} {extra_params}"

            # Pipe to samtools for BAM conversion and sorting
            full_cmd = f"{cmd} | samtools view -bS - | samtools sort -@ {threads} -o {output_bam}"

        elif aligner == "bwa":
            # BWA-MEM alignment
            if fastq_r2:  # Paired-end
                cmd = f"bwa mem -t {threads} {extra_params} {index_path} {fastq_r1} {fastq_r2}"
            else:  # Single-end
                cmd = f"bwa mem -t {threads} {extra_params} {index_path} {fastq_r1}"

            full_cmd = f"{cmd} | samtools view -bS - | samtools sort -@ {threads} -o {output_bam}"

        # Run the command
        process = subprocess.run(
            full_cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=7200,  # 2 hour timeout
        )

        result["stdout"] = process.stdout
        result["stderr"] = process.stderr
        result["success"] = process.returncode == 0

        # Index the BAM file
        if result["success"]:
            subprocess.run(f"samtools index {output_bam}", shell=True)

    except subprocess.TimeoutExpired:
        result["stderr"] = "Alignment timed out after 2 hours"
    except Exception as e:
        result["stderr"] = str(e)

    return result


def get_alignment_stats(bam_file: str) -> dict:
    """Get alignment statistics from BAM file using samtools."""
    stats = {}

    try:
        # Get flagstat
        result = subprocess.run(f"samtools flagstat {bam_file}", shell=True, capture_output=True, text=True)

        if result.returncode == 0:
            lines = result.stdout.strip().split("\n")
            for line in lines:
                if "total" in line and "QC-passed" in line:
                    stats["total_reads"] = int(line.split()[0])
                elif "mapped" in line and "mapped (" in line:
                    parts = line.split()
                    stats["mapped_reads"] = int(parts[0])
                    # Extract percentage
                    pct = line.split("(")[1].split("%")[0]
                    stats["mapping_rate"] = float(pct)
                elif "properly paired" in line:
                    stats["properly_paired"] = int(line.split()[0])
                elif "duplicates" in line:
                    stats["duplicates"] = int(line.split()[0])

    except Exception as e:
        stats["error"] = str(e)

    return stats


def main():
    st.title("🧬 Read Alignment")
    st.markdown("""
    Align sequencing reads (FASTQ) to a reference genome to produce BAM files.
    Supports **Bowtie2** (recommended for ChIP-seq/CUT&Tag) and **BWA-MEM** (for ATAC-seq).
    """)

    # Check available tools
    tools = get_available_tools()

    # Tool status
    st.markdown("### 🔧 Tool Status")
    col1, col2, col3 = st.columns(3)

    with col1:
        if tools["bowtie2"]:
            st.success("✅ Bowtie2 installed")
        else:
            st.error("❌ Bowtie2 not found")

    with col2:
        if tools["bwa"]:
            st.success("✅ BWA installed")
        else:
            st.error("❌ BWA not found")

    with col3:
        if tools["samtools"]:
            st.success("✅ Samtools installed")
        else:
            st.error("❌ Samtools not found")

    if not tools["samtools"]:
        st.error("⚠️ Samtools is required for BAM processing. Please install it first.")
        st.code("conda install -c bioconda samtools")
        return

    if not (tools["bowtie2"] or tools["bwa"]):
        st.error("⚠️ No aligner available. Please install Bowtie2 or BWA.")
        st.code("conda install -c bioconda bowtie2 bwa")
        return

    st.markdown("---")

    # Main interface
    tab1, tab2, tab3 = st.tabs(["📤 New Alignment", "📁 Batch Processing", "📊 View Results"])

    with tab1:
        render_single_alignment(tools)

    with tab2:
        render_batch_alignment(tools)

    with tab3:
        render_alignment_results()


def render_single_alignment(tools):
    """Single sample alignment interface."""
    st.header("Single Sample Alignment")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("📥 Input Files")

        # Input method
        input_method = st.radio("Input method", ["Upload FASTQ", "Specify file path"], horizontal=True)

        if input_method == "Upload FASTQ":
            fastq_r1 = st.file_uploader("Read 1 (R1) FASTQ", type=["fastq", "fq", "gz"], key="fastq_r1")

            is_paired = st.checkbox("Paired-end sequencing", value=True)

            if is_paired:
                fastq_r2 = st.file_uploader("Read 2 (R2) FASTQ", type=["fastq", "fq", "gz"], key="fastq_r2")
            else:
                fastq_r2 = None

        else:
            fastq_r1_path = st.text_input("Path to Read 1 FASTQ", placeholder="/path/to/sample_R1.fastq.gz")

            is_paired = st.checkbox("Paired-end sequencing", value=True)

            if is_paired:
                fastq_r2_path = st.text_input("Path to Read 2 FASTQ", placeholder="/path/to/sample_R2.fastq.gz")
            else:
                fastq_r2_path = None

        # Sample name
        sample_name = st.text_input("Sample name", value="sample_01")

    with col2:
        st.subheader("⚙️ Alignment Settings")

        # Aligner selection
        available_aligners = []
        if tools["bowtie2"]:
            available_aligners.append("bowtie2")
        if tools["bwa"]:
            available_aligners.append("bwa")

        aligner = st.selectbox(
            "Aligner", available_aligners, help="Bowtie2 recommended for ChIP-seq/CUT&Tag, BWA for ATAC-seq"
        )

        # Reference genome
        genome = st.selectbox("Reference genome", ["hg38", "hg19", "mm10", "mm39", "dm6"], key="align_genome")

        # Index path
        default_index = (
            f"~/genomes/{genome}/bowtie2_index/{genome}"
            if aligner == "bowtie2"
            else f"~/genomes/{genome}/bwa_index/{genome}"
        )
        index_path = st.text_input(
            "Index path", value=default_index, help="Path to aligner index (without file extension)"
        )

        # Advanced options
        with st.expander("Advanced Options"):
            threads = st.slider("CPU threads", 1, 16, 4)

            if aligner == "bowtie2":
                preset = st.selectbox("Preset", ["--sensitive", "--very-sensitive", "--fast", "--very-fast"], index=1)
                extra_params = st.text_input("Extra parameters", value=preset)
            else:
                extra_params = st.text_input("Extra parameters", value="-M")

        # Output directory
        output_dir = st.text_input("Output directory", value="./aligned", help="Directory to save BAM files")

    # Run alignment button
    st.markdown("---")

    if st.button("🚀 Run Alignment", type="primary", use_container_width=True):
        # Validate inputs
        if input_method == "Upload FASTQ":
            if fastq_r1 is None:
                st.error("Please upload Read 1 FASTQ file")
                return
            if is_paired and fastq_r2 is None:
                st.error("Please upload Read 2 FASTQ file")
                return
        else:
            if not fastq_r1_path:
                st.error("Please specify path to Read 1 FASTQ")
                return
            if is_paired and not fastq_r2_path:
                st.error("Please specify path to Read 2 FASTQ")
                return

        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        output_bam = os.path.join(output_dir, f"{sample_name}.sorted.bam")

        # Run alignment
        with st.spinner(f"Aligning reads with {aligner}... This may take a while."):
            progress_bar = st.progress(0)
            status_text = st.empty()

            status_text.text("Preparing input files...")
            progress_bar.progress(10)

            # Handle uploaded files
            if input_method == "Upload FASTQ":
                with tempfile.TemporaryDirectory() as tmpdir:
                    # Save uploaded files
                    r1_path = os.path.join(tmpdir, "R1.fastq.gz" if fastq_r1.name.endswith(".gz") else "R1.fastq")
                    with open(r1_path, "wb") as f:
                        f.write(fastq_r1.getvalue())

                    r2_path = None
                    if is_paired and fastq_r2:
                        r2_path = os.path.join(tmpdir, "R2.fastq.gz" if fastq_r2.name.endswith(".gz") else "R2.fastq")
                        with open(r2_path, "wb") as f:
                            f.write(fastq_r2.getvalue())

                    status_text.text(f"Running {aligner} alignment...")
                    progress_bar.progress(30)

                    result = run_alignment(
                        r1_path, r2_path, output_bam, aligner, os.path.expanduser(index_path), threads, extra_params
                    )
            else:
                status_text.text(f"Running {aligner} alignment...")
                progress_bar.progress(30)

                result = run_alignment(
                    fastq_r1_path,
                    fastq_r2_path if is_paired else None,
                    output_bam,
                    aligner,
                    os.path.expanduser(index_path),
                    threads,
                    extra_params,
                )

            progress_bar.progress(80)

            if result["success"]:
                status_text.text("Getting alignment statistics...")
                stats = get_alignment_stats(output_bam)
                progress_bar.progress(100)

                st.success(f"✅ Alignment complete! Output: `{output_bam}`")

                # Show stats
                if stats and "total_reads" in stats:
                    col1, col2, col3, col4 = st.columns(4)
                    col1.metric("Total Reads", f"{stats.get('total_reads', 0):,}")
                    col2.metric("Mapped Reads", f"{stats.get('mapped_reads', 0):,}")
                    col3.metric("Mapping Rate", f"{stats.get('mapping_rate', 0):.1f}%")
                    col4.metric("Properly Paired", f"{stats.get('properly_paired', 0):,}")

                # Save to session state
                if "aligned_bams" not in st.session_state:
                    st.session_state.aligned_bams = []

                st.session_state.aligned_bams.append(
                    {"sample": sample_name, "bam": output_bam, "stats": stats, "aligner": aligner, "genome": genome}
                )

                # Record workflow step
                if HAS_WORKFLOW_MANAGER:
                    WorkflowManager.record_step(
                        step_type="alignment",
                        parameters={
                            "aligner": aligner,
                            "genome": genome,
                            "threads": threads,
                            "paired_end": is_paired,
                        },
                        inputs={
                            "fastq_r1": fastq_r1.name if input_method == "Upload FASTQ" else fastq_r1_path,
                            "fastq_r2": (fastq_r2.name if fastq_r2 else None)
                            if input_method == "Upload FASTQ"
                            else fastq_r2_path,
                        },
                        outputs={"bam": output_bam},
                        tool_versions={aligner.lower(): WorkflowManager.get_tool_version(aligner)},
                        output_metadata={
                            "total_reads": stats.get("total_reads", 0),
                            "mapped_reads": stats.get("mapped_reads", 0),
                            "mapping_rate": stats.get("mapping_rate", 0),
                        },
                        command_line=result.get("command", ""),
                    )

            else:
                st.error("❌ Alignment failed!")
                st.code(result["stderr"])


def render_batch_alignment(tools):
    """Batch alignment interface."""
    st.header("Batch Processing")

    st.markdown("""
    Process multiple samples at once. Upload a sample sheet or specify a directory containing FASTQ files.
    """)

    # Sample sheet format
    st.subheader("📋 Sample Sheet")

    sample_sheet_method = st.radio(
        "Sample sheet input", ["Upload CSV", "Auto-detect from directory", "Manual entry"], horizontal=True
    )

    if sample_sheet_method == "Upload CSV":
        st.markdown("""
        Upload a CSV with columns: `sample_name`, `fastq_r1`, `fastq_r2` (optional for SE)
        """)

        sample_sheet_file = st.file_uploader("Upload sample sheet", type=["csv", "tsv"])

        if sample_sheet_file:
            sample_sheet = pd.read_csv(sample_sheet_file)
            st.dataframe(sample_sheet, use_container_width=True)

    elif sample_sheet_method == "Auto-detect from directory":
        fastq_dir = st.text_input("FASTQ directory", placeholder="/path/to/fastq_files/")

        if fastq_dir and os.path.exists(fastq_dir):
            # Find FASTQ files
            fastq_files = list(Path(fastq_dir).glob("*.fastq*")) + list(Path(fastq_dir).glob("*.fq*"))

            if fastq_files:
                st.success(f"Found {len(fastq_files)} FASTQ files")

                # Try to pair them
                samples = {}
                for f in fastq_files:
                    name = f.stem.replace(".fastq", "").replace(".fq", "").replace(".gz", "")
                    # Try to extract sample name
                    if "_R1" in name or "_1" in name:
                        sample = name.replace("_R1", "").replace("_1", "")
                        if sample not in samples:
                            samples[sample] = {}
                        samples[sample]["r1"] = str(f)
                    elif "_R2" in name or "_2" in name:
                        sample = name.replace("_R2", "").replace("_2", "")
                        if sample not in samples:
                            samples[sample] = {}
                        samples[sample]["r2"] = str(f)
                    else:
                        samples[name] = {"r1": str(f)}

                sample_sheet = pd.DataFrame(
                    [
                        {"sample_name": k, "fastq_r1": v.get("r1", ""), "fastq_r2": v.get("r2", "")}
                        for k, v in samples.items()
                    ]
                )

                st.dataframe(sample_sheet, use_container_width=True)
            else:
                st.warning("No FASTQ files found in directory")

    else:
        st.markdown("Enter samples manually:")

        num_samples = st.number_input("Number of samples", 1, 100, 3)

        samples_data = []
        for i in range(num_samples):
            col1, col2, col3 = st.columns(3)
            with col1:
                name = st.text_input(f"Sample {i + 1} name", key=f"batch_name_{i}")
            with col2:
                r1 = st.text_input("R1 path", key=f"batch_r1_{i}")
            with col3:
                r2 = st.text_input("R2 path (optional)", key=f"batch_r2_{i}")
            samples_data.append({"sample_name": name, "fastq_r1": r1, "fastq_r2": r2})

        sample_sheet = pd.DataFrame(samples_data)

    # Alignment settings
    st.subheader("⚙️ Settings")

    col1, col2 = st.columns(2)

    with col1:
        available_aligners = []
        if tools["bowtie2"]:
            available_aligners.append("bowtie2")
        if tools["bwa"]:
            available_aligners.append("bwa")

        batch_aligner = st.selectbox("Aligner", available_aligners, key="batch_aligner")
        batch_genome = st.selectbox("Genome", ["hg38", "hg19", "mm10", "mm39"], key="batch_genome")

    with col2:
        batch_index = st.text_input("Index path", key="batch_index")
        batch_threads = st.slider("Threads per sample", 1, 8, 4, key="batch_threads")
        batch_output = st.text_input("Output directory", value="./aligned", key="batch_output")

    # Run batch
    if st.button("🚀 Run Batch Alignment", type="primary"):
        st.info("Batch alignment would start here...")
        # TODO: Implement batch processing with progress tracking


def render_alignment_results():
    """View alignment results."""
    st.header("Alignment Results")

    if "aligned_bams" not in st.session_state or not st.session_state.aligned_bams:
        st.info("No alignments completed yet. Run an alignment first!")
        return

    # Summary table
    results_df = pd.DataFrame(st.session_state.aligned_bams)

    # Extract stats into columns
    if "stats" in results_df.columns:
        for stat in ["total_reads", "mapped_reads", "mapping_rate"]:
            results_df[stat] = results_df["stats"].apply(lambda x: x.get(stat, 0) if isinstance(x, dict) else 0)

    st.dataframe(
        results_df[["sample", "bam", "aligner", "genome", "mapping_rate"]], use_container_width=True, hide_index=True
    )

    # Visualization
    if len(results_df) > 1:
        st.subheader("Mapping Rates")

        import plotly.express as px

        fig = px.bar(
            results_df,
            x="sample",
            y="mapping_rate",
            color="mapping_rate",
            color_continuous_scale="RdYlGn",
            title="Mapping Rate by Sample",
        )
        fig.add_hline(y=80, line_dash="dash", line_color="red", annotation_text="80% threshold")
        fig.update_layout(height=400)
        st.plotly_chart(fig, use_container_width=True)

    # Pass to peak calling
    st.markdown("---")
    st.markdown("### Next Step")

    if st.button("📊 Continue to Peak Calling →", use_container_width=True):
        st.switch_page("pages/23_peak_calling.py")


if __name__ == "__main__":
    main()
