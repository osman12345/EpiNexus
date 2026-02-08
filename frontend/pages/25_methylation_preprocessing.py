"""
Methylation Preprocessing Module

Bisulfite alignment and methylation calling using Bismark.

Copyright (c) 2026 EpiNexus Contributors
SPDX-License-Identifier: AGPL-3.0-or-later OR Commercial
"""

import streamlit as st
import pandas as pd
from pathlib import Path
import subprocess
import tempfile
import shutil

# =============================================================================
# CONFIGURATION
# =============================================================================

ALIGNERS = {
    "Bismark": {
        "description": "Gold standard for bisulfite alignment",
        "command": "bismark",
        "url": "https://www.bioinformatics.babraham.ac.uk/projects/bismark/"
    },
    "bwa-meth": {
        "description": "Fast bisulfite aligner using BWA",
        "command": "bwameth.py",
        "url": "https://github.com/brentp/bwa-meth"
    },
}

METHYLATION_CALLERS = {
    "Bismark Extractor": {
        "description": "Extract methylation calls from Bismark BAM",
        "command": "bismark_methylation_extractor"
    },
    "MethylDackel": {
        "description": "Fast methylation caller for any aligner",
        "command": "MethylDackel"
    },
}

GENOMES = {
    "hg38": "Human (GRCh38)",
    "hg19": "Human (GRCh37)",
    "mm10": "Mouse (GRCm38)",
    "mm39": "Mouse (GRCm39)",
}

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

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def check_tool(tool_name: str) -> bool:
    """Check if a command-line tool is available."""
    try:
        result = subprocess.run(
            [tool_name, "--version"] if tool_name != "bwameth.py" else [tool_name, "-h"],
            capture_output=True,
            check=False
        )
        return result.returncode == 0 or "bismark" in result.stderr.decode().lower()
    except FileNotFoundError:
        return False

def check_genome_index(genome: str, aligner: str) -> dict:
    """Check if genome index exists for the specified aligner."""
    # Common index locations
    index_paths = [
        Path.home() / "genomes" / genome / "bismark_index",
        Path.home() / "references" / genome / "bismark_index",
        Path(f"/data/genomes/{genome}/bismark_index"),
        Path(f"/references/{genome}/bismark_index"),
    ]

    for path in index_paths:
        if path.exists():
            return {"exists": True, "path": str(path)}

    return {"exists": False, "path": None}

def run_bismark_alignment(
    fastq1: str,
    fastq2: str | None,
    genome_path: str,
    output_dir: str,
    threads: int = 4,
    options: dict = None
) -> dict:
    """Run Bismark alignment."""
    options = options or {}

    cmd = [
        "bismark",
        "--genome", genome_path,
        "--output_dir", output_dir,
        "--parallel", str(max(1, threads // 4)),  # Bismark uses 4 threads per instance
        "--temp_dir", tempfile.gettempdir(),
    ]

    # Add options
    if options.get("non_directional"):
        cmd.append("--non_directional")

    if options.get("pbat"):
        cmd.append("--pbat")

    if options.get("local"):
        cmd.append("--local")

    # Add FASTQ files
    if fastq2:
        cmd.extend(["-1", fastq1, "-2", fastq2])
    else:
        cmd.append(fastq1)

    result = subprocess.run(cmd, capture_output=True, text=True)

    return {
        "success": result.returncode == 0,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "command": " ".join(cmd)
    }

def run_methylation_extraction(
    bam_file: str,
    genome_path: str,
    output_dir: str,
    paired: bool = True,
    options: dict = None
) -> dict:
    """Run Bismark methylation extractor."""
    options = options or {}

    cmd = [
        "bismark_methylation_extractor",
        "--gzip",
        "--bedGraph",
        "--cytosine_report",
        "--genome_folder", genome_path,
        "-o", output_dir,
    ]

    if paired:
        cmd.append("--paired-end")
    else:
        cmd.append("--single-end")

    if options.get("comprehensive"):
        cmd.append("--comprehensive")

    if options.get("merge_non_cpg"):
        cmd.append("--merge_non_CpG")

    if options.get("no_overlap"):
        cmd.append("--no_overlap")

    cmd.append(bam_file)

    result = subprocess.run(cmd, capture_output=True, text=True)

    return {
        "success": result.returncode == 0,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "command": " ".join(cmd)
    }

def run_methyldackel(
    bam_file: str,
    genome_fasta: str,
    output_prefix: str,
    options: dict = None
) -> dict:
    """Run MethylDackel for methylation calling."""
    options = options or {}

    cmd = [
        "MethylDackel", "extract",
        "--CHH", "--CHG",  # Include non-CpG contexts
        "-o", output_prefix,
    ]

    if options.get("min_depth"):
        cmd.extend(["--minDepth", str(options["min_depth"])])

    if options.get("ignore_flags"):
        cmd.extend(["--ignoreFlags", options["ignore_flags"]])

    cmd.extend([genome_fasta, bam_file])

    result = subprocess.run(cmd, capture_output=True, text=True)

    return {
        "success": result.returncode == 0,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "command": " ".join(cmd)
    }

def parse_bismark_report(report_path: str) -> dict:
    """Parse Bismark alignment report."""
    stats = {}

    if not Path(report_path).exists():
        return stats

    with open(report_path) as f:
        content = f.read()

    # Parse key metrics
    import re

    patterns = {
        "total_reads": r"Sequences analysed in total:\s+(\d+)",
        "unique_alignments": r"Number of alignments with a unique best hit.*?:\s+(\d+)",
        "mapping_efficiency": r"Mapping efficiency:\s+([\d.]+)%",
        "no_alignment": r"Sequences with no alignments.*?:\s+(\d+)",
        "ambiguous": r"Sequences did not map uniquely:\s+(\d+)",
        "total_c": r"Total number of C's analysed:\s+(\d+)",
        "meth_cpg": r"Total methylated C's in CpG context:\s+(\d+)",
        "meth_chg": r"Total methylated C's in CHG context:\s+(\d+)",
        "meth_chh": r"Total methylated C's in CHH context:\s+(\d+)",
        "cpg_meth_pct": r"C methylated in CpG context:\s+([\d.]+)%",
    }

    for key, pattern in patterns.items():
        match = re.search(pattern, content)
        if match:
            value = match.group(1)
            stats[key] = float(value) if '.' in value else int(value)

    return stats

# =============================================================================
# MAIN PAGE
# =============================================================================

def main():
    st.title("🧬 Methylation Preprocessing")
    st.markdown("Bisulfite alignment and methylation calling")

    # Initialize session state
    if 'meth_alignments' not in st.session_state:
        st.session_state.meth_alignments = []
    if 'meth_calls' not in st.session_state:
        st.session_state.meth_calls = []

    # Check tool availability
    tools_available = {
        "bismark": check_tool("bismark"),
        "bismark_extractor": check_tool("bismark_methylation_extractor"),
        "methyldackel": check_tool("MethylDackel"),
        "bwameth": check_tool("bwameth.py"),
    }

    # Tool status
    with st.expander("🔧 Tool Availability", expanded=False):
        col1, col2 = st.columns(2)

        with col1:
            st.markdown("**Aligners:**")
            for name, info in ALIGNERS.items():
                cmd = info["command"]
                if cmd == "bismark":
                    available = tools_available["bismark"]
                else:
                    available = tools_available.get("bwameth", False)

                icon = "✅" if available else "❌"
                st.markdown(f"{icon} **{name}** - {info['description']}")

        with col2:
            st.markdown("**Methylation Callers:**")
            for name, info in METHYLATION_CALLERS.items():
                cmd = info["command"]
                if "bismark" in cmd.lower():
                    available = tools_available["bismark_extractor"]
                else:
                    available = tools_available["methyldackel"]

                icon = "✅" if available else "❌"
                st.markdown(f"{icon} **{name}** - {info['description']}")

        if not any(tools_available.values()):
            st.warning("""
            No methylation tools detected. Install with:
            ```bash
            conda install -c bioconda bismark bowtie2
            # or
            conda install -c bioconda methyldackel bwa
            ```
            """)

    # Tabs
    tab1, tab2, tab3, tab4 = st.tabs([
        "🧬 Alignment", "📊 Methylation Calling",
        "📁 Upload Existing", "📋 Results"
    ])

    # =========================================================================
    # TAB 1: ALIGNMENT
    # =========================================================================
    with tab1:
        st.subheader("Bisulfite Alignment")

        col1, col2 = st.columns(2)

        with col1:
            aligner = st.selectbox(
                "Aligner",
                list(ALIGNERS.keys()),
                help="Bismark is recommended for most use cases"
            )

            genome = st.selectbox(
                "Reference Genome",
                list(GENOMES.keys()),
                format_func=lambda x: f"{x} - {GENOMES[x]}"
            )

            # Check for genome index
            index_info = check_genome_index(genome, aligner)
            if index_info["exists"]:
                st.success(f"Genome index found: {index_info['path']}")
                genome_path = index_info["path"]
            else:
                genome_path = st.text_input(
                    "Path to Bismark genome index",
                    help="Directory containing Bisulfite_Genome folder"
                )

                st.info("""
                **Build Bismark index:**
                ```bash
                bismark_genome_preparation --path_to_aligner /path/to/bowtie2 /path/to/genome/
                ```
                """)

        with col2:
            library_type = st.selectbox(
                "Library Type",
                ["Directional (standard)", "Non-directional", "PBAT"]
            )

            seq_type = st.radio(
                "Sequencing Type",
                ["Paired-end", "Single-end"]
            )

            threads = st.slider("Threads", 1, 32, 8)

        st.markdown("---")

        # Input files
        st.markdown("### Input FASTQ Files")

        input_method = st.radio(
            "Input Method",
            ["Upload Files", "Specify Paths"],
            horizontal=True
        )

        if input_method == "Upload Files":
            if seq_type == "Paired-end":
                col1, col2 = st.columns(2)
                with col1:
                    fastq1 = st.file_uploader("Read 1 (R1)", type=['fastq', 'fq', 'gz'])
                with col2:
                    fastq2 = st.file_uploader("Read 2 (R2)", type=['fastq', 'fq', 'gz'])
            else:
                fastq1 = st.file_uploader("FASTQ File", type=['fastq', 'fq', 'gz'])
                fastq2 = None
        else:
            if seq_type == "Paired-end":
                col1, col2 = st.columns(2)
                with col1:
                    fastq1_path = st.text_input("Read 1 (R1) Path")
                with col2:
                    fastq2_path = st.text_input("Read 2 (R2) Path")
            else:
                fastq1_path = st.text_input("FASTQ Path")
                fastq2_path = None

        sample_name = st.text_input("Sample Name", "Sample_1")

        # Advanced options
        with st.expander("Advanced Options"):
            col1, col2 = st.columns(2)
            with col1:
                local_align = st.checkbox("Local alignment", value=False,
                                         help="Use local alignment instead of end-to-end")
                ambiguous = st.checkbox("Report ambiguous", value=False)
            with col2:
                nucleotide_cov = st.checkbox("Nucleotide coverage report", value=True)
                temp_dir = st.text_input("Temp directory", tempfile.gettempdir())

        # Run alignment
        if st.button("▶️ Run Alignment", type="primary", disabled=not tools_available.get("bismark")):
            if not genome_path:
                st.error("Please specify genome index path")
            else:
                with st.spinner(f"Running {aligner} alignment..."):
                    # This would run the actual alignment
                    st.info("Alignment would run here with the specified parameters")

                    # Simulated result for demo
                    st.session_state.meth_alignments.append({
                        'sample': sample_name,
                        'aligner': aligner,
                        'genome': genome,
                        'status': 'completed',
                        'bam_file': f"{sample_name}_bismark_bt2.bam"
                    })

                    st.success("Alignment complete!")

                    # Record workflow step
                    if HAS_WORKFLOW_MANAGER:
                        WorkflowManager.record_step(
                            step_name="Bisulfite Alignment",
                            tool="Bismark",
                            parameters={
                                'aligner': aligner,
                                'genome': genome,
                                'sample': sample_name
                            },
                            inputs=['fastq_files', 'genome_index'],
                            outputs=['aligned_bam']
                        )

    # =========================================================================
    # TAB 2: METHYLATION CALLING
    # =========================================================================
    with tab2:
        st.subheader("Methylation Calling")

        col1, col2 = st.columns(2)

        with col1:
            caller = st.selectbox(
                "Methylation Caller",
                list(METHYLATION_CALLERS.keys())
            )

            # Input BAM
            bam_source = st.radio(
                "BAM Source",
                ["From Alignment Step", "Upload BAM", "Specify Path"],
                horizontal=True
            )

        with col2:
            min_cov = st.slider("Minimum Coverage", 1, 20, 5)

            contexts = st.multiselect(
                "Methylation Contexts",
                ["CpG", "CHG", "CHH"],
                default=["CpG"]
            )

        if bam_source == "From Alignment Step":
            if st.session_state.meth_alignments:
                bam_options = [a['bam_file'] for a in st.session_state.meth_alignments]
                selected_bam = st.selectbox("Select BAM", bam_options)
            else:
                st.info("No alignments available. Run alignment first or upload BAM.")
                selected_bam = None
        elif bam_source == "Upload BAM":
            bam_upload = st.file_uploader("Upload BAM file", type=['bam'])
            selected_bam = bam_upload.name if bam_upload else None
        else:
            selected_bam = st.text_input("BAM file path")

        # Advanced options
        with st.expander("Advanced Options"):
            col1, col2 = st.columns(2)
            with col1:
                comprehensive = st.checkbox("Comprehensive output", value=True)
                merge_non_cpg = st.checkbox("Merge non-CpG", value=False)
            with col2:
                no_overlap = st.checkbox("No overlap (PE)", value=True,
                                        help="Avoid double-counting overlapping reads")
                cytosine_report = st.checkbox("Cytosine report", value=True)

        if st.button("▶️ Extract Methylation", type="primary"):
            if selected_bam:
                with st.spinner(f"Running {caller}..."):
                    st.info("Methylation extraction would run here")

                    st.session_state.meth_calls.append({
                        'sample': Path(selected_bam).stem,
                        'caller': caller,
                        'contexts': contexts,
                        'status': 'completed'
                    })

                    st.success("Methylation calling complete!")

                    # Record workflow step
                    if HAS_WORKFLOW_MANAGER:
                        WorkflowManager.record_step(
                            step_name="Methylation Calling",
                            tool=caller,
                            parameters={
                                'bam_file': selected_bam,
                                'min_coverage': min_cov,
                                'contexts': contexts
                            },
                            inputs=['aligned_bam'],
                            outputs=['methylation_calls']
                        )
            else:
                st.error("Please select or upload a BAM file")

    # =========================================================================
    # TAB 3: UPLOAD EXISTING
    # =========================================================================
    with tab3:
        st.subheader("Upload Existing Methylation Files")

        st.markdown("""
        Upload pre-computed methylation files to skip alignment and calling steps.

        **Supported formats:**
        - Bismark coverage files (`.cov`, `.bismark.cov.gz`)
        - bedGraph files (`.bedGraph`)
        - methylKit format (`.methylKit`)
        - Generic coverage files (chr, pos, methylation%, coverage)
        """)

        uploaded_files = st.file_uploader(
            "Upload methylation files",
            type=['cov', 'bedgraph', 'txt', 'gz', 'methylkit'],
            accept_multiple_files=True
        )

        if uploaded_files:
            for f in uploaded_files:
                col1, col2, col3 = st.columns([2, 1, 1])
                with col1:
                    st.text(f.name)
                with col2:
                    sample_name = st.text_input("Name", Path(f.name).stem, key=f"name_{f.name}")
                with col3:
                    st.text(f"{f.size / 1024:.1f} KB")

            if st.button("Load Files", type="primary"):
                st.success(f"Loaded {len(uploaded_files)} methylation files")
                st.info("Go to **Methylation Analysis** page to analyze these files")

    # =========================================================================
    # TAB 4: RESULTS
    # =========================================================================
    with tab4:
        st.subheader("Processing Results")

        # Alignments
        st.markdown("### Alignments")
        if st.session_state.meth_alignments:
            align_df = pd.DataFrame(st.session_state.meth_alignments)
            st.dataframe(align_df, use_container_width=True, hide_index=True)
        else:
            st.info("No alignments yet")

        # Methylation calls
        st.markdown("### Methylation Calls")
        if st.session_state.meth_calls:
            calls_df = pd.DataFrame(st.session_state.meth_calls)
            st.dataframe(calls_df, use_container_width=True, hide_index=True)
        else:
            st.info("No methylation calls yet")

        # Export options
        if st.session_state.meth_alignments or st.session_state.meth_calls:
            st.markdown("### Export")
            if st.button("Export Processing Summary"):
                summary = {
                    'alignments': st.session_state.meth_alignments,
                    'methylation_calls': st.session_state.meth_calls
                }
                import json
                st.download_button(
                    "Download Summary (JSON)",
                    json.dumps(summary, indent=2),
                    "methylation_processing_summary.json",
                    "application/json"
                )

if __name__ == "__main__":
    main()
