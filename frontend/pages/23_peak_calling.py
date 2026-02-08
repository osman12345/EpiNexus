"""
Peak Calling Page

Call peaks from aligned BAM files using MACS2 or SEACR.
BAM → Peaks (BED/narrowPeak)
"""

import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import subprocess
import shutil
import os
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

st.set_page_config(
    page_title="Peak Calling - EpiNexus",
    page_icon="📍",
    layout="wide"
)


def check_tool_installed(tool: str) -> bool:
    """Check if a command-line tool is installed."""
    return shutil.which(tool) is not None


def get_available_tools():
    """Get list of available peak calling tools."""
    tools = {}
    tools['macs2'] = check_tool_installed('macs2') or check_tool_installed('macs3')
    tools['seacr'] = check_tool_installed('SEACR_1.3.sh') or os.path.exists(os.path.expanduser('~/SEACR/SEACR_1.3.sh'))
    tools['samtools'] = check_tool_installed('samtools')
    tools['bedtools'] = check_tool_installed('bedtools')
    return tools


def run_macs2(bam_file: str, output_prefix: str, genome: str,
              control_bam: str = None, call_summits: bool = True,
              broad: bool = False, qvalue: float = 0.05,
              extra_params: str = "") -> dict:
    """
    Run MACS2 peak calling.

    Returns dict with status and output files.
    """
    result = {
        'success': False,
        'stdout': '',
        'stderr': '',
        'peaks_file': ''
    }

    # Genome size mapping
    genome_sizes = {
        'hg38': 'hs', 'hg19': 'hs',
        'mm10': 'mm', 'mm39': 'mm',
        'dm6': 'dm'
    }
    gsize = genome_sizes.get(genome, 'hs')

    try:
        # Build MACS2 command
        cmd = f"macs2 callpeak -t {bam_file} -g {gsize} -n {output_prefix} -q {qvalue}"

        if control_bam:
            cmd += f" -c {control_bam}"

        if call_summits and not broad:
            cmd += " --call-summits"

        if broad:
            cmd += " --broad"

        cmd += f" {extra_params}"

        # Run MACS2
        process = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=3600  # 1 hour timeout
        )

        result['stdout'] = process.stdout
        result['stderr'] = process.stderr
        result['success'] = process.returncode == 0

        # Find output peaks file
        if broad:
            result['peaks_file'] = f"{output_prefix}_peaks.broadPeak"
        else:
            result['peaks_file'] = f"{output_prefix}_peaks.narrowPeak"

    except subprocess.TimeoutExpired:
        result['stderr'] = "Peak calling timed out after 1 hour"
    except Exception as e:
        result['stderr'] = str(e)

    return result


def run_seacr(bedgraph_file: str, output_prefix: str,
              control_bedgraph: str = None, threshold: float = 0.01,
              norm: str = "norm", mode: str = "relaxed") -> dict:
    """
    Run SEACR peak calling for CUT&Tag/CUT&RUN data.

    Returns dict with status and output files.
    """
    result = {
        'success': False,
        'stdout': '',
        'stderr': '',
        'peaks_file': ''
    }

    try:
        # Find SEACR script
        seacr_path = shutil.which('SEACR_1.3.sh')
        if not seacr_path:
            seacr_path = os.path.expanduser('~/SEACR/SEACR_1.3.sh')

        if not os.path.exists(seacr_path):
            result['stderr'] = "SEACR not found"
            return result

        # Build SEACR command
        if control_bedgraph:
            cmd = f"bash {seacr_path} {bedgraph_file} {control_bedgraph} {norm} {mode} {output_prefix}"
        else:
            cmd = f"bash {seacr_path} {bedgraph_file} {threshold} {norm} {mode} {output_prefix}"

        # Run SEACR
        process = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=1800
        )

        result['stdout'] = process.stdout
        result['stderr'] = process.stderr
        result['success'] = process.returncode == 0

        result['peaks_file'] = f"{output_prefix}.{mode}.bed"

    except subprocess.TimeoutExpired:
        result['stderr'] = "Peak calling timed out"
    except Exception as e:
        result['stderr'] = str(e)

    return result


def bam_to_bedgraph(bam_file: str, output_bg: str, genome: str) -> bool:
    """Convert BAM to bedGraph for SEACR."""
    try:
        # Get chromosome sizes
        chrom_sizes = os.path.expanduser(f"~/genomes/{genome}/{genome}.chrom.sizes")

        # BAM to bedGraph
        cmd = f"bedtools genomecov -bg -ibam {bam_file} > {output_bg}"

        process = subprocess.run(cmd, shell=True, capture_output=True)
        return process.returncode == 0

    except Exception:
        return False


def get_peak_stats(peaks_file: str) -> dict:
    """Get statistics from peak file."""
    stats = {}

    try:
        if not os.path.exists(peaks_file):
            return stats

        # Read peaks
        peaks = pd.read_csv(peaks_file, sep='\t', header=None)

        stats['n_peaks'] = len(peaks)

        # For narrowPeak/broadPeak format
        if len(peaks.columns) >= 3:
            peaks.columns = ['chr', 'start', 'end'] + [f'col{i}' for i in range(3, len(peaks.columns))]
            widths = peaks['end'] - peaks['start']
            stats['mean_width'] = int(widths.mean())
            stats['median_width'] = int(widths.median())
            stats['total_bp'] = int(widths.sum())

            # Signal if available (column 7 for narrowPeak)
            if len(peaks.columns) >= 7:
                stats['mean_signal'] = float(peaks.iloc[:, 6].mean())

    except Exception as e:
        stats['error'] = str(e)

    return stats


def main():
    st.title("📍 Peak Calling")
    st.markdown("""
    Call peaks from aligned BAM files to identify regions of enrichment.
    Supports **MACS2** (ChIP-seq, ATAC-seq) and **SEACR** (CUT&Tag, CUT&RUN).
    """)

    # Check available tools
    tools = get_available_tools()

    # Tool status
    st.markdown("### 🔧 Tool Status")
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        if tools['macs2']:
            st.success("✅ MACS2/3 installed")
        else:
            st.warning("⚠️ MACS2 not found")

    with col2:
        if tools['seacr']:
            st.success("✅ SEACR installed")
        else:
            st.warning("⚠️ SEACR not found")

    with col3:
        if tools['samtools']:
            st.success("✅ Samtools installed")
        else:
            st.error("❌ Samtools required")

    with col4:
        if tools['bedtools']:
            st.success("✅ Bedtools installed")
        else:
            st.warning("⚠️ Bedtools not found")

    if not (tools['macs2'] or tools['seacr']):
        st.error("⚠️ No peak caller available. Please install MACS2 or SEACR.")
        st.code("conda install -c bioconda macs2")
        return

    st.markdown("---")

    # Main tabs
    tab1, tab2, tab3, tab4 = st.tabs([
        "📤 Upload BAM",
        "🏔️ Call Peaks",
        "📁 Upload Existing Peaks",
        "📊 View Results"
    ])

    with tab1:
        render_bam_upload()

    with tab2:
        render_peak_calling(tools)

    with tab3:
        render_peak_upload()

    with tab4:
        render_peak_results()


def render_bam_upload():
    """Upload or specify BAM files."""
    st.header("Input BAM Files")

    st.markdown("""
    Upload aligned BAM files or specify paths to existing files.
    These will be used for peak calling.
    """)

    # Check for BAMs from alignment step
    if 'aligned_bams' in st.session_state and st.session_state.aligned_bams:
        st.success(f"Found {len(st.session_state.aligned_bams)} BAM files from alignment step!")

        bams_df = pd.DataFrame(st.session_state.aligned_bams)
        st.dataframe(bams_df[['sample', 'bam', 'genome']], use_container_width=True, hide_index=True)

        if st.button("Use these BAM files", type="primary"):
            st.session_state.peak_calling_bams = st.session_state.aligned_bams.copy()
            st.success("BAM files loaded for peak calling!")
            st.rerun()

    st.markdown("---")

    # Manual input
    st.subheader("Add BAM Files Manually")

    input_method = st.radio(
        "Input method",
        ["Specify file path", "Upload BAM"],
        horizontal=True
    )

    if input_method == "Specify file path":
        col1, col2, col3 = st.columns(3)

        with col1:
            sample_name = st.text_input("Sample name", key="bam_sample")
        with col2:
            bam_path = st.text_input("BAM file path", key="bam_path")
        with col3:
            genome = st.selectbox("Genome", ["hg38", "hg19", "mm10", "mm39"], key="bam_genome")

        control_bam = st.text_input(
            "Control/Input BAM (optional)",
            help="IgG or Input control for background subtraction"
        )

        if st.button("Add BAM"):
            if sample_name and bam_path:
                if 'peak_calling_bams' not in st.session_state:
                    st.session_state.peak_calling_bams = []

                st.session_state.peak_calling_bams.append({
                    'sample': sample_name,
                    'bam': bam_path,
                    'control': control_bam if control_bam else None,
                    'genome': genome
                })
                st.success(f"Added {sample_name}")
                st.rerun()

    else:
        uploaded_bam = st.file_uploader("Upload BAM file", type=['bam'])
        if uploaded_bam:
            st.info("BAM file uploaded. For large files, consider specifying the path instead.")

    # Show current BAM list
    if 'peak_calling_bams' in st.session_state and st.session_state.peak_calling_bams:
        st.markdown("---")
        st.subheader("Current BAM Files")

        bams_df = pd.DataFrame(st.session_state.peak_calling_bams)
        st.dataframe(bams_df, use_container_width=True, hide_index=True)

        if st.button("Clear all BAMs"):
            st.session_state.peak_calling_bams = []
            st.rerun()


def render_peak_calling(tools):
    """Peak calling interface."""
    st.header("Call Peaks")

    # Check for BAMs
    if 'peak_calling_bams' not in st.session_state or not st.session_state.peak_calling_bams:
        st.warning("No BAM files loaded. Go to 'Upload BAM' tab first.")
        return

    bams = st.session_state.peak_calling_bams

    st.success(f"{len(bams)} BAM file(s) ready for peak calling")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("🎯 Peak Caller Settings")

        # Select peak caller
        available_callers = []
        if tools['macs2']:
            available_callers.append('MACS2')
        if tools['seacr']:
            available_callers.append('SEACR')

        peak_caller = st.selectbox(
            "Peak caller",
            available_callers,
            help="MACS2 for ChIP-seq/ATAC-seq, SEACR for CUT&Tag/CUT&RUN"
        )

        # Assay type
        assay_type = st.selectbox(
            "Assay type",
            ["ChIP-seq (Histone)", "ChIP-seq (TF)", "CUT&Tag", "CUT&RUN", "ATAC-seq"]
        )

        # Output directory
        output_dir = st.text_input("Output directory", value="./peaks")

    with col2:
        st.subheader("⚙️ Parameters")

        if peak_caller == "MACS2":
            # MACS2 parameters
            if assay_type in ["ChIP-seq (Histone)", "CUT&Tag", "CUT&RUN"]:
                broad_peaks = st.checkbox("Call broad peaks", value=True,
                                         help="For histone marks like H3K27me3, H3K36me3")
            else:
                broad_peaks = st.checkbox("Call broad peaks", value=False)

            call_summits = st.checkbox("Call summits", value=not broad_peaks,
                                      help="Identify peak summits (narrow peaks only)")

            qvalue = st.number_input("Q-value threshold", 0.001, 0.1, 0.05, 0.01)

            with st.expander("Advanced MACS2 options"):
                extra_params = st.text_input("Extra parameters", value="")
                st.markdown("""
                Common options:
                - `--nomodel`: Don't build shifting model
                - `--shift -100 --extsize 200`: For ATAC-seq
                - `--keep-dup all`: Keep duplicate reads
                """)

        else:  # SEACR
            # SEACR parameters
            threshold = st.number_input(
                "FDR threshold (if no control)",
                0.001, 0.1, 0.01, 0.001,
                help="Used when no control BAM is provided"
            )

            seacr_mode = st.selectbox(
                "Peak mode",
                ["relaxed", "stringent"],
                help="Stringent = fewer, higher confidence peaks"
            )

            seacr_norm = st.selectbox(
                "Normalization",
                ["norm", "non"],
                help="norm = normalize signal to control"
            )

    # Sample selection
    st.markdown("---")
    st.subheader("📋 Samples to Process")

    sample_selection = st.multiselect(
        "Select samples",
        [b['sample'] for b in bams],
        default=[b['sample'] for b in bams]
    )

    # Run peak calling
    st.markdown("---")

    if st.button("🚀 Run Peak Calling", type="primary", use_container_width=True):
        os.makedirs(output_dir, exist_ok=True)

        progress = st.progress(0)
        status = st.empty()

        results = []

        for i, sample_name in enumerate(sample_selection):
            # Find BAM info
            bam_info = next((b for b in bams if b['sample'] == sample_name), None)
            if not bam_info:
                continue

            status.text(f"Processing {sample_name} ({i+1}/{len(sample_selection)})...")
            progress.progress((i + 0.5) / len(sample_selection))

            output_prefix = os.path.join(output_dir, sample_name)

            if peak_caller == "MACS2":
                result = run_macs2(
                    bam_info['bam'],
                    output_prefix,
                    bam_info.get('genome', 'hg38'),
                    control_bam=bam_info.get('control'),
                    call_summits=call_summits,
                    broad=broad_peaks,
                    qvalue=qvalue,
                    extra_params=extra_params if 'extra_params' in dir() else ""
                )
            else:  # SEACR
                # First convert BAM to bedGraph
                bg_file = f"{output_prefix}.bedgraph"
                if bam_to_bedgraph(bam_info['bam'], bg_file, bam_info.get('genome', 'hg38')):
                    control_bg = None
                    if bam_info.get('control'):
                        control_bg = f"{output_prefix}_control.bedgraph"
                        bam_to_bedgraph(bam_info['control'], control_bg, bam_info.get('genome', 'hg38'))

                    result = run_seacr(
                        bg_file, output_prefix,
                        control_bedgraph=control_bg,
                        threshold=threshold,
                        norm=seacr_norm,
                        mode=seacr_mode
                    )
                else:
                    result = {'success': False, 'stderr': 'Failed to create bedGraph'}

            # Get stats
            if result['success'] and os.path.exists(result['peaks_file']):
                stats = get_peak_stats(result['peaks_file'])
                result['stats'] = stats

            results.append({
                'sample': sample_name,
                'success': result['success'],
                'peaks_file': result.get('peaks_file', ''),
                'n_peaks': result.get('stats', {}).get('n_peaks', 0),
                'error': result.get('stderr', '') if not result['success'] else ''
            })

            progress.progress((i + 1) / len(sample_selection))

        status.text("Complete!")

        # Show results
        st.markdown("### Results")

        results_df = pd.DataFrame(results)

        # Color code success/failure
        st.dataframe(results_df, use_container_width=True, hide_index=True)

        # Summary
        success_count = sum(r['success'] for r in results)
        st.info(f"✅ {success_count}/{len(results)} samples completed successfully")

        # Save to session state
        if 'called_peaks' not in st.session_state:
            st.session_state.called_peaks = []

        for r in results:
            if r['success']:
                st.session_state.called_peaks.append({
                    'sample': r['sample'],
                    'peaks_file': r['peaks_file'],
                    'n_peaks': r['n_peaks'],
                    'caller': peak_caller
                })


def render_peak_upload():
    """Upload existing peak files."""
    st.header("Upload Existing Peak Files")

    st.markdown("""
    If you already have peak files (from MACS2, SEACR, or other callers),
    you can upload them directly to skip the peak calling step.
    """)

    # File upload
    uploaded_files = st.file_uploader(
        "Upload peak files",
        type=['bed', 'narrowPeak', 'broadPeak', 'csv', 'tsv'],
        accept_multiple_files=True
    )

    if uploaded_files:
        for f in uploaded_files:
            st.markdown(f"**{f.name}**")

            # Preview
            try:
                df = pd.read_csv(f, sep='\t', header=None, nrows=5)
                st.dataframe(df, use_container_width=True)
            except:
                st.warning("Could not preview file")

        # Sample name mapping
        st.subheader("Sample Names")
        sample_names = {}

        for f in uploaded_files:
            default_name = f.name.replace('.bed', '').replace('.narrowPeak', '').replace('.broadPeak', '')
            sample_names[f.name] = st.text_input(f"Name for {f.name}", value=default_name, key=f"name_{f.name}")

        if st.button("Load Peak Files", type="primary"):
            if 'uploaded_peaks' not in st.session_state:
                st.session_state.uploaded_peaks = []

            for f in uploaded_files:
                # Save file
                peaks_dir = "./peaks"
                os.makedirs(peaks_dir, exist_ok=True)
                output_path = os.path.join(peaks_dir, f.name)

                with open(output_path, 'wb') as out:
                    out.write(f.getvalue())

                # Get stats
                stats = get_peak_stats(output_path)

                st.session_state.uploaded_peaks.append({
                    'sample': sample_names[f.name],
                    'peaks_file': output_path,
                    'n_peaks': stats.get('n_peaks', 0)
                })

            st.success(f"Loaded {len(uploaded_files)} peak file(s)!")
            st.rerun()


def render_peak_results():
    """View peak calling results."""
    st.header("Peak Calling Results")

    # Combine called and uploaded peaks
    all_peaks = []

    if 'called_peaks' in st.session_state:
        all_peaks.extend(st.session_state.called_peaks)

    if 'uploaded_peaks' in st.session_state:
        all_peaks.extend(st.session_state.uploaded_peaks)

    if not all_peaks:
        st.info("No peaks available. Run peak calling or upload existing peak files.")
        return

    # Results table
    peaks_df = pd.DataFrame(all_peaks)
    st.dataframe(peaks_df, use_container_width=True, hide_index=True)

    # Visualization
    if len(peaks_df) > 1 and 'n_peaks' in peaks_df.columns:
        st.subheader("Peak Counts by Sample")

        import plotly.express as px
        fig = px.bar(
            peaks_df,
            x='sample',
            y='n_peaks',
            color='n_peaks',
            color_continuous_scale='Viridis',
            title="Number of Peaks per Sample"
        )
        fig.update_layout(height=400)
        st.plotly_chart(fig, use_container_width=True)

    # Next steps
    st.markdown("---")
    st.markdown("### Next Steps")

    col1, col2, col3 = st.columns(3)

    with col1:
        if st.button("📊 Quality Control →", use_container_width=True):
            st.switch_page("pages/02_quality_control.py")

    with col2:
        if st.button("📈 Differential Analysis →", use_container_width=True):
            st.switch_page("pages/03_differential.py")

    with col3:
        if st.button("🏷️ Annotate Peaks →", use_container_width=True):
            st.switch_page("pages/05_annotation.py")


if __name__ == "__main__":
    main()
