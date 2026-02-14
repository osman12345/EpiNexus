"""
Quality Control Page
QC metrics computed from actual uploaded data.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

st.set_page_config(page_title="Quality Control - EpiNexus", page_icon="✅", layout="wide")

# Try to import data manager
try:
    from frontend.components.data_manager import DataManager
    from frontend.components.empty_states import render_empty_state

    HAS_DATA_MANAGER = True
except ImportError:
    HAS_DATA_MANAGER = False

    def render_empty_state(title, icon, message, requirements):
        st.warning(f"{icon} {title}: {message}")


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


def check_data_loaded():
    """Check if we have data to analyze."""
    if HAS_DATA_MANAGER:
        peaks = DataManager.get_data("peaks")
        if peaks is not None and len(peaks) > 0:
            return True
    return len(st.session_state.get("samples", [])) > 0


def main():
    st.title("✅ Quality Control")
    st.markdown("Assess the quality of your epigenomics data.")

    # Check if data is loaded
    if not check_data_loaded():
        render_empty_state(
            title="No Data Loaded",
            icon="✅",
            message="Upload your data in the Data & Project page first.",
            requirements=[
                "Peak files (BED/narrowPeak format)",
                "At least one sample defined",
                "Optional: BAM files for detailed QC",
            ],
        )

        if st.button("Go to Data & Project →"):
            st.switch_page("pages/01_data_project.py")
        return

    # Get actual data
    peaks_df = None
    if HAS_DATA_MANAGER:
        peaks_df = DataManager.get_data("peaks")

    samples = st.session_state.get("samples", [])

    # Show data source
    is_demo = st.session_state.get("using_demo_data", True)
    if is_demo:
        st.info("📌 **Demo Mode** - Showing metrics for demonstration data")
    else:
        st.success("✅ **Your Data** - Metrics computed from your uploaded peaks")

    # Tabs
    tab1, tab2, tab3, tab4 = st.tabs(["📊 Sample Summary", "📈 Peak Metrics", "🔗 Correlation", "📋 QC Report"])

    with tab1:
        render_sample_summary(peaks_df, samples)
    with tab2:
        render_peak_metrics(peaks_df, samples)
    with tab3:
        render_correlation_analysis(peaks_df, samples)
    with tab4:
        render_qc_report(peaks_df, samples)


def compute_sample_qc(peaks_df: pd.DataFrame, samples: list) -> pd.DataFrame:
    """Compute QC metrics from actual peak data."""
    if peaks_df is None or len(peaks_df) == 0:
        return pd.DataFrame()

    qc_data = []

    # Group peaks by source file if available
    if "source_file" in peaks_df.columns:
        grouped = peaks_df.groupby("source_file")
    else:
        # Create a single group for all peaks
        grouped = [("All Peaks", peaks_df)]

    for source, group in grouped:
        # Extract sample name from filename
        if isinstance(source, str):
            sample_name = Path(source).stem
        else:
            sample_name = str(source)

        # Compute actual metrics from peak data
        n_peaks = len(group)

        # Peak width statistics
        if "start" in group.columns and "end" in group.columns:
            widths = group["end"] - group["start"]
            median_width = widths.median()
            mean_width = widths.mean()
        else:
            median_width = 0
            mean_width = 0

        # Signal strength (if available)
        if "signalValue" in group.columns:
            mean_signal = group["signalValue"].mean()
            max_signal = group["signalValue"].max()
        elif "score" in group.columns:
            mean_signal = group["score"].mean()
            max_signal = group["score"].max()
        else:
            mean_signal = 0
            max_signal = 0

        # Chromosomal distribution
        if "chr" in group.columns:
            n_chroms = group["chr"].nunique()
            # Check for unusual chromosome distribution
            chr_counts = group["chr"].value_counts()
            chr_balance = chr_counts.std() / chr_counts.mean() if chr_counts.mean() > 0 else 0
        else:
            n_chroms = 0
            chr_balance = 0

        # Quality score estimation based on metrics
        # Higher peaks, better signal, good distribution = better quality
        quality_score = min(
            100,
            int(
                (min(n_peaks / 1000, 30))  # Peak count contribution (max 30)
                + (min(mean_signal / 10, 30) if mean_signal > 0 else 15)  # Signal contribution (max 30)
                + (40 - min(chr_balance * 20, 20))  # Balance contribution (max 40)
            ),
        )

        # Determine status
        if quality_score >= 70:
            status = "✅ Pass"
        elif quality_score >= 50:
            status = "⚠️ Warning"
        else:
            status = "❌ Fail"

        qc_data.append(
            {
                "Sample": sample_name,
                "Peaks": n_peaks,
                "Median Width": int(median_width),
                "Mean Width": int(mean_width),
                "Mean Signal": round(mean_signal, 1),
                "Max Signal": round(max_signal, 1),
                "Chromosomes": n_chroms,
                "Quality Score": quality_score,
                "Status": status,
            }
        )

    return pd.DataFrame(qc_data)


def render_sample_summary(peaks_df, samples):
    """Sample quality summary computed from real data."""
    st.header("Sample Quality Summary")

    # Compute QC from actual data
    qc_df = compute_sample_qc(peaks_df, samples)

    if len(qc_df) == 0:
        st.warning("No peak data available for QC analysis.")

        # Show sample list if available
        if samples:
            st.markdown("**Defined samples:**")
            st.dataframe(pd.DataFrame(samples), use_container_width=True)
        return

    # Summary metrics
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Total Samples", len(qc_df))
    col2.metric("Passing QC", (qc_df["Status"] == "✅ Pass").sum())
    col3.metric("Total Peaks", f"{qc_df['Peaks'].sum():,}")
    col4.metric("Avg Quality Score", f"{qc_df['Quality Score'].mean():.0f}")

    st.markdown("---")

    # Full QC table
    st.dataframe(
        qc_df,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Sample": st.column_config.TextColumn("Sample"),
            "Peaks": st.column_config.NumberColumn("Peaks", format="%d"),
            "Median Width": st.column_config.NumberColumn("Median Width (bp)"),
            "Mean Signal": st.column_config.NumberColumn("Mean Signal"),
            "Quality Score": st.column_config.ProgressColumn("Quality Score", min_value=0, max_value=100, format="%d"),
            "Status": st.column_config.TextColumn("Status"),
        },
    )

    # QC thresholds explanation
    with st.expander("ℹ️ Quality Score Calculation"):
        st.markdown("""
        **Quality Score (0-100) is computed from:**
        - **Peak count** (max 30 pts): More peaks generally indicate better enrichment
        - **Signal strength** (max 30 pts): Higher signal values indicate cleaner data
        - **Chromosomal balance** (max 40 pts): Even distribution across chromosomes

        **Status thresholds:**
        - ✅ **Pass**: Score ≥ 70
        - ⚠️ **Warning**: Score 50-69
        - ❌ **Fail**: Score < 50

        *For CUT&Tag/CUT&RUN data without IgG controls, these metrics help assess data quality
        without traditional FRiP calculations.*
        """)


def render_peak_metrics(peaks_df, samples):
    """Peak-level quality metrics from real data."""
    st.header("Peak Quality Metrics")

    if peaks_df is None or len(peaks_df) == 0:
        st.warning("No peak data available.")
        return

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Peak Width Distribution")

        if "start" in peaks_df.columns and "end" in peaks_df.columns:
            widths = peaks_df["end"] - peaks_df["start"]

            # Filter reasonable widths for visualization
            widths_filtered = widths[(widths > 0) & (widths < 10000)]

            fig = go.Figure(data=go.Histogram(x=widths_filtered, nbinsx=50, marker_color="#3498db"))

            fig.update_layout(xaxis_title="Peak Width (bp)", yaxis_title="Count", height=400)

            # Add median line
            median_width = widths_filtered.median()
            fig.add_vline(
                x=median_width, line_dash="dash", line_color="red", annotation_text=f"Median: {median_width:.0f} bp"
            )

            st.plotly_chart(fig, use_container_width=True)

            # Statistics
            st.markdown(f"""
            **Peak width statistics:**
            - Median: {widths.median():.0f} bp
            - Mean: {widths.mean():.0f} bp
            - Min: {widths.min():.0f} bp
            - Max: {widths.max():.0f} bp
            """)
        else:
            st.info("Peak coordinates (start/end) not available.")

    with col2:
        st.subheader("Chromosomal Distribution")

        if "chr" in peaks_df.columns:
            chr_counts = peaks_df["chr"].value_counts()

            # Sort chromosomes naturally
            def chr_sort_key(x):
                x = str(x).replace("chr", "")
                if x.isdigit():
                    return (0, int(x))
                return (1, x)

            sorted_chrs = sorted(chr_counts.index, key=chr_sort_key)
            chr_counts = chr_counts[sorted_chrs]

            fig = px.bar(
                x=chr_counts.index[:24],  # Limit to main chromosomes
                y=chr_counts.values[:24],
                labels={"x": "Chromosome", "y": "Peak Count"},
            )
            fig.update_layout(height=400)
            st.plotly_chart(fig, use_container_width=True)

            # Check for chromosome balance
            main_chrs = [c for c in chr_counts.index if str(c).replace("chr", "").isdigit()]
            if main_chrs:
                main_counts = chr_counts[main_chrs]
                cv = main_counts.std() / main_counts.mean() if main_counts.mean() > 0 else 0
                if cv < 0.5:
                    st.success("✅ Peaks are well-distributed across chromosomes")
                elif cv < 1.0:
                    st.warning("⚠️ Some chromosomal bias detected")
                else:
                    st.error("❌ Strong chromosomal bias - check data quality")
        else:
            st.info("Chromosome information not available in peak data.")

    # Signal distribution
    st.markdown("---")
    st.subheader("Signal Distribution")

    signal_col = None
    if "signalValue" in peaks_df.columns:
        signal_col = "signalValue"
    elif "score" in peaks_df.columns:
        signal_col = "score"

    if signal_col:
        col1, col2 = st.columns(2)

        with col1:
            # Log-scale histogram of signal
            signal_values = peaks_df[signal_col].dropna()
            signal_values = signal_values[signal_values > 0]

            fig = go.Figure(data=go.Histogram(x=np.log10(signal_values + 1), nbinsx=50, marker_color="#2ecc71"))
            fig.update_layout(xaxis_title=f"Log10({signal_col} + 1)", yaxis_title="Count", height=350)
            st.plotly_chart(fig, use_container_width=True)

        with col2:
            # Signal statistics
            st.markdown(f"""
            **Signal statistics ({signal_col}):**
            - Mean: {signal_values.mean():.2f}
            - Median: {signal_values.median():.2f}
            - Std Dev: {signal_values.std():.2f}
            - Max: {signal_values.max():.2f}

            **Dynamic range:** {np.log10(signal_values.max() / (signal_values.median() + 1)):.1f} log10 units
            """)

            # Quality assessment based on signal
            if signal_values.mean() > 10 and signal_values.max() > 100:
                st.success("✅ Good signal dynamic range")
            elif signal_values.mean() > 5:
                st.info("ℹ️ Moderate signal - typical for some histone marks")
            else:
                st.warning("⚠️ Low signal - may indicate weak enrichment")
    else:
        st.info("Signal values not available in peak data. Consider using narrowPeak format for signal information.")


def render_correlation_analysis(peaks_df, samples):
    """Sample correlation based on peak overlap."""
    st.header("Sample Correlation")

    if peaks_df is None or "source_file" not in peaks_df.columns:
        st.info("Correlation analysis requires multiple samples with peak data.")
        st.markdown("""
        **To enable correlation analysis:**
        1. Upload peak files for multiple samples
        2. Each file should contain peaks from one sample
        """)
        return

    # Get unique sources
    sources = peaks_df["source_file"].unique()

    if len(sources) < 2:
        st.info("Need at least 2 samples for correlation analysis.")
        return

    st.subheader("Peak Overlap Matrix")

    # Compute overlap between samples
    overlap_matrix = compute_peak_overlap_matrix(peaks_df, sources)

    if overlap_matrix is not None:
        # Plot heatmap
        sample_names = [Path(s).stem for s in sources]

        fig = px.imshow(
            overlap_matrix,
            x=sample_names,
            y=sample_names,
            color_continuous_scale="RdBu_r",
            zmin=0,
            zmax=1,
            aspect="auto",
            labels={"color": "Jaccard Index"},
        )
        fig.update_layout(height=500)
        st.plotly_chart(fig, use_container_width=True)

        # Interpretation
        avg_overlap = np.mean(overlap_matrix[np.triu_indices(len(sources), k=1)])
        st.markdown(f"**Average pairwise overlap (Jaccard):** {avg_overlap:.2%}")

        if avg_overlap > 0.5:
            st.success("✅ High reproducibility between samples")
        elif avg_overlap > 0.3:
            st.info("ℹ️ Moderate overlap - typical for biological replicates")
        else:
            st.warning("⚠️ Low overlap - samples may be quite different or have quality issues")


@st.cache_data(show_spinner="Computing peak overlap matrix...")
def compute_peak_overlap_matrix(peaks_df, sources):
    """Compute Jaccard overlap between samples."""
    n = len(sources)
    matrix = np.eye(n)

    for i, src_i in enumerate(sources):
        peaks_i = peaks_df[peaks_df["source_file"] == src_i]
        if "chr" not in peaks_i.columns or "start" not in peaks_i.columns:
            continue

        set_i = set(zip(peaks_i["chr"], peaks_i["start"] // 500))  # Bin by 500bp windows

        for j, src_j in enumerate(sources):
            if j <= i:
                continue

            peaks_j = peaks_df[peaks_df["source_file"] == src_j]
            set_j = set(zip(peaks_j["chr"], peaks_j["start"] // 500))

            # Jaccard index
            intersection = len(set_i & set_j)
            union = len(set_i | set_j)
            jaccard = intersection / union if union > 0 else 0

            matrix[i, j] = jaccard
            matrix[j, i] = jaccard

    return matrix


def render_qc_report(peaks_df, samples):
    """Generate downloadable QC report."""
    st.header("QC Report")

    qc_df = compute_sample_qc(peaks_df, samples)

    if len(qc_df) == 0:
        st.warning("No data available for QC report.")
        return

    # Summary
    st.subheader("Summary")

    total_samples = len(qc_df)
    passing = (qc_df["Status"] == "✅ Pass").sum()
    warning = (qc_df["Status"] == "⚠️ Warning").sum()
    failing = (qc_df["Status"] == "❌ Fail").sum()

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Total Samples", total_samples)
    col2.metric("Passing", passing, delta=None if passing == total_samples else f"-{total_samples - passing}")
    col3.metric("Warnings", warning)
    col4.metric("Failing", failing)

    # Recommendations
    st.subheader("Recommendations")

    if failing > 0:
        st.error(f"""
        ❌ **{failing} sample(s) failed QC**

        Consider:
        - Checking peak calling parameters
        - Reviewing raw data quality
        - Re-running with different settings
        """)

    if warning > 0:
        st.warning(f"""
        ⚠️ **{warning} sample(s) have warnings**

        These samples may still be usable but review:
        - Peak count compared to expected
        - Signal strength
        - Chromosomal distribution
        """)

    if passing == total_samples:
        st.success("✅ **All samples passed QC** - Ready for downstream analysis!")

    # Record workflow step
    if HAS_WORKFLOW_MANAGER:
        WorkflowManager.record_step(
            step_name="Quality Control",
            tool="EpiNexus QC Module",
            parameters={
                "total_samples": total_samples,
                "passing": passing,
                "warnings": warning,
                "failing": failing,
                "avg_quality_score": float(qc_df["Quality Score"].mean()),
            },
            inputs=["peaks"],
            outputs=["qc_report"],
        )

    # Download report
    st.markdown("---")
    st.subheader("Download Report")

    # Generate report content
    report_content = f"""# EpiNexus QC Report

## Project: {st.session_state.get("project_name", "Unknown")}
## Date: {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M")}

## Summary
- Total Samples: {total_samples}
- Passing: {passing}
- Warnings: {warning}
- Failing: {failing}

## Sample Details

{qc_df.to_markdown(index=False)}

## Quality Metrics Explanation

- **Peaks**: Total number of called peaks
- **Median Width**: Median peak width in base pairs
- **Mean Signal**: Average signal value across peaks
- **Quality Score**: Composite score (0-100) based on multiple metrics
- **Status**: Pass (≥70), Warning (50-69), Fail (<50)

## Assay: {st.session_state.get("assay_type", "Unknown")}
## Genome: {st.session_state.get("selected_genome", "Unknown")}
"""

    col1, col2 = st.columns(2)

    with col1:
        st.download_button(
            "📄 Download Report (Markdown)", report_content, "qc_report.md", "text/markdown", use_container_width=True
        )

    with col2:
        st.download_button(
            "📊 Download QC Table (CSV)",
            qc_df.to_csv(index=False),
            "qc_metrics.csv",
            "text/csv",
            use_container_width=True,
        )


if __name__ == "__main__":
    main()
