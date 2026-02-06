"""
ENCODE Data Integration Page

Browse, download, and integrate public ENCODE datasets:
- ChIP-seq experiments
- Reference epigenomes
- Signal tracks (BigWig)
- Peak files
"""

import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

st.set_page_config(
    page_title="ENCODE Integration - EpiNexus",
    page_icon="🏛️",
    layout="wide"
)

# ENCODE API base URL
ENCODE_API = "https://www.encodeproject.org"


def main():
    st.title("🏛️ ENCODE Data Integration")
    st.markdown("""
    Browse and download public ChIP-seq datasets from ENCODE.
    Compare your data with reference epigenomes.
    """)

    tab1, tab2, tab3, tab4 = st.tabs([
        "🔍 Search ENCODE",
        "📥 Download Data",
        "🔗 Compare with Reference",
        "📊 Public Datasets"
    ])

    with tab1:
        render_encode_search()

    with tab2:
        render_download_manager()

    with tab3:
        render_comparison()

    with tab4:
        render_public_datasets()


def render_encode_search():
    """Search ENCODE database."""
    st.header("Search ENCODE")

    col1, col2, col3 = st.columns(3)

    with col1:
        assay = st.selectbox(
            "Assay Type",
            ["ChIP-seq", "ATAC-seq", "DNase-seq", "RNA-seq", "Hi-C"]
        )

    with col2:
        target = st.selectbox(
            "Target / Histone Mark",
            ["H3K27ac", "H3K4me3", "H3K4me1", "H3K27me3", "H3K36me3",
             "H3K9me3", "CTCF", "POLR2A", "EP300", "Other"]
        )

    with col3:
        organism = st.selectbox(
            "Organism",
            ["Homo sapiens", "Mus musculus", "Drosophila melanogaster"]
        )

    col1, col2 = st.columns(2)

    with col1:
        biosample = st.text_input(
            "Biosample / Cell Type",
            placeholder="e.g., K562, HepG2, GM12878"
        )

    with col2:
        genome = st.selectbox(
            "Genome Assembly",
            ["GRCh38", "hg19", "mm10", "mm39"]
        )

    if st.button("🔍 Search ENCODE", type="primary"):
        with st.spinner("Searching ENCODE database..."):
            results = search_encode_mock(assay, target, organism, biosample, genome)
            st.session_state.encode_results = results

    # Display results
    if 'encode_results' in st.session_state:
        results = st.session_state.encode_results

        st.subheader(f"Found {len(results)} experiments")

        # Filter options
        col1, col2 = st.columns(2)
        with col1:
            status_filter = st.multiselect(
                "Status",
                ["released", "archived"],
                default=["released"]
            )
        with col2:
            lab_filter = st.text_input("Filter by lab")

        # Display results table
        filtered = results[results['status'].isin(status_filter)]
        if lab_filter:
            filtered = filtered[filtered['lab'].str.contains(lab_filter, case=False)]

        st.dataframe(filtered, use_container_width=True, hide_index=True)

        # Select for download
        selected = st.multiselect(
            "Select experiments to download",
            filtered['accession'].tolist()
        )

        if selected:
            st.session_state.selected_encode = selected
            st.success(f"Selected {len(selected)} experiments. Go to 'Download Data' tab.")


def render_download_manager():
    """Manage ENCODE downloads."""
    st.header("Download Manager")

    selected = st.session_state.get('selected_encode', [])

    if not selected:
        st.info("No experiments selected. Use the Search tab to find and select experiments.")

        # Quick download options
        st.subheader("Quick Download - Popular Datasets")

        popular = [
            {"name": "H3K27ac K562", "accession": "ENCSR000AKP", "size": "1.2 GB"},
            {"name": "H3K4me3 K562", "accession": "ENCSR000AKQ", "size": "890 MB"},
            {"name": "H3K27me3 K562", "accession": "ENCSR000AKR", "size": "750 MB"},
            {"name": "CTCF K562", "accession": "ENCSR000AKS", "size": "1.1 GB"},
            {"name": "H3K27ac GM12878", "accession": "ENCSR000AKT", "size": "980 MB"},
        ]

        for ds in popular:
            col1, col2, col3 = st.columns([3, 1, 1])
            with col1:
                st.markdown(f"**{ds['name']}** ({ds['accession']})")
            with col2:
                st.caption(ds['size'])
            with col3:
                if st.button("Add", key=f"add_{ds['accession']}"):
                    if 'selected_encode' not in st.session_state:
                        st.session_state.selected_encode = []
                    st.session_state.selected_encode.append(ds['accession'])
                    st.rerun()

        return

    st.subheader(f"Selected Experiments ({len(selected)})")

    # File type selection
    file_types = st.multiselect(
        "File types to download",
        ["bigWig (signal)", "bed narrowPeak", "bed broadPeak", "bam (alignments)"],
        default=["bigWig (signal)", "bed narrowPeak"]
    )

    # Download location
    download_dir = st.text_input("Download directory", "encode_data/")

    # Show files to download
    st.subheader("Files to Download")

    files_to_download = []
    for acc in selected:
        files = get_experiment_files_mock(acc, file_types)
        files_to_download.extend(files)

    files_df = pd.DataFrame(files_to_download)
    st.dataframe(files_df, use_container_width=True)

    total_size = files_df['size_mb'].sum()
    st.info(f"Total download size: {total_size:.0f} MB ({total_size/1024:.1f} GB)")

    col1, col2 = st.columns(2)

    with col1:
        if st.button("📥 Start Download", type="primary"):
            with st.spinner("Downloading files..."):
                progress = st.progress(0)
                for i, f in enumerate(files_to_download):
                    # Simulate download
                    progress.progress((i + 1) / len(files_to_download),
                                    text=f"Downloading {f['filename']}...")

                st.success(f"Downloaded {len(files_to_download)} files to {download_dir}")

    with col2:
        if st.button("Clear Selection"):
            st.session_state.selected_encode = []
            st.rerun()


def render_comparison():
    """Compare user data with ENCODE reference."""
    st.header("Compare with Reference Epigenomes")

    st.markdown("""
    Upload your peak files and compare overlap with ENCODE reference datasets.
    Useful for validating ChIP-seq quality and identifying cell-type specific signals.
    """)

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Your Data")

        uploaded = st.file_uploader(
            "Upload your peak file",
            type=['bed', 'narrowPeak', 'broadPeak']
        )

        if uploaded:
            st.success(f"Loaded: {uploaded.name}")

    with col2:
        st.subheader("Reference Dataset")

        ref_cell = st.selectbox(
            "Reference cell type",
            ["K562", "GM12878", "HepG2", "H1-hESC", "HUVEC", "IMR90", "NHEK"]
        )

        ref_mark = st.selectbox(
            "Histone mark",
            ["H3K27ac", "H3K4me3", "H3K4me1", "H3K27me3"]
        )

    if st.button("🔗 Compare", type="primary"):
        with st.spinner("Computing overlap..."):
            # Mock comparison results
            results = compute_overlap_mock()

            st.subheader("Overlap Results")

            col1, col2, col3, col4 = st.columns(4)
            col1.metric("Your Peaks", f"{results['user_peaks']:,}")
            col2.metric("Reference Peaks", f"{results['ref_peaks']:,}")
            col3.metric("Overlapping", f"{results['overlap']:,}")
            col4.metric("Jaccard Index", f"{results['jaccard']:.3f}")

            # Overlap visualization
            st.subheader("Overlap Visualization")

            import plotly.graph_objects as go

            fig = go.Figure(go.Sunburst(
                labels=["Total", "Your Peaks Only", "Overlap", "Reference Only"],
                parents=["", "Total", "Total", "Total"],
                values=[results['user_peaks'] + results['ref_peaks'] - results['overlap'],
                       results['user_peaks'] - results['overlap'],
                       results['overlap'],
                       results['ref_peaks'] - results['overlap']],
                branchvalues="total"
            ))
            fig.update_layout(height=400)
            st.plotly_chart(fig, use_container_width=True)

            # Interpretation
            st.subheader("Interpretation")

            if results['jaccard'] > 0.3:
                st.success(f"""
                ✅ **High similarity** with {ref_cell} {ref_mark}

                Your peaks show strong overlap with the reference dataset,
                suggesting your ChIP-seq is of good quality and the cell type
                matches the reference.
                """)
            elif results['jaccard'] > 0.1:
                st.warning(f"""
                ⚠️ **Moderate similarity** with {ref_cell} {ref_mark}

                Some overlap detected. This could indicate:
                - Different cell state or treatment condition
                - Biological differences between samples
                - Technical variation
                """)
            else:
                st.error(f"""
                ❌ **Low similarity** with {ref_cell} {ref_mark}

                Limited overlap with reference. Consider:
                - Verifying cell type identity
                - Checking ChIP-seq quality metrics
                - Comparing with other reference cell types
                """)


def render_public_datasets():
    """Browse curated public datasets."""
    st.header("Curated Public Datasets")

    st.markdown("""
    Pre-processed datasets ready for comparison and integration.
    """)

    # Categories
    category = st.selectbox(
        "Category",
        ["Reference Epigenomes", "Cancer Cell Lines", "Primary Tissues",
         "Developmental Series", "Drug Treatment"]
    )

    if category == "Reference Epigenomes":
        datasets = [
            {"name": "Roadmap Epigenomics - 127 Epigenomes", "samples": 127,
             "marks": "H3K4me1, H3K4me3, H3K27ac, H3K27me3, H3K36me3, H3K9me3",
             "description": "NIH Roadmap comprehensive human epigenome reference"},
            {"name": "ENCODE Tier 1 Cell Lines", "samples": 6,
             "marks": "Full histone panel + TFs",
             "description": "K562, GM12878, H1-hESC, HepG2, HeLa-S3, HUVEC"},
            {"name": "Blueprint Epigenome", "samples": 300,
             "marks": "H3K4me1, H3K4me3, H3K27ac, H3K27me3, H3K9me3",
             "description": "European reference epigenomes for hematopoietic cells"}
        ]

    elif category == "Cancer Cell Lines":
        datasets = [
            {"name": "ENCODE Cancer Cell Lines", "samples": 15,
             "marks": "H3K27ac, H3K4me3, CTCF",
             "description": "Common cancer cell lines from ENCODE"},
            {"name": "CCLE Epigenomics", "samples": 60,
             "marks": "H3K27ac, ATAC-seq",
             "description": "Cancer Cell Line Encyclopedia epigenomic data"}
        ]
    else:
        datasets = [
            {"name": "Sample Dataset", "samples": 10,
             "marks": "H3K27ac", "description": "Example dataset"}
        ]

    for ds in datasets:
        with st.container():
            col1, col2 = st.columns([4, 1])

            with col1:
                st.markdown(f"### {ds['name']}")
                st.markdown(f"**Samples:** {ds['samples']} | **Marks:** {ds['marks']}")
                st.caption(ds['description'])

            with col2:
                if st.button("Load", key=f"load_{ds['name'][:10]}"):
                    st.success(f"Loading {ds['name']}...")

            st.markdown("---")


# Mock functions (replace with real ENCODE API calls)
def search_encode_mock(assay, target, organism, biosample, genome):
    """Mock ENCODE search results."""
    np.random.seed(42)

    n_results = np.random.randint(20, 50)

    results = pd.DataFrame({
        'accession': [f'ENCSR{np.random.randint(100000, 999999)}' for _ in range(n_results)],
        'description': [f'{target} ChIP-seq in {biosample or "various cell types"}' for _ in range(n_results)],
        'biosample': np.random.choice(['K562', 'GM12878', 'HepG2', 'H1-hESC'], n_results),
        'lab': np.random.choice(['Bernstein', 'Snyder', 'Ren', 'Stamatoyannopoulos'], n_results),
        'status': np.random.choice(['released', 'archived'], n_results, p=[0.9, 0.1]),
        'date': pd.date_range('2020-01-01', periods=n_results, freq='W').strftime('%Y-%m-%d')
    })

    return results


def get_experiment_files_mock(accession, file_types):
    """Mock file list for an experiment."""
    files = []

    if "bigWig" in str(file_types):
        files.append({
            'accession': accession,
            'filename': f'{accession}_signal.bigWig',
            'type': 'bigWig',
            'size_mb': np.random.randint(200, 1500)
        })

    if "narrowPeak" in str(file_types):
        files.append({
            'accession': accession,
            'filename': f'{accession}_peaks.narrowPeak',
            'type': 'narrowPeak',
            'size_mb': np.random.randint(1, 50)
        })

    return files


def compute_overlap_mock():
    """Mock overlap computation."""
    return {
        'user_peaks': 35000,
        'ref_peaks': 42000,
        'overlap': 18500,
        'jaccard': 0.32
    }


if __name__ == "__main__":
    main()
