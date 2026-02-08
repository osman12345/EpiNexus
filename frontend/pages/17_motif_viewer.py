"""
Motif Logo Visualization Page

Interactive visualization of transcription factor binding motifs.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from app.core.motif_logo import MotifLogoGenerator, MotifMatrix, KNOWN_MOTIFS, list_available_motifs

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

st.set_page_config(
    page_title="Motif Viewer - EpiNexus",
    page_icon="🔤",
    layout="wide"
)


def main():
    st.title("🔤 Motif Logo Viewer")
    st.markdown("""
    Visualize transcription factor binding motifs as sequence logos.
    Information content at each position indicates conservation.
    """)

    tab1, tab2, tab3 = st.tabs([
        "📊 View Motifs",
        "🔀 Compare Motifs",
        "📤 Custom Motif"
    ])

    with tab1:
        render_motif_viewer()

    with tab2:
        render_motif_comparison()

    with tab3:
        render_custom_motif()


def render_motif_viewer():
    """View pre-defined motifs."""
    st.header("Motif Logo Viewer")

    col1, col2 = st.columns([1, 3])

    with col1:
        st.subheader("Select Motif")

        motif_name = st.selectbox(
            "Transcription Factor",
            list_available_motifs()
        )

        color_scheme = st.selectbox(
            "Color Scheme",
            ['weblogo', 'classic', 'colorblind', 'grayscale']
        )

        st.markdown("---")

        # Motif info
        motif = KNOWN_MOTIFS.get(motif_name)
        if motif:
            st.markdown(f"**Consensus:** `{motif.consensus}`")
            st.markdown(f"**Length:** {motif.length} bp")
            st.markdown(f"**Source:** {motif.source}")

    with col2:
        if motif:
            generator = MotifLogoGenerator(color_scheme)
            fig = generator.generate_plotly_logo(motif, width=700, height=300)
            st.plotly_chart(fig, use_container_width=True)

            # Show matrix
            with st.expander("View Position Weight Matrix"):
                matrix_df = pd.DataFrame(
                    motif.matrix,
                    index=['A', 'C', 'G', 'T'],
                    columns=[f"Pos {i+1}" for i in range(motif.length)]
                )
                st.dataframe(matrix_df.round(3), use_container_width=True)

    # Multiple motif gallery
    st.markdown("---")
    st.subheader("Motif Gallery")

    cols = st.columns(2)
    generator = MotifLogoGenerator(color_scheme)

    for i, (name, motif) in enumerate(KNOWN_MOTIFS.items()):
        with cols[i % 2]:
            fig = generator.generate_plotly_logo(motif, width=400, height=200)
            fig.update_layout(title=f"{name} ({motif.consensus})")
            st.plotly_chart(fig, use_container_width=True)


def render_motif_comparison():
    """Compare two motifs."""
    st.header("Motif Comparison")

    col1, col2 = st.columns(2)

    with col1:
        motif1_name = st.selectbox(
            "Motif 1",
            list_available_motifs(),
            key='compare_m1'
        )

    with col2:
        motif2_name = st.selectbox(
            "Motif 2",
            list_available_motifs(),
            index=1,
            key='compare_m2'
        )

    motif1 = KNOWN_MOTIFS.get(motif1_name)
    motif2 = KNOWN_MOTIFS.get(motif2_name)

    if motif1 and motif2:
        generator = MotifLogoGenerator()

        # Show both logos
        col1, col2 = st.columns(2)

        with col1:
            fig1 = generator.generate_plotly_logo(motif1, width=400, height=250)
            st.plotly_chart(fig1, use_container_width=True)

        with col2:
            fig2 = generator.generate_plotly_logo(motif2, width=400, height=250)
            st.plotly_chart(fig2, use_container_width=True)

        # Comparison metrics
        st.subheader("Similarity Metrics")

        comparison = generator.compare_motifs(motif1, motif2)

        col1, col2, col3, col4 = st.columns(4)

        with col1:
            corr = comparison['pearson_correlation']
            st.metric(
                "Pearson Correlation",
                f"{corr:.3f}",
                delta="Similar" if corr > 0.7 else "Different"
            )

        with col2:
            st.metric(
                "Euclidean Distance",
                f"{comparison['euclidean_distance']:.3f}"
            )

        with col3:
            st.metric(
                "KL Divergence",
                f"{comparison['kl_divergence']:.3f}"
            )

        with col4:
            st.metric(
                "Length Difference",
                f"{comparison['length_diff']} bp"
            )

        # Interpretation
        st.markdown("---")
        if comparison['pearson_correlation'] > 0.8:
            st.success("✅ **High similarity** - These motifs likely recognize similar DNA sequences.")
        elif comparison['pearson_correlation'] > 0.5:
            st.warning("⚠️ **Moderate similarity** - Some overlap in binding preferences.")
        else:
            st.info("ℹ️ **Low similarity** - These TFs likely recognize distinct sequences.")


def render_custom_motif():
    """Upload or create custom motif."""
    st.header("Custom Motif")

    input_method = st.radio(
        "Input Method",
        ["Paste Matrix", "Upload File", "From Sequences"],
        horizontal=True
    )

    if input_method == "Paste Matrix":
        st.markdown("""
        Enter a Position Frequency Matrix (4 rows for A, C, G, T):
        """)

        matrix_text = st.text_area(
            "Matrix (tab or comma separated)",
            value="0.1, 0.7, 0.1, 0.1, 0.1, 0.7\n0.7, 0.1, 0.8, 0.8, 0.1, 0.1\n0.1, 0.1, 0.1, 0.1, 0.8, 0.1\n0.1, 0.1, 0.0, 0.0, 0.0, 0.1",
            height=120
        )

        motif_name = st.text_input("Motif Name", "Custom Motif")

        if st.button("Generate Logo"):
            try:
                # Parse matrix
                lines = matrix_text.strip().split('\n')
                matrix = []
                for line in lines:
                    values = [float(x.strip()) for x in line.replace('\t', ',').split(',')]
                    matrix.append(values)

                matrix = np.array(matrix)

                if matrix.shape[0] != 4:
                    st.error("Matrix must have exactly 4 rows (A, C, G, T)")
                else:
                    custom_motif = MotifMatrix(
                        name=motif_name,
                        matrix=matrix,
                        source="User input"
                    )

                    generator = MotifLogoGenerator()
                    fig = generator.generate_plotly_logo(custom_motif, width=600, height=300)
                    st.plotly_chart(fig, use_container_width=True)

                    st.success(f"Consensus: {custom_motif.consensus}")

            except Exception as e:
                st.error(f"Error parsing matrix: {e}")

    elif input_method == "Upload File":
        st.markdown("Upload a JASPAR, MEME, or TRANSFAC format file.")

        uploaded = st.file_uploader(
            "Upload motif file",
            type=['txt', 'jaspar', 'meme', 'pfm']
        )

        if uploaded:
            st.info(f"Uploaded: {uploaded.name}")
            st.warning("File parsing not yet implemented. Use 'Paste Matrix' for now.")

    elif input_method == "From Sequences":
        st.markdown("Generate a motif from aligned sequences.")

        sequences = st.text_area(
            "Enter aligned sequences (one per line)",
            value="CACGTG\nCACGTG\nCACGTG\nCATGTG\nCAGGTG",
            height=150
        )

        if st.button("Build Motif"):
            seqs = [s.strip().upper() for s in sequences.strip().split('\n') if s.strip()]

            if len(seqs) < 2:
                st.error("Need at least 2 sequences")
            elif len(set(len(s) for s in seqs)) > 1:
                st.error("All sequences must be the same length")
            else:
                # Build PFM
                length = len(seqs[0])
                pfm = np.zeros((4, length))
                base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

                for seq in seqs:
                    for pos, base in enumerate(seq):
                        if base in base_to_idx:
                            pfm[base_to_idx[base], pos] += 1

                custom_motif = MotifMatrix(
                    name="From Sequences",
                    matrix=pfm,
                    source=f"Built from {len(seqs)} sequences"
                )

                generator = MotifLogoGenerator()
                fig = generator.generate_plotly_logo(custom_motif, width=600, height=300)
                st.plotly_chart(fig, use_container_width=True)

                st.success(f"Consensus: {custom_motif.consensus}")

                # Show PFM
                with st.expander("View Position Frequency Matrix"):
                    matrix_df = pd.DataFrame(
                        pfm,
                        index=['A', 'C', 'G', 'T'],
                        columns=[f"Pos {i+1}" for i in range(length)]
                    )
                    st.dataframe(matrix_df.astype(int), use_container_width=True)


if __name__ == "__main__":
    main()
