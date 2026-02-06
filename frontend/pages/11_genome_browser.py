"""
Interactive Genome Browser - IGV.js Integration

Embedded genome browser for visualizing:
- ChIP-seq signal tracks
- Peak regions
- Gene annotations
- Custom BED files
"""

import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
from pathlib import Path
import json
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

st.set_page_config(
    page_title="Genome Browser - EpiNexus",
    page_icon="🧬",
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
            <div style="font-size: 3rem; margin-bottom: 1rem;">🧬</div>
            <h2 style="color: #6c757d; margin-bottom: 0.5rem;">No Data Loaded</h2>
            <p style="color: #6c757d; font-size: 1.1rem;">
                Upload your data to visualize in the genome browser.
            </p>
        </div>
        """, unsafe_allow_html=True)
        st.markdown("")
        if st.button("📁 Go to Data & Project", type="primary", use_container_width=True):
            st.switch_page("pages/01_data_project.py")
        st.markdown("")
        st.markdown("**Browser features:**")
        st.markdown("- IGV.js interactive visualization")
        st.markdown("- Multiple signal tracks")
        st.markdown("- Peak annotations")
        st.markdown("- Gene search & navigation")


def main():
    st.title("🧬 Interactive Genome Browser")
    st.markdown("Visualize your data in an interactive IGV.js genome browser.")

    # Check if data is loaded
    if not has_data():
        render_empty_state()
        return

    # Sidebar controls
    with st.sidebar:
        st.header("Navigation")

        genome = st.selectbox(
            "Reference Genome",
            ["hg38", "hg19", "mm10", "mm39", "dm6"],
            index=0
        )

        # Gene/region search
        st.subheader("Go to Region")

        search_type = st.radio("Search by", ["Gene", "Coordinates"], horizontal=True)

        if search_type == "Gene":
            gene = st.text_input("Gene symbol", value="MYC")
            if st.button("Go to Gene"):
                st.session_state.igv_locus = gene
        else:
            col1, col2 = st.columns(2)
            with col1:
                chrom = st.selectbox("Chr", [f"chr{i}" for i in list(range(1, 23)) + ['X', 'Y']])
            with col2:
                start = st.number_input("Start", value=127735434, step=1000)
                end = st.number_input("End", value=127742951, step=1000)

            if st.button("Go to Region"):
                st.session_state.igv_locus = f"{chrom}:{start}-{end}"

        st.markdown("---")

        # Track management
        st.subheader("Add Tracks")

        track_type = st.selectbox(
            "Track type",
            ["BigWig Signal", "BED Peaks", "VCF Variants", "BAM Alignments"]
        )

        track_url = st.text_input("Track URL or file path")
        track_name = st.text_input("Track name")
        track_color = st.color_picker("Color", "#1f77b4")

        if st.button("Add Track"):
            if 'custom_tracks' not in st.session_state:
                st.session_state.custom_tracks = []

            st.session_state.custom_tracks.append({
                "name": track_name,
                "url": track_url,
                "type": track_type,
                "color": track_color
            })
            st.success(f"Added track: {track_name}")

    # Main browser area
    render_igv_browser(genome)

    # Quick navigation buttons
    st.markdown("---")
    render_quick_navigation()

    # Track list
    render_track_list()


def render_igv_browser(genome: str):
    """Render IGV.js browser component."""

    locus = st.session_state.get('igv_locus', 'MYC')
    custom_tracks = st.session_state.get('custom_tracks', [])

    # Build tracks configuration
    tracks_config = []

    # Add demo tracks
    demo_tracks = get_demo_tracks(genome)
    tracks_config.extend(demo_tracks)

    # Add custom tracks
    for track in custom_tracks:
        track_conf = {
            "name": track["name"],
            "url": track["url"],
            "color": track["color"]
        }

        if "BigWig" in track["type"]:
            track_conf["type"] = "wig"
        elif "BED" in track["type"]:
            track_conf["type"] = "annotation"
            track_conf["format"] = "bed"
        elif "VCF" in track["type"]:
            track_conf["type"] = "variant"
        elif "BAM" in track["type"]:
            track_conf["type"] = "alignment"

        tracks_config.append(track_conf)

    # IGV.js HTML component
    igv_html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <script src="https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js"></script>
        <style>
            body {{ margin: 0; padding: 0; }}
            #igv-div {{ width: 100%; height: 600px; }}
        </style>
    </head>
    <body>
        <div id="igv-div"></div>
        <script>
            var igvDiv = document.getElementById("igv-div");

            var options = {{
                genome: "{genome}",
                locus: "{locus}",
                tracks: {json.dumps(tracks_config)}
            }};

            igv.createBrowser(igvDiv, options)
                .then(function (browser) {{
                    console.log("IGV browser created");
                }});
        </script>
    </body>
    </html>
    """

    components.html(igv_html, height=650, scrolling=False)


def get_demo_tracks(genome: str):
    """Get demo tracks from ENCODE for the selected genome."""

    if genome == "hg38":
        return [
            {
                "name": "H3K27ac (ENCODE K562)",
                "type": "wig",
                "url": "https://www.encodeproject.org/files/ENCFF469RDC/@@download/ENCFF469RDC.bigWig",
                "color": "#f39c12"
            },
            {
                "name": "H3K4me3 (ENCODE K562)",
                "type": "wig",
                "url": "https://www.encodeproject.org/files/ENCFF422UFO/@@download/ENCFF422UFO.bigWig",
                "color": "#e74c3c"
            },
            {
                "name": "CTCF (ENCODE K562)",
                "type": "wig",
                "url": "https://www.encodeproject.org/files/ENCFF368TYI/@@download/ENCFF368TYI.bigWig",
                "color": "#3498db"
            }
        ]
    elif genome == "hg19":
        return [
            {
                "name": "H3K27ac (ENCODE)",
                "type": "wig",
                "url": "https://www.encodeproject.org/files/ENCFF001WIO/@@download/ENCFF001WIO.bigWig",
                "color": "#f39c12"
            }
        ]
    elif genome in ["mm10", "mm39"]:
        return [
            {
                "name": "H3K27ac (ENCODE)",
                "type": "wig",
                "url": "https://www.encodeproject.org/files/ENCFF029NCA/@@download/ENCFF029NCA.bigWig",
                "color": "#f39c12"
            }
        ]
    else:
        return []


def render_quick_navigation():
    """Quick navigation to common loci."""
    st.subheader("Quick Navigation")

    col1, col2, col3, col4, col5 = st.columns(5)

    genes = {
        "MYC": "chr8:127,735,434-127,742,951",
        "TP53": "chr17:7,668,421-7,687,490",
        "BRCA1": "chr17:43,044,295-43,170,245",
        "EGFR": "chr7:55,019,017-55,211,628",
        "GAPDH": "chr12:6,534,517-6,538,371"
    }

    for col, (gene, coords) in zip([col1, col2, col3, col4, col5], genes.items()):
        with col:
            if st.button(f"🧬 {gene}"):
                st.session_state.igv_locus = gene
                st.rerun()

    # Also show enhancer/promoter regions
    st.markdown("**Example regulatory regions:**")

    col1, col2, col3 = st.columns(3)

    with col1:
        if st.button("Super-enhancer (MYC)"):
            st.session_state.igv_locus = "chr8:128,700,000-128,900,000"
            st.rerun()

    with col2:
        if st.button("Bivalent promoter"):
            st.session_state.igv_locus = "chr12:54,350,000-54,380,000"
            st.rerun()

    with col3:
        if st.button("Polycomb domain"):
            st.session_state.igv_locus = "chr7:27,100,000-27,300,000"
            st.rerun()


def render_track_list():
    """Display and manage loaded tracks."""
    st.subheader("Loaded Tracks")

    custom_tracks = st.session_state.get('custom_tracks', [])

    if not custom_tracks:
        st.info("No custom tracks loaded. Add tracks from the sidebar or use demo tracks from ENCODE.")
        return

    for i, track in enumerate(custom_tracks):
        col1, col2, col3 = st.columns([3, 1, 1])

        with col1:
            st.markdown(f"**{track['name']}** ({track['type']})")
            st.caption(track['url'][:60] + "..." if len(track['url']) > 60 else track['url'])

        with col2:
            st.color_picker("", track['color'], key=f"color_{i}", disabled=True)

        with col3:
            if st.button("🗑️", key=f"remove_{i}"):
                st.session_state.custom_tracks.pop(i)
                st.rerun()


# Additional utility: Bookmark regions
def render_bookmarks():
    """Manage bookmarked regions."""
    st.markdown("---")
    st.subheader("📌 Bookmarked Regions")

    if 'bookmarks' not in st.session_state:
        st.session_state.bookmarks = []

    # Add bookmark
    with st.expander("Add Bookmark"):
        bm_name = st.text_input("Region name")
        bm_locus = st.text_input("Coordinates", placeholder="chr1:1000000-2000000")
        bm_notes = st.text_area("Notes")

        if st.button("Save Bookmark"):
            st.session_state.bookmarks.append({
                "name": bm_name,
                "locus": bm_locus,
                "notes": bm_notes
            })
            st.success("Bookmark saved!")

    # Display bookmarks
    for i, bm in enumerate(st.session_state.bookmarks):
        col1, col2 = st.columns([4, 1])
        with col1:
            if st.button(f"📍 {bm['name']}", key=f"bm_{i}"):
                st.session_state.igv_locus = bm['locus']
                st.rerun()
            st.caption(f"{bm['locus']} - {bm['notes'][:50]}")
        with col2:
            if st.button("🗑️", key=f"rm_bm_{i}"):
                st.session_state.bookmarks.pop(i)
                st.rerun()


if __name__ == "__main__":
    main()
    render_bookmarks()
