"""
Publication-Ready Figure Export

Provides:
- High-resolution PNG/PDF/SVG export
- Journal-specific formatting presets
- Figure legends and annotations
- Multi-panel figure assembly
"""

from typing import Dict, List, Optional, Tuple
import plotly.graph_objects as go
from dataclasses import dataclass


@dataclass
class JournalPreset:
    """Journal-specific figure requirements."""

    name: str
    max_width_mm: float  # Single column width
    max_width_full_mm: float  # Full page width
    min_dpi: int
    font_family: str
    min_font_size: int
    formats: List[str]
    color_mode: str  # 'rgb' or 'cmyk'


# Journal presets
JOURNAL_PRESETS = {
    "nature": JournalPreset(
        name="Nature",
        max_width_mm=89,
        max_width_full_mm=183,
        min_dpi=300,
        font_family="Arial",
        min_font_size=5,
        formats=["pdf", "eps", "tiff"],
        color_mode="cmyk",
    ),
    "cell": JournalPreset(
        name="Cell",
        max_width_mm=85,
        max_width_full_mm=174,
        min_dpi=300,
        font_family="Arial",
        min_font_size=6,
        formats=["pdf", "eps", "tiff"],
        color_mode="cmyk",
    ),
    "science": JournalPreset(
        name="Science",
        max_width_mm=57,
        max_width_full_mm=120,
        min_dpi=300,
        font_family="Helvetica",
        min_font_size=6,
        formats=["pdf", "eps"],
        color_mode="rgb",
    ),
    "plos": JournalPreset(
        name="PLOS",
        max_width_mm=83,
        max_width_full_mm=173,
        min_dpi=300,
        font_family="Arial",
        min_font_size=8,
        formats=["tiff", "eps", "pdf"],
        color_mode="rgb",
    ),
    "biorxiv": JournalPreset(
        name="bioRxiv (Preprint)",
        max_width_mm=180,
        max_width_full_mm=180,
        min_dpi=150,
        font_family="Arial",
        min_font_size=8,
        formats=["png", "pdf", "svg"],
        color_mode="rgb",
    ),
    "presentation": JournalPreset(
        name="Presentation",
        max_width_mm=254,
        max_width_full_mm=254,
        min_dpi=150,
        font_family="Arial",
        min_font_size=14,
        formats=["png", "svg"],
        color_mode="rgb",
    ),
}


class FigureExporter:
    """Export publication-ready figures."""

    def __init__(self, preset: str = "nature"):
        self.preset = JOURNAL_PRESETS.get(preset, JOURNAL_PRESETS["nature"])

    def set_preset(self, preset_name: str):
        """Change journal preset."""
        if preset_name in JOURNAL_PRESETS:
            self.preset = JOURNAL_PRESETS[preset_name]

    def format_figure(
        self,
        fig: go.Figure,
        width: str = "single",  # 'single' or 'full'
        height_mm: Optional[float] = None,
        title: str = None,
        add_panel_label: str = None,  # e.g., 'A', 'B', 'C'
    ) -> go.Figure:
        """Format figure for publication."""

        # Calculate dimensions
        width_mm = self.preset.max_width_mm if width == "single" else self.preset.max_width_full_mm
        width_px = int(width_mm * self.preset.min_dpi / 25.4)

        if height_mm:
            height_px = int(height_mm * self.preset.min_dpi / 25.4)
        else:
            # Default aspect ratio
            height_px = int(width_px * 0.75)

        # Update layout
        fig.update_layout(
            width=width_px,
            height=height_px,
            font=dict(family=self.preset.font_family, size=max(self.preset.min_font_size, 10)),
            title=dict(text=title, font=dict(size=max(self.preset.min_font_size + 2, 12))) if title else None,
            margin=dict(l=50, r=20, t=40 if title else 20, b=50),
        )

        # Add panel label
        if add_panel_label:
            fig.add_annotation(
                text=f"<b>{add_panel_label}</b>",
                xref="paper",
                yref="paper",
                x=-0.05,
                y=1.05,
                showarrow=False,
                font=dict(size=14, family=self.preset.font_family),
            )

        return fig

    def export_png(self, fig: go.Figure, scale: float = 2.0) -> bytes:
        """Export figure as PNG."""
        return fig.to_image(format="png", scale=scale)

    def export_svg(self, fig: go.Figure) -> str:
        """Export figure as SVG."""
        return fig.to_image(format="svg").decode("utf-8")

    def export_pdf(self, fig: go.Figure) -> bytes:
        """Export figure as PDF."""
        return fig.to_image(format="pdf")

    def get_download_link(
        self, fig: go.Figure, filename: str, format: str = "png", scale: float = 2.0
    ) -> Tuple[bytes, str]:
        """Get downloadable figure data."""

        if format == "png":
            data = self.export_png(fig, scale)
            mime = "image/png"
        elif format == "svg":
            data = self.export_svg(fig).encode("utf-8")
            mime = "image/svg+xml"
        elif format == "pdf":
            data = self.export_pdf(fig)
            mime = "application/pdf"
        else:
            raise ValueError(f"Unsupported format: {format}")

        return data, mime

    def create_figure_legend(self, items: List[Dict[str, str]]) -> str:
        """Create a figure legend text.

        Args:
            items: List of dicts with 'label' and 'description'

        Returns:
            Formatted legend text
        """
        legend_parts = []

        for item in items:
            label = item.get("label", "")
            desc = item.get("description", "")
            legend_parts.append(f"**{label}**: {desc}")

        return " ".join(legend_parts)


class MultiPanelFigure:
    """Create multi-panel figures (A, B, C, etc.)."""

    def __init__(self, preset: str = "nature"):
        self.exporter = FigureExporter(preset)
        self.panels: List[Tuple[str, go.Figure]] = []

    def add_panel(self, label: str, fig: go.Figure):
        """Add a panel to the figure."""
        self.panels.append((label, fig))

    def create_layout(
        self,
        layout: str = "2x2",  # '1x2', '2x1', '2x2', '3x1', etc.
    ) -> go.Figure:
        """Create combined multi-panel figure."""
        from plotly.subplots import make_subplots

        # Parse layout
        rows, cols = map(int, layout.split("x"))

        # Create subplots
        fig = make_subplots(rows=rows, cols=cols, subplot_titles=[label for label, _ in self.panels[: rows * cols]])

        # Add each panel
        for i, (label, panel_fig) in enumerate(self.panels[: rows * cols]):
            row = i // cols + 1
            col = i % cols + 1

            for trace in panel_fig.data:
                fig.add_trace(trace, row=row, col=col)

        return fig

    def export_all(self, format: str = "png", scale: float = 2.0) -> List[Tuple[str, bytes]]:
        """Export all panels as separate files."""
        results = []

        for label, fig in self.panels:
            formatted = self.exporter.format_figure(fig, add_panel_label=label)
            data, _ = self.exporter.get_download_link(formatted, f"panel_{label}", format, scale)
            results.append((label, data))

        return results


def render_export_ui(fig: go.Figure, default_filename: str = "figure"):
    """Render Streamlit export UI for a figure."""
    import streamlit as st

    with st.expander("📥 Export Figure"):
        col1, col2, col3 = st.columns(3)

        with col1:
            preset = st.selectbox(
                "Journal/Format", list(JOURNAL_PRESETS.keys()), format_func=lambda x: JOURNAL_PRESETS[x].name
            )

        with col2:
            format = st.selectbox("File Format", ["png", "svg", "pdf"])

        with col3:
            width = st.selectbox(
                "Width", ["single", "full"], format_func=lambda x: "Single Column" if x == "single" else "Full Page"
            )

        # Scale for PNG
        if format == "png":
            scale = st.slider("Resolution Scale", 1.0, 4.0, 2.0)
        else:
            scale = 1.0

        filename = st.text_input("Filename", default_filename)

        if st.button("📥 Download"):
            exporter = FigureExporter(preset)
            formatted_fig = exporter.format_figure(fig, width=width)

            data, mime = exporter.get_download_link(formatted_fig, filename, format, scale)

            st.download_button(f"Download {format.upper()}", data, f"{filename}.{format}", mime)
