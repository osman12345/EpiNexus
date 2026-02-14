"""
Theme Manager - Dark/Light Mode Support

Provides:
- Theme switching (dark/light/auto)
- Custom color palettes
- Consistent styling across app
- Semantic color maps for genomics visualisations
"""

import streamlit as st
from typing import Dict, Any


# ============================================================================
# Semantic color constants for genomics plots
# ============================================================================
# Import these directly: ``from frontend.components.theme import COLORS``


class COLORS:
    """Centralized color constants for consistent plot styling.

    Usage::

        from frontend.components.theme import COLORS
        fig = px.scatter(..., color_discrete_map=COLORS.DIRECTION_MAP)
    """

    # Directional change (differential analysis, volcano, MA)
    UP = "#e74c3c"
    DOWN = "#3498db"
    NOT_SIG = "#95a5a6"
    DIRECTION_MAP = {"Up": UP, "Down": DOWN, "Not Significant": NOT_SIG}

    # Accessibility (ATAC-seq)
    MORE_ACCESSIBLE = "#27ae60"
    LESS_ACCESSIBLE = "#e74c3c"
    ACCESS_MAP = {
        "More Accessible": MORE_ACCESSIBLE,
        "Less Accessible": LESS_ACCESSIBLE,
        "Not Significant": NOT_SIG,
    }

    # Expression integration
    CONCORDANT_ACT = "#4DAF4A"
    CONCORDANT_REP = "#E41A1C"
    DISCORDANT = "#984EA3"
    NO_CHANGE = "#999999"
    EXPRESSION_MAP = {
        "Concordant Activation": CONCORDANT_ACT,
        "Concordant Repression": CONCORDANT_REP,
        "Discordant": DISCORDANT,
        "No Change": NO_CHANGE,
    }

    # Methylation
    HYPER = "#e74c3c"
    HYPO = "#3498db"
    METH_MAP = {"hyper": HYPER, "hypo": HYPO}

    # Super-enhancers
    SUPER_ENHANCER = "#e74c3c"
    TYPICAL_ENHANCER = "#3498db"
    SE_MAP = {
        "Super-Enhancer": SUPER_ENHANCER,
        "Typical Enhancer": TYPICAL_ENHANCER,
    }

    # Significance
    SIG_HIGH = "#E41A1C"
    SIG_MEDIUM = "#FF7F00"
    SIG_LOW = "#999999"
    SIG_MAP = {"High": SIG_HIGH, "Medium": SIG_MEDIUM, "Low": SIG_LOW}

    # General qualitative palette (matches light theme default)
    QUALITATIVE = [
        "#3498db",
        "#e74c3c",
        "#27ae60",
        "#f39c12",
        "#9b59b6",
        "#2ecc71",
        "#e67e22",
        "#1abc9c",
        "#34495e",
        "#7f8c8d",
    ]

    # Continuous scales (strings accepted by Plotly)
    SCALE_DIVERGING = "RdBu_r"
    SCALE_SEQUENTIAL = "Viridis"
    SCALE_HEAT = "Reds"


# Theme definitions
THEMES = {
    "light": {
        "name": "Light",
        "primaryColor": "#1f77b4",
        "backgroundColor": "#ffffff",
        "secondaryBackgroundColor": "#f8f9fa",
        "textColor": "#262730",
        "font": "sans-serif",
        # Plot colors
        "plot_bg": "#ffffff",
        "plot_paper_bg": "#ffffff",
        "plot_grid": "#e5e5e5",
        "plot_text": "#262730",
        # Chart palette
        "palette": [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
        ],
        # Semantic colors
        "success": "#28a745",
        "warning": "#ffc107",
        "error": "#dc3545",
        "info": "#17a2b8",
    },
    "dark": {
        "name": "Dark",
        "primaryColor": "#4dabf7",
        "backgroundColor": "#0e1117",
        "secondaryBackgroundColor": "#1a1d24",
        "textColor": "#fafafa",
        "font": "sans-serif",
        # Plot colors
        "plot_bg": "#0e1117",
        "plot_paper_bg": "#0e1117",
        "plot_grid": "#2d3139",
        "plot_text": "#fafafa",
        # Chart palette (brighter for dark mode)
        "palette": [
            "#4dabf7",
            "#ffa94d",
            "#69db7c",
            "#ff6b6b",
            "#cc5de8",
            "#f783ac",
            "#a9e34b",
            "#adb5bd",
            "#ffe066",
            "#66d9e8",
        ],
        # Semantic colors
        "success": "#51cf66",
        "warning": "#fcc419",
        "error": "#ff6b6b",
        "info": "#4dabf7",
    },
    "midnight": {
        "name": "Midnight Blue",
        "primaryColor": "#7c3aed",
        "backgroundColor": "#0f0f23",
        "secondaryBackgroundColor": "#1a1a2e",
        "textColor": "#e0e0e0",
        "font": "sans-serif",
        "plot_bg": "#0f0f23",
        "plot_paper_bg": "#0f0f23",
        "plot_grid": "#2a2a4a",
        "plot_text": "#e0e0e0",
        "palette": [
            "#7c3aed",
            "#06b6d4",
            "#10b981",
            "#f59e0b",
            "#ef4444",
            "#ec4899",
            "#8b5cf6",
            "#6366f1",
            "#14b8a6",
            "#f97316",
        ],
        "success": "#10b981",
        "warning": "#f59e0b",
        "error": "#ef4444",
        "info": "#06b6d4",
    },
}


def init_theme():
    """Initialize theme in session state."""
    if "theme" not in st.session_state:
        st.session_state.theme = "light"


def get_current_theme() -> Dict[str, Any]:
    """Get current theme configuration."""
    init_theme()
    return THEMES.get(st.session_state.theme, THEMES["light"])


def set_theme(theme_name: str):
    """Set the current theme."""
    if theme_name in THEMES:
        st.session_state.theme = theme_name


def get_theme_css() -> str:
    """Generate CSS for current theme."""
    theme = get_current_theme()

    return f"""
    <style>
        /* Main app styling */
        .stApp {{
            background-color: {theme["backgroundColor"]};
            color: {theme["textColor"]};
        }}

        /* Sidebar */
        [data-testid="stSidebar"] {{
            background-color: {theme["secondaryBackgroundColor"]};
        }}

        /* Headers */
        h1, h2, h3, h4, h5, h6 {{
            color: {theme["textColor"]} !important;
        }}

        /* Cards and containers */
        .stMetric {{
            background-color: {theme["secondaryBackgroundColor"]};
            border-radius: 8px;
            padding: 10px;
        }}

        /* Dataframes */
        .stDataFrame {{
            background-color: {theme["secondaryBackgroundColor"]};
        }}

        /* Expanders */
        .streamlit-expanderHeader {{
            background-color: {theme["secondaryBackgroundColor"]};
            color: {theme["textColor"]};
        }}

        /* Tabs */
        .stTabs [data-baseweb="tab"] {{
            color: {theme["textColor"]};
        }}

        .stTabs [data-baseweb="tab"][aria-selected="true"] {{
            color: {theme["primaryColor"]};
            border-bottom-color: {theme["primaryColor"]};
        }}

        /* Buttons */
        .stButton > button {{
            background-color: {theme["primaryColor"]};
            color: white;
            border: none;
        }}

        .stButton > button:hover {{
            background-color: {theme["primaryColor"]}dd;
        }}

        /* Success/Warning/Error messages */
        .stSuccess {{
            background-color: {theme["success"]}22;
            border-left-color: {theme["success"]};
        }}

        .stWarning {{
            background-color: {theme["warning"]}22;
            border-left-color: {theme["warning"]};
        }}

        .stError {{
            background-color: {theme["error"]}22;
            border-left-color: {theme["error"]};
        }}

        /* Code blocks */
        code {{
            background-color: {theme["secondaryBackgroundColor"]};
            color: {theme["primaryColor"]};
        }}

        /* Links */
        a {{
            color: {theme["primaryColor"]};
        }}
    </style>
    """


def apply_plotly_theme(fig):
    """Apply current theme to a Plotly figure."""
    theme = get_current_theme()

    fig.update_layout(
        paper_bgcolor=theme["plot_paper_bg"],
        plot_bgcolor=theme["plot_bg"],
        font_color=theme["plot_text"],
        title_font_color=theme["plot_text"],
        legend_font_color=theme["plot_text"],
        xaxis=dict(gridcolor=theme["plot_grid"], linecolor=theme["plot_grid"], tickfont=dict(color=theme["plot_text"])),
        yaxis=dict(gridcolor=theme["plot_grid"], linecolor=theme["plot_grid"], tickfont=dict(color=theme["plot_text"])),
    )

    return fig


def get_plot_colors(n: int = 10) -> list:
    """Get color palette for plots."""
    theme = get_current_theme()
    palette = theme["palette"]

    if n <= len(palette):
        return palette[:n]
    else:
        # Repeat palette if needed
        return (palette * (n // len(palette) + 1))[:n]


def render_theme_selector():
    """Render theme selector widget."""
    init_theme()

    col1, col2 = st.columns([3, 1])

    with col2:
        theme_options = list(THEMES.keys())
        current_idx = theme_options.index(st.session_state.theme)

        new_theme = st.selectbox(
            "🎨 Theme", theme_options, index=current_idx, format_func=lambda x: THEMES[x]["name"], key="theme_selector"
        )

        if new_theme != st.session_state.theme:
            set_theme(new_theme)
            st.rerun()


def render_theme_preview():
    """Render a preview of available themes."""
    st.subheader("Theme Preview")

    cols = st.columns(len(THEMES))

    for col, (theme_name, theme) in zip(cols, THEMES.items()):
        with col:
            st.markdown(
                f"""
            <div style="
                background-color: {theme["backgroundColor"]};
                border: 2px solid {theme["primaryColor"]};
                border-radius: 10px;
                padding: 15px;
                text-align: center;
            ">
                <div style="color: {theme["textColor"]}; font-weight: bold;">
                    {theme["name"]}
                </div>
                <div style="
                    display: flex;
                    justify-content: center;
                    gap: 5px;
                    margin-top: 10px;
                ">
                    {"".join([f'<div style="width:20px;height:20px;background:{c};border-radius:3px;"></div>' for c in theme["palette"][:5]])}
                </div>
            </div>
            """,
                unsafe_allow_html=True,
            )

            if st.button("Apply", key=f"apply_{theme_name}"):
                set_theme(theme_name)
                st.rerun()
