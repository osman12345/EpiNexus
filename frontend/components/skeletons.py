"""
Reusable skeleton / placeholder loading screens for Streamlit pages.

Replace bare ``st.spinner(...)`` with these helpers to give users a
richer visual cue while long operations run.  Each skeleton renders
lightweight HTML/CSS that mimics the shape of the final content.

Usage::

    from frontend.components.skeletons import (
        skeleton_metrics_row,
        skeleton_chart,
        skeleton_table,
        skeleton_card,
    )

    with skeleton_metrics_row():
        results = expensive_computation()  # blocks here
    # after the ``with`` block, the skeleton vanishes and real content shows

Copyright (c) 2026 EpiNexus Contributors
SPDX-License-Identifier: AGPL-3.0-or-later OR Commercial
"""

from contextlib import contextmanager
from typing import Optional

import streamlit as st

# ---------------------------------------------------------------------------
# CSS shared across all skeletons
# ---------------------------------------------------------------------------
_SKELETON_CSS = """
<style>
@keyframes epinexus-shimmer {
    0%   { background-position: -200px 0; }
    100% { background-position: calc(200px + 100%) 0; }
}
.epinexus-skeleton {
    background: linear-gradient(90deg, #f0f0f0 25%, #e0e0e0 50%, #f0f0f0 75%);
    background-size: 200px 100%;
    animation: epinexus-shimmer 1.5s infinite;
    border-radius: 6px;
}
.epinexus-skeleton-dark {
    background: linear-gradient(90deg, #2a2a2a 25%, #333 50%, #2a2a2a 75%);
    background-size: 200px 100%;
    animation: epinexus-shimmer 1.5s infinite;
    border-radius: 6px;
}
</style>
"""

_css_injected = False


def _inject_css() -> None:
    """Inject skeleton CSS once per session."""
    global _css_injected
    if not _css_injected:
        st.markdown(_SKELETON_CSS, unsafe_allow_html=True)
        _css_injected = True


# ---------------------------------------------------------------------------
# Individual skeleton primitives
# ---------------------------------------------------------------------------

def skeleton_metrics_row(n_metrics: int = 4, height: int = 80) -> None:
    """Render a row of metric-card placeholders.

    Args:
        n_metrics: Number of metric cards to show.
        height: Height of each card in pixels.
    """
    _inject_css()
    cols = st.columns(n_metrics)
    for col in cols:
        col.markdown(
            f'<div class="epinexus-skeleton" '
            f'style="height:{height}px;margin-bottom:8px;"></div>',
            unsafe_allow_html=True,
        )


def skeleton_chart(height: int = 350, label: Optional[str] = None) -> None:
    """Render a chart-shaped placeholder.

    Args:
        height: Height of the chart area in pixels.
        label: Optional text shown above the skeleton.
    """
    _inject_css()
    if label:
        st.caption(label)
    st.markdown(
        f'<div class="epinexus-skeleton" '
        f'style="height:{height}px;width:100%;margin-bottom:12px;"></div>',
        unsafe_allow_html=True,
    )


def skeleton_table(n_rows: int = 5, n_cols: int = 4) -> None:
    """Render a table-shaped placeholder with shimmering rows.

    Args:
        n_rows: Number of rows to display.
        n_cols: Number of columns.
    """
    _inject_css()
    # Header row
    header = "".join(
        f'<div class="epinexus-skeleton" '
        f'style="flex:1;height:18px;margin:0 4px;"></div>'
        for _ in range(n_cols)
    )
    rows = ""
    for _ in range(n_rows):
        cells = "".join(
            f'<div class="epinexus-skeleton" '
            f'style="flex:1;height:14px;margin:0 4px;"></div>'
            for _ in range(n_cols)
        )
        rows += f'<div style="display:flex;gap:4px;margin-top:6px;">{cells}</div>'

    st.markdown(
        f'<div style="padding:8px 0;">'
        f'<div style="display:flex;gap:4px;margin-bottom:10px;">{header}</div>'
        f'{rows}</div>',
        unsafe_allow_html=True,
    )


def skeleton_card(height: int = 120, label: Optional[str] = None) -> None:
    """Render a generic card placeholder.

    Args:
        height: Card height in pixels.
        label: Optional caption.
    """
    _inject_css()
    if label:
        st.caption(label)
    st.markdown(
        f'<div class="epinexus-skeleton" '
        f'style="height:{height}px;width:100%;border-radius:12px;'
        f'margin-bottom:12px;"></div>',
        unsafe_allow_html=True,
    )


# ---------------------------------------------------------------------------
# Context-manager wrappers (show skeleton → run code → replace with result)
# ---------------------------------------------------------------------------

@contextmanager
def loading_metrics(n_metrics: int = 4, message: str = "Loading..."):
    """Context manager: show metric skeletons while code runs.

    Usage::

        with loading_metrics(4, "Computing QC..."):
            results = run_qc()
        # skeletons disappear; render real metrics below
    """
    placeholder = st.empty()
    with placeholder.container():
        st.caption(message)
        skeleton_metrics_row(n_metrics)
    try:
        yield
    finally:
        placeholder.empty()


@contextmanager
def loading_chart(height: int = 350, message: str = "Generating plot..."):
    """Context manager: show chart skeleton while code runs."""
    placeholder = st.empty()
    with placeholder.container():
        skeleton_chart(height, label=message)
    try:
        yield
    finally:
        placeholder.empty()


@contextmanager
def loading_table(n_rows: int = 5, message: str = "Fetching data..."):
    """Context manager: show table skeleton while code runs."""
    placeholder = st.empty()
    with placeholder.container():
        st.caption(message)
        skeleton_table(n_rows)
    try:
        yield
    finally:
        placeholder.empty()
