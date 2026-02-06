"""
Empty State Components for EpiNexus

Provides consistent empty state UI when no data is loaded.
"""

import streamlit as st


def render_empty_state(
    title: str = "No Data Loaded",
    icon: str = "📊",
    message: str = "Upload your data to begin analysis.",
    requirements: list = None,
    help_topic: str = None
):
    """
    Render a consistent empty state when no data is loaded.

    Args:
        title: Main title to display
        icon: Emoji icon
        message: Description message
        requirements: List of requirements for this analysis
        help_topic: Help page topic to link to
    """
    st.markdown("---")

    col1, col2, col3 = st.columns([1, 2, 1])

    with col2:
        st.markdown(f"""
        <div style="text-align: center; padding: 3rem; background: #f8f9fa;
                    border-radius: 12px; border: 2px dashed #dee2e6;">
            <div style="font-size: 3rem; margin-bottom: 1rem;">{icon}</div>
            <h2 style="color: #6c757d; margin-bottom: 0.5rem;">{title}</h2>
            <p style="color: #6c757d; font-size: 1.1rem;">{message}</p>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("")

        if st.button("📁 Go to Data & Project", type="primary", use_container_width=True):
            st.switch_page("pages/01_data_project.py")

        if requirements:
            st.markdown("")
            st.markdown("**What you need:**")
            for req in requirements:
                st.markdown(f"- {req}")

        if help_topic:
            st.markdown("")
            st.markdown(f"**Need help?** Check the [Documentation]({help_topic}) for guidance.")


def render_no_results_state(
    title: str = "No Results Yet",
    icon: str = "📋",
    message: str = "Run the analysis to see results.",
    action_text: str = None
):
    """
    Render a state when analysis hasn't been run yet.

    Args:
        title: Main title
        icon: Emoji icon
        message: Description
        action_text: Additional action instructions
    """
    st.info(f"{icon} **{title}**: {message}")

    if action_text:
        st.caption(action_text)


def check_data_loaded() -> bool:
    """Check if user has loaded any data."""
    # Check data manager
    try:
        from frontend.components.data_manager import DataManager
        peaks = DataManager.get_data('peaks')
        if peaks is not None and len(peaks) > 0:
            return True
    except ImportError:
        pass

    # Check session state
    if len(st.session_state.get('samples', [])) > 0:
        return True

    return False


def check_peaks_loaded() -> bool:
    """Check if peak data is loaded."""
    try:
        from frontend.components.data_manager import DataManager
        peaks = DataManager.get_data('peaks')
        return peaks is not None and len(peaks) > 0
    except ImportError:
        return False


def check_expression_loaded() -> bool:
    """Check if expression data is loaded."""
    try:
        from frontend.components.data_manager import DataManager
        expr = DataManager.get_data('expression')
        return expr is not None and len(expr) > 0
    except ImportError:
        return False
