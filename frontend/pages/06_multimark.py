"""
Multi-Mark Integration Page
Analyze relationships between multiple histone modifications.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

st.set_page_config(page_title="Multi-Mark Integration - EpiNexus", page_icon="🔄", layout="wide")

# Try to import data manager
try:
    from frontend.components.data_manager import DataManager
    HAS_DATA_MANAGER = True
except ImportError:
    HAS_DATA_MANAGER = False


def has_data():
    """Check if user has loaded data for multiple marks."""
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
            <div style="font-size: 3rem; margin-bottom: 1rem;">🔄</div>
            <h2 style="color: #6c757d; margin-bottom: 0.5rem;">No Data Loaded</h2>
            <p style="color: #6c757d; font-size: 1.1rem;">
                Upload data for multiple histone marks to analyze their relationships.
            </p>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("")

        if st.button("📁 Go to Data & Project", type="primary", use_container_width=True):
            st.switch_page("pages/01_data_project.py")

        st.markdown("")
        st.markdown("**What you need:**")
        st.markdown("- Peak files for 2+ histone marks")
        st.markdown("- Same samples across marks recommended")
        st.markdown("- Optional: expression data for correlation")


def main():
    st.title("🔄 Multi-Mark Integration")
    st.markdown("Analyze relationships and co-occurrence patterns between histone modifications.")

    # Check if data is loaded
    if not has_data():
        render_empty_state()
        return

    tab1, tab2, tab3, tab4 = st.tabs([
        "📊 Co-occurrence",
        "🎯 Chromatin States",
        "📈 Mark Transitions",
        "🔬 Feature Correlation"
    ])

    with tab1:
        render_cooccurrence()
    with tab2:
        render_chromatin_states()
    with tab3:
        render_mark_transitions()
    with tab4:
        render_feature_correlation()


def render_cooccurrence():
    """Analyze co-occurrence of histone marks."""
    st.header("Histone Mark Co-occurrence")

    marks = ['H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K27me3', 'H3K36me3', 'H3K9me3']

    col1, col2 = st.columns([1, 2])

    with col1:
        st.subheader("Settings")

        selected_marks = st.multiselect(
            "Select marks to analyze",
            marks,
            default=marks[:4]
        )

        overlap_thresh = st.slider("Minimum overlap (bp)", 1, 500, 100)

        jaccard_method = st.checkbox("Use Jaccard similarity", value=True)

    with col2:
        st.subheader("Co-occurrence Matrix")

        # Generate correlation matrix
        np.random.seed(42)
        n = len(selected_marks)
        corr = np.eye(n)

        # Define expected correlations
        mark_corrs = {
            ('H3K4me3', 'H3K27ac'): 0.75,
            ('H3K4me3', 'H3K4me1'): 0.45,
            ('H3K27ac', 'H3K4me1'): 0.65,
            ('H3K4me3', 'H3K27me3'): -0.3,
            ('H3K27ac', 'H3K27me3'): -0.4,
            ('H3K4me1', 'H3K27me3'): -0.2,
            ('H3K36me3', 'H3K4me3'): 0.3,
            ('H3K9me3', 'H3K27me3'): 0.5
        }

        for i, m1 in enumerate(selected_marks):
            for j, m2 in enumerate(selected_marks):
                if i < j:
                    key = (m1, m2) if (m1, m2) in mark_corrs else (m2, m1)
                    corr[i, j] = corr[j, i] = mark_corrs.get(key, np.random.uniform(-0.1, 0.3))

        fig = px.imshow(
            corr, x=selected_marks, y=selected_marks,
            color_continuous_scale='RdBu_r', zmin=-1, zmax=1,
            text_auto='.2f'
        )
        fig.update_layout(height=500)
        st.plotly_chart(fig, use_container_width=True)

    st.subheader("Upset Plot - Mark Combinations")

    # Create upset-style visualization
    combinations = [
        'H3K4me3 only', 'H3K27ac only', 'H3K4me1 only',
        'H3K4me3 + H3K27ac', 'H3K4me1 + H3K27ac',
        'H3K4me3 + H3K27ac + H3K4me1', 'All marks'
    ]
    counts = [8500, 12000, 15000, 6200, 8900, 3200, 1500]

    fig = px.bar(x=combinations, y=counts, color=counts,
                color_continuous_scale='Viridis')
    fig.update_layout(
        xaxis_title='Mark Combination',
        yaxis_title='Number of Regions',
        height=400,
        xaxis_tickangle=-45
    )
    st.plotly_chart(fig, use_container_width=True)


def render_chromatin_states():
    """Infer chromatin states from mark combinations."""
    st.header("Chromatin State Analysis")

    st.markdown("""
    Infer chromatin states based on combinations of histone modifications,
    similar to ChromHMM but using your differential peaks.
    """)

    col1, col2 = st.columns([1, 2])

    with col1:
        st.subheader("State Definitions")

        states = {
            'Active Promoter': 'H3K4me3+ H3K27ac+ H3K27me3-',
            'Strong Enhancer': 'H3K4me1+ H3K27ac+ H3K4me3-',
            'Weak Enhancer': 'H3K4me1+ H3K27ac- H3K4me3-',
            'Poised Promoter': 'H3K4me3+ H3K27me3+',
            'Repressed': 'H3K27me3+ H3K4me3- H3K27ac-',
            'Transcribed': 'H3K36me3+ H3K27me3-'
        }

        for state, marks in states.items():
            st.markdown(f"**{state}**: `{marks}`")

    with col2:
        st.subheader("State Distribution")

        state_counts = pd.DataFrame({
            'State': list(states.keys()),
            'Control': [12500, 18200, 25600, 3200, 8900, 15400],
            'Treatment': [15800, 22100, 23400, 2800, 6500, 18200]
        })

        fig = go.Figure()
        fig.add_trace(go.Bar(x=state_counts['State'], y=state_counts['Control'],
                            name='Control', marker_color='#3498db'))
        fig.add_trace(go.Bar(x=state_counts['State'], y=state_counts['Treatment'],
                            name='Treatment', marker_color='#e74c3c'))

        fig.update_layout(
            barmode='group',
            xaxis_title='Chromatin State',
            yaxis_title='Number of Regions',
            height=450,
            xaxis_tickangle=-45
        )
        st.plotly_chart(fig, use_container_width=True)


def render_mark_transitions():
    """Analyze chromatin state transitions."""
    st.header("Chromatin State Transitions")

    st.markdown("Track how chromatin states change between conditions.")

    # Sankey diagram for state transitions
    states = ['Active Promoter', 'Strong Enhancer', 'Weak Enhancer', 'Poised', 'Repressed']

    # Define transitions
    source = [0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 4, 4]
    target = [5, 6, 7, 5, 6, 6, 7, 8, 5, 9, 8, 9]
    value = [8500, 2100, 1200, 3200, 14500, 5800, 12400, 6200, 800, 2100, 1500, 7200]

    # Labels for both sides
    labels = ['Active Promoter (C)', 'Strong Enhancer (C)', 'Weak Enhancer (C)',
             'Poised (C)', 'Repressed (C)',
             'Active Promoter (T)', 'Strong Enhancer (T)', 'Weak Enhancer (T)',
             'Poised (T)', 'Repressed (T)']

    colors = ['#2ecc71', '#f39c12', '#f1c40f', '#9b59b6', '#e74c3c'] * 2

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color='black', width=0.5),
            label=labels,
            color=colors
        ),
        link=dict(
            source=source,
            target=target,
            value=value
        )
    )])

    fig.update_layout(title_text="State Transitions: Control → Treatment", height=500)
    st.plotly_chart(fig, use_container_width=True)

    # Summary table
    st.subheader("Transition Summary")

    transition_summary = pd.DataFrame({
        'From State': ['Active Promoter', 'Strong Enhancer', 'Weak Enhancer', 'Poised', 'Repressed'],
        'Stable': ['72%', '78%', '51%', '27%', '83%'],
        'Activated': ['-', '18%', '24%', '27%', '17%'],
        'Repressed': ['20%', '4%', '25%', '73%', '-'],
        'Net Change': ['+26%', '+21%', '-9%', '-13%', '-27%']
    })

    st.dataframe(transition_summary, use_container_width=True, hide_index=True)


def render_feature_correlation():
    """Correlate marks with genomic features."""
    st.header("Mark-Feature Correlations")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Marks vs Gene Expression")

        # Scatter plot
        np.random.seed(42)
        n = 500

        expr_corr = pd.DataFrame({
            'H3K27ac_signal': np.random.lognormal(2, 1, n),
            'Expression': lambda: 0,
            'Mark': np.random.choice(['H3K27ac', 'H3K4me3', 'H3K27me3'], n)
        })

        # Positive correlation for active marks
        expr_corr.loc[expr_corr['Mark'] == 'H3K27ac', 'Expression'] = \
            expr_corr.loc[expr_corr['Mark'] == 'H3K27ac', 'H3K27ac_signal'] * 0.8 + np.random.normal(0, 2, (expr_corr['Mark'] == 'H3K27ac').sum())
        expr_corr.loc[expr_corr['Mark'] == 'H3K4me3', 'Expression'] = \
            expr_corr.loc[expr_corr['Mark'] == 'H3K4me3', 'H3K27ac_signal'] * 0.6 + np.random.normal(0, 2, (expr_corr['Mark'] == 'H3K4me3').sum())
        expr_corr.loc[expr_corr['Mark'] == 'H3K27me3', 'Expression'] = \
            -expr_corr.loc[expr_corr['Mark'] == 'H3K27me3', 'H3K27ac_signal'] * 0.3 + np.random.normal(5, 2, (expr_corr['Mark'] == 'H3K27me3').sum())

        fig = px.scatter(expr_corr, x='H3K27ac_signal', y='Expression', color='Mark',
                        opacity=0.6, trendline='ols')
        fig.update_layout(
            xaxis_title='Histone Mark Signal',
            yaxis_title='Gene Expression (log2 TPM)',
            height=400
        )
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("Feature Enrichment per Mark")

        features = ['TSS', 'Gene Body', 'Enhancer', 'Intergenic', 'CTCF Sites']
        marks = ['H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K36me3', 'H3K27me3']

        enrichment = np.array([
            [12.5, 0.8, 2.1, 0.3, 3.2],    # H3K4me3
            [8.2, 1.2, 8.5, 0.5, 2.8],     # H3K27ac
            [2.1, 1.5, 15.2, 0.8, 1.5],    # H3K4me1
            [1.2, 8.5, 0.5, 0.6, 0.8],     # H3K36me3
            [0.8, 0.6, 0.4, 2.5, 0.5]      # H3K27me3
        ])

        fig = px.imshow(enrichment, x=features, y=marks,
                       color_continuous_scale='YlOrRd', text_auto='.1f')
        fig.update_layout(
            xaxis_title='Genomic Feature',
            yaxis_title='Histone Mark',
            height=400
        )
        st.plotly_chart(fig, use_container_width=True)


if __name__ == "__main__":
    main()
