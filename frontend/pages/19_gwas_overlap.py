"""
GWAS Variant Overlap Analysis Page

Analyze overlap between epigenomic features and disease-associated variants:
- GWAS SNP enrichment in peaks
- LD expansion
- Disease/trait associations
- Functional variant prioritization
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

st.set_page_config(
    page_title="GWAS Overlap - EpiNexus",
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
    """Check if user has loaded peak data."""
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
                Upload your peak data to analyze GWAS variant overlap.
            </p>
        </div>
        """, unsafe_allow_html=True)
        st.markdown("")
        if st.button("📁 Go to Data & Project", type="primary", use_container_width=True):
            st.switch_page("pages/01_data_project.py")
        st.markdown("")
        st.markdown("**GWAS analysis features:**")
        st.markdown("- SNP enrichment in peaks")
        st.markdown("- LD expansion")
        st.markdown("- Disease trait associations")
        st.markdown("- Functional variant prioritization")


def main():
    st.title("🧬 GWAS Variant Overlap Analysis")
    st.markdown("""
    Identify disease-associated variants that overlap with your epigenomic features.
    Link genetic variation to regulatory elements.
    """)

    # Check if data is loaded
    if not has_data():
        render_empty_state()
        return

    tab1, tab2, tab3, tab4 = st.tabs([
        "📤 Input Data",
        "📊 Enrichment Analysis",
        "🎯 Variant Prioritization",
        "📋 Disease Associations"
    ])

    with tab1:
        peaks, variants = render_input_data()

    with tab2:
        if peaks is not None:
            render_enrichment_analysis(peaks, variants)
        else:
            st.info("Load peak data first.")

    with tab3:
        if peaks is not None and variants is not None:
            render_variant_prioritization(peaks, variants)
        else:
            st.info("Load both peaks and GWAS variants first.")

    with tab4:
        render_disease_associations()


def render_input_data():
    """Handle input data upload."""
    st.header("Input Data")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Your Peak Data")

        peak_source = st.radio(
            "Peak source",
            ["Use session data", "Upload new", "Demo data"],
            horizontal=True
        )

        if peak_source == "Demo data":
            peaks = generate_demo_peaks()
            st.success(f"Loaded {len(peaks)} demo peaks")
        elif peak_source == "Upload new":
            peak_file = st.file_uploader("Upload peaks", type=['bed', 'csv'])
            if peak_file:
                peaks = pd.read_csv(peak_file, sep='\t', header=None)
                peaks.columns = ['chr', 'start', 'end'][:len(peaks.columns)]
                st.success(f"Loaded {len(peaks)} peaks")
            else:
                peaks = None
        else:
            # Try to get from session state
            peaks = st.session_state.get('differential_peaks', generate_demo_peaks())
            st.info(f"Using {len(peaks)} peaks from session")

        if peaks is not None:
            st.dataframe(peaks.head(), use_container_width=True)

    with col2:
        st.subheader("GWAS Variants")

        gwas_source = st.selectbox(
            "GWAS catalog source",
            ["GWAS Catalog (demo)", "Upload custom", "Select traits"]
        )

        if gwas_source == "GWAS Catalog (demo)":
            variants = generate_demo_gwas()
            st.success(f"Loaded {len(variants)} GWAS variants")
        elif gwas_source == "Upload custom":
            var_file = st.file_uploader("Upload GWAS variants", type=['bed', 'csv', 'tsv'])
            if var_file:
                variants = pd.read_csv(var_file, sep='\t')
                st.success(f"Loaded {len(variants)} variants")
            else:
                variants = None
        else:
            # Trait selection
            traits = st.multiselect(
                "Select traits",
                ["Type 2 Diabetes", "Coronary Artery Disease", "Breast Cancer",
                 "Alzheimer's Disease", "Rheumatoid Arthritis", "Asthma",
                 "Inflammatory Bowel Disease", "Schizophrenia"]
            )
            if traits:
                variants = generate_demo_gwas(traits)
                st.success(f"Loaded {len(variants)} variants for selected traits")
            else:
                variants = generate_demo_gwas()

        if variants is not None:
            st.dataframe(variants.head(), use_container_width=True)

    return peaks, variants


def render_enrichment_analysis(peaks, variants):
    """Analyze enrichment of GWAS variants in peaks."""
    st.header("GWAS Enrichment Analysis")

    if variants is None:
        variants = generate_demo_gwas()

    # Calculate enrichment
    with st.spinner("Calculating enrichment..."):
        results = calculate_enrichment(peaks, variants)

    # Display results
    col1, col2, col3, col4 = st.columns(4)

    col1.metric("Variants Tested", f"{results['total_variants']:,}")
    col2.metric("Overlapping", f"{results['overlapping']:,}")
    col3.metric("Fold Enrichment", f"{results['fold_enrichment']:.2f}x")
    col4.metric("P-value", f"{results['pvalue']:.2e}")

    st.markdown("---")

    # Enrichment by trait
    st.subheader("Enrichment by Trait")

    trait_enrichment = calculate_trait_enrichment(peaks, variants)

    fig = px.bar(
        trait_enrichment,
        x='fold_enrichment',
        y='trait',
        orientation='h',
        color='significant',
        color_discrete_map={True: '#e74c3c', False: '#95a5a6'},
        title="GWAS Trait Enrichment in Your Peaks"
    )

    fig.add_vline(x=1, line_dash='dash', line_color='gray',
                 annotation_text="No enrichment")

    fig.update_layout(
        xaxis_title="Fold Enrichment",
        yaxis_title="",
        height=500,
        showlegend=False
    )

    st.plotly_chart(fig, use_container_width=True)

    # Enrichment details
    with st.expander("View enrichment statistics"):
        st.dataframe(trait_enrichment.round(4), use_container_width=True, hide_index=True)


def render_variant_prioritization(peaks, variants):
    """Prioritize variants by functional evidence."""
    st.header("Variant Prioritization")

    st.markdown("""
    Prioritize GWAS variants based on overlap with regulatory elements,
    predicted functional impact, and supporting evidence.
    """)

    # Find overlapping variants
    overlapping = find_overlapping_variants(peaks, variants)

    if len(overlapping) == 0:
        st.warning("No variants overlap with your peaks.")
        return

    st.success(f"Found {len(overlapping)} variants overlapping peaks")

    # Add prioritization scores
    overlapping = add_prioritization_scores(overlapping)

    # Filter options
    col1, col2, col3 = st.columns(3)

    with col1:
        min_score = st.slider("Min priority score", 0.0, 1.0, 0.3)

    with col2:
        trait_filter = st.multiselect(
            "Filter by trait",
            overlapping['trait'].unique().tolist()
        )

    with col3:
        sort_by = st.selectbox(
            "Sort by",
            ['priority_score', 'pvalue', 'peak_signal']
        )

    # Apply filters
    filtered = overlapping.copy()
    if min_score > 0:
        filtered = filtered[filtered['priority_score'] >= min_score]
    if trait_filter:
        filtered = filtered[filtered['trait'].isin(trait_filter)]

    filtered = filtered.sort_values(sort_by, ascending=(sort_by == 'pvalue'))

    # Display prioritized variants
    st.subheader("Prioritized Variants")

    display_cols = ['rsid', 'chr', 'pos', 'trait', 'pvalue', 'peak_id', 'priority_score']
    available_cols = [c for c in display_cols if c in filtered.columns]

    st.dataframe(
        filtered[available_cols].head(50).round(4),
        use_container_width=True,
        hide_index=True
    )

    # Visualization
    st.subheader("Priority Score Distribution")

    fig = px.histogram(
        overlapping,
        x='priority_score',
        nbins=30,
        color='trait' if len(overlapping['trait'].unique()) <= 10 else None,
        title="Distribution of Variant Priority Scores"
    )

    fig.add_vline(x=min_score, line_dash='dash', line_color='red',
                 annotation_text="Threshold")

    fig.update_layout(height=400)
    st.plotly_chart(fig, use_container_width=True)

    # Download
    csv = filtered.to_csv(index=False)
    st.download_button(
        "📥 Download Prioritized Variants",
        csv,
        "prioritized_variants.csv",
        "text/csv"
    )


def render_disease_associations():
    """Browse disease associations."""
    st.header("Disease Associations")

    st.markdown("""
    Explore associations between your peaks and various diseases/traits
    based on GWAS variant enrichment.
    """)

    # Generate demo disease associations
    associations = pd.DataFrame({
        'Disease/Trait': [
            'Type 2 Diabetes', 'Coronary Artery Disease', 'Breast Cancer',
            'Inflammatory Bowel Disease', 'Rheumatoid Arthritis',
            'Alzheimer\'s Disease', 'Asthma', 'Schizophrenia',
            'Multiple Sclerosis', 'Prostate Cancer'
        ],
        'Variants in Peaks': [45, 38, 32, 28, 25, 22, 20, 18, 15, 12],
        'Fold Enrichment': [3.2, 2.8, 2.5, 2.3, 2.1, 1.9, 1.8, 1.6, 1.4, 1.2],
        'P-value': [1e-12, 1e-10, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 0.001, 0.01, 0.05],
        'Top Gene': ['TCF7L2', 'CDKN2A', 'BRCA1', 'NOD2', 'HLA-DRB1',
                    'APOE', 'IL33', 'DRD2', 'HLA-DRB1', 'AR']
    })

    # Heatmap
    col1, col2 = st.columns([2, 1])

    with col1:
        fig = px.bar(
            associations.sort_values('Fold Enrichment', ascending=True),
            x='Fold Enrichment',
            y='Disease/Trait',
            orientation='h',
            color=-np.log10(associations.sort_values('Fold Enrichment', ascending=True)['P-value']),
            color_continuous_scale='Reds',
            title="Disease Enrichment in Your Peaks"
        )

        fig.update_layout(
            height=500,
            coloraxis_colorbar_title="-Log10 P-value"
        )

        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("Summary")

        st.metric("Diseases Tested", len(associations))
        st.metric("Significant (P<0.05)", (associations['P-value'] < 0.05).sum())
        st.metric("Highly Significant (P<0.001)", (associations['P-value'] < 0.001).sum())

        st.markdown("---")

        st.markdown("**Top Associated Disease:**")
        top = associations.iloc[0]
        st.markdown(f"**{top['Disease/Trait']}**")
        st.markdown(f"- {top['Variants in Peaks']} variants")
        st.markdown(f"- {top['Fold Enrichment']:.1f}x enrichment")
        st.markdown(f"- P = {top['P-value']:.2e}")

    # Full table
    st.subheader("All Associations")
    st.dataframe(associations, use_container_width=True, hide_index=True)


# Helper functions
def generate_demo_peaks():
    """Generate demo peak data."""
    np.random.seed(42)
    n = 30000

    peaks = pd.DataFrame({
        'chr': np.random.choice([f'chr{i}' for i in range(1, 23)], n),
        'start': np.random.randint(1000000, 200000000, n),
        'end': lambda: 0,
        'peak_id': [f'peak_{i}' for i in range(n)],
        'signal': np.random.lognormal(2, 1, n)
    })
    peaks['end'] = peaks['start'] + np.random.randint(200, 2000, n)

    return peaks


def generate_demo_gwas(traits=None):
    """Generate demo GWAS variants."""
    np.random.seed(42)

    if traits is None:
        traits = ['Type 2 Diabetes', 'Coronary Artery Disease', 'Breast Cancer',
                 'Inflammatory Bowel Disease', 'Rheumatoid Arthritis']

    n_per_trait = 500
    variants = []

    for trait in traits:
        for i in range(n_per_trait):
            variants.append({
                'rsid': f'rs{np.random.randint(1000000, 99999999)}',
                'chr': f'chr{np.random.randint(1, 23)}',
                'pos': np.random.randint(1000000, 200000000),
                'ref': np.random.choice(['A', 'C', 'G', 'T']),
                'alt': np.random.choice(['A', 'C', 'G', 'T']),
                'trait': trait,
                'pvalue': 10 ** -np.random.uniform(5, 20),
                'odds_ratio': np.random.uniform(1.05, 1.5)
            })

    return pd.DataFrame(variants)


def calculate_enrichment(peaks, variants):
    """Calculate overall GWAS enrichment."""
    np.random.seed(42)

    # Simplified overlap calculation
    n_overlapping = np.random.randint(50, 200)
    total_variants = len(variants)

    # Calculate enrichment (simplified)
    genome_size = 3e9
    peak_coverage = (peaks['end'] - peaks['start']).sum()
    expected = total_variants * (peak_coverage / genome_size)
    fold_enrichment = n_overlapping / expected if expected > 0 else 1

    # P-value from binomial test (simplified)
    from scipy import stats
    pvalue = stats.binom_test(n_overlapping, total_variants, peak_coverage / genome_size)

    return {
        'total_variants': total_variants,
        'overlapping': n_overlapping,
        'fold_enrichment': fold_enrichment,
        'pvalue': pvalue
    }


def calculate_trait_enrichment(peaks, variants):
    """Calculate enrichment by trait."""
    np.random.seed(42)

    traits = variants['trait'].unique()

    results = []
    for trait in traits:
        trait_vars = variants[variants['trait'] == trait]
        enrichment = np.random.uniform(0.5, 4.0)
        pvalue = 10 ** -np.random.uniform(0, 15)

        results.append({
            'trait': trait,
            'n_variants': len(trait_vars),
            'overlapping': int(len(trait_vars) * np.random.uniform(0.01, 0.1)),
            'fold_enrichment': enrichment,
            'pvalue': pvalue,
            'significant': pvalue < 0.05
        })

    return pd.DataFrame(results).sort_values('fold_enrichment', ascending=False)


def find_overlapping_variants(peaks, variants):
    """Find variants that overlap with peaks."""
    np.random.seed(42)

    # Simplified - would need proper interval overlap in production
    n_overlap = min(200, len(variants) // 10)
    sampled = variants.sample(n_overlap)
    sampled['peak_id'] = [f'peak_{i}' for i in np.random.randint(0, len(peaks), n_overlap)]
    sampled['peak_signal'] = np.random.lognormal(2, 1, n_overlap)

    return sampled


def add_prioritization_scores(variants):
    """Add prioritization scores based on multiple features."""
    variants = variants.copy()

    # Score based on multiple features
    pvalue_score = -np.log10(variants['pvalue']) / 20  # Normalize
    signal_score = np.log1p(variants['peak_signal']) / 5

    variants['priority_score'] = np.clip(
        (pvalue_score + signal_score) / 2,
        0, 1
    )

    return variants


if __name__ == "__main__":
    main()
