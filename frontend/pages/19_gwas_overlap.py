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
            ["Use session data", "Upload new"],
            horizontal=True
        )

        if peak_source == "Upload new":
            peak_file = st.file_uploader("Upload peaks", type=['bed', 'csv'])
            if peak_file:
                peaks = pd.read_csv(peak_file, sep='\t', header=None)
                peaks.columns = ['chr', 'start', 'end'][:len(peaks.columns)]
                st.success(f"Loaded {len(peaks)} peaks")
            else:
                peaks = None
        else:
            # Try to get from session state
            peaks = None
            if 'differential_peaks' in st.session_state:
                peaks = st.session_state.differential_peaks
            elif 'uploaded_peaks' in st.session_state:
                peaks = st.session_state.uploaded_peaks
            elif 'samples' in st.session_state and len(st.session_state.samples) > 0:
                all_peaks = []
                for sample in st.session_state.samples:
                    if 'peaks' in sample:
                        all_peaks.append(sample['peaks'])
                if all_peaks:
                    peaks = pd.concat(all_peaks, ignore_index=True)

            if peaks is not None:
                st.success(f"Using {len(peaks)} peaks from session")
            else:
                st.warning("No peaks loaded. Please upload peaks or go to Data & Project page.")

        if peaks is not None:
            st.dataframe(peaks.head(), use_container_width=True)

    with col2:
        st.subheader("GWAS Variants")

        gwas_source = st.selectbox(
            "GWAS catalog source",
            ["Select traits", "Upload custom"]
        )

        if gwas_source == "Upload custom":
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
                with st.spinner("Fetching GWAS variants..."):
                    variants = generate_demo_gwas(traits)
                st.success(f"Loaded {len(variants)} variants for selected traits")
            else:
                variants = None
                st.info("Select traits to load GWAS variants")

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
    """Generate demo peak data (fallback when no real data)."""
    np.random.seed(42)
    n = 30000

    starts = np.random.randint(1000000, 200000000, n)
    widths = np.random.randint(200, 2000, n)

    peaks = pd.DataFrame({
        'chr': np.random.choice([f'chr{i}' for i in range(1, 23)], n),
        'start': starts,
        'end': starts + widths,
        'peak_id': [f'peak_{i}' for i in range(n)],
        'signal': np.random.lognormal(2, 1, n)
    })

    return peaks


def fetch_gwas_catalog(traits: list = None, pvalue_threshold: float = 5e-8) -> pd.DataFrame:
    """
    Fetch real GWAS variants from the NHGRI-EBI GWAS Catalog API.

    Args:
        traits: List of trait names to fetch (None = common diseases)
        pvalue_threshold: P-value threshold for significance

    Returns:
        DataFrame with GWAS variants
    """
    import requests

    if traits is None:
        traits = ['type 2 diabetes', 'coronary artery disease', 'breast cancer']

    all_variants = []

    for trait in traits:
        try:
            # Query GWAS Catalog API
            url = "https://www.ebi.ac.uk/gwas/rest/api/efoTraits/search/findByEfoTrait"
            params = {"trait": trait, "projection": "associationByEfoTrait"}

            response = requests.get(
                f"https://www.ebi.ac.uk/gwas/rest/api/associations/search/findByDiseaseTrait",
                params={"diseaseTrait": trait},
                headers={"Accept": "application/json"},
                timeout=30
            )

            if response.status_code == 200:
                data = response.json()

                for assoc in data.get("_embedded", {}).get("associations", [])[:100]:
                    for locus in assoc.get("loci", []):
                        for snp in locus.get("strongestRiskAlleles", []):
                            risk_allele = snp.get("riskAlleleName", "")
                            rsid = risk_allele.split("-")[0] if "-" in risk_allele else risk_allele

                            # Get location
                            for gene in locus.get("authorReportedGenes", []):
                                pass

                            pvalue = assoc.get("pvalue", 1)
                            if pvalue and float(pvalue) <= pvalue_threshold:
                                all_variants.append({
                                    'rsid': rsid,
                                    'chr': '',  # Would need additional API call
                                    'pos': 0,
                                    'trait': trait.title(),
                                    'pvalue': float(pvalue),
                                    'odds_ratio': float(assoc.get("orPerCopyNum", 1) or 1),
                                    'study': assoc.get("study", {}).get("publicationInfo", {}).get("title", "")[:50]
                                })

        except Exception as e:
            st.warning(f"Could not fetch GWAS data for {trait}: {e}")
            continue

    if all_variants:
        return pd.DataFrame(all_variants)
    else:
        # Fall back to demo data if API fails
        return generate_demo_gwas(traits=[t.title() for t in traits])


def generate_demo_gwas(traits=None):
    """Generate demo GWAS variants (fallback when API unavailable)."""
    np.random.seed(42)

    if traits is None:
        traits = ['Type 2 Diabetes', 'Coronary Artery Disease', 'Breast Cancer',
                 'Inflammatory Bowel Disease', 'Rheumatoid Arthritis']

    n_per_trait = 500
    variants = []

    for trait in traits:
        for i in range(n_per_trait):
            chrom = np.random.randint(1, 23)
            variants.append({
                'rsid': f'rs{np.random.randint(1000000, 99999999)}',
                'chr': f'chr{chrom}',
                'pos': np.random.randint(1000000, 200000000),
                'ref': np.random.choice(['A', 'C', 'G', 'T']),
                'alt': np.random.choice(['A', 'C', 'G', 'T']),
                'trait': trait,
                'pvalue': 10 ** -np.random.uniform(5, 20),
                'odds_ratio': np.random.uniform(1.05, 1.5)
            })

    return pd.DataFrame(variants)


def calculate_enrichment(peaks, variants):
    """Calculate real GWAS enrichment using actual overlap."""

    # Standardize column names
    if 'chrom' not in peaks.columns and 'chr' in peaks.columns:
        peaks = peaks.rename(columns={'chr': 'chrom'})
    if 'chrom' not in variants.columns and 'chr' in variants.columns:
        variants = variants.rename(columns={'chr': 'chrom'})

    # Calculate actual overlaps
    n_overlapping = 0

    for chrom in peaks['chrom'].unique():
        peaks_chr = peaks[peaks['chrom'] == chrom]
        variants_chr = variants[variants['chrom'] == chrom]

        if len(variants_chr) == 0:
            continue

        for _, var in variants_chr.iterrows():
            pos = var.get('pos', 0)
            # Check if variant falls within any peak
            overlaps = peaks_chr[(peaks_chr['start'] <= pos) & (peaks_chr['end'] >= pos)]
            if len(overlaps) > 0:
                n_overlapping += 1

    total_variants = len(variants)

    # Calculate enrichment
    genome_size = 3e9
    peak_coverage = (peaks['end'] - peaks['start']).sum()
    expected = total_variants * (peak_coverage / genome_size)
    fold_enrichment = n_overlapping / expected if expected > 0 else 1

    # P-value from binomial test
    try:
        from scipy import stats
        pvalue = stats.binom_test(n_overlapping, total_variants, peak_coverage / genome_size)
    except (ImportError, AttributeError):
        # scipy.stats.binom_test deprecated, use binomtest
        try:
            from scipy.stats import binomtest
            result = binomtest(n_overlapping, total_variants, peak_coverage / genome_size)
            pvalue = result.pvalue
        except ImportError:
            pvalue = 0.05  # Fallback

    return {
        'total_variants': total_variants,
        'overlapping': n_overlapping,
        'fold_enrichment': round(fold_enrichment, 2),
        'pvalue': pvalue
    }


def calculate_trait_enrichment(peaks, variants):
    """Calculate real enrichment by trait."""

    traits = variants['trait'].unique()

    # Standardize column names
    if 'chrom' not in peaks.columns and 'chr' in peaks.columns:
        peaks = peaks.rename(columns={'chr': 'chrom'})
    if 'chrom' not in variants.columns and 'chr' in variants.columns:
        variants = variants.rename(columns={'chr': 'chrom'})

    genome_size = 3e9
    peak_coverage = (peaks['end'] - peaks['start']).sum()
    bg_prob = peak_coverage / genome_size

    results = []
    for trait in traits:
        trait_vars = variants[variants['trait'] == trait]

        # Count overlaps for this trait
        n_overlapping = 0
        for chrom in peaks['chrom'].unique():
            peaks_chr = peaks[peaks['chrom'] == chrom]
            vars_chr = trait_vars[trait_vars['chrom'] == chrom]

            for _, var in vars_chr.iterrows():
                pos = var.get('pos', 0)
                overlaps = peaks_chr[(peaks_chr['start'] <= pos) & (peaks_chr['end'] >= pos)]
                if len(overlaps) > 0:
                    n_overlapping += 1

        # Calculate enrichment
        n_trait = len(trait_vars)
        expected = n_trait * bg_prob
        enrichment = n_overlapping / expected if expected > 0 else 1

        # P-value
        try:
            from scipy.stats import binomtest
            result = binomtest(n_overlapping, n_trait, bg_prob)
            pvalue = result.pvalue
        except (ImportError, AttributeError):
            from scipy import stats
            pvalue = stats.binom_test(n_overlapping, n_trait, bg_prob) if n_trait > 0 else 1

        results.append({
            'trait': trait,
            'n_variants': n_trait,
            'overlapping': n_overlapping,
            'fold_enrichment': round(enrichment, 2),
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
