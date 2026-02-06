"""
Unified Data Management System for EpiNexus

Provides:
- Central data store for all loaded data
- Clear demo vs real data indicators
- Session state management
- Data loading/saving utilities
"""

import streamlit as st
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
import json
from datetime import datetime


class DataSource(Enum):
    """Source of the data."""
    DEMO = "demo"
    USER_UPLOAD = "user_upload"
    PREPROCESSING = "preprocessing"
    ENCODE = "encode"
    EXTERNAL = "external"


@dataclass
class DatasetInfo:
    """Information about a loaded dataset."""
    name: str
    source: DataSource
    data_type: str  # 'peaks', 'expression', 'annotation', etc.
    n_records: int
    loaded_at: str
    description: str = ""
    file_path: Optional[str] = None
    metadata: Dict = field(default_factory=dict)


class DataManager:
    """
    Central data management for EpiNexus.

    Manages all data in session state with clear tracking of
    data sources (demo vs user data).
    """

    # Session state keys
    _STATE_KEYS = {
        'peaks': 'epinexus_peaks',
        'expression': 'epinexus_expression',
        'annotation': 'epinexus_annotation',
        'samples': 'epinexus_samples',
        'comparisons': 'epinexus_comparisons',
        'results': 'epinexus_results',
        'genes': 'epinexus_genes',
        'datasets_info': 'epinexus_datasets_info',
        'using_demo': 'using_demo_data',
    }

    @classmethod
    def init(cls):
        """Initialize session state for data management."""
        if 'epinexus_initialized' not in st.session_state:
            for key in cls._STATE_KEYS.values():
                if key not in st.session_state:
                    st.session_state[key] = None
            st.session_state['epinexus_datasets_info'] = {}
            st.session_state['using_demo_data'] = True
            st.session_state['epinexus_initialized'] = True

    @classmethod
    def load_data(
        cls,
        data_type: str,
        data: Any,
        source: DataSource,
        name: str = "",
        description: str = "",
        file_path: str = None,
        metadata: Dict = None
    ):
        """
        Load data into the session state.

        Args:
            data_type: Type of data ('peaks', 'expression', etc.)
            data: The actual data (DataFrame, dict, etc.)
            source: Where the data came from
            name: Human-readable name
            description: Description of the data
            file_path: Original file path if applicable
            metadata: Additional metadata
        """
        cls.init()

        key = cls._STATE_KEYS.get(data_type)
        if key:
            st.session_state[key] = data

            # Track dataset info
            n_records = len(data) if hasattr(data, '__len__') else 0

            info = DatasetInfo(
                name=name or data_type,
                source=source,
                data_type=data_type,
                n_records=n_records,
                loaded_at=datetime.now().isoformat(),
                description=description,
                file_path=file_path,
                metadata=metadata or {}
            )

            st.session_state['epinexus_datasets_info'][data_type] = info

            # Update demo status
            if source != DataSource.DEMO:
                st.session_state['using_demo_data'] = False

    @classmethod
    def get_data(cls, data_type: str) -> Optional[Any]:
        """Get data from session state."""
        cls.init()
        key = cls._STATE_KEYS.get(data_type)
        return st.session_state.get(key) if key else None

    @classmethod
    def get_info(cls, data_type: str) -> Optional[DatasetInfo]:
        """Get information about a loaded dataset."""
        cls.init()
        info_dict = st.session_state.get('epinexus_datasets_info', {})
        return info_dict.get(data_type)

    @classmethod
    def is_demo_data(cls, data_type: str = None) -> bool:
        """
        Check if current data is demo data.

        Args:
            data_type: Check specific data type, or None for overall status
        """
        cls.init()

        if data_type:
            info = cls.get_info(data_type)
            return info is None or info.source == DataSource.DEMO

        return st.session_state.get('using_demo_data', True)

    @classmethod
    def clear_data(cls, data_type: str = None):
        """Clear data from session state."""
        cls.init()

        if data_type:
            key = cls._STATE_KEYS.get(data_type)
            if key:
                st.session_state[key] = None
                info_dict = st.session_state.get('epinexus_datasets_info', {})
                info_dict.pop(data_type, None)
        else:
            # Clear all data
            for key in cls._STATE_KEYS.values():
                if key in st.session_state:
                    st.session_state[key] = None
            st.session_state['epinexus_datasets_info'] = {}
            st.session_state['using_demo_data'] = True

    @classmethod
    def get_all_datasets(cls) -> Dict[str, DatasetInfo]:
        """Get info about all loaded datasets."""
        cls.init()
        return st.session_state.get('epinexus_datasets_info', {})

    @classmethod
    def export_session(cls) -> Dict:
        """Export current session data for saving."""
        cls.init()

        export = {
            'metadata': {
                'exported_at': datetime.now().isoformat(),
                'version': '1.0'
            },
            'datasets': {},
            'info': {}
        }

        for data_type, key in cls._STATE_KEYS.items():
            if data_type in ['datasets_info', 'using_demo']:
                continue

            data = st.session_state.get(key)
            if data is not None:
                if isinstance(data, pd.DataFrame):
                    export['datasets'][data_type] = data.to_dict('records')
                elif isinstance(data, (list, dict)):
                    export['datasets'][data_type] = data

        export['info'] = {
            k: v.__dict__ if hasattr(v, '__dict__') else str(v)
            for k, v in st.session_state.get('epinexus_datasets_info', {}).items()
        }

        return export


def render_data_status_indicator():
    """
    Render a clear indicator of data source status.

    Call this at the top of each analysis page.
    """
    DataManager.init()

    if st.session_state.get('using_demo_data', True):
        st.markdown("""
        <div style="background: linear-gradient(90deg, #fff3cd, #ffeeba);
                    border: 1px solid #ffc107;
                    border-radius: 8px;
                    padding: 0.75rem 1rem;
                    margin-bottom: 1rem;
                    display: flex;
                    align-items: center;
                    gap: 0.5rem;">
            <span style="font-size: 1.2rem;">📌</span>
            <div>
                <strong style="color: #856404;">Demo Mode</strong>
                <span style="color: #856404; font-size: 0.9rem;">
                    — Viewing simulated data. Upload your own data or run preprocessing for real analysis.
                </span>
            </div>
        </div>
        """, unsafe_allow_html=True)
    else:
        # Show what data is loaded
        datasets = DataManager.get_all_datasets()
        user_datasets = [
            info for info in datasets.values()
            if info.source != DataSource.DEMO
        ]

        if user_datasets:
            dataset_names = ", ".join([d.name for d in user_datasets[:3]])
            if len(user_datasets) > 3:
                dataset_names += f" + {len(user_datasets) - 3} more"

            st.markdown(f"""
            <div style="background: linear-gradient(90deg, #d4edda, #c3e6cb);
                        border: 1px solid #28a745;
                        border-radius: 8px;
                        padding: 0.75rem 1rem;
                        margin-bottom: 1rem;
                        display: flex;
                        align-items: center;
                        gap: 0.5rem;">
                <span style="font-size: 1.2rem;">✅</span>
                <div>
                    <strong style="color: #155724;">Your Data Loaded</strong>
                    <span style="color: #155724; font-size: 0.9rem;">
                        — {dataset_names}
                    </span>
                </div>
            </div>
            """, unsafe_allow_html=True)


def render_data_source_badge(data_type: str):
    """
    Render a small badge indicating data source.

    Args:
        data_type: The type of data to check
    """
    DataManager.init()

    info = DataManager.get_info(data_type)

    if info is None or info.source == DataSource.DEMO:
        st.markdown("""
        <span style="background: #fff3cd; color: #856404;
                     padding: 0.2rem 0.5rem; border-radius: 4px;
                     font-size: 0.75rem; font-weight: 600;">
            📌 DEMO DATA
        </span>
        """, unsafe_allow_html=True)
    else:
        source_labels = {
            DataSource.USER_UPLOAD: "📤 UPLOADED",
            DataSource.PREPROCESSING: "⚙️ PROCESSED",
            DataSource.ENCODE: "🏛️ ENCODE",
            DataSource.EXTERNAL: "🔗 EXTERNAL"
        }
        label = source_labels.get(info.source, "✅ YOUR DATA")

        st.markdown(f"""
        <span style="background: #d4edda; color: #155724;
                     padding: 0.2rem 0.5rem; border-radius: 4px;
                     font-size: 0.75rem; font-weight: 600;">
            {label}
        </span>
        """, unsafe_allow_html=True)


def render_data_summary_sidebar():
    """Render a data summary in the sidebar."""
    DataManager.init()

    with st.sidebar:
        st.markdown("#### 📊 Loaded Data")

        datasets = DataManager.get_all_datasets()

        if not datasets:
            st.caption("No data loaded yet")
            return

        for data_type, info in datasets.items():
            icon = "📌" if info.source == DataSource.DEMO else "✅"
            st.markdown(f"""
            <div style="display: flex; justify-content: space-between;
                        padding: 0.25rem 0; font-size: 0.85rem;">
                <span>{icon} {info.name}</span>
                <span style="color: #666;">{info.n_records:,}</span>
            </div>
            """, unsafe_allow_html=True)


def render_load_demo_button(data_type: str, demo_loader_func, button_text: str = "Load Demo Data"):
    """
    Render a button to load demo data with clear labeling.

    Args:
        data_type: Type of data to load
        demo_loader_func: Function that returns demo data
        button_text: Text for the button
    """
    col1, col2 = st.columns([1, 3])

    with col1:
        if st.button(f"📌 {button_text}", key=f"load_demo_{data_type}"):
            demo_data = demo_loader_func()

            DataManager.load_data(
                data_type=data_type,
                data=demo_data,
                source=DataSource.DEMO,
                name=f"Demo {data_type.title()}",
                description="Simulated demonstration data"
            )

            st.success("Demo data loaded!")
            st.rerun()

    with col2:
        st.caption("⚠️ This loads simulated data for demonstration purposes only.")


def render_upload_data_section(
    data_type: str,
    accepted_types: List[str],
    parser_func,
    label: str = "Upload Data"
):
    """
    Render a standardized data upload section.

    Args:
        data_type: Type of data being uploaded
        accepted_types: List of accepted file extensions
        parser_func: Function to parse uploaded file
        label: Label for the upload section
    """
    st.markdown(f"#### 📤 {label}")

    uploaded = st.file_uploader(
        f"Choose file ({', '.join(accepted_types)})",
        type=accepted_types,
        key=f"upload_{data_type}"
    )

    if uploaded:
        try:
            data = parser_func(uploaded)

            DataManager.load_data(
                data_type=data_type,
                data=data,
                source=DataSource.USER_UPLOAD,
                name=uploaded.name,
                description=f"Uploaded from {uploaded.name}",
                file_path=uploaded.name
            )

            n_records = len(data) if hasattr(data, '__len__') else "N/A"
            st.success(f"✅ Loaded {n_records} records from {uploaded.name}")

            return data

        except Exception as e:
            st.error(f"Error parsing file: {str(e)}")
            return None

    return None


# Demo data generators
def generate_demo_peaks(
    n_peaks: int = 5000,
    assay_type: str = "ChIP-seq",
    mark_or_tf: str = "H3K27ac"
) -> pd.DataFrame:
    """Generate realistic demo peak data."""
    np.random.seed(42)

    chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX']
    chr_weights = np.array([8, 7, 6, 5, 5, 5, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 3])
    chr_weights = chr_weights / chr_weights.sum()

    # Adjust parameters based on assay type
    if assay_type == "ATAC-seq":
        width_mean, width_std = 150, 50
        signal_mean = 3.5
    elif "TF" in mark_or_tf or mark_or_tf in ["MYC", "P53", "CTCF"]:
        width_mean, width_std = 200, 60
        signal_mean = 4.0
    else:  # Histone
        width_mean, width_std = 500, 200
        signal_mean = 3.0

    peaks = pd.DataFrame({
        'chr': np.random.choice(chromosomes, n_peaks, p=chr_weights),
        'start': np.random.randint(1000000, 200000000, n_peaks),
        'peak_id': [f'peak_{i}' for i in range(n_peaks)],
        'signal': np.random.lognormal(signal_mean, 1, n_peaks),
        'pvalue': 10 ** -np.random.uniform(2, 20, n_peaks),
        'qvalue': 10 ** -np.random.uniform(1, 15, n_peaks)
    })

    widths = np.maximum(50, np.random.normal(width_mean, width_std, n_peaks).astype(int))
    peaks['end'] = peaks['start'] + widths
    peaks['summit'] = peaks['start'] + widths // 2

    return peaks[['chr', 'start', 'end', 'peak_id', 'signal', 'summit', 'pvalue', 'qvalue']]


def generate_demo_expression(n_genes: int = 500) -> pd.DataFrame:
    """Generate realistic demo expression data."""
    np.random.seed(42)

    gene_symbols = [f'GENE{i}' for i in range(n_genes)]

    # Add some known genes
    known = ['MYC', 'TP53', 'BRCA1', 'EGFR', 'KRAS', 'PTEN', 'BCL2', 'GAPDH', 'ACTB']
    gene_symbols[:len(known)] = known

    expression = pd.DataFrame({
        'gene_id': [f'ENSG{i:011d}' for i in range(n_genes)],
        'gene_symbol': gene_symbols,
        'baseMean': np.random.lognormal(5, 2, n_genes),
        'log2FoldChange': np.random.normal(0, 1.5, n_genes),
        'lfcSE': np.random.uniform(0.1, 0.5, n_genes),
        'pvalue': 10 ** -np.random.uniform(0, 10, n_genes),
        'padj': 10 ** -np.random.uniform(0, 8, n_genes)
    })

    return expression


def generate_demo_samples(n_samples: int = 8) -> pd.DataFrame:
    """Generate demo sample sheet."""
    conditions = ['Control', 'Treatment']
    marks = ['H3K27ac', 'H3K4me3']

    samples = []
    for mark in marks:
        for condition in conditions:
            for rep in range(1, n_samples // 4 + 1):
                samples.append({
                    'SampleID': f'{condition}_{mark}_rep{rep}',
                    'Condition': condition,
                    'Factor': mark,
                    'Replicate': rep,
                    'bamReads': f'/data/{condition}_{mark}_rep{rep}.bam',
                    'Peaks': f'/data/{condition}_{mark}_rep{rep}_peaks.bed'
                })

    return pd.DataFrame(samples)
