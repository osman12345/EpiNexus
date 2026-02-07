# EpiNexus

**A comprehensive epigenomics analysis platform for ChIP-seq, CUT&Tag, CUT&RUN, and ATAC-seq data.**

Analyze histone modifications, transcription factor binding, and chromatin accessibility with an intuitive web interface — all in Python, no R required.

---

## Table of Contents

- [Features](#features)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Usage Guide](#usage-guide)
- [Supported Assays](#supported-assays)
- [Analysis Workflows](#analysis-workflows)
- [Configuration](#configuration)
- [Project Structure](#project-structure)
- [API Reference](#api-reference)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)

---

## Features

### Data Input & Processing
- **Multiple techniques**: ChIP-seq, CUT&Tag, CUT&RUN, ATAC-seq
- **Multiple targets**: Histone modifications and Transcription Factors
- **Flexible input**: Start from FASTQ, BAM, or peak files
- **Multi-species support**: Human (hg38, hg19), Mouse (mm10, mm39), Rat (rn6, rn7)

### Core Analysis
- **Differential Analysis**: PyDESeq2-based (equivalent to DiffBind/DESeq2)
- **Quality Control**: Comprehensive QC metrics (FRiP, NSC, RSC, fragment analysis)
- **Peak Annotation**: PyRanges-based gene annotation and overlap analysis
- **Visualization**: Heatmaps, profile plots, volcano plots, genome browser (IGV.js)

### Advanced Features
- **Multi-mark Integration**: Chromatin state analysis from multiple histone marks
- **Super-Enhancer Detection**: ROSE algorithm implementation
- **TF Analysis**: Motif enrichment, binding site visualization
- **Expression Integration**: Correlate epigenetic changes with gene expression
- **GWAS Overlap**: Link epigenomic features to disease variants
- **Peak-Gene Linking**: ABC model for enhancer-gene predictions
- **ENCODE Integration**: Browse and compare with public reference data

### Why EpiNexus?

| Feature | EpiNexus | Traditional R Pipeline |
|---------|----------|------------------------|
| Language | Python only | R + Python |
| Installation | `pip install` | Complex R/Bioconductor setup |
| Results | Equivalent (PyDESeq2) | DESeq2/DiffBind |
| Interface | Web UI + API | Command line / RStudio |
| Learning curve | Minimal | Requires R knowledge |

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/osman12345/epinexus.git
cd epinexus
```

### 2. Set up environment (choose one)

**Option A: Conda (Recommended)**
```bash
conda env create -f environment.yml
conda activate epinexus
```

**Option B: pip**
```bash
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

**Option C: Docker**
```bash
cd docker && docker-compose up -d
```

### 3. Launch the application

```bash
streamlit run frontend/EpiNexus.py
```

Open your browser to **http://localhost:8501**

### 4. Start analyzing!

1. Go to **Data & Project** → Create a new project
2. Upload your peak files (BED/narrowPeak format)
3. Navigate to **Quality Control** to assess sample quality
4. Run **Differential Analysis** to find significantly changed regions
5. Explore results with **Visualization** and **Annotation**

---

## Installation

### System Requirements

- **Python**: 3.10 or higher
- **Memory**: 8GB RAM minimum (16GB recommended for large datasets)
- **Storage**: 1GB for installation + space for your data
- **OS**: Linux, macOS, or Windows (WSL recommended)

### Detailed Installation

#### Using Conda (Recommended)

Conda handles all dependencies including bioinformatics tools:

```bash
# Create and activate environment
conda env create -f environment.yml
conda activate epinexus

# Verify installation
python -c "from app.core.differential import DifferentialAnalyzer; print('OK')"
```

The conda environment includes command-line tools (samtools, bedtools, MACS2) needed for FASTQ processing.

#### Using pip

For pip installation, you'll need to install bioinformatics tools separately:

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install Python packages
pip install -r requirements.txt

# Optional: Install CLI tools for FASTQ processing
# Ubuntu/Debian:
sudo apt install samtools bedtools
# macOS:
brew install samtools bedtools
```

#### Using Docker

```bash
cd docker
docker-compose up -d

# Access at http://localhost:8501
# Stop with: docker-compose down
```

### Verifying Installation

Run the test suite to ensure everything is working:

```bash
pytest tests/ -v
```

---

## Usage Guide

### Creating a Project

1. Navigate to **Data & Project** page
2. Click **Create New Project**
3. Fill in project details:
   - **Name**: Descriptive project name
   - **Reference Genome**: hg38, hg19, mm10, mm39, rn6, or rn7
   - **Technique**: ChIP-seq, CUT&Tag, CUT&RUN, or ATAC-seq
   - **Target Type**: Histone Modification or Transcription Factor

### Uploading Data

EpiNexus accepts multiple input formats:

| Format | Extension | Description |
|--------|-----------|-------------|
| narrowPeak | `.narrowPeak` | ENCODE peak format with summit info |
| BED | `.bed` | Standard genomic intervals |
| broadPeak | `.broadPeak` | For broad histone marks (H3K27me3) |
| BAM | `.bam` | Aligned reads (optional, for read counting) |

**Minimum requirements:**
- At least 2 samples per condition for differential analysis
- Peak files are required; BAM files are optional but improve quantification

### Sample Metadata

Define your experimental design:

```
Sample_ID,Condition,Replicate,Factor
H3K27ac_Treat_1,Treatment,1,H3K27ac
H3K27ac_Treat_2,Treatment,2,H3K27ac
H3K27ac_Ctrl_1,Control,1,H3K27ac
H3K27ac_Ctrl_2,Control,2,H3K27ac
```

### Running Analysis

#### Quality Control
- **Sample Summary**: Read counts, mapping rates, duplication
- **Peak Metrics**: Peak counts, FRiP scores, width distributions
- **Correlation Analysis**: Sample clustering and PCA
- **Fragment Analysis**: Insert size distributions (ATAC-seq)

#### Differential Analysis
1. Select comparison: Treatment vs Control
2. Set thresholds: FDR (default 0.05), Log2FC (default 1.0)
3. Choose method: DESeq2 (recommended) or edgeR
4. Run analysis and explore results

#### Visualization
- **Volcano Plot**: Overview of differential peaks
- **Heatmaps**: Signal intensity across samples
- **Profile Plots**: Average signal at genomic features
- **Genome Browser**: IGV.js integration for detailed viewing

### Exporting Results

Export options available from each analysis page:
- **Tables**: CSV, Excel, TSV
- **Figures**: PNG, SVG, PDF
- **Reports**: Comprehensive HTML/PDF reports
- **Methods**: Auto-generated methods section for publications

---

## Supported Assays

### Techniques

| Technique | Description | Best For |
|-----------|-------------|----------|
| **ChIP-seq** | Chromatin Immunoprecipitation + Sequencing | Gold standard, high signal |
| **CUT&Tag** | Cleavage Under Targets and Tagmentation | Low input, low background |
| **CUT&RUN** | Cleavage Under Targets and Release Using Nuclease | Low input, sharp peaks |
| **ATAC-seq** | Assay for Transposase-Accessible Chromatin | Chromatin accessibility |

### Target Types

| Target | Examples | Peak Type |
|--------|----------|-----------|
| **Active Histone Marks** | H3K27ac, H3K4me3, H3K4me1 | Narrow |
| **Repressive Marks** | H3K27me3, H3K9me3 | Broad |
| **Transcription Elongation** | H3K36me3 | Gene body |
| **Transcription Factors** | MYC, P53, CTCF, STAT3 | Narrow |

---

## Analysis Workflows

### Workflow 1: Basic Differential Analysis

```
Peak Files → Quality Control → Differential Analysis → Annotation → Export
```

Best for: Comparing two conditions with histone ChIP-seq

### Workflow 2: Multi-mark Chromatin States

```
Multiple Marks → Multi-mark Page → State Learning → Differential States
```

Best for: Comprehensive chromatin characterization (e.g., H3K4me3 + H3K27ac + H3K27me3)

### Workflow 3: ATAC-seq Analysis

```
ATAC Peaks → ATAC-seq Page → Accessibility QC → Differential Accessibility → TF Footprinting
```

Best for: Chromatin accessibility studies

### Workflow 4: Expression Integration

```
Peaks + RNA-seq → Expression Page → Correlation → Target Gene Identification
```

Best for: Linking epigenetic changes to gene expression

---

## Configuration

### Environment Variables

Create a `.env` file in the project root:

```bash
# Server settings
EPINEXUS_PORT=8501
EPINEXUS_HOST=0.0.0.0

# Data paths
EPINEXUS_DATA_DIR=/path/to/data
EPINEXUS_TEMP_DIR=/tmp/epinexus

# Analysis defaults
EPINEXUS_FDR_THRESHOLD=0.05
EPINEXUS_LFC_THRESHOLD=1.0
EPINEXUS_NORMALIZE_METHOD=RLE

# Resource limits
EPINEXUS_MAX_WORKERS=4
EPINEXUS_MEMORY_LIMIT=8G
```

### Streamlit Configuration

Edit `.streamlit/config.toml`:

```toml
[server]
port = 8501
maxUploadSize = 500  # MB

[theme]
primaryColor = "#1E88E5"
backgroundColor = "#FFFFFF"
secondaryBackgroundColor = "#F5F5F5"
```

---

## Project Structure

```
epinexus/
├── app/
│   ├── core/                    # Core analysis modules
│   │   ├── differential.py      # PyDESeq2-based differential analysis
│   │   ├── diffbind.py          # DiffBind-compatible interface
│   │   ├── atacseq.py           # ATAC-seq specific analysis
│   │   ├── super_enhancers.py   # ROSE algorithm
│   │   ├── tf_analysis.py       # TF motif analysis
│   │   ├── batch_processor.py   # Batch job management
│   │   └── autosave.py          # State management
│   └── models/                  # Data models
├── frontend/
│   ├── EpiNexus.py              # Main Streamlit application
│   ├── pages/                   # Analysis pages (01-21)
│   └── components/              # Reusable UI components
│       ├── data_manager.py      # Data handling
│       └── empty_states.py      # Empty state UI
├── tests/                       # Unit tests
│   ├── test_differential.py
│   ├── test_atacseq.py
│   └── test_diffbind.py
├── docker/                      # Docker configuration
├── environment.yml              # Conda environment
├── requirements.txt             # pip requirements
└── README.md                    # This file
```

---

## API Reference

### Core Classes

#### DifferentialAnalyzer

```python
from app.core.differential import DifferentialAnalyzer, Sample, DifferentialConfig

# Create samples
samples = [
    Sample("treat_1", "Treatment", "H3K27ac", replicate=1, peak_file="peaks.bed"),
    Sample("treat_2", "Treatment", "H3K27ac", replicate=2, peak_file="peaks2.bed"),
    Sample("ctrl_1", "Control", "H3K27ac", replicate=1, peak_file="ctrl_peaks.bed"),
    Sample("ctrl_2", "Control", "H3K27ac", replicate=2, peak_file="ctrl_peaks2.bed"),
]

# Configure analysis
config = DifferentialConfig(
    comparison_name="Treatment_vs_Control",
    group1="Treatment",
    group2="Control",
    fdr_threshold=0.05,
    lfc_threshold=1.0
)

# Run analysis
analyzer = DifferentialAnalyzer()
results = analyzer.run(samples, config)

# Access results
print(f"Total peaks: {results.total_peaks}")
print(f"Significant: {results.significant_peaks}")
print(results.significant.head())
```

#### ATACSeqAnalyzer

```python
from app.core.atacseq import ATACSeqAnalyzer, generate_demo_atac_data

analyzer = ATACSeqAnalyzer()

# Calculate QC metrics
metrics = analyzer.calculate_qc_metrics(bam_stats, peaks, tss_scores)

# Analyze fragment sizes
dist = analyzer.analyze_fragment_sizes(fragment_sizes)

# Differential accessibility
results = analyzer.differential_accessibility(
    count_matrix,
    condition1=["treat_1", "treat_2"],
    condition2=["ctrl_1", "ctrl_2"]
)
```

### Convenience Functions

```python
from app.core.differential import run_differential_analysis

# Quick analysis with dictionaries
results = run_differential_analysis(
    samples=[
        {"sample_id": "s1", "condition": "Treatment", "histone_mark": "H3K27ac", "peak_file": "s1.bed"},
        {"sample_id": "s2", "condition": "Control", "histone_mark": "H3K27ac", "peak_file": "s2.bed"},
    ],
    comparison_name="Test",
    group1="Treatment",
    group2="Control",
    output_dir="results/"
)
```

---

## Troubleshooting

### Common Issues

#### "No module named 'pydeseq2'"

```bash
pip install pydeseq2
# or
conda install -c conda-forge pydeseq2
```

#### "pysam not available"

```bash
# Conda (recommended)
conda install -c bioconda pysam

# pip (may require compiler)
pip install pysam
```

#### Memory errors with large datasets

- Increase available memory or use batch processing
- Process chromosomes separately
- Use the batch processing page for large sample sets

#### Slow performance

- Enable multiprocessing: Set `EPINEXUS_MAX_WORKERS=8`
- Use SSD storage for data files
- Consider subsetting data for initial exploration

### Getting Help

1. Check the **Help** page in the application (page 21)
2. Search existing [GitHub Issues](https://github.com/osman12345/epinexus/issues)
3. Open a new issue with:
   - Python version (`python --version`)
   - Package versions (`pip freeze`)
   - Error message and stack trace
   - Sample data (if possible)

---

## Contributing

We welcome contributions! Here's how to get started:

### Development Setup

```bash
# Clone and install in development mode
git clone https://github.com/osman12345/epinexus.git
cd epinexus/histone_analyzer

# Create environment
conda env create -f environment.yml
conda activate epinexus

# Install development dependencies
pip install pytest black flake8

# Run tests
pytest tests/ -v
```

### Code Style

- Use [Black](https://black.readthedocs.io/) for formatting
- Follow PEP 8 guidelines
- Add docstrings to all public functions
- Include type hints

### Pull Request Process

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Run tests (`pytest tests/ -v`)
5. Commit (`git commit -m 'Add amazing feature'`)
6. Push (`git push origin feature/amazing-feature`)
7. Open a Pull Request

---

## Citation

If you use EpiNexus in your research, please cite:

```bibtex
@software{epinexus2026,
  title = {EpiNexus: A comprehensive Python-based epigenomics analysis platform},
  author = {Osman, A.},
  year = {2026},
  url = {https://github.com/osman12345/epinexus},
  note = {Version 1.0.0}
}
```

Also consider citing the underlying tools:
- **PyDESeq2**: Muzellec et al. (2023) - Python implementation of DESeq2
- **Streamlit**: Streamlit Inc. - Web framework

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- The [DESeq2](https://bioconductor.org/packages/DESeq2/) team for the statistical methods
- The [PyDESeq2](https://github.com/owkin/PyDESeq2) developers for the Python implementation
- The [scverse](https://scverse.org/) community for bioinformatics tools
- The [Streamlit](https://streamlit.io/) team for the web framework

---

**Happy analyzing!** If you find EpiNexus useful, please consider starring the repository on GitHub.
