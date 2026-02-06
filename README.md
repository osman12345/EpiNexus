# EpiNexus

A comprehensive **Python-only** epigenomics analysis platform for ChIP-seq, CUT&Tag, CUT&RUN, and ATAC-seq data. Analyze histone modifications, transcription factor binding, and chromatin accessibility with an intuitive web interface.

**No R required** - Uses PyDESeq2 for differential analysis (equivalent results to DESeq2/DiffBind).

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

## Quick Start

### Option 1: Conda (Recommended)

```bash
# Create environment
conda env create -f environment.yml
conda activate epinexus

# Run the app
./run.sh
```

Access at: **http://localhost:8501**

### Option 2: pip

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # or: venv\Scripts\activate on Windows

# Install dependencies
pip install -r requirements.txt

# Run
streamlit run frontend/EpiNexus.py
```

### Option 3: Docker

```bash
cd docker
docker-compose up -d
```

## Supported Assays

### Techniques (Methods)
| Technique | Description |
|-----------|-------------|
| ChIP-seq | Chromatin Immunoprecipitation Sequencing |
| CUT&Tag | Cleavage Under Targets and Tagmentation |
| CUT&RUN | Cleavage Under Targets and Release Using Nuclease |
| ATAC-seq | Assay for Transposase-Accessible Chromatin |

### Targets
| Target Type | Examples |
|-------------|----------|
| Histone Modifications | H3K27ac, H3K4me3, H3K4me1, H3K27me3, H3K9me3, H3K36me3 |
| Transcription Factors | Any TF with available antibody (MYC, P53, CTCF, etc.) |

## Python Dependencies

**Core:**
- Python ≥3.10
- Streamlit (web interface)
- Pandas, NumPy, SciPy (data processing)
- Plotly, Matplotlib, Seaborn (visualization)

**Bioinformatics:**
- PyDESeq2 (differential analysis - replaces R DiffBind/DESeq2)
- pysam (BAM file handling)
- pyranges (genomic interval operations - replaces R ChIPseeker)
- pybedtools (BED file operations)

**Optional CLI tools** (for FASTQ processing):
- samtools, bedtools, MACS2, bowtie2

## Usage Guide

### 1. Create a Project
Navigate to **Data & Project** and create a new project:
- Project name and description
- Reference genome (hg38, mm10, etc.)
- Technique (ChIP-seq, CUT&Tag, etc.)
- Target type (Histone or TF)

### 2. Upload Data
- **From Peaks**: Upload BED/narrowPeak files + optional BAM
- **From FASTQ**: Raw reads → alignment → peak calling (requires CLI tools)

### 3. Run Analysis
- Quality Control → Differential Analysis → Annotation → Visualization

### 4. Export Results
- CSV/Excel tables
- Publication-ready figures
- Methods section for papers

## Project Structure

```
histone_analyzer/
├── app/
│   └── core/
│       ├── differential.py    # PyDESeq2-based analysis
│       ├── diffbind.py        # DiffBind-compatible interface
│       ├── super_enhancers.py # ROSE algorithm
│       ├── tf_analysis.py     # TF motif analysis
│       └── atacseq.py         # ATAC-seq analysis
├── frontend/
│   ├── EpiNexus.py            # Main application
│   └── pages/                 # Analysis pages
├── environment.yml            # Conda environment (Python-only)
└── requirements.txt           # pip requirements
```

## Differential Analysis: Python vs R

EpiNexus uses **PyDESeq2** instead of R's DESeq2/DiffBind:

| Aspect | PyDESeq2 (EpiNexus) | DESeq2 (R) |
|--------|---------------------|------------|
| Results | Equivalent (<1% difference) | Reference |
| Installation | `pip install pydeseq2` | R + Bioconductor |
| Speed | Comparable | Comparable |
| Maintenance | Active development | Stable |

**Reference**: PyDESeq2 is a peer-reviewed Python implementation published in Bioinformatics.

## Citation

If you use EpiNexus in your research, please cite:

```bibtex
@software{epinexus,
  title = {EpiNexus: A comprehensive epigenomics analysis platform.},
  year = {2026},
  url = {https://github.com/osman12345/epinexus}
}
```



