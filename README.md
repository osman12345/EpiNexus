# EpiNexus

A comprehensive epigenomics analysis platform for ChIP-seq, CUT&Tag, CUT&RUN, and ATAC-seq data. Analyze histone modifications, transcription factor binding, and chromatin accessibility with an intuitive web interface.

## Features

### Data Input & Processing
- **Multiple techniques**: ChIP-seq, CUT&Tag, CUT&RUN, ATAC-seq
- **Multiple targets**: Histone modifications and Transcription Factors
- **Flexible input**: Start from FASTQ, BAM, or peak files
- **Multi-species support**: Human (hg38, hg19), Mouse (mm10, mm39), Rat (rn6, rn7), custom genomes

### Core Analysis
- **Differential Analysis**: DiffBind integration with DESeq2 for robust statistical testing
- **Quality Control**: Comprehensive QC metrics (FRiP, NSC, RSC, fragment analysis)
- **Peak Annotation**: Gene annotation, pathway enrichment, regulatory element overlap
- **Visualization**: Heatmaps, profile plots, volcano plots, genome browser (IGV.js)

### Advanced Features
- **Multi-mark Integration**: Chromatin state analysis from multiple histone marks
- **Super-Enhancer Detection**: ROSE algorithm for identifying super-enhancers
- **TF Analysis**: Motif enrichment, binding site visualization, target gene prediction
- **Expression Integration**: Correlate epigenetic changes with gene expression
- **GWAS Overlap**: Link epigenomic features to disease-associated variants
- **Peak-Gene Linking**: ABC model for enhancer-gene predictions
- **ENCODE Integration**: Browse and compare with public reference data

## Supported Assays

### Techniques (Methods)
| Technique | Description | Use Cases |
|-----------|-------------|-----------|
| ChIP-seq | Chromatin Immunoprecipitation Sequencing | Gold standard for histone and TF profiling |
| CUT&Tag | Cleavage Under Targets and Tagmentation | Low-input, high signal-to-noise |
| CUT&RUN | Cleavage Under Targets and Release Using Nuclease | Low background, single-cell compatible |
| ATAC-seq | Assay for Transposase-Accessible Chromatin | Chromatin accessibility profiling |

### Targets
| Target Type | Examples |
|-------------|----------|
| Histone Modifications | H3K27ac, H3K4me3, H3K4me1, H3K27me3, H3K9me3, H3K36me3 |
| Transcription Factors | Any TF with available antibody (MYC, P53, CTCF, etc.) |

## Quick Start

### Option 1: Conda (Recommended)

```bash
# Create environment
conda env create -f environment.yml
conda activate epinexus

# Install R packages
Rscript -e 'BiocManager::install(c("DiffBind", "ChIPseeker", "DESeq2"))'

# Run the app
./run.sh
```

Access at: http://localhost:8501

### Option 2: Docker

```bash
cd docker
docker-compose up -d
```

Access the application:
- Web UI: http://localhost:8501
- API: http://localhost:8000/docs

### Option 3: Manual Installation

1. **System dependencies**:
```bash
# Ubuntu/Debian
sudo apt install samtools bedtools r-base bowtie2

# macOS
brew install samtools bedtools r bowtie2
```

2. **Python dependencies**:
```bash
pip install -r requirements.txt
```

3. **R packages**:
```r
BiocManager::install(c("DiffBind", "ChIPseeker", "DESeq2", "csaw"))
```

4. **Run**:
```bash
streamlit run frontend/EpiNexus.py
```

## Usage Guide

### 1. Create a Project
Navigate to **Data & Project** → Create a new project with:
- Project name and description
- Reference genome (hg38, mm10, etc.)
- Experimental technique (ChIP-seq, CUT&Tag, etc.)
- Target type (Histone Modification or Transcription Factor)

### 2. Upload Data

**From FASTQ files:**
- Upload raw sequencing reads
- Configure alignment parameters
- Automatic peak calling with MACS2

**From Peak files:**
- Upload pre-called peaks (BED, narrowPeak)
- Optional: BAM files for signal quantification
- Sample metadata CSV

### 3. Run Analysis

- **Quality Control**: Assess sample quality metrics
- **Differential Analysis**: Compare conditions using DiffBind
- **Annotation**: Annotate peaks with genes and pathways
- **Visualization**: Generate publication-ready figures

### 4. Export Results

- Download differential peak tables (CSV, Excel)
- Export publication-ready figures
- Generate methods section for papers

## Project Structure

```
histone_analyzer/
├── app/                          # Backend modules
│   ├── main.py                   # FastAPI application
│   ├── config.py                 # Configuration
│   └── core/                     # Analysis engines
│       ├── diffbind.py           # Differential binding
│       ├── super_enhancers.py    # ROSE algorithm
│       ├── tf_analysis.py        # TF analysis
│       ├── atacseq.py            # ATAC-seq analysis
│       ├── peak_gene_linking.py  # ABC model
│       └── ...
├── frontend/                     # Streamlit UI
│   ├── EpiNexus.py              # Main application
│   ├── pages/                    # Application pages
│   │   ├── 01_data_project.py    # Data & Project management
│   │   ├── 02_quality_control.py # QC metrics
│   │   ├── 03_differential.py    # Differential analysis
│   │   ├── 04_visualization.py   # Visualizations
│   │   ├── 05_annotation.py      # Peak annotation
│   │   ├── 06_multimark.py       # Multi-mark integration
│   │   ├── 08_tf_analysis.py     # TF motif analysis
│   │   ├── 12_super_enhancers.py # Super-enhancer detection
│   │   ├── 16_atacseq.py         # ATAC-seq analysis
│   │   ├── 19_gwas_overlap.py    # GWAS variant analysis
│   │   └── 21_help.py            # Help & documentation
│   └── components/               # Reusable UI components
├── scripts/                      # R and shell scripts
├── docker/                       # Docker configuration
├── data/                         # Sample data
└── tests/                        # Unit tests
```

## Analysis Modules

### Differential Analysis
Uses DiffBind with DESeq2 for:
- Consensus peak generation
- Read counting and normalization
- Statistical testing with multiple testing correction
- Volcano and MA plot visualization

### Multi-Mark Integration
Identifies chromatin states based on mark combinations:
- **Active Promoter**: H3K4me3+ H3K27ac+
- **Strong Enhancer**: H3K4me1+ H3K27ac+
- **Poised Enhancer**: H3K4me1+ H3K27ac-
- **Bivalent/Poised**: H3K4me3+ H3K27me3+
- **Repressed**: H3K27me3+

### Super-Enhancer Analysis
ROSE algorithm implementation:
- Stitch nearby H3K27ac peaks
- Rank by total signal
- Hockey-stick plot for cutoff determination
- Gene association and pathway enrichment

### TF ChIP-seq Analysis
Dedicated workflow for transcription factors:
- Peak characterization and genomic distribution
- Motif enrichment (known and de novo)
- Target gene identification
- Co-binding analysis with other TFs

## Configuration

Environment variables (optional):
```bash
EPINEXUS_GENOME=hg38          # Default genome
EPINEXUS_THREADS=4            # Parallel threads
EPINEXUS_DATA_DIR=/data       # Data directory
```

## API Reference

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/health` | GET | Health check |
| `/api/samples` | GET/POST | Sample management |
| `/api/analysis/diffbind` | POST | Run differential analysis |
| `/api/analysis/super-enhancers` | POST | Run SE detection |
| `/api/results/{job_id}` | GET | Get analysis results |

## Citation

If you use EpiNexus in your research, please cite:

```bibtex
@software{epinexus2024,
  title = {EpiNexus: Comprehensive Epigenomics Analysis Platform},
  year = {2024},
  url = {https://github.com/your-repo/epinexus}
}
```

Also cite the underlying tools:
- **DiffBind**: Ross-Innes et al. (2012) Nature
- **ChIPseeker**: Yu et al. (2015) Bioinformatics
- **MACS2**: Zhang et al. (2008) Genome Biology
- **DESeq2**: Love et al. (2014) Genome Biology

## Contributing

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/new-feature`
3. Make changes and add tests
4. Run tests: `pytest tests/`
5. Submit a pull request

## License

MIT License - see LICENSE file for details.

## Support

- **Documentation**: See the Help page in the application
- **Issues**: Report bugs via GitHub Issues
- **Contact**: your-email@example.com
