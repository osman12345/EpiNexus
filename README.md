# Histone Mark Analyzer

A comprehensive Python-based web application for analyzing histone modifications from CUT&Tag and ChIP-seq data.

## Features

- **Multi-input support**: Start from FASTQ, BAM, or peak files
- **Multi-species**: Mouse (mm10, mm39), Human (hg38, hg19), Rat (rn6, rn7), and custom genomes
- **DiffBind integration**: Differential peak analysis with DESeq2
- **Multi-mark analysis**: Two-mark and three-mark chromatin state integration
- **Interactive visualization**: Volcano plots, MA plots, heatmaps, genome browser
- **Background job queue**: Handle long-running analyses asynchronously

## Supported Histone Marks

| Mark | Type | Description |
|------|------|-------------|
| H3K27ac | Active | Active enhancers and promoters |
| H3K27me3 | Repressive | Polycomb repression |
| H3K4me1 | Priming | Enhancer priming/poised |
| H3K4me3 | Active | Active promoters |
| H3K9me3 | Repressive | Heterochromatin |
| H3K36me3 | Active | Transcribed gene bodies |

## Quick Start

### Option 1: Docker (Recommended)

```bash
cd docker
docker-compose up -d
```

Access the application:
- Web UI: http://localhost:8501
- API: http://localhost:8000
- API Docs: http://localhost:8000/docs

### Option 2: Local Installation

1. **Install system dependencies**:
```bash
# Ubuntu/Debian
sudo apt install samtools bedtools r-base bowtie2

# macOS
brew install samtools bedtools r bowtie2
```

2. **Install R packages**:
```r
install.packages("BiocManager")
BiocManager::install(c("DiffBind", "ChIPseeker", "DESeq2", "csaw"))
BiocManager::install(c("TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db"))
```

3. **Install Python dependencies**:
```bash
pip install -r requirements.txt
```

4. **Start the application**:
```bash
# Start API server
uvicorn app.main:app --reload

# In another terminal, start Streamlit
streamlit run frontend/app.py
```

## Usage

### 1. Upload Sample Data

Create a sample sheet CSV with the following columns:
- `SampleID`: Unique sample identifier
- `Condition`: Treatment group (e.g., "Treatment", "Control")
- `Factor`: Histone mark (e.g., "H3K27ac")
- `Replicate`: Replicate number
- `bamReads`: Path to BAM file
- `Peaks`: Path to peak file (BED format)

Example:
```csv
SampleID,Condition,Factor,Replicate,bamReads,Peaks
Sample1_Trt,Treatment,H3K27ac,1,/data/sample1.bam,/data/sample1.bed
Sample2_Trt,Treatment,H3K27ac,2,/data/sample2.bam,/data/sample2.bed
Sample3_Ctrl,Control,H3K27ac,1,/data/sample3.bam,/data/sample3.bed
Sample4_Ctrl,Control,H3K27ac,2,/data/sample4.bam,/data/sample4.bed
```

### 2. Configure Comparison

Define your comparison:
- Group 1 (Treatment): Samples to compare
- Group 2 (Control): Reference samples
- Histone mark
- FDR threshold (default: 0.1)
- Log2 fold change threshold (default: 0.5)

### 3. Run Analysis

Click "Run DiffBind Analysis" and monitor progress in the Jobs tab.

### 4. View Results

- **Volcano Plot**: Visualize significant differential peaks
- **MA Plot**: Mean-difference plot
- **Peak Table**: Browse and filter peaks
- **Download**: Export CSV or PDF results

## Multi-Mark Integration

### Two-Mark Integration (H3K27ac vs H3K27me3)

Identifies antagonistic chromatin changes:
- **Strong Activation**: Gained H3K27ac + Lost H3K27me3
- **Strong Repression**: Lost H3K27ac + Gained H3K27me3

### Three-Mark Integration

Classifies chromatin states:
- **Active Enhancer**: H3K27ac + H3K4me1
- **Poised Enhancer**: H3K4me1 only
- **Active Promoter**: H3K27ac only
- **Bivalent**: H3K27me3 + H3K4me1
- **Repressed**: H3K27me3 only

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/health` | GET | Health check |
| `/samples` | GET/POST | List/create samples |
| `/comparisons` | GET/POST | List/create comparisons |
| `/jobs` | GET/POST | List/submit jobs |
| `/analysis/diffbind` | POST | Run DiffBind analysis |
| `/analysis/integration` | POST | Run multi-mark integration |
| `/results/{job_id}` | GET | Get analysis results |

## Project Structure

```
histone_analyzer/
├── app/                    # Backend API
│   ├── main.py            # FastAPI application
│   ├── config.py          # Configuration
│   ├── models/            # Database models
│   └── core/              # Analysis modules
├── frontend/              # Streamlit UI
│   ├── app.py            # Main application
│   └── components/       # Reusable components
├── scripts/              # R and bash scripts
├── docker/               # Docker configuration
└── tests/                # Unit tests
```

## Configuration

Environment variables:
- `DATABASE_URL`: Database connection string
- `REDIS_URL`: Redis URL for job queue
- `DEFAULT_GENOME`: Default reference genome (mm10)

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Run tests: `pytest`
5. Submit a pull request

## License

MIT License

## Citation

If you use this tool in your research, please cite:
- DiffBind: Ross-Innes et al. (2012) Nature
- ChIPseeker: Yu et al. (2015) Bioinformatics
