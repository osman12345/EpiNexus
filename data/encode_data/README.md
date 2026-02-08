# ENCODE Data Directory

Store downloaded ENCODE data here.

## Directory Structure

```
encode_data/
├── peaks/          # Peak files (.bed, .narrowPeak, .broadPeak)
├── bigwig/         # Signal tracks (.bigWig, .bw)
└── metadata/       # Experiment metadata (.json, .tsv)
```

## Downloading from ENCODE

### Option 1: Using the EpiNexus Interface
1. Go to **ENCODE Data** page in EpiNexus
2. Search for experiments
3. Click "Download" - files will be saved here automatically

### Option 2: Manual Download from ENCODE Portal
1. Go to https://www.encodeproject.org
2. Search for your experiment
3. Download files and place them in the appropriate subdirectory

### Option 3: Command Line (xargs + wget)
```bash
# Download all H3K27ac ChIP-seq peaks for K562
curl -s "https://www.encodeproject.org/search/?type=Experiment&assay_title=Histone+ChIP-seq&target.label=H3K27ac&biosample_ontology.term_name=K562&files.file_type=bed+narrowPeak&format=json&limit=all" | \
  jq -r '.["@graph"][].files[] | select(.file_type=="bed narrowPeak") | .href' | \
  xargs -I {} wget -P peaks/ "https://www.encodeproject.org{}"
```

## Naming Convention

Recommended file naming:
```
{ExperimentID}_{Target}_{Biosample}_{FileType}.{ext}

Examples:
ENCFF123ABC_H3K27ac_K562_peaks.narrowPeak
ENCFF456DEF_H3K4me3_HepG2_signal.bigWig
```

## Supported File Types

| Type | Extensions | Description |
|------|------------|-------------|
| Peaks | `.bed`, `.narrowPeak`, `.broadPeak` | Called peaks |
| Signal | `.bigWig`, `.bw` | Normalized signal tracks |
| Metadata | `.json`, `.tsv` | Experiment info |
