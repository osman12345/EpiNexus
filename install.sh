#!/bin/bash
# EpiNexus Installation Script

echo "🧬 EpiNexus Installer"
echo "====================="
echo ""

# Check for conda
if command -v conda &> /dev/null; then
    echo "✓ Conda detected"
    read -p "Install with conda? (recommended) [Y/n]: " use_conda
    use_conda=${use_conda:-Y}

    if [[ $use_conda =~ ^[Yy]$ ]]; then
        echo ""
        echo "Creating conda environment 'epinexus'..."
        conda env create -f environment.yml
        echo ""
        echo "✅ Installation complete!"
        echo ""
        echo "The conda environment includes:"
        echo "  • Python packages for analysis"
        echo "  • Bowtie2 & BWA for alignment"
        echo "  • MACS2 for peak calling"
        echo "  • Samtools & Bedtools for BAM/BED processing"
        echo ""
        echo "To activate: conda activate epinexus"
        echo "To run:      streamlit run frontend/EpiNexus.py"
        exit 0
    fi
fi

# Fallback to pip
echo ""
echo "Installing with pip..."

# Check for venv
if [[ ! -d "venv" ]]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
fi

source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt

echo ""
echo "✅ Python packages installed!"
echo ""
echo "⚠️  Note: For preprocessing pipeline (FASTQ → Peaks), you need to install:"
echo ""
echo "  • Alignment tools:"
echo "    conda install -c bioconda bowtie2 bwa"
echo "    # or: brew install bowtie2 bwa (macOS)"
echo ""
echo "  • Peak calling:"
echo "    conda install -c bioconda macs2"
echo "    # or: pip install MACS2"
echo ""
echo "  • BAM/BED processing:"
echo "    conda install -c bioconda samtools bedtools"
echo "    # or: brew install samtools bedtools (macOS)"
echo ""
echo "To activate: source venv/bin/activate"
echo "To run:      streamlit run frontend/EpiNexus.py"
