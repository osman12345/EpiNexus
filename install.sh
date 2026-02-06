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
        echo "To activate: conda activate epinexus"
        echo "To run:      python -m streamlit run frontend/app.py"
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
echo "✅ Installation complete!"
echo ""
echo "To activate: source venv/bin/activate"
echo "To run:      python -m streamlit run frontend/app.py"
