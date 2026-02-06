#!/bin/bash
# Launch EpiNexus

cd "$(dirname "$0")"

echo "🧬 Starting EpiNexus..."
echo ""
echo "The app will open at: http://localhost:8501"
echo "Press Ctrl+C to stop the server"
echo ""

python3 -m streamlit run frontend/EpiNexus.py --server.headless true
