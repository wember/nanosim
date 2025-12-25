#!/bin/bash
# Setup script for Nanosim project

set -e  # Exit on error

echo "🔧 Setting up Nanosim environment..."
echo

# Check if Python 3 is available
if ! command -v python3 &> /dev/null; then
    echo "❌ Error: Python 3 is not installed"
    echo "Please install Python 3.11 or higher"
    exit 1
fi

# Display Python version
PYTHON_VERSION=$(python3 --version)
echo "✓ Found $PYTHON_VERSION"
echo

# Create virtual environment if it doesn't exist
if [ -d "venv" ]; then
    echo "⚠️  Virtual environment already exists at ./venv"
    read -p "Do you want to recreate it? (y/N) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "🗑️  Removing existing virtual environment..."
        rm -rf venv
    else
        echo "Skipping virtual environment creation"
        exit 0
    fi
fi

echo "📦 Creating virtual environment..."
python3 -m venv venv
echo "✓ Virtual environment created"
echo

# Upgrade pip
echo "⬆️  Upgrading pip..."
./venv/bin/pip install --upgrade pip > /dev/null
echo "✓ pip upgraded"
echo

# Install dependencies
echo "📥 Installing dependencies from requirements.txt..."
./venv/bin/pip install -r requirements.txt
echo "✓ Dependencies installed"
echo

# Verify installation
echo "🔍 Verifying installation..."
./venv/bin/python -c "import numpy, scipy, pandas, plotly" && echo "✓ All packages imported successfully"
echo

echo "✅ Setup complete!"
echo
echo "To activate the virtual environment, run:"
echo "    source venv/bin/activate"
echo
echo "To run simulations:"
echo "    python creutz-sim/sim.py"
echo "    python creutz-sim/irr_sim.py"
