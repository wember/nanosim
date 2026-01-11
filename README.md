# nanosim

Nanoscale simulation tools for statistical mechanics.

## Introduction

Monte Carlo simulations for studying statistical mechanics at the nanoscale. Built with Python and NumPy/SciPy for number crunching, with Plotly for visualization.

## Environment

### Prerequisites

- Python 3.11 or higher
- Make (optional, but recommended)

### Setup

#### Using Make (Recommended)

```bash
# Create virtual environment and install dependencies
make setup

# Activate the virtual environment
source venv/bin/activate
```

#### Manual Setup

```bash
# Create virtual environment
python3 -m venv venv

# Activate the virtual environment
source venv/bin/activate

# Install dependencies
pip install --upgrade pip
pip install -r requirements.txt
```

### Cleanup

```bash
# Remove the virtual environment
make clean-venv

# Deactivate when done
deactivate
```

## Simulation

### Running Simulations

```bash
# Run the Monte Carlo simulation
make sim

# Or manually:
source venv/bin/activate
cd creutz-sim
python sim.py
```

### Generating Plots

```bash
# Generate comparison plots from simulation data
make plot

# Or manually:
source venv/bin/activate
cd creutz-sim
python Sk_comparison.py
```

### Cleaning Data

```bash
# Remove simulation output directories
make clean
```
