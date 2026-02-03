# nanosim

Nanoscale simulation tools for statistical mechanics.

## Introduction

Monte Carlo simulations for studying statistical mechanics at the nanoscale. Built with Python and NumPy/SciPy for number crunching, with Plotly for visualization.

## Official Thesis Data

The complete simulation data and interactive visualizations from the original thesis research are publicly available at:

**[https://plots.myember.org](https://plots.myember.org)**

This hosted instance provides:

- Full dataset of completed simulation runs
- Interactive Plotly visualizations with detailed hover information
- Parameter comparisons across different lattice sizes, radii, and sweep counts
- Combined reversible and irreversible dynamics analysis
- Downloadable simulation archives for further analysis

The hosted data represents the final validated results used in the thesis. You can use this repository to reproduce the simulations or run new parameter explorations.

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

### Viewing Results

#### Web Interface (Recommended)

The simulation browser provides a comprehensive web interface to view, analyze, and manage all your simulation runs:

```bash
# Launch the web browser interface
make browse

# Access at http://localhost:8000
```

Features:

- **Interactive plots** - View Plotly visualizations with hover tooltips showing radius, sweep, and entropy data
- **Archive management** - Browse all simulation runs with status, parameters, and completion times
- **Combine simulations** - Merge reversible and irreversible runs into combined archives
- **Add notes** - Document observations and findings for each simulation
- **Export/Import** - Share simulation archives as validated zip files
- **File browser** - Inspect all data files in each archive
- **Responsive design** - Works on desktop and mobile with dark/light themes

The browser automatically discovers all simulations in the `data/` directory and provides a clean interface for exploring results without manually running plot scripts.

#### Direct Plot Generation

For quick single-run visualization or automated workflows:

```bash
# Generate comparison plots from simulation data
make plot

# Or manually:
source venv/bin/activate
cd creutz-sim
python Sk_comparison.py
```

This creates a standalone HTML file with Plotly plots comparing reversible and irreversible dynamics.

### Cleaning Data

```bash
# Remove simulation output directories
make clean
```
