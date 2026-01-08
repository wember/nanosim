# nanosim

Nanoscale simulation tools for statistical mechanics.

## What's This?

Monte Carlo simulations for studying statistical mechanics at the nanoscale. Built with Python and NumPy/SciPy for number crunching, with Plotly for visualization.

## Setup

### Prerequisites

- Python 3.11 or higher
- Make (optional, but recommended)

### Quick Start

#### Using Make (Recommended)

```bash
# Create virtual environment and install dependencies
make setup

# Activate the virtual environment
source venv/bin/activate
```

#### Manual Setup

Use these steps if you don't have `make` or want to go old school

```bash
# Create virtual environment
python3 -m venv venv

# Activate the virtual environment
source venv/bin/activate

# Install dependencies
pip install --upgrade pip
pip install -r requirements.txt
```

## Usage

```bash
# Activate the virtual environment
source venv/bin/activate

# Run simulations
cd creutz-sim
python sim.py
```

## Deactivating the Environment

```bash
deactivate
```

## Cleaning Up

To remove the virtual environment:

```bash
make clean
# or manually:
rm -rf venv
```
