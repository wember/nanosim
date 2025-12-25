# Nanosim - Microcanonical Monte Carlo Simulation

A Python implementation of microcanonical ensemble Monte Carlo simulation for a 1D Ising lattice using **Creutz's demon algorithm** (Creutz, 1983). This project explores thermodynamic irreversibility by comparing reversible and irreversible dynamics.

> **Note**: The Creutz demon algorithm is a well-established computational physics method. This implementation focuses on the novel comparison of reversible vs. irreversible dynamics using the algorithm.

## Documentation

- **[ARCHITECTURE.md](ARCHITECTURE.md)** - Comprehensive technical documentation covering physics concepts, implementation details, and algorithm specifics
- **[BEST_PRACTICES.md](BEST_PRACTICES.md)** - Summary of implemented improvements and development practices
- **[CHANGELOG.md](CHANGELOG.md)** - Version history and changelog

## Installation

### Prerequisites

- Python 3.11 or higher
- pip

### Quick Setup (Recommended)

Use the Makefile for automated setup:

```bash
make setup    # Create virtual environment and install dependencies
make test     # Verify installation
```

This will:

- Create a virtual environment in `venv/`
- Upgrade pip
- Install all dependencies from `requirements.txt`
- Verify the installation

**Alternative:** Run `./setup.sh` directly if Make is not available.

### Manual Setup

If you prefer manual setup:

1. Clone the repository:

```bash
git clone <repository-url>
cd nanosim
```

2. Create and activate a virtual environment:

```bash
# Create virtual environment
python3 -m venv venv

# Activate on macOS/Linux
source venv/bin/activate

# Activate on Windows
# venv\Scripts\activate
```

3. Install dependencies:

```bash
pip install -r requirements.txt
```

### Verify Installation

```bash
make test
```

## Usage

### Quick Start

Use the Makefile for all common tasks:

```bash
make help              # Show all available commands
make test              # Verify environment setup
make run-sim-test      # Quick test (n=100, s=10, ~seconds)
make run-sim-small     # Small test (n=1000, s=100, ~minutes)
make run-sim           # Full simulation (n=1000000, s=5000, ~hours)
make run-irr-sim       # Irreversible simulation (full parameters)
make run-tests         # Run unit tests
make run-examples      # Run all example scripts
make plot              # Run plotting script
make clean             # Remove virtual environment and cache files
```

### Command-Line Arguments

Both simulation scripts now accept command-line arguments:

```bash
# Default parameters (n=1000000, s=5000, r=11, m=5)
python creutz-sim/sim.py

# Custom parameters
python creutz-sim/sim.py --n 1000 --s 100 --r 3 --m 2

# Get help
python creutz-sim/sim.py --help
```

**Parameters:**

- `--n`: Lattice size (number of spins)
- `--s`: Number of sweeps per phase (forward/reverse)
- `--r`: Max demon-coupling radius (tests R=1 to r-1)
- `--m`: Number of independent runs for statistics

### Examples

The `examples/` directory contains three demonstration scripts:

1. **quick_test.py** - Minimal working example

   ```bash
   python examples/quick_test.py
   ```

2. **custom_parameters.py** - Parameter customization

   ```bash
   python examples/custom_parameters.py
   ```

3. **analysis_pipeline.py** - Complete workflow
   ```bash
   python examples/analysis_pipeline.py
   ```

Or run all at once:

```bash
make run-examples
```

### Manual Execution (Advanced)

If you need to run commands manually, first activate the virtual environment:

```bash
source venv/bin/activate
python creutz-sim/sim.py --n 1000 --s 100     # Reversible simulation
python creutz-sim/irr_sim.py --n 1000 --s 100 # Irreversible simulation
```

### HPC Cluster (SLURM)

**What is SLURM/sbatch?**

SLURM (Simple Linux Utility for Resource Management) is a workload manager used on HPC clusters. The `sbatch` command submits batch jobs to the cluster's job queue, allowing you to run long simulations in the background with allocated resources.

**When to use it:**

- Running simulations that take hours or days
- Need more computational resources than your laptop provides
- Working on a university/research HPC cluster (e.g., ASU's Agave or Sol)

**Important:** SLURM (the workload manager with `sbatch`) is **not available** for local installation on macOS/Windows via Homebrew or other package managers. The batch scripts are designed to run on HPC clusters only. For local development and testing, use `make run-sim` instead.

**Getting HPC Access at ASU:**

ASU students have free access to high-performance computing resources:

1. **Apply for an account** at https://cores.research.asu.edu/research-computing
2. **Available clusters:**
   - **Sol** - General purpose HPC cluster
   - **Agave** - High-memory computing cluster
3. **Requirements:**
   - ASU affiliation (students, faculty, staff)
   - Faculty advisor approval (for thesis work)
4. **Support:** Email rc-help@asu.edu for assistance

Once you have access, you can SSH to the cluster and use the batch scripts in this repository.

**Submitting jobs to a cluster:**

```bash
# Using Make (recommended - run from project root)
make sbatch-sim          # Submit reversible simulation
make sbatch-irr-sim      # Submit irreversible simulation

# Manual submission (advanced)
cd creutz-sim/batch_jobs
sbatch sim_sbatch.sh      # Reversible simulation
sbatch irr_sim_sbatch.sh  # Irreversible simulation
```

**Monitoring your jobs:**

```bash
squeue -u $USER          # Check job status
scancel <job_id>         # Cancel a job
cat slurm.<job_id>.out   # View output logs
```

### Visualization

```bash
make plot    # Run plotting script
```

For specific plots, use:

```bash
source venv/bin/activate
python creutz-sim/sim_plot.py         # Single simulation plot
python creutz-sim/Sk_comparison.py    # Compare entropy across radii
```

## Project Structure

```
nanosim/
├── creutz-sim/           # Main simulation code
│   ├── inferno.py        # Reversible simulation class
│   ├── irr_inferno.py    # Irreversible simulation class
│   ├── sim.py            # Reversible simulation runner
│   ├── irr_sim.py        # Irreversible simulation runner
│   ├── sim_plot.py       # Single simulation plotter
│   ├── sim_plot_r.py     # Radius comparison plotter
│   ├── Sk_comparison.py  # Entropy comparison plotter
│   └── batch_jobs/       # SLURM batch scripts
├── tests/                # Unit tests
│   ├── test_inferno.py   # Tests for reversible simulation
│   └── test_irr_inferno.py # Tests for irreversible simulation
├── examples/             # Example scripts
│   ├── quick_test.py     # Minimal working example
│   ├── custom_parameters.py # Parameter customization
│   └── analysis_pipeline.py # Complete workflow
├── data/                 # Output data (generated, in .gitignore)
├── logs/                 # Simulation logs (generated, in .gitignore)
├── requirements.txt      # Python dependencies
├── Makefile             # Task automation
├── CHANGELOG.md         # Version history
├── ARCHITECTURE.md      # Technical documentation
└── README.md            # This file
```

## Testing

Run the unit tests to verify everything works:

```bash
make run-tests
```

Or manually:

```bash
source venv/bin/activate
pytest tests/ -v
```

## Features

- ✅ **Command-line configuration** - No need to edit source files
- ✅ **Progress bars** - Visual feedback during long simulations
- ✅ **Structured logging** - Logs saved to `logs/` directory
- ✅ **Metadata output** - Each run saves JSON metadata alongside CSV data
- ✅ **Type hints** - Better IDE support and code clarity
- ✅ **Unit tests** - Comprehensive test suite with pytest
- ✅ **Examples** - Three demonstration scripts to get started
- ✅ **Portable paths** - Works on any machine without configuration

## Key Parameters

Configure via command-line arguments (see `--help` for details):

- `--n`: Lattice size (default: 1000000)
- `--s`: Number of sweeps (default: 5000)
- `--r`: Maximum demon-coupling radius (default: 11, testing R=1 to R=10)
- `--m`: Number of simulation runs for averaging (default: 5)

## Runtime Estimates

- **Test** (n=100, s=10): ~1 second
- **Small** (n=1000, s=100): ~10 seconds
- **Medium** (n=10000, s=1000): ~5 minutes
- **Full** (n=1000000, s=5000): ~30-60 minutes per radius

## Output

### CSV Data

Simulation data is saved in `data/r{R}/` directories with columns:

- `t`: Sweep number
- `K`: Average demon energy per site
- `U`: Lattice energy per site
- `N0`: Broken bonds per site
- `Nx`: Anti-aligned neighbor pairs per site
- `S/nk`: Total entropy per site
- `n`: Lattice size

### Metadata

JSON files alongside each CSV contain:

- Simulation parameters (n, s, r, m)
- Timestamp
- Simulation type (reversible/irreversible)

### Logs

Detailed execution logs saved in `logs/` directory with timestamps.

## License

[Add your license here]

## Citation

This project implements the Creutz demon algorithm for microcanonical Monte Carlo simulation:

**Primary Reference:**

> Creutz, M. (1983). "Microcanonical Monte Carlo Simulation." _Physical Review Letters_, 50(19), 1411-1414.
> DOI: [10.1103/PhysRevLett.50.1411](https://doi.org/10.1103/PhysRevLett.50.1411)

**Additional Reading:**

- Creutz, M., & Freedman, B. (1981). "A statistical approach to quantum mechanics." _Annals of Physics_, 132(2), 427-462.
- Newman, M. E. J., & Barkema, G. T. (1999). _Monte Carlo Methods in Statistical Physics_. Oxford University Press. (Chapter 3: Microcanonical Methods)

**This Implementation:**
The code in this repository is an original implementation developed for studying thermodynamic irreversibility. While the Creutz demon algorithm is a standard technique, the specific focus on comparing reversible vs. irreversible dynamics and the software implementation are original contributions.
