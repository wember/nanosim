# Nanosim - Microcanonical Monte Carlo Simulation

A Python implementation of microcanonical ensemble Monte Carlo simulation for a 1D Ising lattice using Creutz's demon algorithm. This project explores thermodynamic irreversibility by comparing reversible and irreversible dynamics.

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
make help           # Show all available commands
make test           # Verify environment setup
make run-sim        # Run reversible simulation
make run-irr-sim    # Run irreversible simulation
make plot           # Run plotting script
make clean          # Remove virtual environment and cache files
```

### Manual Execution (Advanced)

If you need to run commands manually, first activate the virtual environment:

```bash
source venv/bin/activate
python creutz-sim/sim.py        # Reversible simulation
python creutz-sim/irr_sim.py    # Irreversible simulation
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
├── my_venv/              # Virtual environment (deprecated - use venv/)
├── requirements.txt      # Python dependencies
└── README.md            # This file
```

## Configuration

The simulation uses hostname-based path switching for portability between local and HPC environments. Modify paths in `sim.py` and `irr_sim.py` as needed.

## Key Parameters

Modify these in `sim.py` or `irr_sim.py`:

- `n`: Lattice size (default: 10000 for reversible, 1000000 for irreversible)
- `s`: Number of sweeps (default: 10000)
- `r`: Maximum demon-coupling radius (default: 11, testing R=1 to R=10)
- `m`: Number of simulation runs for averaging (default: 5)

## Output

Simulation data is saved as CSV files with columns:

- `t`: Sweep number
- `K`: Average demon energy per site
- `U`: Lattice energy per site
- `N0`: Broken bonds per site
- `Nx`: Anti-aligned neighbor pairs per site
- `S/nk`: Total entropy per site
- `n`: Lattice size

## License

[Add your license here]

## Citation

[Add citation information if applicable]
