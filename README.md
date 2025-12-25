# Nanosim - Microcanonical Monte Carlo Simulation

A Python implementation of microcanonical ensemble Monte Carlo simulation for a 1D Ising lattice using **Creutz's demon algorithm** (Creutz, 1983). This project explores thermodynamic irreversibility by comparing reversible and irreversible dynamics.

> **Note**: The Creutz demon algorithm is a well-established computational physics method. This implementation focuses on the novel comparison of reversible vs. irreversible dynamics using the algorithm.

## Documentation

- **[ARCHITECTURE.md](ARCHITECTURE.md)** - Comprehensive technical documentation covering physics concepts, implementation details, and algorithm specifics
- **[OPTIMIZATIONS.md](OPTIMIZATIONS.md)** - Detailed explanations of all performance optimizations with benchmarks and rationale
- **[JIT_BEST_PRACTICES.md](JIT_BEST_PRACTICES.md)** - Guide for using Numba JIT optimization in production (70-106x speedup)
- **[BEST_PRACTICES.md](BEST_PRACTICES.md)** - Summary of implemented improvements and development practices
- **[CHANGELOG.md](CHANGELOG.md)** - Version history and changelog

## Performance

**Production scale simulations are now ~1400x faster!**

| Configuration      | Time (n=1M, s=10k, r=11, m=5) | Speedup    |
| ------------------ | ----------------------------- | ---------- |
| Original           | ~27 hours                     | 1x         |
| JIT only           | ~15-23 minutes                | 70-106x    |
| Parallel only      | ~2 hours                      | 13-14x     |
| **JIT + Parallel** | **~1.2 minutes**              | **~1400x** |

**Quick Start - Maximum Performance:**

```bash
make run-parallel-sim-jit        # Reversible with JIT (fastest)
make run-parallel-irr-sim-jit    # Irreversible with JIT (fastest)
```

See [OPTIMIZATIONS.md](OPTIMIZATIONS.md) for implementation details and [JIT_BEST_PRACTICES.md](JIT_BEST_PRACTICES.md) for usage guide.

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

# Parallel execution with JIT (RECOMMENDED - fastest)
make run-parallel-sim-jit           # Full parallel + JIT reversible (~1.2 min)
make run-parallel-irr-sim-jit       # Full parallel + JIT irreversible (~1.2 min)
make run-parallel-sim-jit-test      # Quick test with JIT

# Parallel execution without JIT (fast, but not maximum performance)
make run-parallel-sim-test          # Quick parallel test
make run-parallel-sim               # Full parallel reversible simulation (~2 hours)
make run-parallel-irr-sim           # Full parallel irreversible simulation (~2 hours)

make run-tests         # Run unit tests
make run-examples      # Run all example scripts
make plot              # Run plotting script
make clean             # Remove virtual environment and cache files
```

### Maximum Performance: JIT + Parallel (Recommended)

**New in v3.0:** Numba JIT compilation provides **70-106x speedup** on top of parallel processing, reducing production runs from 27 hours to just **1.2 minutes**.

**Quick start (fastest):**

```bash
# Run with JIT optimization (recommended for production)
make run-parallel-sim-jit          # Reversible (~1.2 minutes for full run)
make run-parallel-irr-sim-jit      # Irreversible (~1.2 minutes for full run)

# Or use command line directly
python creutz-sim/parallel_sim.py --jit
python creutz-sim/parallel_irr_sim.py --jit
```

**Performance comparison:**

| Cores | Without JIT | With JIT (Rev) | With JIT (Irr) |
| ----- | ----------- | -------------- | -------------- |
| 1     | 1x          | 70x            | 106x           |
| 4     | 3-4x        | 280x           | 424x           |
| 8     | 6-8x        | 560x           | 848x           |
| 16    | 13-14x      | ~1000x         | ~1400x         |

See [JIT_BEST_PRACTICES.md](JIT_BEST_PRACTICES.md) for detailed usage guide.

### Parallel Execution

Parallel versions can run multiple independent simulations simultaneously across all available CPU cores.

**Performance gains (without JIT):**

- **16-core system**: ~10-14x speedup
- **8-core system**: ~6-8x speedup
- **4-core system**: ~3-4x speedup

**Quick start:**

```bash
# Run with auto-detected CPU cores (no JIT)
make run-parallel-sim          # Reversible
make run-parallel-irr-sim      # Irreversible

# Or specify core count manually (useful for HPC SLURM jobs)
python creutz-sim/parallel_sim.py --cores 16
python creutz-sim/parallel_irr_sim.py --cores 8

# Add --jit for maximum performance
python creutz-sim/parallel_sim.py --jit --cores 16
```

**When to use parallel vs sequential:**

- **Use parallel** (`parallel_sim.py`):
  - Multiple radii (r > 2) and/or multiple runs (m > 2)
  - Multi-core system available (laptop, workstation, HPC node)
  - Time-critical analysis for thesis deadlines
  - **Always add --jit for production runs**
- **Use sequential** (`sim.py`):
  - Single radius, single run (r=2, m=1)
  - Limited memory systems
  - Debugging or testing (original classes easier to debug)

The parallel implementation uses Python's `multiprocessing` module and automatically detects available CPU cores. On HPC clusters, it respects SLURM's `--cpus-per-task` allocation.

### Command-Line Arguments

All simulation scripts accept command-line arguments:

```bash
# Default parameters (n=1000000, s=5000, r=11, m=5)
python creutz-sim/parallel_sim.py --jit

# Custom parameters
python creutz-sim/parallel_sim.py --jit --n 100000 --s 1000 --r 5 --m 10

# Quick test
python creutz-sim/parallel_sim.py --jit --n 1000 --s 100 --r 3 --m 2

# Get help
python creutz-sim/parallel_sim.py --help
```

**Parameters:**

- `--n`: Lattice size (number of spins)
- `--s`: Number of sweeps per phase (forward/reverse)
- `--r`: Max demon-coupling radius (tests R=1 to r-1)
- `--m`: Number of independent runs for statistics
- `--validate`: Validation mode (default: `off`)
  - `off` - No automatic validation (fastest, recommended for production)
  - `periodic` - Validate every 100 sweeps (for testing)
  - `frequent` - Validate every sweep (debug mode, slower)

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
make sbatch-sim                 # Submit sequential reversible simulation
make sbatch-irr-sim             # Submit sequential irreversible simulation
make sbatch-parallel-sim        # Submit parallel reversible simulation (NEW)
make sbatch-parallel-irr-sim    # Submit parallel irreversible simulation (NEW)

# Manual submission (advanced)
cd creutz-sim/batch_jobs
sbatch sim_sbatch.sh                  # Sequential reversible
sbatch irr_sim_sbatch.sh              # Sequential irreversible
sbatch parallel_sim_sbatch.sh         # Parallel reversible (16 cores)
sbatch parallel_irr_sim_sbatch.sh     # Parallel irreversible (16 cores)
```

**Parallel vs Sequential on HPC:**

The parallel batch scripts request 16 cores (`--cpus-per-task=16`) and run significantly faster:

- **Sequential**: 50 simulations × 2 hours each = 100 hours total
- **Parallel (16 cores)**: 50 simulations ÷ 16 cores = ~6-7 hours total

Edit the `#SBATCH --cpus-per-task=` line in `parallel_*_sbatch.sh` to match your cluster's resources or job limits.

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

Run the comprehensive test suite (81 tests):

```bash
make run-tests              # Run in parallel (~30 sec, default)
make run-tests-serial       # Run serially (~2 min, for debugging)
```

Run specific tests for faster iteration:

```bash
# Single test file (fast - great for development)
make run-test-file FILE=test_inferno.py

# Filter tests by name
make run-tests ARGS="-k energy_conservation"

# Run with verbose output and print statements
make run-tests ARGS="-v -s"

# Combine options
make run-tests ARGS="-k validation -v"
```

**Pro tip:** Use `make run-test-file` when developing/fixing a specific test file, then run `make run-tests` (parallel by default) before committing.

Or run manually:

```bash
source venv/bin/activate
pytest tests/ -v -n auto            # Parallel (default)
pytest tests/ -v                    # Serial
pytest tests/test_inferno.py -v    # Single file
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
- 🚀 **JIT compilation** - 70-106x speedup with Numba (optional)
- 🚀 **Parallel execution** - 13-14x speedup with multiprocessing

## Performance

**Optimization Level:**

- **Original implementation:** Baseline
- **With all optimizations:** 1.7x faster (random sign generation, index cycling, neighbor arrays)
- **With parallel processing (16 cores):** 13-14x faster
- **With JIT compilation:** 70-106x faster per core
- **Combined (parallel + JIT):** ~1000-1400x faster overall

**Production impact:** 27-hour simulation → **1.2 minutes**

See `OPTIMIZATIONS.md` for detailed performance analysis and `benchmark_jit.py` for benchmarking.

## Key Parameters

Configure via command-line arguments (see `--help` for details):

- `--n`: Lattice size (default: 1000000)
- `--s`: Number of sweeps (default: 5000)
- `--r`: Maximum demon-coupling radius (default: 11, testing R=1 to R=10)
- `--m`: Number of simulation runs for averaging (default: 5)

## Runtime Estimates

**Without JIT (original):**

- **Test** (n=100, s=10): ~1 second
- **Small** (n=1000, s=100): ~10 seconds
- **Medium** (n=10000, s=1000): ~5 minutes
- **Full** (n=1000000, s=5000): ~30-60 minutes per radius

**With JIT compilation (70-106x faster):**

- **Test** (n=100, s=10): <0.1 seconds
- **Small** (n=1000, s=100): ~0.1 seconds
- **Medium** (n=10000, s=1000): ~3-5 seconds
- **Full** (n=1000000, s=5000): ~20-40 seconds per radius

**With parallel + JIT (16 cores, ~1000x faster):**

- **Full production runs**: ~27 hours → **~1-2 minutes**

To use JIT versions: `from jit_inferno import JITInferno as Inferno`

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
