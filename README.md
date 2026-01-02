# Nanosim - Microcanonical Monte Carlo Simulation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

A high-performance Python implementation of microcanonical ensemble Monte Carlo simulation for a 1D Ising lattice using **Creutz's demon algorithm** (Creutz, 1983). This project explores thermodynamic irreversibility by comparing reversible and irreversible dynamics.

**Key Features:**

- 🚀 **~8,900x speedup** with parallel JIT compilation (Numba + multiprocessing)
- 🔄 **Reversible vs irreversible dynamics** for studying entropy production
- 📊 **Interactive visualizations** with Plotly
- ✅ **337 unit tests** with 84% coverage ensuring correctness
- 🎯 **Production-ready** with comprehensive documentation

## What Does This Simulate?

This code simulates a **1D chain of magnetic spins** (an Ising model) that can flip up or down while conserving total energy. Think of it like a row of tiny magnets that influence their neighbors.

**The Physics:**

- **Microcanonical ensemble**: Total energy is fixed (isolated system)
- **Creutz's demon**: An energy "bookkeeper" that allows spins to flip by borrowing/returning energy
- **Reversible vs Irreversible**: We test whether running dynamics backward produces the same entropy evolution

**Why It Matters:**

- Explores fundamental questions about thermodynamic irreversibility
- Tests whether microscopic reversibility affects macroscopic entropy production
- Relevant to statistical mechanics, thermodynamics, and computational physics

> **Note**: The Creutz demon algorithm is a well-established computational physics method (Creutz, 1983). This implementation focuses on the novel comparison of reversible vs. irreversible dynamics using the algorithm.

## Quick Start

Get running in under 2 minutes!

### 1. Installation

Clone and set up the environment:

```bash
git clone <repository-url>
cd nanosim
make setup    # Automated setup (creates venv, installs dependencies, compiles JIT)
make test-env # Verify installation works
```

This will install all dependencies and compile the JIT functions (~15-20 seconds one-time compilation). All subsequent runs will be fast.

**Requirements:** Python 3.11+, pip

**Alternative:** Run `./setup.sh` then `make compile` if Make is unavailable, or see [Manual Installation](#manual-installation) below.

### 2. Run Your First Simulation

Test with a small, fast simulation (completes in ~5 seconds):

```bash
make run-sim-small
```

This runs a reversible simulation with n=100 spins, s=10 sweeps. You should see:

- Progress bars showing forward and reverse phases
- Output: `data/r*/sim_data_*.csv` files with entropy and energy data

> **Note:** Since JIT functions were compiled during `make setup`, this runs at full speed immediately!

### 3. View Results

Generate plots of entropy evolution:

```bash
make plot
```

Opens interactive plots showing how entropy changes during forward vs reverse dynamics.

### 4. Try Different Parameters

Run with custom settings:

```bash
# Medium-sized simulation (30 seconds)
make run-sim ARGS="--n 10000 --s 1000"

# Full production run (1-2 minutes with JIT)
# Uses: n=1000000, s=5000, r=11 (tests R=1-10), m=5 runs
make run-sim
```

**Common parameters:**

- `--n`: Number of spins (default: 1000000)
- `--s`: Sweeps per phase (default: 5000)
- `--r`: Max radius to test (default: 11, tests R=1-10)
- `--m`: Independent runs for averaging (default: 5)

## Performance

**Full production runs** (`make run-sim`: n=1M, s=5k, r=11, m=5 = 55 simulations):

| Configuration                   | Total Time     | Speedup  |
| ------------------------------- | -------------- | -------- |
| Original (no optimization)      | ~87 hours      | 1x       |
| **Production (JIT + Parallel)** | **~1.8 hours** | **~48x** |

> **Benchmark System:** Apple M3 Max (16 cores). Performance will vary on different hardware.

**All production commands use parallel JIT by default** (`make run-sim`, `make run-sim-irr`).

> **First-Time Setup Note:** Running `make setup` compiles all JIT functions once (~15-20 seconds). After that, all simulations run at full speed immediately. If you skip `make setup` and run a simulation directly, the first run will include compilation time.

**Why Not GPU?** We evaluated GPU acceleration but determined it would only add 1.5-2x speedup at significant complexity cost. The sequential nature of energy-conserving Monte Carlo limits GPU effectiveness. See [OPTIMIZATIONS.md](docs/OPTIMIZATIONS.md#9-gpucuda-acceleration-not-implemented) for analysis.

## Documentation

### Core Documentation

Start here for understanding the physics and implementation:

- **[ARCHITECTURE.md](docs/ARCHITECTURE.md)** - Physics concepts, algorithm details, implementation specifics
- **[JIT_BEST_PRACTICES.md](docs/JIT_BEST_PRACTICES.md)** - Using Numba JIT optimization (~3,881x speedup single-core)
- **[BEST_PRACTICES.md](docs/BEST_PRACTICES.md)** - Development practices and code patterns
- **[CHANGELOG.md](CHANGELOG.md)** - Version history

### Reference Documents

Historical context explaining design decisions:

- **[OPTIMIZATIONS.md](docs/OPTIMIZATIONS.md)** - Detailed performance optimization journey with benchmarks
- **[JIT_OPTIMIZATION_SUMMARY.md](docs/JIT_OPTIMIZATION_SUMMARY.md)** - Benchmark results proving JIT effectiveness
- **[JIT_TEST_SUMMARY.md](docs/JIT_TEST_SUMMARY.md)** - Test methodology ensuring correctness
- **[PRODUCTION_REFINEMENTS.md](docs/PRODUCTION_REFINEMENTS.md)** - Workflow improvements
- **[PROFILE_RESULTS.md](docs/PROFILE_RESULTS.md)** - Profiling data

## Advanced Usage

### Command-Line Options

All simulation scripts accept these arguments:

```bash
# See all available options
python creutz-sim/parallel_sim.py --help

# Common patterns
make run-sim ARGS="--n 50000 --s 2000"           # Custom size
make run-sim ARGS="--cores 8"                    # Specify cores
make run-sim ARGS="--n 10000 --s 1000 --r 5"     # Test fewer radii
```

**Parameters:**

- `--n`: Lattice size (number of spins, default: 1000000)
- `--s`: Sweeps per phase (forward/reverse, default: 5000)
- `--r`: Max demon-coupling radius (tests R=1 to r-1, default: 11)
- `--m`: Independent runs for statistics (default: 5)
- `--cores`: CPU cores to use (default: auto-detect)
- `--no-jit`: Disable JIT compilation (default: JIT enabled)
- `--validate`: Validation mode - `off` (default), `periodic`, `frequent`

### Parallel Execution Details

The parallel implementation automatically distributes work across CPU cores:

```bash
# Auto-detect cores (recommended)
make run-sim

# Specify core count (useful for HPC)
make run-sim ARGS="--cores 16"

# Single-core (for debugging)
make run-sim ARGS="--cores 1"
```

**Progress tracking:**

During execution, you'll see real-time progress with estimated completion time:

- **ETA is adaptive**: Based on actual completion times of finished simulations on _your_ hardware
- Shows active parallel workers and current simulation details
- Updates every 10% of sweeps within each simulation

**Performance scaling (measured on Apple M3 Max):**

| Cores | Speedup vs Legacy      |
| ----- | ---------------------- |
| 1     | ~3,881x (JIT)          |
| 16    | ~8,902x (JIT+Parallel) |

**When to use what:**

- **Parallel + JIT** (default): Multiple radii/runs, production work, thesis deadlines
- **Single-core + JIT**: Single radius/run, limited memory, simple testing
- **Legacy (no JIT)**: Educational reference, debugging core algorithm

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
make sbatch-sim-irr             # Submit sequential irreversible simulation
make sbatch-sim                 # Submit reversible (parallel JIT)
make sbatch-sim-irr             # Submit irreversible (parallel JIT)

# Manual submission (advanced)
cd creutz-sim/batch_jobs
sbatch sim_sbatch.sh                  # Sequential reversible
sbatch sim_sbatch_irr.sh              # Sequential irreversible
sbatch sim_sbatch.sh                  # Reversible with parallel JIT (16 cores)
sbatch sim_sbatch_irr.sh              # Irreversible with parallel JIT (16 cores)
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

Generate plots from simulation data:

```bash
make plot    # Run default plotting script
```

For specific plots:

```bash
source venv/bin/activate
python creutz-sim/sim_plot.py         # Single simulation plot
python creutz-sim/Sk_comparison.py    # Compare entropy across radii
python creutz-sim/sim_plot_radii.py   # Radius comparison
```

### Manual Environment Activation

If you need to run commands manually:

```bash
source venv/bin/activate   # Activate environment

# Then run any Python script
python creutz-sim/parallel_sim.py
pytest tests/ -v
```

## Project Structure

```
nanosim/
├── creutz-sim/              # Main simulation code
│   ├── inferno.py           # Core reversible simulation class
│   ├── inferno_irr.py       # Core irreversible simulation class
│   ├── jit_inferno.py       # JIT-compiled reversible wrapper
│   ├── jit_inferno_irr.py   # JIT-compiled irreversible wrapper
│   ├── parallel_sim.py      # Production reversible (parallel + JIT)
│   ├── parallel_sim_irr.py  # Production irreversible (parallel + JIT)
│   ├── sim_plot.py          # Single simulation plotter
│   ├── sim_plot_radii.py    # Radius comparison plotter
│   ├── Sk_comparison.py     # Entropy comparison plotter
│   ├── legacy/              # Original single-core implementations
│   └── batch_jobs/          # SLURM HPC batch scripts
├── tests/                   # Comprehensive test suite (337 tests, 84% coverage)
├── tools/                   # Benchmarking and profiling tools
├── docs/                    # Detailed documentation
├── data/                    # Output data (auto-generated)
├── logs/                    # Execution logs (auto-generated)
├── Makefile                 # Task automation
└── requirements.txt         # Python dependencies
```

## Manual Installation

If `make setup` doesn't work, install manually:

````bash
# Clone repository
git clone <repository-url>
cd nanosim

# Create virtual environment
python3 -m venv venv

# Activate (macOS/Linux)
source venv/bin/activate

# Activate (Windows)
# venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Verify
python -c "import numpy, scipy, pandas, numba; print('✓ All dependencies installed')"
```ingle test file (fast)
make run-test-file FILE=test_inferno.py

# Filter by test name
make run-tests ARGS="-k energy_conservation"

# Verbose output with print statements
make run-tests ARGS="-v -s"
````

**Pro tip:** Use `make run-test-file` during development, then `make run-tests` before committing.

### Runtime Estimates

**Quick reference for planning simulation runs:**

| Configuration | n      | s     | Time (1 core + JIT) | Time (16 cores + JIT) |
| ------------- | ------ | ----- | ------------------- | --------------------- |
| Test          | 100    | 10    | <0.1 sec            | <0.1 sec              |
| Small         | 1,000  | 100   | ~0.1 sec            | ~0.1 sec              |
| Medium        | 10,000 | 1,000 | ~3-5 sec            | ~1 sec                |
| Production    | 1M     | 5,000 | ~20-40 sec/radius   | **~1-2 min total**    |

**Without JIT (legacy):** ~60-70x slower (test/educational use only)

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

## Output Files

Simulations generate three types of output:

### 1. CSV Data Files

Saved in `data/r{R}/sim_data_*.csv` with columns:

| Column | Description                              |
| ------ | ---------------------------------------- |
| `t`    | Sweep number (0 to 2s-1)                 |
| `K`    | Average demon energy per site            |
| `U`    | Lattice energy per site                  |
| `N0`   | Broken bonds per site                    |
| `Nx`   | Anti-aligned neighbor pairs per site     |
| `S/nk` | Total entropy per site (in units of k_B) |
| `n`    | Lattice size                             |

### 2. Metadata (JSON)

Each CSV has a companion `*_metadata.json` file with:

- Simulation parameters (n, s, r, m)
- Timestamp
- Simulation type (reversible/irreversible)

## Common Commands Reference

Quick reference for daily use:

```bash
make help              # Show all available commands
make setup             # Initial installation (includes JIT compilation)
make compile           # Recompile JIT functions if needed
make test-env          # Verify environment
make run-sim-small     # Quick test (~5 sec)
make run-sim           # Production run (~1-2 min)
make run-sim-irr       # Irreversible run (~1-2 min)
make plot              # Generate plots
make run-tests         # Run test suite (~30 sec)
make clean             # Clean up generated files
```

### Features

- ✅ **Command-line configuration** - No source file editing needed
- ✅ **Progress bars** - Visual feedback during long runs
- ✅ **Structured logging** - Detailed logs in `logs/` directory
- ✅ **Metadata tracking** - JSON files document each run
- ✅ **Type hints** - Better IDE support and code clarity
- ✅ **Comprehensive tests** - 150+ unit tests with pytest
- ✅ **Example scripts** - Three demos to get started
- ✅ **Portable paths** - Works anywhere without configuration
- 🚀 **JIT compilation** - ~3,881x speedup with Numba (single-core)
- 🚀 **Parallel execution** - ~2.3x additional speedup with multiprocessing (16 cores)

## Citation

If you use this software in your research, please cite it:

```bibtex
@software{ember2025nanosim,
  author = {Ember, Winry},
  title = {Nanosim: Microcanonical Monte Carlo Simulation of 1D Ising Lattice},
  year = {2025},
  url = {https://github.com/wember/nanosim},
  version = {0.1.0}
}
```

See [CITATION.cff](CITATION.cff) for more details.

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on:

- Reporting bugs
- Suggesting enhancements
- Submitting pull requests
- Development setup
- Testing procedures
- Code style guidelines

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Development, Testing, and Coverage

### Development Environment Setup

For development (including code formatting, linting, and coverage), run:

```bash
make setup-dev
```

This installs all base and development dependencies (see `requirements-dev.txt`).

### Adding and Running Tests

- Add new tests in the `tests/` directory. Name files as `test_*.py` and functions as `def test_*`.
- Run all tests:
  ```bash
  make run-tests
  ```
- Run a single test file:
  ```bash
  make run-test-file FILE=test_inferno.py
  ```
- Filter by test name:
  ```bash
  make run-tests ARGS="-k <test_name>"
  ```

### Coverage Reporting

To check code coverage:

```bash
make coverage
```

This runs the full test suite and generates a coverage report for the `creutz-sim` directory. For an HTML report:

```bash
make coverage-html
```

Open `htmlcov/index.html` in your browser to view detailed coverage results.

### Code Quality Tools

- Format code:
  ```bash
  make format
  ```
- Lint code:
  ```bash
  make lint
  ```

See [CONTRIBUTING.md](CONTRIBUTING.md) for more on development workflow and best practices.

## Development, Testing, and Coverage

### Development Environment Setup

For development (including code formatting, linting, and coverage), run:

```bash
make setup-dev
```

This installs all base and development dependencies (see `requirements-dev.txt`).

### Adding and Running Tests

- Add new tests in the `tests/` directory. Name files as `test_*.py` and functions as `def test_*`.
- Run all tests:
  ```bash
  make run-tests
  ```
- Run a single test file:
  ```bash
  make run-test-file FILE=test_inferno.py
  ```
- Filter by test name:
  ```bash
  make run-tests ARGS="-k <test_name>"
  ```

### Coverage Reporting

To check code coverage:

```bash
make coverage
```

This runs the full test suite and generates a coverage report for the `creutz-sim` directory. For an HTML report:

```bash
make coverage-html
```

Open `htmlcov/index.html` in your browser to view detailed coverage results.

### Code Quality Tools

- Format code:
  ```bash
  make format
  ```
- Lint code:
  ```bash
  make lint
  ```

See [CONTRIBUTING.md](CONTRIBUTING.md) for more on development workflow and best practices.

## Contact

- **Authors**: Winry Ember, Wayne Mock
- **Repository**: [github.com/wember/nanosim](https://github.com/wember/nanosim)
- **Issues**: [github.com/wember/nanosim/issues](https://github.com/wember/nanosim/issues)

## Acknowledgments

Based on Creutz's microcanonical demon algorithm:

- Creutz, M. (1983). "Microcanonical Monte Carlo Simulation". _Physical Review Letters_, 50(19), 1411-1414.
