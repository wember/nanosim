# Nanosim - Microcanonical Monte Carlo Simulation

## Project Overview

This project implements a **microcanonical ensemble Monte Carlo simulation** of a 1D Ising lattice with Creutz's demon algorithm. The simulation explores thermodynamic irreversibility by comparing reversible and irreversible dynamics.

## Architecture

### Core Simulation Classes

- **`inferno.py`**: `Inferno` class - reversible simulation with pre-computed random walk patterns
  - Uses fixed arrays `radius_spin`/`radius_bond` and their reversed counterparts
  - `demon_move()` and `demon_reverse()` traverse lattice in opposite orders for time reversibility
- **`irr_inferno.py`**: `irrInferno` class - irreversible simulation with truly random dynamics
  - Generates new random radii on each call to `demon_move()`/`demon_reverse()`
  - Otherwise identical structure to `Inferno`

### Simulation Runners

- **Production**: `parallel_sim.py` and `parallel_irr_sim.py` - Parallel JIT-compiled runners (~1400x speedup)
- **Legacy**: `legacy/sim.py` and `legacy/irr_sim.py` - Original single-core implementations (educational reference only)
- All runners output CSV data per radius/iteration to parameterized folder paths
- Production runners used by all Makefile targets (`make run-sim`, `make run-irr-sim`)

### Key Physics Concepts

- **Lattice**: 1D Ising model with nearest-neighbor interactions and periodic boundary conditions
- **Energy**: Total energy (2\*N) conserved between lattice energy `E_lattice` and demon energy `E_demon` (array of N oscillators)
- **Entropy**: Calculated as `(Sk + Su)/n` where:
  - `Sk`: Kinetic entropy from demon energy distribution (uses `loggamma`)
  - `Su`: Configurational entropy from bond states (broken bonds N0, anti-aligned spins Nx)
  - Uses `N0_exp = max(N0, 1)` to avoid log(0) in entropy calculation

### Data Flow

1. Initialize lattice (half +1 spins, half -1 spins) with random demon energy distribution
2. Each **sweep** = N attempted moves (spin flips + bond changes)
3. Forward phase (s sweeps) → data written to CSV → Reverse phase (s sweeps)
4. Multiple runs (M) per radius (R) generate statistical ensembles

## Environment & Dependencies

### Python Environment Setup

**Quick Setup** (recommended):

```bash
./setup.sh  # Automated setup script
make test-env # Verify installation
```

**Manual Setup**:

```bash
# Create new virtual environment
python3 -m venv venv

# Activate virtual environment
source venv/bin/activate  # macOS/Linux
# venv\Scripts\activate   # Windows

# Install dependencies
pip install -r requirements.txt
```

**Common Commands**:

```bash
make help        # Show all available make targets (organized by category)
make test-env    # Verify environment works
make run-sim     # Run reversible simulation (parallel JIT, fastest)
make run-irr-sim # Run irreversible simulation (parallel JIT, fastest)
make clean       # Remove venv and cache files
```

- **Dependencies**: Managed via `requirements.txt` at project root
- Core packages: `numpy>=2.0.0`, `scipy>=1.14.0`, `pandas>=2.0.0`, `plotly>=5.0.0`
- Note: Old `my_venv/` directory is deprecated; use `venv/` instead

### File Paths

The code uses **relative paths from the project root** for portability:

```python
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
folder = os.path.join(project_root, 'data') + os.sep
```

Data is organized in `data/` subdirectories relative to the project root, making the code portable across different machines and environments.

## Running Simulations

### Local Execution

```bash
source venv/bin/activate

# Production (recommended - parallel JIT)
python creutz-sim/parallel_sim.py --jit
python creutz-sim/parallel_irr_sim.py --jit

# Legacy single-core (educational only)
python creutz-sim/legacy/sim.py
python creutz-sim/legacy/irr_sim.py
```

### HPC Cluster (SLURM)

Submit batch jobs from `creutz-sim/batch_jobs/`:

```bash
sbatch batch_jobs/sim_sbatch.sh      # Reversible
sbatch batch_jobs/irr_sim_sbatch.sh  # Irreversible
```

Job specs: 1 node, 1 core, 1GB memory, 7-day limit

## Visualization & Analysis

### Plotting Scripts

- **`sim_plot.py`**: Plots single simulation CSV (demon energy, lattice temp, entropy vs sweeps)
- **`sim_plot_r.py`**: Likely plots results across radii variations
- **`Sk_comparison.py`**: Compares entropy between reversible/irreversible runs
  - Averages multiple runs per radius (finds all CSVs in `data/r{R}/` and `data/irr/r{R}/`)
  - Applies rolling average smoothing (default `bin_size=10`)
  - Uses consistent color scheme across radii (11 colors for R=0-10)

### Data Structure

CSV columns: `['t', 'K', 'U', 'N0', 'Nx', 'S/nk', 'n']`

- `t`: sweep number (forward: 0 to s-1, reverse: s to 2s-1)
- `K`: avg demon energy per site
- `U`: lattice energy per site
- `N0`: broken bonds per site
- `Nx`: anti-aligned neighbor pairs per site
- `S/nk`: total entropy per site (in units of k_B)

## Code Conventions

### Naming Patterns

- Classes: CamelCase (`Inferno`, `irrInferno`)
- Physics variables: single capital letters matching equations (`N`, `K`, `U`, `R`)
- Lambda functions for formulas: `Sk`, `Su` (entropy equations directly from theory)

### Common Modifications

- **Lattice size**: Change `n=10000` (or `n=1000000`) in sim files
- **Sweep count**: Change `s=10000` in sim files
- **Radius range**: Change `r=11` (tests R=1 to R=10)
- **Number of runs**: Change `m=5` for statistical averaging

### Output Management

- `sim.py` writes `sim_status.csv` with timestamped progress updates
- Data organized by radius: `data/r{R}/` and `data/irr/r{R}/` subdirectories
- Each run appends `_{M}.csv` suffix (e.g., `sim_data_r3_2.csv`)

## Important Notes

- **Energy conservation**: Always verify `E_total = E_lattice + sum(E_demon)` remains constant (2\*N)
- **Reversed dynamics**: The `rev_order` and reversed radius arrays must mirror `order` exactly for reversibility tests
- **Bond updates**: Bond state updates happen after spin flips in `spin_flip()` method
- **Entropy edge case**: N0=0 handled specially (uses 2^(N0+1) instead of 2^N0 in Su calculation)
- **Testing**: ALWAYS use `make run-tests` for running tests (uses pytest-xdist for ~3x faster parallel execution), never run `pytest` directly for full test suite
