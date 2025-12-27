# Nanosim Architecture Documentation

## Overview

Nanosim is a Monte Carlo simulation implementing **Creutz's demon algorithm** to study thermodynamic irreversibility in a 1D Ising lattice. The project compares reversible vs. irreversible dynamics by tracking entropy changes during forward and reverse simulation phases.

## Physics Background

### The 1D Ising Model

- **Lattice**: Linear chain of N spins, each can be +1 or -1
- **Periodic boundary conditions**: The chain wraps around (spin N connects to spin 1)
- **Nearest-neighbor interactions**: Energy depends on adjacent spin alignment
- **Bonds**: Connection between adjacent spins (aligned=-1, anti-aligned=+1, broken=0)

### Creutz's Demon Algorithm

Instead of traditional Monte Carlo methods using random thermal fluctuations, this uses a "demon":

- **Demon**: Array of N Einstein oscillators (energy carriers)
- **Energy conservation**: Total energy E_total = E_lattice + E_demon (constant at 2\*N)
- **Microcanonical ensemble**: Fixed total energy, exploring all microstates with that energy
- **Demon coupling**: Each demon oscillator can exchange energy with spins within radius R

### Key Physics Quantities

- **E_lattice**: Lattice potential energy (sum of bond energies)
- **E_demon**: Kinetic energy distributed among N oscillators
- **S_k** (Kinetic entropy): From demon energy distribution
- **S_u** (Configurational entropy): From lattice bond configurations
- **Total entropy**: S = (S_k + S_u) / N (per site, in units of k_B)

## System Architecture

```
nanosim/
├── Core Simulation Classes
│   ├── inferno.py              # Reversible simulation
│   ├── inferno_irr.py          # Irreversible simulation
│   ├── jit_inferno.py          # JIT-compiled reversible
│   └── jit_inferno_irr.py      # JIT-compiled irreversible
├── Production Runners
│   ├── parallel_sim.py         # Reversible (parallel + JIT)
│   └── parallel_sim_irr.py     # Irreversible (parallel + JIT)
├── Legacy (Single-Core)
│   ├── legacy/README.md        # Documentation for legacy code
│   ├── legacy/sim.py           # Original single-core reversible
│   └── legacy/irr_sim.py       # Original single-core irreversible
├── Visualization
│   ├── sim_plot.py             # Plot single simulation
│   ├── sim_plot_radii.py       # Plot across radii
│   └── Sk_comparison.py        # Compare rev vs. irr entropy
├── HPC Integration
│   └── batch_jobs/
│       ├── sim_sbatch.sh       # SLURM job: reversible
│       └── sim_sbatch_irr.sh   # SLURM job: irreversible
└── Infrastructure
    ├── venv/                   # Python virtual environment
    ├── requirements.txt        # Dependencies
    ├── Makefile               # Task automation
    └── setup.sh               # Environment setup
```

## Core Components

### 1. Inferno Class (`inferno.py`)

**Purpose**: Implements reversible dynamics using pre-computed random patterns.

**Key Design:**

```python
class Inferno:
    def __init__(self, N, R):
        # N = lattice size
        # R = demon coupling radius

        # Random walk patterns (fixed for entire simulation)
        self.order = random_permutation(0...N-1)
        self.rev_order = reversed(order)
        self.radius_spin = random_offsets(size=N, range=[-R, R])
        self.radius_bond = reversed(radius_spin)

        # Initial state
        self.lattice = [+1, +1, ..., -1, -1]  # Half up, half down
        self.bonds = [-1, -1, ..., 1, 1]      # Aligned/anti-aligned
        self.E_demon = random_energy_distribution(total=2*N)
```

**Critical Methods:**

1. **`demon_move()`** - Forward time evolution

   - Visits spins in `order` sequence
   - For each spin: attempt flip using `radius_spin` offset
   - For each bond: attempt change using `radius_bond` offset
   - Updates bond counts for entropy calculation

2. **`demon_reverse()`** - Backward time evolution

   - Visits spins in `rev_order` (reversed sequence)
   - Uses reversed radius arrays
   - **Order of operations reversed**: bond change, then spin flip
   - Should return to initial state if dynamics are reversible

3. **`spin_flip(a, i)`** - Attempt to flip spin at position `a`

   - Calculate energy cost using neighbor spins
   - If favorable (cost < 0): flip and give energy to demon[i]
   - If unfavorable (cost > 0): only flip if demon[i] has enough energy
   - Update adjacent bonds after flip

4. **`bond_change(a, i)`** - Attempt to break/remake bond at position `a`
   - If spins aligned: cost = -1 to break, +1 to make
   - If spins anti-aligned: cost = +1 to break, -1 to make
   - Exchange energy with demon[i]

**Reversibility Guarantee:**

- Fixed random patterns ensure deterministic dynamics
- Reversed order + reversed radii = exact time reversal
- Used to test fundamental reversibility of microcanonical dynamics

### 2. irrInferno Class (`inferno_irr.py`)

**Purpose**: Implements irreversible dynamics with true randomness.

**Key Difference:**

```python
class irrInferno:
    def __init__(self, N, R):
        # Same initialization as Inferno
        # BUT: No fixed radius arrays
        self.R = R  # Store radius for random generation

    def demon_move(self):
        # Generate NEW random radii each call
        r_spin = random.randint(0, R) * random.choice([-1, 1])
        r_bond = random.randint(0, R) * random.choice([-1, 1])
        # Use for current move, then discard
```

**Why This Matters:**

- Cannot reverse precisely - randomness breaks time-reversal symmetry
- Used to study entropy production in irreversible processes
- Compare with `Inferno` to quantify irreversibility

### 3. Simulation Runners

#### Production: `parallel_sim.py` and `parallel_sim_irr.py`

The production runners use multiprocessing and JIT compilation for maximum performance (~8,900x speedup).

**Key Features:**

- Parallel execution across multiple CPU cores
- JIT compilation enabled by default (use `--no-jit` to disable)
- Same algorithm and output format as legacy versions
- Command-line control: `--cores N`, `--n`, `--s`, `--r`, `--m`

#### Legacy: `legacy/sim.py` and `legacy/irr_sim.py`

Original single-core implementations preserved for educational purposes. See `creutz-sim/legacy/README.md` for details.

**Core Algorithm Structure (shared by all runners):**

```python
# Parameters
n = 10000      # Lattice size
s = 10000      # Sweeps per phase (forward/reverse)
r = 11         # Test radii R=1 to R=10
m = 5          # Number of independent runs

# Entropy formulas (from statistical mechanics)
Sk = lambda N, K: loggamma(K+N) - loggamma(K+1) - loggamma(N)
Su = lambda N, N0, Nx, N0_exp: loggamma(N+1) + log(2**N0_exp) -
                                (loggamma(N-N0-Nx+1) + loggamma(N0+1) + loggamma(Nx+1))

# Main loop
for M in range(m):                    # Multiple runs
    for R in range(r):                # Each radius
        x = Inferno(n, R+1)           # Create simulation

        # Forward phase
        for i in range(s):            # Each sweep
            for j in range(n):        # Each site
                x.demon_move()
                # Calculate entropy
                total_entropy = (Sk(n, sum(x.E_demon)) +
                                Su(n, x.bond_count[1], x.bond_count[2], N0_exp)) / n
            # Save data to CSV

        # Reverse phase
        for i in range(s):
            for j in range(n):
                x.demon_reverse()     # Run backward
                # Calculate entropy
            # Save data to CSV
```

**Data Output:**

- CSV per run: `data/r{R}/sim_data_r{R}_{M}.csv`
- Columns: `['t', 'K', 'U', 'N0', 'Nx', 'S/nk', 'n']`
  - `t`: Sweep number (0 to 2s-1, transition at s)
  - `K`: Average demon energy per site
  - `U`: Lattice energy per site
  - `N0`: Broken bonds per site
  - `Nx`: Anti-aligned spins per site
  - `S/nk`: Total entropy per site
  - `n`: Lattice size (for reference)

**Path Configuration:**

```python
# Use relative path from project root for portability
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
folder = os.path.join(project_root, 'data') + os.sep
```

## Data Flow

```
1. Initialization
   ├─> Create lattice (half +1, half -1 spins)
   ├─> Distribute energy 2*N between lattice and demon
   └─> Pre-compute random patterns (Inferno only)

2. Forward Phase (s sweeps)
   ├─> Each sweep: N demon moves
   ├─> Each move: attempt spin flip + bond change
   ├─> Update: lattice, bonds, demon energies
   ├─> Calculate: entropy (Sk + Su)
   └─> Write: average quantities to CSV

3. Reverse Phase (s sweeps)
   ├─> Each sweep: N demon_reverse moves
   ├─> Reversed order + operations
   ├─> Calculate: entropy
   └─> Write: continued CSV (t = s to 2s-1)

4. Analysis
   ├─> Multiple runs (M) per radius (R)
   ├─> Average across runs
   └─> Compare reversible vs irreversible
```

## Key Algorithms

### Entropy Calculation

**Kinetic Entropy (Sk):**

```python
# Multiplicity: ways to distribute K quanta among N oscillators
Ω_k = (K + N - 1)! / (K! * (N - 1)!)
S_k = log(Ω_k) = loggamma(K + N) - loggamma(K + 1) - loggamma(N)
```

**Configurational Entropy (Su):**

```python
# Multiplicity: ways to arrange N0 broken bonds and Nx anti-aligned pairs
Ω_u = N! * 2^N0 / ((N - N0 - Nx)! * N0! * Nx!)
S_u = log(Ω_u) = loggamma(N+1) + N0*log(2) - sum(loggamma(...))

# Special case: N0 = 0
# Implementation uses N0_exp = max(N0, 1), so when N0 = 0 we compute
# N0_exp * log(2) = 1 * log(2) = log(2^1) instead of log(2^0) = log(1) = 0.
```

### Energy Exchange Rules

**Spin Flip:**

```python
# Cost to flip spin[a]
cost = 2 * spin[a] * (neighbor_left + neighbor_right)

if cost < 0:
    # Energetically favorable - always flip
    spin[a] *= -1
    E_demon[i] -= cost  # Gain energy
    E_lattice += cost   # Lose energy
elif cost <= E_demon[i]:
    # Demon has enough energy
    spin[a] *= -1
    E_demon[i] -= cost  # Spend energy
    E_lattice += cost   # Increase energy
# else: no change
```

**Bond Change:**

```python
# Cost to break/make bond[a]
if spin[a] == spin[a+1]:
    cost = -1  # Aligned: costs -1 to break, +1 to make
else:
    cost = +1  # Anti-aligned: costs +1 to break, -1 to make

if bond[a] == 0 and E_demon[i] >= cost:
    # Remake broken bond
    bond[a] = cost
    E_demon[i] -= cost
    E_lattice += cost
elif bond[a] != 0 and E_demon[i] >= -cost:
    # Break existing bond
    bond[a] = 0
    E_demon[i] += cost
    E_lattice -= cost
```

## Performance Considerations

### Local Development (n=10000, s=10000)

- Single run: ~10-30 seconds
- Full suite (r=11, m=5): ~30 minutes
- Memory: ~100 MB

### HPC Production (n=1000000, s=10000)

- Single run: ~10-60 minutes
- Full suite: ~5-10 hours
- Memory: ~1 GB
- Benefits from cluster: parallel runs across radii

### Optimization Opportunities

1. **NumPy vectorization**: Currently uses Python loops
2. **Parallel runs**: Different radii/runs independent
3. **Checkpoint/restart**: Save state for long simulations
4. **Compiled extensions**: Cython for inner loops

## Testing & Validation

### Energy Conservation Test

```python
assert x.E_total == x.E_lattice + sum(x.E_demon)
assert x.E_total == 2 * N  # Should always be true
```

### Reversibility Test

```python
# Run forward then reverse
initial_state = copy(x.lattice), copy(x.E_demon)
for i in range(s):
    x.demon_move()
for i in range(s):
    x.demon_reverse()
final_state = x.lattice, x.E_demon

# For Inferno: final_state should equal initial_state
# For irrInferno: final_state will differ (irreversible)
```

### Entropy Sanity Checks

1. Entropy should be non-negative
2. For isolated system: entropy ≥ initial entropy
3. For reversible dynamics: entropy returns to initial value

## Visualization

### `sim_plot.py` - Single Run Analysis

- Time series of demon energy, lattice energy, entropy
- Red line marks forward/reverse transition (t = s)
- Used to debug individual runs

### `Sk_comparison.py` - Statistical Analysis

- Averages multiple runs per radius
- Rolling average smoothing (bin_size=10)
- Compares reversible vs. irreversible entropy evolution
- Color-coded by radius (11 colors for R=0-10)
- Zoomed view around transition point

## Extension Points

### Adding New Observables

1. Calculate in main loop (after `demon_move()`)
2. Add to `data` accumulator
3. Add column to CSV output
4. Update plotting scripts

### Different Lattice Geometries

- Modify neighbor calculation in `spin_flip()`
- Update bond structure
- Adjust periodic boundary conditions

### Alternative Dynamics

- Implement new `demon_move()` variant
- Keep same energy conservation rules
- Compare entropy production

### Parallel Processing

The project includes parallel processing capabilities with JIT compilation:

- **Implementation**: `parallel_sim.py` and `parallel_sim_irr.py`
- **Default mode**: Parallel JIT enabled (provides ~8,900x speedup)
- **Auto-detection**: Uses all available CPU cores by default
- **Manual control**: `--cores N` flag or `make run-sim ARGS="--cores N"`
- **SLURM integration**: Batch scripts configure core allocation
- **File management**: Each run generates independent CSV files
- **Aggregation**: Results combined for analysis across radii

See [JIT_BEST_PRACTICES.md](JIT_BEST_PRACTICES.md) and [OPTIMIZATIONS.md](OPTIMIZATIONS.md) for detailed performance information.

## References

### Primary Algorithm

**Creutz's Demon Algorithm:**

- Creutz, M. (1983). "Microcanonical Monte Carlo Simulation." _Physical Review Letters_, 50(19), 1411-1414. DOI: [10.1103/PhysRevLett.50.1411](https://doi.org/10.1103/PhysRevLett.50.1411)

  - Original paper introducing the demon algorithm for microcanonical simulation
  - Key innovation: Using an auxiliary "demon" particle to enable Monte Carlo moves while conserving total energy

- Creutz, M., & Freedman, B. (1981). "A statistical approach to quantum mechanics." _Annals of Physics_, 132(2), 427-462.
  - Earlier work developing the theoretical foundation

### Statistical Physics Background

**Ising Model:**

- Ising, E. (1925). "Beitrag zur Theorie des Ferromagnetismus." _Zeitschrift für Physik_, 31(1), 253-258.
  - Original formulation of the Ising model for magnetic systems

**Microcanonical Ensemble:**

- Huang, K. (2008). _Statistical Mechanics_ (2nd ed.). Wiley.
  - Comprehensive treatment of microcanonical ensembles (Chapter 7)

**Monte Carlo Methods:**

- Newman, M. E. J., & Barkema, G. T. (1999). _Monte Carlo Methods in Statistical Physics_. Oxford University Press.

  - Chapter 3: Microcanonical methods including Creutz demon
  - Practical implementation details

- Landau, D. P., & Binder, K. (2014). _A Guide to Monte Carlo Simulations in Statistical Physics_ (4th ed.). Cambridge University Press.
  - Section 2.4: Microcanonical Monte Carlo

### Entropy Calculations

**Statistical Mechanics Formulas:**

- Boltzmann entropy: S = k_B ln(Ω)
- Multiplicity calculations for microcanonical ensembles
- Configurational and kinetic entropy decomposition

### Implementation Notes

This codebase is an **original implementation** of the Creutz demon algorithm, developed to study:

- Comparison of reversible vs. irreversible dynamics
- Entropy production in microcanonical ensembles
- Time-reversal symmetry in Monte Carlo simulations

While the algorithm itself (Creutz, 1983) is standard, the specific research focus and all code are original contributions.
