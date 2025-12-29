# Legacy Single-Core Simulation Runners

This directory contains the original single-core implementations of the Creutz demon simulation runners. These files are preserved for educational and reference purposes but are **not used in production**.

## Files

- **`sim.py`** - Original single-core reversible simulation runner
- **`irr_sim.py`** - Original single-core irreversible simulation runner
- **`inferno.py`** - Original `Inferno` class implementation (pre-optimization)

## Why These Are Legacy

These files represent the **original implementation from December 24, 2025**, before the optimization work that introduced:

- Parallel execution with multiprocessing
- JIT compilation with Numba
- Pre-allocated arrays and batch CSV writing
- Enhanced validation and error handling
- Refactored class structure

The legacy versions are **single-core, interpreted Python** with simpler implementations. **These files have been restored to their exact Dec 24 state** to preserve the original codebase before the optimization sprint.

They were superseded by the parallel JIT-compiled versions:

- `parallel_sim.py` - Production reversible runner (uses multiprocessing + JIT)
- `parallel_sim_irr.py` - Production irreversible runner (uses multiprocessing + JIT)

**Performance comparison:**

- Legacy (single-core): 2.48 hours for benchmark run (n=10k, s=100, r=3, m=2)
- Production (parallel JIT): 1.0 second for same run (~8,900x faster)

## When to Use Legacy Versions

These legacy scripts may be useful for:

1. **Educational purposes** - Simpler code without parallel/JIT complexity
2. **Understanding the algorithm** - Easier to read and follow
3. **Debugging** - Isolate issues without multiprocessing complexity
4. **Historical reference** - See how the code evolved

## Using Legacy Scripts

Both scripts support command-line arguments:

```bash
# Reversible simulation
python creutz-sim/legacy/sim.py --n 1000 --s 100 --r 3 --m 2

# Irreversible simulation
python creutz-sim/legacy/irr_sim.py --n 1000 --s 100 --r 3 --m 2
```

**Arguments:**

- `--n`: Lattice size (default: 1000000)
- `--s`: Sweeps per phase (default: 5000)
- `--r`: Max demon-coupling radius (default: 11, tests R=1 to R=10)
- `--m`: Number of independent runs (default: 5)
- `--validate`: Validation mode - `off`, `periodic`, or `frequent` (default: off)

## Migration to Production

For production use, always use the parallel JIT versions:

```bash
# Use Makefile targets (recommended)
make run-sim           # Reversible (parallel JIT)
make run-sim-irr       # Irreversible (parallel JIT)
make run-sim-small     # Quick test

# Or directly with Python
python creutz-sim/parallel_sim.py
python creutz-sim/parallel_sim_irr.py

# Single-core mode (if needed)
python creutz-sim/parallel_sim.py --cores 1
```

## Technical Notes

- Both legacy scripts import from the same core classes (`Inferno`, `irrInferno`)
- Output format is identical to production versions (same CSV columns)
- Path configuration uses project root for portability (same as production)
- All unit tests validate both legacy and production code paths

## History

- **Pre-Dec 25, 2025**: Original single-core implementation
- **Dec 25, 2025**: Parallel processing + JIT compilation added
- **Dec 26, 2025**: Moved to legacy/ directory for preservation
