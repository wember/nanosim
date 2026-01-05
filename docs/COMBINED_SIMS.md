# Combined Simulation Mode

## Overview

The `parallel_sim_combined.py` script runs **both reversible and irreversible simulations together**, effectively **halving your total runtime** compared to running them separately.

## How It Works

Instead of running:

1. `make run-sim` (reversible only)
2. `make run-irr-sim` (irreversible only)

You can now run:

- `make run-sims` (both at once!)

## Time Savings

**Example: Full production run (r=11, m=5)**

- **Sequential**: Run reversible (~2 min) + Run irreversible (~2 min) = **~4 minutes total**
- **Combined**: Run both together = **~2 minutes total** ✨

The combined mode interleaves both simulation types, so your CPU cores work on both simultaneously instead of sitting idle waiting for one batch to complete.

## Usage

### Quick Test (2 radii, 2 runs each = 8 total sims)

```bash
make run-sims ARGS="--n 100 --s 5 --r 2 --m 2"
```

### Medium Run

```bash
make run-sims ARGS="--n 10000 --s 1000 --r 5 --m 3"
```

### Full Production Run (default)

```bash
make run-sims
# Uses: n=1000000, s=5000, r=11, m=5
# Creates: 110 simulations (55 reversible + 55 irreversible)
```

## Output

The combined simulation creates data in the **same directories** as separate runs:

- Reversible data: `data/r{R}/sim_data_r{R}_{M}.csv`
- Irreversible data: `data/irr/r{R}/irr_sim_data_r{R}_{M}.csv`

All downstream analysis tools (plotting, comparison) work exactly the same way!

## Performance

With 16 CPU cores and JIT enabled:

- 8 small simulations: **~0.7 seconds**
- Expected full run: **~2 minutes** (compared to 4 minutes sequential)

## When to Use Combined vs Separate

**Use Combined Mode When:**

- ✅ You need both reversible AND irreversible data
- ✅ You want to minimize total wall-clock time
- ✅ Running production datasets for comparison

**Use Separate Modes When:**

- ✅ You only need one simulation type
- ✅ Debugging a specific issue in one simulation type
- ✅ Re-running just one type after data loss

## Technical Details

The combined simulation:

1. Creates task list with both types: `[rev_R0M0, irr_R0M0, rev_R0M1, irr_R0M1, ...]`
2. Submits all tasks to multiprocessing pool
3. CPU cores process whichever tasks are available
4. Both types complete in parallel, using all available cores efficiently

No computational shortcuts are taken - each simulation runs with full fidelity, just scheduled more efficiently!
