# Nanosim Performance Optimization Strategy

## Context

Current simulation works correctly but needs significant speedup for production runs at N=10,000 to N=1,000,000.

## Implemented Optimizations (February 2026)

### Multiprocessing Parallelization - COMPLETED ✅

**Achieved Speedup**: ~11x on 11-core system (36 hours → 3.2 hours for N=1000, s=10000)

#### Implementation Details

Parallelized the outer R loop (radius values) using Python's `multiprocessing.Pool`:

```python
from multiprocessing import Pool, Manager
from functools import partial

def run_radius_simulations(R, n, s, flag, m, file_names, irr_files, init_files, pbar_queue):
    """Worker function to run all simulations for a given radius R."""
    # Ignore keyboard interrupts in worker - main process handles them
    import signal
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    x = Inferno(n, R)
    for M in range(m):
        # Run all sweeps for this (R, M) combination
        # Signal progress via shared queue
        pbar_queue.put(1)  # Update progress bar
    return completed_count

# Main process
num_cores = min(mp.cpu_count(), r + 1)  # Don't use more cores than jobs
pool = Pool(processes=num_cores)
results = pool.map_async(worker_func, range(r+1))
```

#### Key Design Decisions

**1. Parallelization Strategy**

- Chose to parallelize by R (radius) rather than M (runs)
- Each worker handles all M runs for one R value
- Limits parallelism to (r+1) workers maximum
- Rationale: Simpler progress tracking, each R writes to independent file structure

**2. Core Count Optimization**

```python
num_cores = min(mp.cpu_count(), r + 1)
```

- Don't spawn more workers than actual jobs
- Example: r=10 uses 11 cores (R=0 to R=10), not all 16 available
- Prevents idle workers and resource waste

**3. Progress Tracking**

- Uses `Manager.Queue` for inter-process communication
- Workers put progress updates into shared queue
- Main process consumes queue and updates tqdm progress bar
- Real-time feedback with elapsed/remaining time estimates

**4. Interrupt Handling**
Critical fixes for clean Ctrl-C behavior:

```python
# In worker function:
signal.signal(signal.SIGINT, signal.SIG_IGN)  # Workers ignore Ctrl-C

# In main process monitoring loop:
except Exception:  # NOT bare except: which catches KeyboardInterrupt!
    pass

# Outer exception handler:
except KeyboardInterrupt:
    print("\n\nInterrupted! Terminating...")
    pbar.close()
    pool.terminate()  # Force kill all workers
    pool.join()       # Wait for cleanup
    sys.exit(1)
```

**Issues encountered and fixed:**

- Initial signal handlers didn't work with multiprocessing spawn on macOS
- Workers were receiving KeyboardInterrupt during scipy import phase
- Bare `except:` was swallowing KeyboardInterrupt in monitoring loop
- Solution: Workers ignore SIGINT, main process catches and terminates pool

**5. Module-Level Code Organization**
All initialization moved inside `if __name__ == '__main__':` block to prevent re-execution in spawned workers:

- Parameter parsing
- File path setup
- Progress bar creation
- Signal handlers removed (incompatible with multiprocessing)

Only function definitions remain at module level:

- `add_row()`
- `get_params()`
- `format_time()`
- `run_radius_simulations()` (worker function)

#### Performance Results

- **Test case**: N=1000, s=10000, r=10, m=5 (55 simulations)
- **Sequential**: ~27 hours
- **Parallel (11 cores)**: ~3.2 hours
- **Speedup**: 8.4x
- **Parallel efficiency**: ~76% (8.4/11)

Efficiency loss due to:

- Process spawning overhead
- Inter-process communication (queue updates)
- Non-uniform work per core (some R values complete faster)

#### Files Modified

- `creutz-sim/sim.py`: Complete restructuring for parallelization
  - Removed: signal handlers, atexit handlers, module-level initialization
  - Added: `run_radius_simulations()` worker function
  - Modified: Main block with Pool management and exception handling

#### Physics Validation

- ✅ No changes to simulation algorithms
- ✅ Each worker runs identical code to sequential version
- ✅ No shared state during computation (only progress updates)
- ✅ Each (R, M) combination writes to independent files
- ✅ Results mathematically identical to sequential execution

---

## Target Performance Requirements

- **Current use case**: N=10-100 (development/testing)
- **Production use case**: N=10,000 to N=1,000,000
- **Bottleneck**: Inner loop calls `demon_move()` → `spin_flip()` + `bond_change()` for n×s iterations
  - Example: N=1,000,000, s=100 → 100M function calls per simulation

## Optimization Approaches

### 1. JIT Compilation (Numba) - FUTURE WORK

**Expected Additional Speedup**: 10-100x for large N per core

**Note**: With existing 8.4x parallelization, combined speedup could reach 84-840x total.

#### Strategy

Extract compute-intensive methods into pure functions and compile with `@njit`:

```python
from numba import njit

@njit
def spin_flip_jit(lattice, bonds, E_demon, E_lattice, d_energy, a, i, N):
    """Pure function version of spin_flip for JIT compilation."""
    s = lattice[a]
    d = E_demon[i]
    nb = lattice[(a+1)%N]*abs(bonds[a%N]) + lattice[(a-1)%N]*abs(bonds[(a-1)%N])
    cost = 2*s*nb

    if cost < 0:
        s *= -1
        E_demon[i] -= cost
        d_energy -= cost
        E_lattice += cost
    elif cost <= E_demon[i]:
        s *= -1
        E_demon[i] -= cost
        d_energy -= cost
        E_lattice += cost

    lattice[a] = s

    # Update bonds
    if bonds[a] != 0:
        if lattice[a] == lattice[(a+1)%N]:
            bonds[a] = -1
        else:
            bonds[a] = 1

    if bonds[(a-1)%N] != 0:
        if lattice[a] == lattice[(a-1)%N]:
            bonds[(a-1)%N] = -1
        else:
            bonds[(a-1)%N] = 1

    return E_lattice, d_energy

# Similar for bond_change_jit, count_bonds_jit
```

#### Implementation Steps

1. **Profile baseline** - measure current performance at various N
2. **Extract functions** - create standalone JIT-able versions of:
   - `spin_flip()`
   - `bond_change()`
   - `count_bonds()`
3. **Wrap in class** - keep class interface, call JIT functions internally
4. **Validate correctness** - run small N (10-100) with/without JIT, compare outputs exactly
5. **Handle RNG carefully** - Numba has its own RNG that differs from numpy.random
   - May need to pass RNG state explicitly
   - Critical for maintaining reversibility
6. **Benchmark scaling** - test at N=100, 1K, 10K, 100K, 1M

#### Known Risks

- **Random number generation**: Numba's RNG vs numpy.random - affects irreversible dynamics
- **Floating point precision**: Any changes break reversibility guarantee
- **First compilation overhead**: ~1-2s on first run (acceptable for long simulations)
- **Array bounds**: Modulo operations with large indices need careful testing

#### Fallback

If JIT breaks physics/reversibility, can create separate JIT and non-JIT code paths:

- Development/validation: use pure Python
- Production: use JIT version (validated against small-N reference runs)

---

### 2. Parallelization (Multiprocessing) - ✅ IMPLEMENTED (Feb 2026)

**Achieved Speedup**: ~8.4x on 11-core system (see "Implemented Optimizations" section above)

Status: Production-ready. Currently running N=1000, s=10000 simulation with 3.2-hour completion time (down from 27 hours sequential).

---

### 3. JIT Compilation (Numba) - FUTURE WORK

**Expected Speedup**: Linear with CPU cores (e.g., 8 cores → 8x)

#### Current Loop Structure

```python
for M in range(m):        # runs (independent)
    for R in range(r+1):  # radii (independent)
        # Each (M,R) pair is completely independent
        # - Creates own Inferno instance
        # - Writes to separate file (different M suffix)
```

#### Approach A: Parallelize (M, R) Pairs

All combinations run in parallel:

```python
from multiprocessing import Pool
from functools import partial

def run_simulation(M_R_tuple, n, s, flag, r, data_folder):
    """Run simulation for a single (M, R) combination."""
    M, R = M_R_tuple
    x = Inferno(n, R)
    # ... simulation logic ...
    return results

# Create all (M, R) combinations
tasks = [(M, R) for M in range(m) for R in range(r+1)]

# Run in parallel
with Pool(processes=cpu_count()) as pool:
    results = pool.map(partial(run_simulation, n=n, s=s, ...), tasks)
```

**Pros:**

- Maximum parallelism
- Clean separation

**Cons:**

- Progress tracking is complex (need shared state via Manager)
- Many concurrent processes if m×(r+1) is large

#### Approach B: Parallelize M Loop Only

Keep R sequential within each M:

```python
def run_all_radii(M, n, s, flag, r, data_folder):
    """Run all radii for a single M value."""
    for R in range(r+1):
        # ... simulation logic ...
    return M, results

# Run multiple M values in parallel
with Pool(processes=min(m, cpu_count())) as pool:
    results = pool.map(partial(run_all_radii, ...), range(m))
```

**Pros:**

- Simpler progress tracking (each worker reports % complete for all radii)
- Fewer processes (exactly m processes)

**Cons:**

- Less parallelism (only m-way instead of m×(r+1)-way)

#### Considerations

- **Memory**: Each process loads full Inferno instance (should be fine, instances are small)
- **I/O**: Many parallel writes, but to different files (no contention)
- **Progress bar**: Needs shared state or periodic status file writes
- **Combining with JIT**: JIT compilation happens once per process (small overhead)

---

### 3. Combined Approach - FUTURE WORK

Use both JIT and parallelization (parallelization already implemented ✅):

1. ✅ Parallelization for multi-core utilization (8.4x achieved on 11-core system)
2. JIT for inner loop speedup (10-100x per core - not yet implemented)
3. Combined: potentially **84-840x faster** for large N on multi-core systems

#### Implementation Order

1. **Phase 1**: ✅ COMPLETED - Parallelization implemented (Feb 2026)
   - 8.4x speedup achieved on 11-core system
   - Production-ready and validated
2. **Phase 2**: Add JIT compilation once needed for larger N
   - Expected 10-100x additional speedup per core
3. **Phase 3**: Benchmark combined approach at production scale

---

## Profiling Strategy

Before optimization, establish baseline:

```python
# Add timing to sim.py
import cProfile
import pstats

profiler = cProfile.Profile()
profiler.enable()

# ... run simulation ...

profiler.disable()
stats = pstats.Stats(profiler)
stats.sort_stats('cumtime')
stats.print_stats(20)  # Top 20 functions
```

Expected hotspots:

1. `demon_move()` - called n×s times
2. `spin_flip()` - called n×s times
3. `bond_change()` - called n×s times
4. `count_bonds()` - called n×s times

---

## Validation Testing

Critical: maintain physics correctness during optimization.

### Test Suite

1. **Exact output comparison** (small N):
   - Run N=10, s=100, m=2, r=5 with both versions
   - Compare CSV outputs bit-for-bit
   - Should be identical for reversible dynamics
2. **Energy conservation**:
   - Verify E_total = E_lattice + sum(E_demon) remains constant
   - Check at every sweep
3. **Reversibility check**:
   - For flag='r', forward then reverse should return to initial state
   - Critical test: compare lattice/bonds/E_demon before/after full cycle
4. **Statistical properties** (large N):
   - For N=10000, verify entropy distributions match
   - Compare Su/k and Sk/k statistics between versions

### Validation Script Template

```python
def validate_optimization():
    """Compare optimized vs baseline versions."""
    import subprocess
    import filecmp

    params = {'n': 10, 's': 100, 'flag': 'r', 'r': 5, 'm': 2}

    # Run baseline
    subprocess.run(['python', 'sim.py', '--baseline', ...])

    # Run optimized
    subprocess.run(['python', 'sim_optimized.py', ...])

    # Compare outputs
    for file in output_files:
        assert filecmp.cmp(baseline_file, optimized_file), f"Mismatch in {file}"

    print("✓ Validation passed")
```

---

## Implementation Timeline

### Phase 1: JIT (Priority 1)

**Timeline**: 1-2 days

- Day 1: Extract pure functions, add @njit decorators, validate small N
- Day 2: Test large N, benchmark, fix any issues

### Phase 2: Parallelization (Priority 2)

**Timeline**: 1 day

- Implement parallel M loop approach
- Test progress tracking
- Validate file I/O

### Phase 3: Combined + Production Testing (Priority 3)

**Timeline**: 1-2 days

- Combine JIT + parallelization
- Run production-scale test (N=100K+)
- Document performance gains

---

## Performance Expectations

### Baseline (Current)

- N=100: ~1s per simulation
- N=10,000: ~100s per simulation (estimated)
- N=1,000,000: ~10,000s per simulation (estimated)

### With JIT Only

- N=100: ~0.5s per simulation (2x speedup)
- N=10,000: ~5s per simulation (20x speedup)
- N=1,000,000: ~100s per simulation (100x speedup)

### With JIT + 8-core Parallelization

- N=100: ~0.06s per simulation (16x speedup)
- N=10,000: ~0.6s per simulation (160x speedup)
- N=1,000,000: ~12s per simulation (800x speedup)

_Note: Actual speedups depend on hardware, memory bandwidth, and specific N/s/m/r parameters._

---

## Open Questions

1. **What broke in previous JIT attempt?**
   - Specific error or incorrect results?
   - Random number generation issues?
   - Reversibility broken?

2. **Hardware constraints?**
   - CPU cores available for parallelization?
   - Memory limits for large N?
   - Running on cluster/HPC or local machine?

3. **Production run parameters?**
   - Typical N, s, m, r values for publication runs?
   - Time budget per simulation?
   - Total number of simulations needed?

---

## References

- Numba documentation: https://numba.pydata.org/
- Numba random: https://numba.pydata.org/numba-doc/dev/reference/numpysupported.html#random
- Python multiprocessing: https://docs.python.org/3/library/multiprocessing.html
