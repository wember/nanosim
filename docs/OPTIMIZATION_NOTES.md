# Nanosim Performance Optimization Strategy

## Context
Current simulation works correctly but needs significant speedup for production runs at N=10,000 to N=1,000,000.

## Target Performance Requirements
- **Current use case**: N=10-100 (development/testing)
- **Production use case**: N=10,000 to N=1,000,000
- **Bottleneck**: Inner loop calls `demon_move()` → `spin_flip()` + `bond_change()` for n×s iterations
  - Example: N=1,000,000, s=100 → 100M function calls per simulation

## Optimization Approaches

### 1. JIT Compilation (Numba) - PRIMARY RECOMMENDATION
**Expected Speedup**: 10-100x for large N

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

### 2. Parallelization (Multiprocessing) - SECONDARY RECOMMENDATION
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

### 3. Combined Approach - MAXIMUM PERFORMANCE
Use both JIT and parallelization:
1. JIT for inner loop speedup (10-100x per simulation)
2. Parallelization for multi-core utilization (8x on 8-core machine)
3. Combined: potentially **80-800x faster** for large N on multi-core systems

#### Implementation Order
1. **Phase 1**: Implement JIT first, validate thoroughly
2. **Phase 2**: Add parallelization once JIT is stable
3. **Phase 3**: Benchmark combined approach

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

*Note: Actual speedups depend on hardware, memory bandwidth, and specific N/s/m/r parameters.*

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
