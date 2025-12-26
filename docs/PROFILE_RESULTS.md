# Profiling Results Summary

**Date:** December 25, 2025  
**Configuration:** N=10,000, sweeps=100, R=5 (2,000,000 total moves)

---

## Executive Summary

### Key Findings

1. **Inferno is 3.3x faster than irrInferno** (8.1s vs 26.9s)
2. **irrInferno's bottleneck:** Random number generation via `np.random.choice()` - takes **28% of total runtime**
3. **Core physics operations dominate:**
   - `bond_change()`: 31-32% of runtime
   - `spin_flip()`: 28-30% of runtime
   - `update_bonds_incremental()`: 16-17% of runtime
4. **Validation overhead negligible** with `validate_mode='off'`

---

## Detailed Breakdown

### Inferno Performance (8.151 seconds)

**Time distribution:**

| Function                        | Time (s) | % Total | Calls     | Purpose                     |
| ------------------------------- | -------- | ------- | --------- | --------------------------- |
| `bond_change()`                 | 2.569    | 31.5%   | 2,000,000 | Bond creation/destruction   |
| `spin_flip()`                   | 2.250    | 27.6%   | 2,000,000 | Spin flip attempts          |
| `update_bonds_incremental()`    | 1.332    | 16.3%   | 1,427,476 | Update bonds after flip     |
| `demon_move()` (overhead)       | 0.767    | 9.4%    | 1,000,000 | Loop overhead + indexing    |
| `demon_reverse()` (overhead)    | 0.760    | 9.3%    | 1,000,000 | Loop overhead + indexing    |
| `abs()` calls                   | 0.189    | 2.3%    | 4,000,000 | Bond sign checks            |
| `perform_periodic_validation()` | 0.083    | 1.0%    | 2,000,000 | Validation check (disabled) |
| Other                           | 0.201    | 2.5%    | -         | Misc operations             |

**Key observations:**

- Core physics (bond_change + spin_flip + update_bonds) = **75.4%** of runtime
- Loop overhead (demon_move + demon_reverse) = **18.7%**
- Validation negligible when disabled ✅
- Integer arithmetic is very efficient

---

### irrInferno Performance (26.869 seconds)

**Time distribution:**

| Function                          | Time (s) | % Total  | Calls         | Purpose                     |
| --------------------------------- | -------- | -------- | ------------- | --------------------------- |
| `demon_move()` (total)            | 5.929    | 22.1%    | 1,000,000     | Main loop + RNG overhead    |
| `demon_reverse()` (total)         | 5.916    | 22.0%    | 1,000,000     | Reverse loop + RNG overhead |
| **`np.random.choice()` overhead** | **~7.6** | **~28%** | **4,000,001** | **Random ±1 selection**     |
| `bond_change()`                   | 2.771    | 10.3%    | 2,000,000     | Bond operations             |
| `spin_flip()`                     | 2.431    | 9.0%     | 2,000,000     | Spin flip attempts          |
| `update_bonds_incremental()`      | 1.422    | 5.3%     | 1,452,175     | Bond updates                |
| Other NumPy overhead              | 2.8-3.0  | 10-11%   | -             | Array operations            |

**Key observations:**

- **Random number generation dominates**: ~28% of total runtime!
- `np.random.choice([-1, 1])` called 4 million times for generating ±1
- This is the primary difference between Inferno (3.3x slower)
- Core physics still efficient, but buried under RNG overhead

---

## Critical Bottleneck Identified: Random Number Generation

### The Problem

In `irrInferno.demon_move()` and `demon_reverse()`:

```python
# Called 2 million times each (4 million total)
radius_spin = np.random.randint(0, self.R) * np.random.choice([-1, 1])
radius_bond = np.random.randint(0, self.R) * np.random.choice([-1, 1])
```

**Analysis:**

- `np.random.choice([-1, 1])` is **extremely expensive** for generating single values
- Creates array wrapper, validates input, performs generic selection
- Called 4 million times → 7.6 seconds of overhead (28% of runtime!)
- This is **unnecessary complexity** for a simple ±1 coin flip

### Optimization Opportunity #1: Fast Random Sign

**Replace this:**

```python
np.random.choice([-1, 1])  # ~2µs per call
```

**With this:**

```python
2 * np.random.randint(0, 2) - 1  # ~0.5µs per call (4x faster)
# Or even better:
np.random.randint(0, 2) * 2 - 1  # Generates 0→-1 or 1→1
```

**Expected impact:**

- Save ~5-6 seconds per 2M moves
- **irrInferno speedup: ~1.3-1.5x** (26.9s → 18-20s)
- **Quick win:** Single line change, no logic changes

---

## Secondary Optimization Opportunities

### 2. Numba JIT Compilation (Highest Potential Impact)

**Why Numba is perfect for this code:**

✅ Pure NumPy array operations  
✅ Integer-only arithmetic (no Python objects)  
✅ Tight inner loops called millions of times  
✅ No dynamic types or complex Python features

**Target functions for @jit decoration:**

- `spin_flip()` - 2.25s → 0.1-0.2s (10-20x faster)
- `bond_change()` - 2.57s → 0.2-0.3s (8-10x faster)
- `update_bonds_incremental()` - 1.33s → 0.1-0.15s (8-10x faster)

**Expected combined impact:**

- Inferno: 8.1s → **2-3s** (3-4x speedup)
- irrInferno: 26.9s → **8-12s** (2-3x speedup)

**Implementation complexity:** Low - add decorators and type hints

### 3. Parallel Processing (Multiple Runs)

**Current simulation pattern:**

```python
for radius in range(1, 11):  # R=1 to R=10
    for run in range(5):     # 5 runs per radius
        run_simulation()      # Sequential!
```

**Optimization:**

- Run different (R, run) combinations in parallel
- 50 total simulations × 8 cores = **~6-8x speedup**

**Expected impact:**

- Full simulation suite: hours → minutes
- Simple to implement with `multiprocessing.Pool`

---

## Recommended Optimization Priority

### Phase 1: Quick Wins (1-2 hours work)

**1a. Fix random sign generation in irrInferno**

- Impact: 1.3-1.5x speedup on irrInferno
- Complexity: Trivial (1 line change)
- Risk: Zero

**1b. Add parallel processing for sim.py/irr_sim.py**

- Impact: 6-8x speedup for full simulation runs
- Complexity: Low (add multiprocessing)
- Risk: Low

**Expected result:** Full simulation suite runs 8-10x faster

### Phase 2: Numba JIT (4-6 hours work)

**2. Add @jit decorators to hot functions**

- Impact: 3-4x speedup on individual simulations
- Complexity: Medium (testing, type hints)
- Risk: Medium (must verify physics correctness)

**Expected result:** Individual simulations 3-4x faster

### Phase 3: Advanced (Optional)

**3. Vectorization**

- Impact: 2-5x additional speedup
- Complexity: High (algorithm restructuring)
- Risk: High (correctness verification)

---

## Profiling Commands Reference

```bash
# Profile with custom parameters
python profile_sim.py --mode inferno --n 50000 --s 50
python profile_sim.py --both --n 10000 --s 100

# Interactive analysis
python -m pstats profile_inferno_n10000_s100.stats
>>> sort cumulative
>>> stats 20
>>> sort tottime
>>> stats 20

# Generate call graph (requires gprof2dot + graphviz)
gprof2dot -f pstats profile_inferno_n10000_s100.stats | dot -Tpng -o profile.png
```

---

## Conclusions

1. **Code is already well-optimized** - no obvious algorithmic inefficiencies
2. **Random number generation** is the #1 bottleneck in irrInferno (easy fix!)
3. **Numba JIT** would provide dramatic speedups with minimal code changes
4. **Parallel processing** is essential for running multiple parameter sweeps
5. **Validation overhead is negligible** - optimization worked ✅

**Next action:** Implement Phase 1 optimizations (random sign + parallel processing)
