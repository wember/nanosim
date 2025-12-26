# JIT Optimization - Production Best Practices

This guide covers best practices for using Numba JIT optimization in production nanosim workflows.

## Overview

Numba JIT compilation provides **70-106x speedup** by compiling Python functions to optimized machine code. Combined with parallel processing, this reduces production runs from **27 hours to 1.2 minutes** (~1400x total speedup).

## Quick Start

### Command-Line Usage

```bash
# Enable JIT for parallel runs
python creutz-sim/parallel_sim.py --jit        # Reversible (70x per core)
python creutz-sim/parallel_irr_sim.py --jit    # Irreversible (106x per core)

# Or use Makefile targets
make run-sim                                   # Production scale (parallel JIT default)
make run-irr-sim                               # Irreversible (parallel JIT default)
make run-sim-small                             # Quick test
```

### Python API

```python
from jit_inferno import JITInferno
from jit_irr_inferno import JITirrInferno

# Create JIT-optimized simulation
sim = JITInferno(n=1000000, s=10000, R=5)

# Run full sweeps (JIT does N moves per call)
for _ in range(s):
    sim.demon_move()
```

## When to Use JIT

### ✅ Use JIT For:

1. **Production runs** (n ≥ 10,000)

   - Maximum performance for large-scale simulations
   - Essential for thesis/publication work

2. **Parallel execution**

   - Multiplies per-core speedup across all CPUs
   - Best combined performance

3. **Multiple runs** (m ≥ 5)

   - JIT compilation overhead amortized across runs
   - First run pays ~1-2s compilation cost

4. **Time-critical analysis**
   - Rapid iteration during data analysis phase
   - Quick parameter exploration

### ❌ Don't Use JIT For:

1. **Small test runs** (n < 1000, single iteration)

   - Compilation overhead dominates runtime
   - Original version faster for tiny tests

2. **Development/debugging**

   - JIT errors harder to debug
   - Use original classes with validation mode

3. **Modifying core algorithms**
   - Original code easier to understand and modify
   - Switch to JIT after validating changes

## Performance Expectations

### Speedup by System Size

| Lattice Size (n) | Sweeps (s) | Original | JIT (Rev) | JIT (Irr) |
| ---------------- | ---------- | -------- | --------- | --------- |
| 100              | 10         | 0.01s    | ~0.01s    | ~0.01s    |
| 1,000            | 100        | 1.5s     | 0.05s     | 0.03s     |
| 10,000           | 1,000      | 150s     | 2.1s      | 1.4s      |
| 100,000          | 5,000      | 1.5hr    | 1.3min    | 0.8min    |
| 1,000,000        | 10,000     | 27hr     | 23min     | 15min     |

### Combined with Parallel Processing (16 cores)

| Configuration      | Speedup    | Production Time |
| ------------------ | ---------- | --------------- |
| Original           | 1x         | 27 hours        |
| JIT only           | 70-106x    | 15-23 minutes   |
| Parallel only      | 13-14x     | ~2 hours        |
| **JIT + Parallel** | **~1400x** | **1.2 minutes** |

## First-Run Behavior

### JIT Compilation Overhead

The first call to a JIT function triggers compilation:

```python
sim = JITInferno(n=10000, s=1000, R=5)

# First call: ~1-2 seconds (includes compilation)
sim.demon_move()

# Subsequent calls: ~0.02 seconds (compiled code)
for _ in range(999):
    sim.demon_move()
```

**Best Practice:** Ignore first-call timing when benchmarking.

### Warmup Pattern

For accurate timing:

```python
import time

sim = JITInferno(n=100000, s=5000, R=5)

# Warmup call (triggers compilation)
sim.demon_move()

# Now measure actual performance
start = time.perf_counter()
for _ in range(4999):
    sim.demon_move()
elapsed = time.perf_counter() - start
```

## Integration Patterns

### Pattern 1: Command-Line Flag (Recommended)

```python
# In your script
parser.add_argument('--jit', action='store_true', help='Enable JIT optimization')
args = parser.parse_args()

if args.jit:
    from jit_inferno import JITInferno as Inferno
else:
    from inferno import Inferno

# Use Inferno normally (works with both)
sim = Inferno(n=args.n, s=args.s, R=args.R)
```

### Pattern 2: Environment Variable

```python
import os

USE_JIT = os.getenv('NANOSIM_USE_JIT', '0') == '1'

if USE_JIT:
    from jit_inferno import JITInferno as Inferno
else:
    from inferno import Inferno
```

Usage:

```bash
NANOSIM_USE_JIT=1 python my_script.py
```

### Pattern 3: Configuration File

```python
import json

with open('config.json') as f:
    config = json.load(f)

if config.get('use_jit', False):
    from jit_inferno import JITInferno as Inferno
else:
    from inferno import Inferno
```

## Validation

### Verify JIT Correctness

Always validate JIT results match original implementation:

```python
from inferno import Inferno
from jit_inferno import JITInferno

n, s, R = 10000, 100, 5

# Run both versions
sim_orig = Inferno(n, s, R)
sim_jit = JITInferno(n, s, R)

for _ in range(s):
    sim_orig.demon_move()
    sim_jit.demon_move()

# Compare results
assert abs(sim_orig.K - sim_jit.K) < 1e-10, "Energy mismatch!"
assert abs(sim_orig.U - sim_jit.U) < 1e-10, "Lattice energy mismatch!"
print("✓ JIT results verified correct")
```

### Energy Conservation Check

```python
def verify_energy_conservation(sim, tolerance=1e-10):
    """Verify total energy conserved throughout simulation."""
    E_initial = 2.0  # Total energy always 2*N in reduced units
    E_final = sim.U + sim.K

    assert abs(E_final - E_initial) < tolerance, \
        f"Energy not conserved: {E_initial} → {E_final}"
    print(f"✓ Energy conserved: {E_final:.12f}")
```

## Common Issues

### Issue 1: ImportError for Numba

**Error:** `ModuleNotFoundError: No module named 'numba'`

**Solution:**

```bash
source venv/bin/activate
pip install -r requirements.txt  # Includes numba>=0.63.1
```

### Issue 2: JIT Compilation Warnings

**Warning:** `NumbaTypeSafetyWarning: unsafe cast from int64 to int32`

**Solution:** These are informational and don't affect correctness. Suppress with:

```python
import warnings
from numba.core.errors import NumbaTypeSafetyWarning
warnings.filterwarnings('ignore', category=NumbaTypeSafetyWarning)
```

### Issue 3: Slower Than Expected

**Symptoms:** JIT not showing expected speedup

**Debug checklist:**

1. ✓ Confirm using JIT classes (`JITInferno`, not `Inferno`)
2. ✓ Ignore first call (compilation overhead)
3. ✓ Use `demon_move()` not internal methods
4. ✓ Check system load (other processes competing)
5. ✓ Verify n ≥ 10,000 (overhead dominates for small n)

### Issue 4: Different Results

**Symptoms:** JIT produces different output than original

**Causes:**

- Random number generation order differs
- Floating-point rounding in different order

**Solution:** Results will differ slightly due to floating-point arithmetic order. Verify:

```python
# Energy should match within machine precision
assert abs(E_jit - E_original) < 1e-10

# Statistical properties should match (run multiple times, compare means)
```

## Migration Guide

### Converting Existing Code to JIT

**Before:**

```python
from inferno import Inferno

sim = Inferno(n=1000000, s=10000, R=5)
for _ in range(s):
    sim.demon_move()
```

**After:**

```python
from jit_inferno import JITInferno

sim = JITInferno(n=1000000, s=10000, R=5)
for _ in range(s):
    sim.demon_move()
```

**That's it!** API is 100% compatible.

### Updating Analysis Scripts

Add optional JIT support without breaking existing workflows:

```python
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--jit', action='store_true')
args = parser.parse_args()

# Import appropriate class
if args.jit:
    from jit_inferno import JITInferno as SimClass
else:
    from inferno import Inferno as SimClass

# Rest of code unchanged
sim = SimClass(n, s, R)
```

## Examples

### Quick Demo

```bash
# Run performance comparison
python examples/jit_quick_demo.py
```

Shows side-by-side comparison of original vs JIT.

### Production Template

```bash
# Validate setup
python examples/production_run_jit.py --validate-only

# Run production simulation
python examples/production_run_jit.py

# Custom parameters
python examples/production_run_jit.py --n 500000 --s 5000 --r 6 --m 10
```

Full-featured production template with error handling and monitoring.

## Testing

### Run JIT Tests

```bash
# IMPORTANT: Always use make run-tests for parallel execution (~3x faster)
make run-tests                           # Run all tests in parallel (recommended)

# Run specific test files
make run-test-file FILE=test_jit_implementation.py

# Run tests serially (for debugging race conditions)
make run-tests-serial

# Direct pytest commands (slower, not recommended for full suite)
pytest tests/test_jit_implementation.py -v              # JIT implementation tests
pytest tests/test_jit_integration.py -v                 # JIT integration tests
pytest tests/test_jit_implementation.py -v -m "not slow"  # Quick tests only
pytest tests/test_jit_implementation.py::TestJITPerformance -v  # Performance tests
```

### Benchmarking

```bash
# Compare JIT vs original performance
python benchmark_jit.py

# Profile JIT execution
python -m cProfile -o jit_profile.stats creutz-sim/parallel_sim.py --jit
python -m pstats jit_profile.stats
```

## Summary

**Key Takeaways:**

- ✅ Use `--jit` flag for all production runs (n ≥ 10,000)
- ✅ Combine with parallel processing for maximum speedup (~1400x)
- ✅ Validate results match original implementation
- ✅ Ignore first-call compilation overhead when timing
- ⚠️ Original version better for small tests and debugging
- ⚠️ JIT classes are drop-in replacements (100% API compatible)

**Production Workflow:**

```bash
# 1. Validate environment
make test-env

# 2. Quick test
make run-sim-small

# 3. Production run
make run-sim

# 4. Analyze results
python creutz-sim/Sk_comparison.py
```

For questions or issues, see [test_jit_implementation.py](../tests/test_jit_implementation.py) for comprehensive examples.
