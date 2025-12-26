# Examples Directory

This directory contains example scripts demonstrating various features and usage patterns of nanosim. These examples are designed to help you understand the codebase and serve as templates for your own simulations.

## Quick Start

Run all examples:

```bash
make run-examples
```

Or run individual examples:

```bash
python examples/quick_test.py
python examples/jit_quick_demo.py
```

---

## Example Scripts

### 1. Basic Usage

#### `quick_test.py` - Minimal Working Example

**Purpose:** Simplest possible demonstration of the simulation  
**Audience:** New users, quick verification  
**Runtime:** < 1 second

Demonstrates:

- Creating an `Inferno` lattice (N=100 sites)
- Running forward and reverse Monte Carlo sweeps
- Energy conservation verification
- Basic state inspection

**When to use:**

- First time running nanosim
- Verifying installation works
- Understanding basic API

```bash
python examples/quick_test.py
```

---

#### `custom_parameters.py` - Parameter Exploration

**Purpose:** Shows how different parameters affect simulation behavior  
**Audience:** Users designing their own experiments  
**Runtime:** < 5 seconds

Demonstrates:

- Testing different lattice sizes (N=50 to N=200)
- Varying demon-coupling radius (R=2 to R=10)
- Comparing reversible vs irreversible dynamics
- Entropy calculation for different configurations

**When to use:**

- Planning production runs
- Understanding parameter effects
- Choosing appropriate N and R values

```bash
python examples/custom_parameters.py
```

---

### 2. Performance & Optimization

#### `jit_quick_demo.py` - JIT Performance Comparison ⚡

**Purpose:** Demonstrates the dramatic speedup from Numba JIT compilation  
**Audience:** Users considering JIT optimization  
**Runtime:** ~30 seconds (includes both versions)

Demonstrates:

- Side-by-side comparison: original vs JIT
- Both reversible and irreversible simulations
- **70x speedup** for reversible (Inferno)
- **107x speedup** for irreversible (irrInferno)
- Energy conservation verification for both

**Key Results:**

```
Reversible:  3.23s → 0.046s  (70x faster)
Irreversible: 6.77s → 0.064s  (107x faster)
```

**When to use:**

- Deciding whether to use JIT
- Benchmarking your system
- Verifying JIT installation

```bash
python examples/jit_quick_demo.py
```

---

#### `benchmark_validation.py` - Validation Overhead Analysis

**Purpose:** Measures performance impact of different validation modes  
**Audience:** Performance-conscious users  
**Runtime:** ~30 seconds

Demonstrates:

- Three validation modes: `off`, `periodic`, `frequent`
- Performance across different lattice sizes
- Overhead quantification (1-5% for periodic, 10-20% for frequent)

**When to use:**

- Optimizing production runs
- Debugging vs performance tradeoff
- Understanding validation costs

```bash
python examples/benchmark_validation.py
```

---

#### `benchmark_neighbor_optimization.py` - Neighbor Array Performance

**Purpose:** Benchmarks pre-computed neighbor array optimization  
**Audience:** Developers, optimization enthusiasts  
**Runtime:** ~1 minute

Demonstrates:

- Performance comparison: modulo ops vs array lookups
- Speedup from neighbor pre-computation
- Effect across different lattice sizes

```bash
python examples/benchmark_neighbor_optimization.py
```

---

#### `benchmark_index_cycling.py` - Index Cycling Optimization

**Purpose:** Benchmarks index cycling vs modulo operation  
**Audience:** Developers interested in micro-optimizations  
**Runtime:** ~1 minute

Demonstrates:

- Conditional reset vs modulo performance
- Branch prediction benefits
- Effect on hot loop performance

```bash
python examples/benchmark_index_cycling.py
```

---

### 3. Production Workflows

#### `production_run_jit.py` - Production Template with JIT ⭐

**Purpose:** Best-practice template for production-scale simulations  
**Audience:** Users running thesis/publication simulations  
**Runtime:** Depends on parameters (default: ~2 minutes)

Features:

- **Automatic JIT optimization** (70-106x speedup)
- Parallel processing across all CPU cores
- Progress monitoring with time estimates
- Energy conservation validation
- CSV output compatible with analysis tools
- Comprehensive error handling
- Flexible command-line interface

**Usage:**

```bash
# Default production run (N=1M, s=10k, r=11, m=5)
python examples/production_run_jit.py

# Custom parameters - reversible
python examples/production_run_jit.py --reversible --n 100000 --s 5000 --r 8 --m 3

# Custom parameters - irreversible
python examples/production_run_jit.py --irreversible --n 100000 --s 5000 --r 8 --m 3

# Quick validation run (fast parameters)
python examples/production_run_jit.py --validation
```

**When to use:**

- Generating thesis/publication data
- Maximum performance production runs
- Systematic parameter studies

**Performance:**

- ~1400x faster than original implementation
- ~1.2 minutes for full production run (N=1M, s=10k, r=11, m=5)
- Scales linearly with CPU cores

---

#### `analysis_pipeline.py` - Complete Workflow Example

**Purpose:** End-to-end demonstration from simulation to analysis  
**Audience:** Users new to the analysis workflow  
**Runtime:** < 5 seconds

Demonstrates:

- Running a simulation
- Collecting data during run
- Saving results to CSV
- Loading and analyzing data
- Calculating physical quantities (entropy, temperature)
- Basic plotting suggestions

**When to use:**

- First time analyzing results
- Understanding data structure
- Building custom analysis scripts

```bash
python examples/analysis_pipeline.py
```

---

## Choosing the Right Example

### I want to...

**...understand the basics**
→ Start with `quick_test.py`, then `custom_parameters.py`

**...see if JIT will help me**
→ Run `jit_quick_demo.py` to see the 70-107x speedup

**...run production simulations**
→ Use `production_run_jit.py` as your template

**...analyze results**
→ Check out `analysis_pipeline.py` for the workflow

**...optimize performance**
→ Review benchmark scripts to understand optimization trade-offs

**...understand the code**
→ Read through examples in increasing complexity:

1. `quick_test.py` (simplest)
2. `custom_parameters.py` (parameter exploration)
3. `analysis_pipeline.py` (full workflow)
4. `production_run_jit.py` (production-ready)

---

## Performance Comparison

| Script                    | Lattice Size | Runtime | Purpose                     |
| ------------------------- | ------------ | ------- | --------------------------- |
| `quick_test.py`           | N=100        | < 1s    | Verify installation         |
| `custom_parameters.py`    | N=50-200     | < 5s    | Explore parameters          |
| `jit_quick_demo.py`       | N=10,000     | ~30s    | Benchmark JIT (both modes)  |
| `benchmark_validation.py` | N=1k-50k     | ~30s    | Measure validation overhead |
| `benchmark_neighbor_*.py` | Varies       | ~1m     | Measure optimizations       |
| `analysis_pipeline.py`    | N=100        | < 5s    | Learn analysis workflow     |
| `production_run_jit.py`   | N=1M         | ~1-2m   | Production simulations      |

---

## Common Patterns

### Pattern 1: Basic Simulation

```python
from inferno import Inferno

# Create lattice
sim = Inferno(N=1000, R=5)

# Run forward sweeps
for sweep in range(100):
    for _ in range(sim.N):
        sim.demon_move()

# Check energy conservation
state = sim.get_validated_state()
print(f"Total energy: {state['E_total']}")
```

### Pattern 2: JIT-Optimized Simulation

```python
from jit_inferno import JITInferno

# Create lattice (same API)
sim = JITInferno(N=10000, R=5)

# Run forward sweeps (one call per sweep!)
for sweep in range(1000):
    sim.demon_move()  # Does full sweep internally

# 70x faster than original!
```

### Pattern 3: Production with Parallel Processing

```bash
# Use the production template
python examples/production_run_jit.py \
    --reversible \
    --n 1000000 \
    --s 10000 \
    --r 11 \
    --m 5 \
    --cores 16

# Or use the parallel scripts directly
python creutz-sim/parallel_sim.py \
    --jit \
    --n 1000000 \
    --s 10000 \
    --r 11 \
    --m 5
```

---

## Tips for New Users

1. **Start small:** Run `quick_test.py` first to verify everything works

2. **Use JIT for production:** The `jit_quick_demo.py` shows why—70-107x faster!

3. **Test parameters:** Use `custom_parameters.py` to understand how N and R affect results

4. **Production template:** Copy `production_run_jit.py` as starting point for your own scripts

5. **Validation modes:** Use `--validate off` for production, `--validate frequent` for debugging

6. **Parallel processing:** Production scripts automatically use all CPU cores

---

## Further Reading

- **[docs/JIT_BEST_PRACTICES.md](../docs/JIT_BEST_PRACTICES.md)** - Comprehensive JIT usage guide
- **[docs/OPTIMIZATIONS.md](../docs/OPTIMIZATIONS.md)** - All performance optimizations explained
- **[docs/ARCHITECTURE.md](../docs/ARCHITECTURE.md)** - Physics and algorithm details
- **[README.md](../README.md)** - Main documentation

---

## Contributing Examples

If you develop a useful pattern or workflow, consider contributing it as an example! Good examples:

- Solve a common problem
- Are well-commented
- Run quickly (< 1 minute preferred)
- Include clear documentation
- Follow existing code style

See existing examples as templates.
