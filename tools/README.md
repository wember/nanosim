# Tools

Benchmarking and profiling utilities for Creutz demon simulations.

## Benchmark Scripts

### `benchmark_production.py`

Comprehensive production-scale benchmark comparing optimization strategies.

**Quick Test** (~1 minute):

```bash
make benchmark-production
# or
python tools/benchmark_production.py --config test
```

**Small Test** (~10-15 minutes):

```bash
make benchmark-production-small
# or
python tools/benchmark_production.py --config small
```

**Full Production** (hours - recommended overnight):

```bash
make benchmark-production-full
# or
python tools/benchmark_production.py --config production
```

**Options:**

- `--skip-legacy` - Skip legacy baseline (very slow)
- `--skip-parallel-only` - Skip parallel-without-JIT test
- `--cores N` - Specify number of cores for parallel tests
- `--output FILE` - Save results to custom JSON file

**What it measures:**

1. **Legacy Baseline** - Original single-core implementation (no JIT, no parallelization)
2. **JIT Only** - Single-core with JIT compilation
3. **Parallel Only** - Multi-core without JIT
4. **Production** - JIT + Parallel (what `make run-sim` uses)

**Output:**

- Terminal: Summary table with speedup factors
- `benchmark_results.json`: Detailed timing data, system info, speedup calculations

### `benchmark_jit.py`

Quick JIT compilation speedup test with small configurations.

```bash
make benchmark-jit
# or
python tools/benchmark_jit.py [n] [sweeps] [r]
```

**Default:** n=10,000, sweeps=100, r=5

Compares original vs JIT-compiled versions for both reversible and irreversible simulations.

## Profiling Scripts

### `profile_sim.py`

Profile simulation performance to identify bottlenecks.

```bash
make profile
# or
python tools/profile_sim.py
```

Generates `.stats` files in `profiling/` directory. View with:

```bash
make view-profile
```

## Usage Examples

### Compare All Optimizations

```bash
# Quick comparison (recommended for development)
python tools/benchmark_production.py --config test

# Skip slow tests
python tools/benchmark_production.py --config test --skip-legacy
```

### Measure JIT Speedup

```bash
# Default configuration
make benchmark-jit

# Custom size
python tools/benchmark_jit.py 50000 200 10
```

### Profile for Optimization

```bash
# Profile and view results
make profile
make view-profile
```

### Production Benchmarking

```bash
# Run full benchmark suite (takes hours)
python tools/benchmark_production.py --config production --output prod_benchmark_$(date +%Y%m%d).json
```

## Interpreting Results

### Speedup Factors

**Measured performance gains (Apple M3 Max, 16 cores):**

Benchmark workload: n=10,000 sites, s=100 sweeps, r=3 radii, m=2 runs (12 total simulations)

- JIT alone: ~3,881x (single-core)
- Parallel alone: ~1,147x (16 cores, without JIT)
- Combined: ~8,902x (JIT + 16 cores)

**System dependence:**

- CPU architecture affects JIT speedup (AVX, vectorization)
- Core count directly impacts parallel speedup
- Memory bandwidth affects large lattices (n > 1M)

### Benchmark Configurations

| Config     | n    | s   | r   | m   | Time   | Use Case                    |
| ---------- | ---- | --- | --- | --- | ------ | --------------------------- |
| test       | 10k  | 100 | 3   | 2   | ~1 min | Quick validation            |
| small      | 100k | 1k  | 5   | 3   | ~15min | Medium-scale testing        |
| production | 1M   | 5k  | 11  | 5   | hours  | Full production measurement |

### Understanding Speedup

Speedup is relative to baseline (legacy single-core):

```
Speedup = Time_baseline / Time_optimized
```

- **< 1x**: Slower than baseline (overhead dominates)
- **1-10x**: Modest improvement (small configs, overhead limits gains)
- **10-100x**: Good parallelization or JIT effectiveness
- **100-1000x+**: Combined optimization success (production scale)

### Key Metrics

**Time per simulation:** Total time / number of simulations

- Useful for comparing efficiency across different parameter sets
- Should decrease with optimization strategies

**Parallel efficiency:** Actual speedup / (# cores)

- 100% = perfect scaling
- 70-90% = good scaling (typical for Monte Carlo)
- < 50% = bottleneck (I/O, memory, or sequential portions)

## Tips

1. **First run includes JIT compilation time** (~15-20s). Use `make compile` to pre-compile or ignore first-run timing.

2. **Memory matters** for large lattices. Monitor with:

   ```bash
   python tools/benchmark_production.py --config production &
   watch -n 1 "ps aux | grep python"
   ```

3. **Save results** for comparison:

   ```bash
   python tools/benchmark_production.py --config test --output baseline_$(date +%Y%m%d).json
   ```

4. **Background runs** for long benchmarks:

   ```bash
   nohup python tools/benchmark_production.py --config production > benchmark.log 2>&1 &
   tail -f benchmark.log
   ```

5. **Compare branches** or hardware:

   ```bash
   # On machine A
   python tools/benchmark_production.py --config small --output machine_a.json

   # On machine B
   python tools/benchmark_production.py --config small --output machine_b.json

   # Compare results
   diff <(jq .speedups machine_a.json) <(jq .speedups machine_b.json)
   ```
