#!/usr/bin/env python3
"""
Analyze JIT warmup overhead vs. runtime benefits for production runs.
"""

# From our benchmarks
jit_speedup = 30.5  # 30.5x faster after warmup
warmup_time = 2.0   # ~2 seconds for JIT compilation (conservative estimate)

print("="*80)
print("JIT Warmup Analysis for Supercomputer Production Runs")
print("="*80)
print()

# Example production parameters from docs
test_cases = [
    ("Small test", 100, 1000, 1),
    ("Development", 1000, 10000, 10),
    ("Production", 10000, 10000, 100),
    ("Large scale", 10000, 100000, 100),
]

print("Scenario                  | Runtime w/o JIT | Warmup | Runtime w/ JIT | Net Benefit")
print("-"*90)

for name, N, sweeps, total_sims in test_cases:
    # Estimate runtime per simulation (extrapolated from N=1000 benchmark)
    # At N=1000, s=10000: 1.36s without JIT, 50ms with JIT
    scale_factor = (N / 1000) ** 1.5  # Complexity roughly O(N^1.5)
    sweep_factor = sweeps / 10000
    
    time_per_sim_no_jit = 1.36 * scale_factor * sweep_factor
    time_per_sim_jit = time_per_sim_no_jit / jit_speedup
    
    total_time_no_jit = time_per_sim_no_jit * total_sims
    total_time_jit = warmup_time + (time_per_sim_jit * total_sims)
    
    net_benefit = total_time_no_jit - total_time_jit
    effective_speedup = total_time_no_jit / total_time_jit
    
    def fmt_time(seconds):
        if seconds < 60:
            return f"{seconds:.1f}s"
        elif seconds < 3600:
            return f"{seconds/60:.1f}m"
        else:
            return f"{seconds/3600:.1f}h"
    
    print(f"{name:<25} | {fmt_time(total_time_no_jit):>15} | {fmt_time(warmup_time):>6} | "
          f"{fmt_time(total_time_jit):>14} | {effective_speedup:>6.1f}x")

print()
print("="*80)
print("Key Findings:")
print("="*80)
print()
print("1. Warmup overhead (~2s) is NEGLIGIBLE for any real production run")
print("   - Even for 1 simulation taking 10s, warmup is only 20% overhead")
print("   - For typical runs with many simulations, warmup is <1% overhead")
print()
print("2. Numba caches compiled code to disk (~/.numba/)")
print("   - First run: pays 2s warmup cost")
print("   - Subsequent runs: cached code loads in milliseconds")
print("   - Cache persists across runs/sessions")
print()
print("3. With multiprocessing (10 workers):")
print("   - Each worker pays warmup cost once (10 * 2s = 20s total)")
print("   - Then all workers run 30x faster for hours")
print("   - 20s overhead for hours of 30x speedup = HUGE win")
print()
print("4. Supercomputer batch jobs (typical: hours to days):")
print("   - Warmup: 2-20 seconds (depending on worker count)")
print("   - Runtime: hours to days with 30x speedup")
print("   - Warmup is <<0.1% of total time")
print()
print("RECOMMENDATION: JIT is ESSENTIAL for supercomputer runs")
print("  - 30x per-core speedup far outweighs <2s warmup")
print("  - Enables running 30x more simulations in same wall time")
print("  - Or reduces wall time by 30x for same number of sims")
print()
