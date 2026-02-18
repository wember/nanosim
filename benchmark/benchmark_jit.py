#!/usr/bin/env python3
"""
Benchmark script to measure optimization performance improvements.
Compares three versions:
1. Original (from repo)
2. Optimized (morning's improvements, no JIT)
3. Optimized + JIT (Numba compilation)
"""

import time
import sys
import numpy as np

sys.path.insert(0, '../creutz-sim')

# Import original version (no optimizations)
from inferno_original import Inferno as InfernoOriginal

# Import optimized version (has morning optimizations, no Numba)
from inferno_nojit import Inferno as InfernoOptimized

# Import JIT version (has morning optimizations + Numba @njit)
from inferno import Inferno as InfernoJIT

def benchmark_simulation(InfernoClass, N, sweeps, warmup=False):
    """
    Benchmark a simulation run.
    
    Args:
        InfernoClass: The Inferno class to use (JIT or non-JIT)
        N: Lattice size
        sweeps: Number of sweeps to run
        warmup: If True, run once before timing to warm up JIT
    
    Returns:
        elapsed_time: Time in seconds
    """
    # Warmup run to trigger JIT compilation (only matters for JIT version)
    if warmup:
        x = InfernoClass(N, R=5)
        for i in range(min(100, sweeps)):
            x.demon_move(flag=0, sweep_count=i)
    
    # Create fresh instance for timed run
    x = InfernoClass(N, R=5)
    
    # Timed run
    start = time.time()
    for i in range(sweeps):
        x.demon_move(flag=0, sweep_count=i)
    elapsed = time.time() - start
    
    return elapsed

def format_time(seconds):
    """Format time in human-readable format."""
    if seconds < 1:
        return f"{seconds*1000:.2f}ms"
    elif seconds < 60:
        return f"{seconds:.2f}s"
    else:
        mins = int(seconds // 60)
        secs = seconds % 60
        return f"{mins}m {secs:.1f}s"

def main():
    """Run benchmarks comparing all three versions."""
    print("="*95)
    print("Complete Optimization Journey Benchmark")
    print("="*95)
    print()
    print("Comparing three versions:")
    print("  1. Original     - From repo (no optimizations)")
    print("  2. Optimized    - Morning improvements (np.unique→loop, simplified logic)")
    print("  3. Optimized+JIT - Morning improvements + Numba JIT compilation")
    print()
    
    # Test cases: (N, sweeps, description)
    test_cases = [
        (50, 500, "Tiny (N=50, s=500)"),
        (100, 1000, "Small (N=100, s=1000)"),
        (200, 2000, "Small+ (N=200, s=2000)"),
        (500, 5000, "Medium (N=500, s=5000)"),
        (1000, 10000, "Large (N=1000, s=10000)"),
        (2000, 10000, "XLarge (N=2000, s=10000)"),
    ]
    
    print("Warmup: Triggering JIT compilation...")
    benchmark_simulation(InfernoJIT, 100, 100, warmup=True)
    print("✓ JIT compilation complete\n")
    
    print(f"{'Test Case':<23} {'Original':<11} {'Optimized':<11} {'Opt+JIT':<11} {'vs Orig':<9} {'JIT gain':<9}")
    print("-"*95)
    
    total_original = 0
    total_optimized = 0
    total_jit = 0
    
    for N, sweeps, description in test_cases:
        # Benchmark original version
        elapsed_original = benchmark_simulation(InfernoOriginal, N, sweeps, warmup=False)
        total_original += elapsed_original
        
        # Benchmark optimized (no JIT)
        elapsed_optimized = benchmark_simulation(InfernoOptimized, N, sweeps, warmup=False)
        total_optimized += elapsed_optimized
        
        # Benchmark with JIT
        elapsed_jit = benchmark_simulation(InfernoJIT, N, sweeps, warmup=False)
        total_jit += elapsed_jit
        
        speedup_total = elapsed_original / elapsed_jit
        speedup_jit_only = elapsed_optimized / elapsed_jit
        
        print(f"{description:<23} {format_time(elapsed_original):<11} {format_time(elapsed_optimized):<11} "
              f"{format_time(elapsed_jit):<11} {speedup_total:>6.1f}x   {speedup_jit_only:>6.1f}x")
    
    print("-"*95)
    total_speedup = total_original / total_jit
    optimized_speedup = total_original / total_optimized
    jit_speedup = total_optimized / total_jit
    
    print(f"{'OVERALL':<23} {format_time(total_original):<11} {format_time(total_optimized):<11} "
          f"{format_time(total_jit):<11} {total_speedup:>6.1f}x   {jit_speedup:>6.1f}x")
    
    print()
    print("="*95)
    print("Summary of Optimization Journey")
    print("="*95)
    print(f"  1. Original → Optimized (no JIT):   {optimized_speedup:>5.2f}x {'speedup' if optimized_speedup > 1 else 'SLOWDOWN'}")
    print(f"  2. Optimized → Optimized+JIT:       {jit_speedup:>5.2f}x speedup")
    print(f"  3. Original → Optimized+JIT:        {total_speedup:>5.2f}x TOTAL speedup")
    print()
    if optimized_speedup < 1:
        print("Key Insight:")
        print(f"  The morning optimizations made pure Python {1/optimized_speedup:.1f}x SLOWER, but this was")
        print("  necessary to enable JIT compilation! NumPy's np.unique() is fast C code,")
        print("  but not JIT-compatible. We traded Python performance for JIT compatibility,")
        print(f"  which gave us a {jit_speedup:.1f}x JIT speedup → {total_speedup:.1f}x total gain.")
        print()
    print(f"  Combined with 10x multiprocessing = {total_speedup*10:>5.0f}x total speedup for production runs!")
    print()
    print("Optimizations applied:")
    print("  ✓ Replaced np.unique() with simple loop (enables JIT, but slower in Python)")
    print("  ✓ Removed unused variables")
    print("  ✓ Simplified conditional logic")
    print("  ✓ Added Numba JIT compilation to hot-path functions (THE BIG WIN)")
    print()
    
    print()

if __name__ == "__main__":
    main()
