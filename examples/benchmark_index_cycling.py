#!/usr/bin/env python3
"""
Benchmark: Index Cycling Optimization

Measures the performance impact of replacing modulo operations with
conditional reset for index cycling in the main Monte Carlo loop.

Optimization:
  Old: self.order_idx = (self.order_idx + 1) % self.N
  New: self.order_idx += 1
       if self.order_idx >= self.N:
           self.order_idx = 0

This benchmark runs before/after comparisons to validate the optimization.
"""

import sys
import os
import time

# Add creutz-sim directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'creutz-sim'))

from inferno import Inferno
from inferno_irr import irrInferno


def benchmark_inferno(N, sweeps, R=5, runs=3):
    """Benchmark Inferno with current implementation."""
    times = []
    
    for run in range(runs):
        x = Inferno(N, R, validate_mode='off')
        
        start = time.perf_counter()
        for _ in range(sweeps):
            for _ in range(x.N):
                x.demon_move()
        elapsed = time.perf_counter() - start
        
        times.append(elapsed)
        print(f"  Inferno run {run+1}: {elapsed:.3f}s")
    
    avg = sum(times) / len(times)
    print(f"  Average: {avg:.3f}s")
    return avg


def benchmark_irr_inferno(N, sweeps, R=5, runs=3):
    """Benchmark irrInferno with current implementation."""
    times = []
    
    for run in range(runs):
        x = irrInferno(N, R, validate_mode='off')
        
        start = time.perf_counter()
        for _ in range(sweeps):
            for _ in range(x.N):
                x.demon_move()
        elapsed = time.perf_counter() - start
        
        times.append(elapsed)
        print(f"  irrInferno run {run+1}: {elapsed:.3f}s")
    
    avg = sum(times) / len(times)
    print(f"  Average: {avg:.3f}s")
    return avg


def main():
    print("=" * 70)
    print("INDEX CYCLING OPTIMIZATION BENCHMARK")
    print("=" * 70)
    print()
    print("This benchmark measures the performance impact of replacing")
    print("modulo operations with conditional reset for index cycling.")
    print()
    
    # Small lattice
    print("Small lattice (N=1,000, 1000 sweeps, R=5)")
    print("-" * 70)
    benchmark_inferno(1000, 1000)
    print()
    benchmark_irr_inferno(1000, 1000)
    print("\n")
    
    # Medium lattice
    print("Medium lattice (N=10,000, 100 sweeps, R=5)")
    print("-" * 70)
    benchmark_inferno(10000, 100)
    print()
    benchmark_irr_inferno(10000, 100)
    print("\n")
    
    # Large lattice
    print("Large lattice (N=100,000, 10 sweeps, R=5)")
    print("-" * 70)
    benchmark_inferno(100000, 10)
    print()
    benchmark_irr_inferno(100000, 10)
    print("\n")
    
    print("=" * 70)
    print("BENCHMARK COMPLETE")
    print("=" * 70)
    print()
    print("NOTES:")
    print("- Expected improvement: 2-5% speedup")
    print("- Branch predictors should handle conditional reset efficiently")
    print("- Operations replaced: 6 modulo per demon_move() in Inferno")
    print("- Operations replaced: 2 modulo per demon_move() in irrInferno")
    print("- Conditional resets: idx += 1; if idx >= N: idx = 0")


if __name__ == '__main__':
    main()
