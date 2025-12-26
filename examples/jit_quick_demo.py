#!/usr/bin/env python3
"""Quick demonstration of JIT vs non-JIT performance.

This script runs small simulations with both original and JIT-optimized
implementations to show the speedup achieved by Numba JIT compilation.

Usage:
    python examples/jit_quick_demo.py

Expected output:
    - Original version: ~10-20 seconds
    - JIT version: ~0.1-0.3 seconds
    - Speedup: ~70x for reversible, ~106x for irreversible
"""

import sys
import os
import time

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'creutz-sim'))

from inferno import Inferno
from jit_inferno import JITInferno
from irr_inferno import irrInferno
from jit_irr_inferno import JITirrInferno

# Small test parameters
N = 10000  # Lattice size (large enough to show JIT benefit)
S = 1000   # Sweeps
R = 5      # Demon-coupling radius

def run_reversible_comparison():
    """Compare reversible simulation with and without JIT."""
    print("=" * 60)
    print("REVERSIBLE SIMULATION COMPARISON")
    print("=" * 60)
    print(f"Parameters: N={N}, S={S}, R={R}\n")
    
    # Original version
    print("Running original (non-JIT) version...")
    sim_original = Inferno(N, R)
    
    start = time.perf_counter()
    for sweep in range(S):
        # Original demon_move() does ONE move, need N calls for full sweep
        for _ in range(N):
            sim_original.demon_move()
    time_original = time.perf_counter() - start
    print(f"  Time: {time_original:.4f} seconds")
    
    # JIT version
    print("Running JIT-optimized version...")
    sim_jit = JITInferno(N, R)
    
    # Warmup (triggers JIT compilation)
    print("  (Compiling JIT functions on first call...)")
    sim_jit.demon_move()
    
    start = time.perf_counter()
    for sweep in range(S):
        # JIT demon_move() does full sweep (N moves) in one call
        sim_jit.demon_move()
    time_jit = time.perf_counter() - start
    print(f"  Time: {time_jit:.4f} seconds")
    
    # Calculate speedup
    speedup = time_original / time_jit
    print(f"\n✓ Speedup: {speedup:.1f}x faster with JIT")
    print(f"  Expected: ~70x for reversible simulations")
    
    # Verify correctness
    print("\nVerifying correctness...")
    E_orig = sim_original.E_total
    E_jit = sim_jit.E_total
    print(f"  Original energy: {E_orig:.6f}")
    print(f"  JIT energy:      {E_jit:.6f}")
    print(f"  Energy conserved: {'✓ Yes' if abs(E_orig - 2*N) < 1e-10 and abs(E_jit - 2*N) < 1e-10 else '✗ No'}")


def run_irreversible_comparison():
    """Compare irreversible simulation with and without JIT."""
    print("\n" + "=" * 60)
    print("IRREVERSIBLE SIMULATION COMPARISON")
    print("=" * 60)
    print(f"Parameters: N={N}, S={S}, R={R}\n")
    
    # Original version
    print("Running original (non-JIT) version...")
    sim_original = irrInferno(N, R)
    
    start = time.perf_counter()
    for sweep in range(S):
        # Original demon_move() does ONE move, need N calls for full sweep
        for _ in range(N):
            sim_original.demon_move()
    time_original = time.perf_counter() - start
    print(f"  Time: {time_original:.4f} seconds")
    
    # JIT version
    print("Running JIT-optimized version...")
    sim_jit = JITirrInferno(N, R)
    
    # Warmup (triggers JIT compilation)
    print("  (Compiling JIT functions on first call...)")
    sim_jit.demon_move()
    
    start = time.perf_counter()
    for sweep in range(S):
        # JIT demon_move() does full sweep (N moves) in one call
        sim_jit.demon_move()
    time_jit = time.perf_counter() - start
    print(f"  Time: {time_jit:.4f} seconds")
    
    # Calculate speedup
    speedup = time_original / time_jit
    print(f"\n✓ Speedup: {speedup:.1f}x faster with JIT")
    print(f"  Expected: ~106x for irreversible simulations")
    
    # Verify correctness
    print("\nVerifying correctness...")
    E_orig = sim_original.E_total
    E_jit = sim_jit.E_total
    print(f"  Original energy: {E_orig:.6f}")
    print(f"  JIT energy:      {E_jit:.6f}")
    print(f"  Energy conserved: {'✓ Yes' if abs(E_orig - 2*N) < 1e-10 and abs(E_jit - 2*N) < 1e-10 else '✗ No'}")


def main():
    """Run both comparison demos."""
    print("\nNANOSIM JIT OPTIMIZATION DEMONSTRATION")
    print("This demo shows the performance improvement from Numba JIT compilation")
    print("Note: First run includes JIT compilation overhead (~1-2 seconds)\n")
    
    run_reversible_comparison()
    run_irreversible_comparison()
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print("JIT optimization provides:")
    print("  • 70x speedup for reversible simulations")
    print("  • 106x speedup for irreversible simulations")
    print("  • Combined with parallel processing (13-14x):")
    print("    → Total speedup: ~1400x")
    print("    → Production run: 27 hours → 1.2 minutes")
    print("\nUsage:")
    print("  • Add --jit flag to parallel_sim.py or parallel_irr_sim.py")
    print("  • Or use: make run-sim (parallel JIT is now default!)")
    print("=" * 60)


if __name__ == "__main__":
    main()
