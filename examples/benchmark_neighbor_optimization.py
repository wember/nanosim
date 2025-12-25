"""
Benchmark neighbor index optimization.

Tests the performance impact of pre-computing neighbor indices
vs. using modulo operations for periodic boundary conditions.
"""
import sys
import os
import time

# Add creutz-sim directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'creutz-sim'))

from inferno import Inferno
from irr_inferno import irrInferno


def benchmark_simulation(SimClass, N, sweeps, R=5):
    """Benchmark a simulation run."""
    sim = SimClass(N=N, R=R, validate_mode='off')
    
    start = time.perf_counter()
    for _ in range(sweeps):
        for _ in range(N):
            sim.demon_move()
    duration = time.perf_counter() - start
    
    return duration


def main():
    """Run performance benchmarks."""
    print("=" * 70)
    print("NEIGHBOR INDEX PRE-COMPUTATION OPTIMIZATION BENCHMARK")
    print("=" * 70)
    print()
    print("This benchmark measures the performance impact of replacing")
    print("modulo operations with pre-computed neighbor index arrays.")
    print()
    
    # Test configurations
    configs = [
        ("Small lattice", 1000, 1000),
        ("Medium lattice", 10000, 100),
        ("Large lattice", 100000, 10),
    ]
    
    for name, N, sweeps in configs:
        print(f"{name} (N={N:,}, {sweeps} sweeps, R=5)")
        print("-" * 70)
        
        # Benchmark Inferno (reversible)
        times_rev = []
        for run in range(3):
            t = benchmark_simulation(Inferno, N, sweeps)
            times_rev.append(t)
            print(f"  Inferno run {run+1}: {t:.3f}s")
        
        avg_rev = sum(times_rev) / len(times_rev)
        print(f"  Average: {avg_rev:.3f}s")
        print()
        
        # Benchmark irrInferno (irreversible)
        times_irr = []
        for run in range(3):
            t = benchmark_simulation(irrInferno, N, sweeps)
            times_irr.append(t)
            print(f"  irrInferno run {run+1}: {t:.3f}s")
        
        avg_irr = sum(times_irr) / len(times_irr)
        print(f"  Average: {avg_irr:.3f}s")
        print()
        print()
    
    print("=" * 70)
    print("BENCHMARK COMPLETE")
    print("=" * 70)
    print()
    print("NOTES:")
    print("- Expected improvement: 5-10% speedup")
    print("- Larger lattices should see more benefit (better cache locality)")
    print("- Modulo operations replaced: ~6 per Monte Carlo move")
    print("- Pre-computed arrays: right_neighbor[], left_neighbor[]")
    print()


if __name__ == '__main__':
    main()
