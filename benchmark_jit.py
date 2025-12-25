"""
Benchmark JIT vs non-JIT performance for Inferno simulations.

Compares:
1. Original Inferno vs JIT Inferno (reversible)
2. Original irrInferno vs JIT irrInferno (irreversible)

Reports speedup factors and validates correctness.
"""

import time
import numpy as np
import sys
import os

# Add creutz-sim to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'creutz-sim'))

from inferno import Inferno
from irr_inferno import irrInferno
from jit_inferno import JITInferno
from jit_irr_inferno import JITirrInferno


def benchmark_simulation(sim_class, name, n, r, sweeps, warmup=True):
    """
    Benchmark a simulation class.
    
    Args:
        sim_class: Simulation class to benchmark
        name: Name for reporting
        n: Number of lattice sites
        r: Demon coupling radius
        sweeps: Number of sweeps to run
        warmup: Whether to do warmup run (for JIT compilation)
        
    Returns:
        tuple: (elapsed_time, final_energy, final_entropy)
    """
    print(f"\n{'='*60}")
    print(f"Benchmarking {name}")
    print(f"N={n}, R={r}, sweeps={sweeps}")
    print(f"{'='*60}")
    
    # Warmup run for JIT (compile functions)
    if warmup and 'JIT' in name:
        print("Warming up JIT compiler...")
        sim_warmup = sim_class(n, r, validate_mode='off')
        for _ in range(5):  # Short warmup
            sim_warmup.demon_move()
        print("Warmup complete")
    
    # Initialize simulation
    sim = sim_class(n, r, validate_mode='off')
    
    # Record initial state
    initial_energy = sim.E_total
    initial_lattice_energy = sim.E_lattice
    initial_demon_energy = sim.E_demon_sum
    
    print(f"\nInitial state:")
    print(f"  Total energy: {initial_energy}")
    print(f"  Lattice energy: {initial_lattice_energy}")
    print(f"  Demon energy: {initial_demon_energy}")
    
    # Benchmark
    print(f"\nRunning {sweeps} sweeps...")
    start = time.perf_counter()
    
    # Check if this is a JIT version (does full sweep per call)
    is_jit = 'JIT' in name
    
    for sweep in range(sweeps):
        if is_jit:
            # JIT version: one call = one full sweep
            sim.demon_move()
        else:
            # Original version: N calls = one full sweep
            for j in range(n):
                sim.demon_move()
        
        if (sweep + 1) % 20 == 0:
            print(f"  Sweep {sweep+1}/{sweeps} complete", end='\r')
    
    elapsed = time.perf_counter() - start
    print(f"\n\nCompleted in {elapsed:.4f} seconds")
    print(f"Time per sweep: {elapsed/sweeps*1000:.2f} ms")
    
    # Verify energy conservation
    final_energy = sim.E_lattice + sim.E_demon_sum
    energy_drift = abs(final_energy - initial_energy)
    
    print(f"\nFinal state:")
    print(f"  Total energy: {final_energy}")
    print(f"  Energy drift: {energy_drift}")
    print(f"  Energy conserved: {'✓' if energy_drift == 0 else '✗'}")
    
    # Calculate entropy
    from sim_utils import Sk, Su
    N0 = sim.bond_count[0]
    Nx = sim.bond_count[2]
    K = sim.E_demon_sum  # Total demon energy
    final_entropy = (Sk(n, K) + Su(N0, Nx, n)) / n
    
    print(f"  Final entropy: {final_entropy:.6f}")
    
    return elapsed, final_energy, final_entropy


def compare_reversible(n=10000, r=5, sweeps=100):
    """Compare Inferno vs JITInferno."""
    print("\n" + "="*70)
    print("REVERSIBLE SIMULATION COMPARISON")
    print("="*70)
    
    # Benchmark original
    time_orig, energy_orig, entropy_orig = benchmark_simulation(
        Inferno, "Original Inferno", n, r, sweeps, warmup=False
    )
    
    # Benchmark JIT
    time_jit, energy_jit, entropy_jit = benchmark_simulation(
        JITInferno, "JIT Inferno", n, r, sweeps, warmup=True
    )
    
    # Report speedup
    speedup = time_orig / time_jit
    print("\n" + "="*70)
    print("REVERSIBLE RESULTS")
    print("="*70)
    print(f"Original time:  {time_orig:.4f}s")
    print(f"JIT time:       {time_jit:.4f}s")
    print(f"Speedup:        {speedup:.2f}x")
    print(f"\nEnergy difference: {abs(energy_orig - energy_jit)}")
    print(f"Entropy difference: {abs(entropy_orig - entropy_jit):.6f}")
    
    return speedup


def compare_irreversible(n=10000, r=5, sweeps=100):
    """Compare irrInferno vs JITirrInferno."""
    print("\n" + "="*70)
    print("IRREVERSIBLE SIMULATION COMPARISON")
    print("="*70)
    
    # Benchmark original
    time_orig, energy_orig, entropy_orig = benchmark_simulation(
        irrInferno, "Original irrInferno", n, r, sweeps, warmup=False
    )
    
    # Benchmark JIT
    time_jit, energy_jit, entropy_jit = benchmark_simulation(
        JITirrInferno, "JIT irrInferno", n, r, sweeps, warmup=True
    )
    
    # Report speedup
    speedup = time_orig / time_jit
    print("\n" + "="*70)
    print("IRREVERSIBLE RESULTS")
    print("="*70)
    print(f"Original time:  {time_orig:.4f}s")
    print(f"JIT time:       {time_jit:.4f}s")
    print(f"Speedup:        {speedup:.2f}x")
    print(f"\nEnergy difference: {abs(energy_orig - energy_jit)}")
    print(f"Entropy difference: {abs(entropy_orig - entropy_jit):.6f}")
    
    return speedup


if __name__ == "__main__":
    import sys
    
    # Parse command line arguments
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
    sweeps = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    r = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    
    print(f"""
╔══════════════════════════════════════════════════════════════════╗
║             NUMBA JIT COMPILATION BENCHMARK                      ║
║                                                                  ║
║  Testing performance of JIT-compiled hot functions              ║
║  Expected speedup: 3-10x for single-core simulations            ║
╚══════════════════════════════════════════════════════════════════╝

Configuration:
  Lattice size: N = {n}
  Sweeps:       {sweeps}
  Radius:       R = {r}
    """)
    
    try:
        # Run benchmarks
        speedup_rev = compare_reversible(n, r, sweeps)
        speedup_irr = compare_irreversible(n, r, sweeps)
        
        # Summary
        print("\n" + "="*70)
        print("OVERALL SUMMARY")
        print("="*70)
        print(f"Reversible speedup:   {speedup_rev:.2f}x")
        print(f"Irreversible speedup: {speedup_irr:.2f}x")
        print(f"Average speedup:      {(speedup_rev + speedup_irr)/2:.2f}x")
        
        # Estimate production impact
        print("\n" + "="*70)
        print("PRODUCTION IMPACT ESTIMATE")
        print("="*70)
        print(f"For parallel runs (16 cores):")
        print(f"  Combined speedup: {16 * (speedup_rev + speedup_irr)/2:.1f}x")
        print(f"  Original time: 27 hours → ~{27*60/(16*(speedup_rev+speedup_irr)/2):.1f} minutes")
        
    except KeyboardInterrupt:
        print("\n\nBenchmark interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n\nError during benchmark: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
