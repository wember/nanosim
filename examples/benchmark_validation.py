"""
Benchmark script to measure validation overhead reduction.
Compares performance with different validation modes.
"""

import os
import sys
import time

# Add creutz-sim directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))

from inferno import Inferno


def benchmark_validation_mode(N, sweeps, R, mode):
    """Run benchmark with specific validation mode."""
    x = Inferno(N, R, validate_mode=mode)

    start = time.time()
    for _ in range(sweeps):
        for _ in range(x.N):
            x.demon_move()
    duration = time.time() - start

    return duration


def main():
    print("=" * 70)
    print("Validation Overhead Benchmark")
    print("=" * 70)

    # Test configurations
    configs = [
        ("Small", 1000, 1000, 5),
        ("Medium", 10000, 100, 5),
        ("Large", 50000, 20, 5),
    ]

    modes = ["off", "periodic", "frequent"]

    results = {}

    for name, N, sweeps, R in configs:
        print(f"\n{name} lattice (N={N:,}, {sweeps} sweeps, R={R})")
        print("-" * 70)

        results[name] = {}

        for mode in modes:
            duration = benchmark_validation_mode(N, sweeps, R, mode)
            results[name][mode] = duration

            if mode == "off":
                baseline = duration
                speedup = 1.0
            else:
                speedup = baseline / duration

            print(f"  {mode:12s}: {duration:6.3f}s  (speedup: {speedup:.2f}x)")

    # Summary
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)

    for name in results:
        off_time = results[name]["off"]
        periodic_time = results[name]["periodic"]
        frequent_time = results[name]["frequent"]

        periodic_speedup = frequent_time / periodic_time
        off_speedup = frequent_time / off_time

        print(f"\n{name} lattice:")
        print(f"  'periodic' vs 'frequent': {periodic_speedup:.2f}x faster")
        print(f"  'off' vs 'frequent':      {off_speedup:.2f}x faster")
        print(f"  'off' vs 'periodic':      {off_speedup/periodic_speedup:.2f}x faster")


if __name__ == "__main__":
    main()
