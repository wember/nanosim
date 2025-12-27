#!/usr/bin/env python3
"""
Profile the nanosim simulation to identify performance bottlenecks.

Usage:
    python profile_sim.py
        [--mode {inferno,irr_inferno}]
        [--n N]
        [--s SWEEPS]
        [--r RADIUS]
"""

import argparse
import cProfile
import io
import pstats
import sys
from pathlib import Path

# Add creutz-sim to path (going up one level from tools/)
sys.path.insert(0, str(Path(__file__).parent.parent / "creutz-sim"))

from inferno import Inferno
from inferno_irr import irrInferno


def run_simulation(mode="inferno", n=10000, sweeps=100, radius=5):
    """Run a simulation for profiling."""
    if mode == "inferno":
        x = Inferno(N=n, R=radius, validate_mode="off")
    else:
        x = irrInferno(N=n, R=radius, validate_mode="off")

    # Run forward sweeps
    for _ in range(sweeps):
        for _ in range(n):
            x.demon_move()

    # Run reverse sweeps
    for _ in range(sweeps):
        for _ in range(n):
            x.demon_reverse()

    return x


def profile_simulation(mode="inferno", n=10000, sweeps=100, radius=5):
    """Profile the simulation and print results."""
    print(f"\n{'='*80}")
    print(f"Profiling {mode.upper()} simulation")
    print(f"Parameters: N={n}, sweeps={sweeps}, R={radius}")
    print(f"Total moves: {2 * n * sweeps:,}")
    print(f"{'='*80}\n")

    # Create profiler
    profiler = cProfile.Profile()

    # Run with profiling
    profiler.enable()
    x = run_simulation(mode=mode, n=n, sweeps=sweeps, radius=radius)
    profiler.disable()

    # Verify energy conservation
    state = x.get_validated_state()
    print(f"Final energy check: {state['E_total']} (expected: {2*n})")

    # Print results
    print("\n" + "=" * 80)
    print("TOP 30 FUNCTIONS BY CUMULATIVE TIME")
    print("=" * 80 + "\n")

    s = io.StringIO()
    ps = pstats.Stats(profiler, stream=s)
    ps.sort_stats("cumulative")
    ps.print_stats(30)
    print(s.getvalue())

    print("\n" + "=" * 80)
    print("TOP 30 FUNCTIONS BY TOTAL TIME (excluding subcalls)")
    print("=" * 80 + "\n")

    s = io.StringIO()
    ps = pstats.Stats(profiler, stream=s)
    ps.sort_stats("tottime")
    ps.print_stats(30)
    print(s.getvalue())

    # Save detailed stats to file
    output_file = f"profile_{mode}_n{n}_s{sweeps}.stats"
    profiler.dump_stats(output_file)
    print(f"\n✓ Detailed profile saved to: {output_file}")
    print(f"  View with: python -m pstats {output_file}")

    return profiler


def main():
    parser = argparse.ArgumentParser(description="Profile nanosim simulation")
    parser.add_argument(
        "--mode",
        choices=["inferno", "irr_inferno"],
        default="inferno",
        help="Simulation mode to profile",
    )
    parser.add_argument(
        "--n", type=int, default=10000, help="Lattice size (default: 10000)"
    )
    parser.add_argument(
        "--s", type=int, default=100, help="Number of sweeps (default: 100)"
    )
    parser.add_argument(
        "--r", type=int, default=5, help="Demon coupling radius (default: 5)"
    )
    parser.add_argument(
        "--both", action="store_true", help="Profile both Inferno and irrInferno"
    )

    args = parser.parse_args()

    if args.both:
        # Profile both modes
        print("\n" + "🔬 " * 20)
        print(" " * 20 + "PROFILING BOTH MODES")
        print("🔬 " * 20 + "\n")

        profile_simulation("inferno", args.n, args.s, args.r)
        profile_simulation("irr_inferno", args.n, args.s, args.r)

        print("\n" + "=" * 80)
        print("COMPARISON SUMMARY")
        print("=" * 80)
        print("\nBoth profiles complete. Compare with:")
        print(f"  python -m pstats profile_inferno_n{args.n}_s{args.s}.stats")
        print(f"  python -m pstats profile_irr_inferno_n{args.n}_s{args.s}.stats")
    else:
        profile_simulation(args.mode, args.n, args.s, args.r)


if __name__ == "__main__":
    main()
