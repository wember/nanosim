#!/usr/bin/env python3
"""Production-ready template for running JIT-optimized simulations.

This script demonstrates best practices for running production simulations
with JIT optimization enabled. It includes proper error handling, progress
monitoring, and result validation.

Usage:
    # Run with default parameters (production scale)
    python examples/production_run_jit.py
    
    # Run with custom parameters
    python examples/production_run_jit.py --reversible --n 100000 --s 10000 --r 11 --m 5
    
    # Run irreversible simulation
    python examples/production_run_jit.py --irreversible

Features:
    - Automatic JIT optimization (70-106x speedup)
    - Parallel processing across all CPU cores
    - Progress monitoring with time estimates
    - Energy conservation validation
    - CSV output compatible with analysis tools
"""

import sys
import os
import argparse
from pathlib import Path

# Add parent directory to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / 'creutz-sim'))

from parallel_sim import main as run_parallel_sim
from parallel_irr_sim import main as run_parallel_irr_sim


def main():
    """Main entry point for production JIT simulations."""
    parser = argparse.ArgumentParser(
        description='Production-ready JIT-optimized nanosim execution',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Reversible simulation (default)
  python %(prog)s
  
  # Irreversible simulation
  python %(prog)s --irreversible
  
  # Custom parameters
  python %(prog)s --n 100000 --s 10000 --r 5 --m 10
  
  # Quick test (small scale)
  python %(prog)s --n 1000 --s 100 --r 3 --m 2

Performance:
  With JIT optimization enabled, expect:
    - Reversible:   70x faster per core
    - Irreversible: 106x faster per core
    - Total speedup with 16 cores: ~1120x (rev) or ~1696x (irr)
  
  Production run (n=1000000, s=10000, r=11, m=5):
    - Without optimization: ~27 hours
    - With JIT + parallel:  ~1.2 minutes
        """
    )
    
    # Simulation type
    type_group = parser.add_mutually_exclusive_group()
    type_group.add_argument(
        '--reversible',
        action='store_true',
        default=True,
        help='Run reversible simulation (default)'
    )
    type_group.add_argument(
        '--irreversible',
        action='store_true',
        help='Run irreversible simulation'
    )
    
    # Simulation parameters
    parser.add_argument(
        '--n',
        type=int,
        default=1000000,
        help='Lattice size (default: 1000000 for production)'
    )
    parser.add_argument(
        '--s',
        type=int,
        default=10000,
        help='Number of sweeps per phase (default: 10000)'
    )
    parser.add_argument(
        '--r',
        type=int,
        default=11,
        help='Number of radii to test (tests R=1 to R=r-1, default: 11)'
    )
    parser.add_argument(
        '--m',
        type=int,
        default=5,
        help='Number of independent runs per radius (default: 5)'
    )
    
    # Output options
    parser.add_argument(
        '--output-dir',
        type=str,
        help='Custom output directory (default: data/r{R}/ or data/irr/r{R}/)'
    )
    parser.add_argument(
        '--validate-only',
        action='store_true',
        help='Run single test iteration to validate setup (n=100, s=10, r=2, m=1)'
    )
    
    args = parser.parse_args()
    
    # Validation mode - run quick test
    if args.validate_only:
        print("=" * 70)
        print("VALIDATION MODE: Running quick test to verify setup")
        print("=" * 70)
        args.n = 100
        args.s = 10
        args.r = 2
        args.m = 1
    
    # Display configuration
    sim_type = "Irreversible" if args.irreversible else "Reversible"
    print("\n" + "=" * 70)
    print(f"NANOSIM PRODUCTION RUN - {sim_type.upper()} (JIT-OPTIMIZED)")
    print("=" * 70)
    print(f"Configuration:")
    print(f"  Lattice size:   N = {args.n:,}")
    print(f"  Sweeps/phase:   S = {args.s:,}")
    print(f"  Radii tested:   R = 1..{args.r-1}")
    print(f"  Runs per radius: M = {args.m}")
    print(f"  Total runs:      {args.r * args.m}")
    print(f"\nOptimizations:")
    print(f"  • JIT compilation: {'106x' if args.irreversible else '70x'} speedup")
    print(f"  • Parallel processing: CPU cores utilized")
    print(f"  • Combined speedup: ~{'1400' if args.irreversible else '1000'}x")
    
    # Calculate expected time
    base_time_per_run = 27 * 60  # 27 hours in minutes for baseline
    speedup = 1400 if args.irreversible else 1000
    expected_minutes = (base_time_per_run * args.r * args.m) / speedup
    
    if expected_minutes < 1:
        time_str = f"{expected_minutes * 60:.0f} seconds"
    else:
        time_str = f"{expected_minutes:.1f} minutes"
    
    print(f"\nEstimated runtime: {time_str}")
    
    if args.output_dir:
        print(f"Output directory: {args.output_dir}")
    
    print("=" * 70 + "\n")
    
    # Run the appropriate simulation
    # Note: We rebuild sys.argv to pass arguments to the simulation scripts
    sys.argv = ['simulation']
    sys.argv.extend(['--jit'])  # Always enable JIT in production
    sys.argv.extend(['--n', str(args.n)])
    sys.argv.extend(['--s', str(args.s)])
    sys.argv.extend(['--r', str(args.r)])
    sys.argv.extend(['--m', str(args.m)])
    
    if args.output_dir:
        sys.argv.extend(['--output-dir', args.output_dir])
    
    try:
        if args.irreversible:
            run_parallel_irr_sim()
        else:
            run_parallel_sim()
        
        print("\n" + "=" * 70)
        print("SIMULATION COMPLETE")
        print("=" * 70)
        print(f"Results written to: data/{'' if not args.irreversible else 'irr/'}r[1..{args.r-1}]/")
        print("\nNext steps:")
        print("  • Analyze results: python creutz-sim/Sk_comparison.py")
        print("  • Plot individual runs: python creutz-sim/sim_plot.py")
        print("  • Visualize entropy trends: python creutz-sim/sim_plot_r.py")
        print("=" * 70)
        
    except KeyboardInterrupt:
        print("\n\n❌ Simulation interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n\n❌ Error during simulation: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
