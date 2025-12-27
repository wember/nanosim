"""
Production-scale benchmark for Creutz demon simulations.

Benchmarks actual production configurations to measure real-world performance:
- Full-scale parameters: n=1M, s=5k, r=11, m=5
- Tests all optimization strategies:
  * Legacy (no optimization)
  * JIT only (single-core)
  * Parallel only (no JIT)
  * JIT + Parallel (production)

Outputs timing data, speedup factors, and updates documentation templates.
"""

import argparse
import json
import multiprocessing as mp
import os
import subprocess
import sys
import time
from datetime import datetime

# Add creutz-sim to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))


def format_time(seconds):
    """Format seconds as human-readable string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.1f}m"
    else:
        hours = seconds / 3600
        return f"{hours:.2f}h"


def run_command_with_timing(cmd, description, use_venv=True):
    """Run a command and measure wall-clock time.

    Args:
        cmd: Command to run (list of strings)
        description: Human-readable description
        use_venv: Whether to use virtual environment python

    Returns:
        dict with timing info and success status
    """
    # Use venv python if requested
    if use_venv and cmd[0] == "python":
        project_root = os.path.join(os.path.dirname(__file__), "..")
        venv_python = os.path.join(project_root, "venv", "bin", "python")
        if os.path.exists(venv_python):
            cmd[0] = venv_python

    print(f"\n{'='*70}")
    print(f"{description}")
    print(f"{'='*70}")
    print(f"Command: {' '.join(cmd)}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    start_time = time.perf_counter()

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            cwd=os.path.join(os.path.dirname(__file__), ".."),
        )

        elapsed = time.perf_counter() - start_time

        print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Duration: {format_time(elapsed)}")
        print(f"✓ Success")

        return {
            "description": description,
            "command": " ".join(cmd),
            "elapsed_seconds": elapsed,
            "elapsed_formatted": format_time(elapsed),
            "success": True,
            "timestamp": datetime.now().isoformat(),
            "stdout_lines": len(result.stdout.splitlines()),
            "stderr_lines": len(result.stderr.splitlines()),
        }

    except subprocess.CalledProcessError as e:
        elapsed = time.perf_counter() - start_time
        print(f"✗ Failed after {format_time(elapsed)}")
        print(f"Error: {e}")

        return {
            "description": description,
            "command": " ".join(cmd),
            "elapsed_seconds": elapsed,
            "elapsed_formatted": format_time(elapsed),
            "success": False,
            "error": str(e),
            "timestamp": datetime.now().isoformat(),
        }
    except KeyboardInterrupt:
        elapsed = time.perf_counter() - start_time
        print(f"\n⚠ Interrupted after {format_time(elapsed)}")
        raise


def benchmark_legacy_baseline(n, s, r, m, test_mode=False):
    """Benchmark legacy implementation (no optimization).

    This is the absolute baseline - single-core Python without JIT.
    """
    if test_mode:
        # Use same params as other benchmarks for fair comparison
        n, s, r, m = 10000, 100, 3, 2

    cmd = [
        "python",
        "creutz-sim/legacy/sim.py",
        "--n",
        str(n),
        "--s",
        str(s),
        "--r",
        str(r),
        "--m",
        str(m),
    ]

    description = f"Legacy Baseline (n={n}, s={s}, r={r}, m={m})"
    return run_command_with_timing(cmd, description)


def benchmark_jit_only(n, s, r, m, test_mode=False):
    """Benchmark JIT compilation alone (single-core).

    Uses JIT but no parallelization - runs simulations sequentially.
    """
    if test_mode:
        # Use same params as other benchmarks for fair comparison
        n, s, r, m = 10000, 100, 3, 2

    # Use parallel_sim.py with --jit but cores=1 for serial execution
    cmd = [
        "python",
        "creutz-sim/parallel_sim.py",
        "--jit",
        "--n",
        str(n),
        "--s",
        str(s),
        "--r",
        str(r),
        "--m",
        str(m),
        "--cores",
        "1",  # Single core = serial execution
    ]

    description = f"JIT Only - Single Core (n={n}, s={s}, r={r}, m={m})"
    return run_command_with_timing(cmd, description)


def benchmark_parallel_only(n, s, r, m, cores, test_mode=False):
    """Benchmark parallelization without JIT.

    Uses multiprocessing but no JIT compilation.
    Note: This is slow and mainly for comparison purposes.
    """
    if test_mode:
        # Use same params as other benchmarks for fair comparison
        n, s, r, m = 10000, 100, 3, 2

    cmd = [
        "python",
        "creutz-sim/parallel_sim.py",
        "--n",
        str(n),
        "--s",
        str(s),
        "--r",
        str(r),
        "--m",
        str(m),
        "--cores",
        str(cores),
        # No --jit flag = uses original classes
    ]

    description = f"Parallel Only - No JIT (n={n}, s={s}, r={r}, m={m}, cores={cores})"
    return run_command_with_timing(cmd, description)


def benchmark_production(n, s, r, m, cores, test_mode=False):
    """Benchmark production configuration (JIT + Parallel).

    This is what 'make run-sim' uses.
    """
    if test_mode:
        # Use same params as other benchmarks for fair comparison
        n, s, r, m = 10000, 100, 3, 2

    cmd = [
        "python",
        "creutz-sim/parallel_sim.py",
        "--jit",
        "--n",
        str(n),
        "--s",
        str(s),
        "--r",
        str(r),
        "--m",
        str(m),
        "--cores",
        str(cores),
    ]

    description = (
        f"Production (JIT + Parallel) " f"(n={n}, s={s}, r={r}, m={m}, cores={cores})"
    )
    return run_command_with_timing(cmd, description)


def calculate_speedups(results, baseline_key="legacy"):
    """Calculate speedup factors relative to baseline."""
    if baseline_key not in results or not results[baseline_key]["success"]:
        print(f"Warning: Baseline '{baseline_key}' not available")
        return {}

    baseline_time = results[baseline_key]["elapsed_seconds"]
    speedups = {}

    for key, data in results.items():
        if data["success"]:
            speedup = baseline_time / data["elapsed_seconds"]
            speedups[key] = {
                "time": data["elapsed_seconds"],
                "formatted": data["elapsed_formatted"],
                "speedup": speedup,
                "description": data["description"],
            }

    return speedups


def print_summary(results, speedups):
    """Print formatted benchmark summary."""
    print("\n" + "=" * 70)
    print("BENCHMARK SUMMARY")
    print("=" * 70)

    # Table header
    print(f"\n{'Configuration':<35} {'Time':<15} {'Speedup':>10}")
    print("-" * 70)

    # Order: legacy, jit_only, parallel_only, production
    order = ["legacy", "jit_only", "parallel_only", "production"]

    for key in order:
        if key in speedups:
            data = speedups[key]
            name = key.replace("_", " ").title()
            speedup_str = (
                f"{data['speedup']:.1f}x" if data["speedup"] != 1.0 else "1x (baseline)"
            )
            print(f"{name:<35} {data['formatted']:<15} {speedup_str:>10}")

    print("\n" + "=" * 70)
    print("KEY FINDINGS")
    print("=" * 70)

    if "jit_only" in speedups and "legacy" in speedups:
        jit_speedup = speedups["jit_only"]["speedup"]
        print(f"• JIT compilation alone: {jit_speedup:.1f}x speedup")

    if "parallel_only" in speedups and "legacy" in speedups:
        parallel_speedup = speedups["parallel_only"]["speedup"]
        cores = mp.cpu_count()
        efficiency = (parallel_speedup / cores) * 100
        print(
            f"• Parallelization alone: {parallel_speedup:.1f}x speedup "
            f"({efficiency:.0f}% efficiency on {cores} cores)"
        )

    if "production" in speedups and "legacy" in speedups:
        prod_speedup = speedups["production"]["speedup"]
        print(f"• Combined (JIT + Parallel): {prod_speedup:.1f}x speedup")

        # Compare to theoretical max
        if "jit_only" in speedups and "parallel_only" in speedups:
            theoretical = (
                speedups["jit_only"]["speedup"] * speedups["parallel_only"]["speedup"]
            )
            actual_vs_theoretical = (prod_speedup / theoretical) * 100
            print(
                f"• Actual vs theoretical: {actual_vs_theoretical:.0f}% of "
                f"{theoretical:.1f}x"
            )


def save_results(results, speedups, output_file):
    """Save benchmark results to JSON file with system information."""
    import platform
    import subprocess

    # Gather detailed system information
    system_info = {
        "cpu_count": mp.cpu_count(),
        "python_version": sys.version,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "machine": platform.machine(),
        "system": platform.system(),
    }

    # Try to get more detailed CPU info on macOS
    if platform.system() == "Darwin":
        try:
            result = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.brand_string"],
                capture_output=True,
                text=True,
                check=True,
            )
            system_info["cpu_model"] = result.stdout.strip()
        except Exception:
            pass

        try:
            result = subprocess.run(
                ["sysctl", "-n", "hw.memsize"],
                capture_output=True,
                text=True,
                check=True,
            )
            memory_bytes = int(result.stdout.strip())
            system_info["memory_gb"] = round(memory_bytes / (1024**3), 1)
        except Exception:
            pass

    # Try to get CPU info on Linux
    elif platform.system() == "Linux":
        try:
            with open("/proc/cpuinfo", "r") as f:
                for line in f:
                    if "model name" in line:
                        system_info["cpu_model"] = line.split(":")[1].strip()
                        break
        except Exception:
            pass

        try:
            with open("/proc/meminfo", "r") as f:
                for line in f:
                    if "MemTotal" in line:
                        kb = int(line.split()[1])
                        system_info["memory_gb"] = round(kb / (1024**2), 1)
                        break
        except Exception:
            pass

    output = {
        "timestamp": datetime.now().isoformat(),
        "system_info": system_info,
        "raw_results": results,
        "speedups": speedups,
    }

    with open(output_file, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\n✓ Results saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark production-scale Creutz demon simulations"
    )
    parser.add_argument(
        "--config",
        choices=["test", "small", "production"],
        default="test",
        help="Benchmark configuration (test=quick, small=medium, production=full)",
    )
    parser.add_argument(
        "--skip-legacy", action="store_true", help="Skip legacy baseline (very slow)"
    )
    parser.add_argument(
        "--skip-parallel-only",
        action="store_true",
        help="Skip parallel-without-JIT test (slow)",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=mp.cpu_count(),
        help="Number of cores for parallel tests",
    )
    parser.add_argument(
        "--output",
        default=os.path.join(os.path.dirname(__file__), "benchmark_results.json"),
        help="Output JSON file for results (default: tools/benchmark_results.json)",
    )

    args = parser.parse_args()

    # Configuration presets
    configs = {
        "test": (10000, 100, 3, 2, True),  # Quick test: ~1 minute total
        "small": (100000, 1000, 5, 3, False),  # Small: ~10-15 minutes
        "production": (1000000, 5000, 11, 5, False),  # Full: hours
    }

    n, s, r, m, test_mode = configs[args.config]

    print(
        f"""
╔══════════════════════════════════════════════════════════════════╗
║           PRODUCTION BENCHMARK - CREUTZ DEMON                    ║
╚══════════════════════════════════════════════════════════════════╝

Configuration: {args.config}
  Lattice size: n = {n:,}
  Sweeps:       s = {s:,}
  Radii:        r = 1 to {r}
  Runs:         m = {m}
  Cores:        {args.cores}

Total simulations: {r * m * 2} forward+reverse runs per benchmark

Benchmarks to run:
  {'✓' if not args.skip_legacy else '✗'} Legacy baseline (single-core, no JIT)
  ✓ JIT only (single-core with JIT)
  {'✓' if not args.skip_parallel_only else '✗'} Parallel only (no JIT)
  ✓ Production (JIT + Parallel)
"""
    )

    input("Press Enter to start benchmarks (Ctrl+C to cancel)...")

    results = {}

    try:
        # 1. Legacy baseline (slowest, run first so we can estimate total time)
        if not args.skip_legacy:
            results["legacy"] = benchmark_legacy_baseline(n, s, r, m, test_mode)

        # 2. JIT only
        results["jit_only"] = benchmark_jit_only(n, s, r, m, test_mode)

        # 3. Parallel only (slow without JIT)
        if not args.skip_parallel_only:
            results["parallel_only"] = benchmark_parallel_only(
                n, s, r, m, args.cores, test_mode
            )

        # 4. Production (fastest)
        results["production"] = benchmark_production(n, s, r, m, args.cores, test_mode)

        # Calculate speedups
        baseline = "legacy" if not args.skip_legacy else "jit_only"
        speedups = calculate_speedups(results, baseline_key=baseline)

        # Print summary
        print_summary(results, speedups)

        # Save results
        save_results(results, speedups, args.output)

        print(f"\n{'='*70}")
        print("BENCHMARK COMPLETE")
        print(f"{'='*70}\n")

    except KeyboardInterrupt:
        print("\n\n⚠ Benchmark interrupted by user")
        if results:
            print("\nPartial results:")
            for key, data in results.items():
                if data["success"]:
                    print(f"  {key}: {data['elapsed_formatted']}")
        sys.exit(1)


if __name__ == "__main__":
    main()
