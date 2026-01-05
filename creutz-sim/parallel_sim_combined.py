"""
Copyright (c) 2026 Winry Ember
Licensed under the MIT License.
See LICENSE file in the project root for full license information.

Combined parallel implementation for both reversible and irreversible simulations.

Runs both simulation types simultaneously using the same parameter sets,
effectively halving total execution time by utilizing available cores for both.

JIT compilation is enabled by default for maximum speedup (~3,880x per core).
Use --no-jit to disable JIT for debugging or validation.
"""

import csv
import json
import logging
import multiprocessing as mp
import os
import time
from datetime import datetime
from multiprocessing import Manager
from typing import Dict, Tuple

import numpy as np
from sim_utils import (
    Sk_stable,
    Su_stable,
    create_argument_parser,
    create_data_directories,
    print_simulation_info,
    setup_logging,
)
from tqdm import tqdm


def run_single_simulation(
    args: Tuple[int, int, int, int, str, str, bool, int, str, "mp.Queue"],
) -> Dict:
    """Run a single simulation for given parameters.

    Args:
        args: Tuple of (R, M, n, s, validate_mode, project_root, use_jit, sim_num,
                sim_type, progress_queue)
            R: Demon-coupling radius
            M: Run number
            n: Lattice size
            s: Number of sweeps
            validate_mode: Validation mode ('off', 'periodic', 'frequent')
            project_root: Project root directory
            use_jit: Whether to use JIT-compiled version
            sim_num: Simulation number (1-based)
            sim_type: 'reversible' or 'irreversible'
            progress_queue: Queue for progress updates

    Returns:
        Dictionary with simulation results and metadata
    """
    (
        R,
        M,
        n,
        s,
        validate_mode,
        project_root,
        use_jit,
        sim_num,
        sim_type,
        progress_queue,
    ) = args

    # Send initial status
    progress_queue.put(
        {"type": "start", "sim_num": sim_num, "R": R, "M": M, "sim_type": sim_type}
    )

    # Import appropriate class based on simulation type
    if sim_type == "reversible":
        if use_jit:
            from jit_inferno import JITInferno as SimClass
        else:
            from inferno import Inferno as SimClass
        file_prefix = "sim_data"
    else:  # irreversible
        if use_jit:
            from jit_inferno_irr import JITirrInferno as SimClass
        else:
            from inferno_irr import irrInferno as SimClass
        file_prefix = "irr_sim_data"

    # Create simulation instance
    x = SimClass(n, R + 1, validate_mode=validate_mode)

    # Setup output paths
    if sim_type == "reversible":
        folder = os.path.join(project_root, "data")
    else:
        folder = os.path.join(project_root, "data", "irr")

    if R == 0:
        output_dir = os.path.join(folder, "r0")
        filename = os.path.join(output_dir, f"{file_prefix}_{M}.csv")
    else:
        output_dir = os.path.join(folder, f"r{R}")
        filename = os.path.join(output_dir, f"{file_prefix}_r{R}_{M}.csv")

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Pre-allocate array for all results
    all_results = np.zeros((2 * s, 7), dtype=np.float64)

    # Forward simulation
    for i in range(s):
        # Note: JIT version does full sweep per call, original needs N calls
        if use_jit:
            x.demon_move()  # JIT: full sweep
        else:
            for j in range(n):  # Original: N calls per sweep
                x.demon_move()

        # Send progress update (display refresh is throttled to 1/sec in main loop)
        progress_queue.put(
            {
                "type": "progress",
                "sim_num": sim_num,
                "sim_type": sim_type,
                "phase": "forward",
                "sweep": i + 1,
                "total_sweeps": s,
            }
        )

        # Get validated state after each sweep
        state = x.get_validated_state()

        # Calculate entropy with stable functions
        N0e = max(1, int(state["bond_count"][1]))

        try:
            Sk_val = Sk_stable(n, state["E_demon_sum"])
            Su_val = Su_stable(n, state["bond_count"][1], state["bond_count"][2], N0e)
            total_entropy = (Sk_val + Su_val) / n
        except (ValueError, OverflowError):
            total_entropy = np.nan

        # Store results
        all_results[i] = [
            np.float64(i + 1),
            np.float64(state["E_demon_sum"]) / n,
            np.float64(state["E_lattice"]) / n,
            np.float64(state["bond_count"][1]) / n,
            np.float64(state["bond_count"][2]) / n,
            total_entropy,
            np.float64(n),
        ]

    # Reverse simulation
    for i in range(s):
        # Note: JIT version does full sweep per call, original needs N calls
        if use_jit:
            x.demon_reverse()  # JIT: full sweep
        else:
            for j in range(n):  # Original: N calls per sweep
                x.demon_reverse()

        # Send progress update (display refresh is throttled to 1/sec in main loop)
        progress_queue.put(
            {
                "type": "progress",
                "sim_num": sim_num,
                "sim_type": sim_type,
                "phase": "reverse",
                "sweep": i + 1,
                "total_sweeps": s,
            }
        )

        # Get validated state
        state = x.get_validated_state()

        # Calculate entropy
        N0_exp = max(1, int(state["bond_count"][1]))

        try:
            Sk_val = Sk_stable(n, state["E_demon_sum"])
            Su_val = Su_stable(
                n, state["bond_count"][1], state["bond_count"][2], N0_exp
            )
            total_entropy = (Sk_val + Su_val) / n
        except (ValueError, OverflowError):
            total_entropy = np.nan

        # Store results
        all_results[s + i] = [
            np.float64(s + i + 1),
            np.float64(state["E_demon_sum"]) / n,
            np.float64(state["E_lattice"]) / n,
            np.float64(state["bond_count"][1]) / n,
            np.float64(state["bond_count"][2]) / n,
            total_entropy,
            np.float64(n),
        ]

    # Write results
    with open(filename, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["t", "K", "U", "N0", "Nx", "S/nk", "n"])
        writer.writerows(all_results)

    # Save metadata
    metadata = {
        "lattice_size": n,
        "sweeps": s,
        "radius": R,
        "run": M,
        "timestamp": datetime.now().isoformat(),
        "simulation_type": f"{sim_type}_parallel",
    }
    metadata_file = filename.replace(".csv", "_metadata.json")
    with open(metadata_file, "w") as f:
        json.dump(metadata, f, indent=2)

    # Send completion status
    progress_queue.put(
        {"type": "complete", "sim_num": sim_num, "R": R, "M": M, "sim_type": sim_type}
    )

    return {
        "R": R,
        "M": M,
        "sim_type": sim_type,
        "filename": filename,
        "success": True,
    }


def main():
    """Main execution function."""
    # Parse arguments
    parser = create_argument_parser("combined")
    args = parser.parse_args()

    # Simulation parameters
    n = args.n  # lattice size
    s = args.s  # sweeps
    r = args.r  # max bond-demon couple radius
    m = args.m  # number of sims per type
    validate_mode = args.validate
    use_jit = not args.no_jit

    # Determine number of cores
    if args.cores is None:
        num_cores = mp.cpu_count()
    else:
        num_cores = min(args.cores, mp.cpu_count())

    total_sims = r * m * 2  # Double for both types

    # Print simulation info
    print_simulation_info("combined", n, s, r, m, num_cores, validate_mode, use_jit)

    # Setup logging and project root
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    log_dir = setup_logging(project_root, "combined")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"parallel_sim_combined_{timestamp}.log")

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
    )

    logger = logging.getLogger(__name__)
    logger.info("=" * 80)
    logger.info("COMBINED PARALLEL SIMULATION (Reversible + Irreversible)")
    logger.info("=" * 80)
    logger.info(
        f"Starting combined simulation: n={n}, s={s}, r={r}, m={m}, cores={num_cores}"
    )

    # Create data directories for both types
    create_data_directories(project_root, r, "combined")

    # Build task list - interleave both simulation types
    sim_params = []
    sim_num = 1

    # For each radius and run, create tasks for BOTH reversible and irreversible
    for R in range(r):
        for M in range(m):
            # Reversible simulation (placeholder None queue replaced later)
            sim_params.append(
                (
                    R,
                    M,
                    n,
                    s,
                    validate_mode,
                    project_root,
                    use_jit,
                    sim_num,
                    "reversible",
                    None,
                )
            )
            sim_num += 1

            # Irreversible simulation
            sim_params.append(
                (
                    R,
                    M,
                    n,
                    s,
                    validate_mode,
                    project_root,
                    use_jit,
                    sim_num,
                    "irreversible",
                    None,
                )
            )
            sim_num += 1

    logger.info(
        f"Total simulations to run: {total_sims} "
        f"({total_sims // 2} reversible + {total_sims // 2} irreversible)"
    )
    logger.info(f"Parameters per simulation: R=[0-{r-1}], M=[0-{m-1}], n={n}, s={s}")

    # Start timer
    start_time = datetime.now()

    # Create manager and progress queue
    manager = Manager()
    progress_queue = manager.Queue()

    # Update sim_params with the actual queue
    sim_params = [
        (R, M, n, s, val, root, jit, num, stype, progress_queue)
        for R, M, n, s, val, root, jit, num, stype, _ in sim_params
    ]

    print(f"\nStarting {total_sims} simulations on {num_cores} cores...\n")

    # Track active simulations and completion
    active_sims = {}
    completed = 0
    sim_times = []

    with mp.Pool(processes=num_cores) as pool:
        # Start async processing
        result = pool.map_async(run_single_simulation, sim_params)

        # Monitor progress
        pbar = tqdm(
            total=total_sims,
            unit="sim",
            bar_format="{desc}: {percentage:3.0f}%|{bar}| {n}/{total} [{elapsed}]",
        )
        pbar.set_description("Overall Progress")
        last_refresh = time.time()

        # Process messages from queue
        try:
            while not result.ready() or not progress_queue.empty():
                try:
                    msg = progress_queue.get(timeout=0.1)

                    if msg["type"] == "start":
                        active_sims[msg["sim_num"]] = {
                            "R": msg["R"],
                            "M": msg["M"],
                            "sim_type": msg.get("sim_type", "unknown"),
                            "start_time": time.time(),
                            "phase": "forward",
                            "progress": 0,
                        }

                    elif msg["type"] == "progress":
                        sim_num = msg["sim_num"]
                        if sim_num in active_sims:
                            active_sims[sim_num]["phase"] = msg["phase"]
                            active_sims[sim_num]["progress"] = msg["sweep"]

                        # Update display once per second
                        now = time.time()
                        if now - last_refresh >= 1.0:
                            # Count active by type
                            rev_active = sum(
                                1
                                for s in active_sims.values()
                                if s.get("sim_type") == "reversible"
                            )
                            irr_active = sum(
                                1
                                for s in active_sims.values()
                                if s.get("sim_type") == "irreversible"
                            )
                            active_count = len(active_sims)
                            desc = (
                                f"Progress ({active_count} active: "
                                f"{rev_active} rev, {irr_active} irr)"
                            )
                            pbar.set_description(desc)
                            last_refresh = now

                    elif msg["type"] == "complete":
                        sim_num = msg["sim_num"]
                        if sim_num in active_sims:
                            elapsed_time = (
                                time.time() - active_sims[sim_num]["start_time"]
                            )
                            sim_times.append(elapsed_time)
                            del active_sims[sim_num]
                        completed += 1
                        pbar.update(1)

                except Exception:
                    # Queue timeout - check if we're done
                    if result.ready() and progress_queue.empty():
                        break

        except KeyboardInterrupt:
            logger.warning("Received interrupt signal, terminating...")
            pool.terminate()
            pool.join()
            raise

        # Get results
        results = result.get()
        pbar.close()

    end_time = datetime.now()
    elapsed = (end_time - start_time).total_seconds()

    # Print summary
    print("\n" + "=" * 80)
    print("SIMULATION COMPLETE")
    print("=" * 80)
    print(f"Total simulations: {total_sims}")
    print(f"Time elapsed: {elapsed:.2f} seconds")
    if sim_times:
        avg_time = sum(sim_times) / len(sim_times)
        print(f"Average time per simulation: {avg_time:.2f} seconds")
    success_count = len([r for r in results if r.get("success", False)])
    print(f"Success rate: {success_count} / {total_sims}")

    logger.info("=" * 80)
    logger.info("COMBINED SIMULATION COMPLETE")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
