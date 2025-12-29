"""
Copyright (c) 2026 Winry Ember
Licensed under the MIT License.
See LICENSE file in the project root for full license information.

Parallel implementation of irreversible Creutz demon simulation.

Uses multiprocessing to run multiple independent simulations simultaneously,
significantly reducing total execution time on multi-core systems.

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

# Dynamic import based on JIT (default on, use --no-jit to disable)
import numpy as np
from sim_utils import (
    Sk_stable,
    Su_stable,
    build_simulation_parameters,
    create_argument_parser,
    create_data_directories,
    print_final_results,
    print_simulation_info,
    process_message_queue,
    setup_logging,
)
from tqdm import tqdm


def run_single_simulation(
    args: Tuple[int, int, int, int, str, str, bool, int, "mp.Queue"],
) -> Dict:
    """Run a single simulation for given parameters.

    Args:
        args: Tuple of (R, M, n, s, validate_mode, project_root, use_jit, sim_num,
                progress_queue)
            R: Demon-coupling radius
            M: Run number
            n: Lattice size
            s: Number of sweeps
            validate_mode: Validation mode ('off', 'periodic', 'frequent')
            project_root: Project root directory
            use_jit: Whether to use JIT-compiled version
            sim_num: Simulation number (1-based)
            progress_queue: Queue for progress updates

    Returns:
        Dictionary with simulation results and metadata
    """
    R, M, n, s, validate_mode, project_root, use_jit, sim_num, progress_queue = args

    # Send initial status
    progress_queue.put({"type": "start", "sim_num": sim_num, "R": R, "M": M})

    # Import appropriate irrInferno class
    if use_jit:
        from jit_inferno_irr import JITirrInferno as irrInferno
    else:
        from inferno_irr import irrInferno

    # Create irrInferno instance
    x = irrInferno(n, R + 1, validate_mode=validate_mode)

    # Setup output paths
    folder = os.path.join(project_root, "data", "irr")
    if R == 0:
        output_dir = os.path.join(folder, "r0")
        filename = os.path.join(output_dir, f"irr_sim_data_{M}.csv")
    else:
        output_dir = os.path.join(folder, f"r{R}")
        filename = os.path.join(output_dir, f"irr_sim_data_r{R}_{M}.csv")

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
        "simulation_type": "irreversible_parallel",
    }
    metadata_file = filename.replace(".csv", "_metadata.json")
    with open(metadata_file, "w") as f:
        json.dump(metadata, f, indent=2)

    # Send completion status
    progress_queue.put({"type": "complete", "sim_num": sim_num, "R": R, "M": M})

    return {
        "R": R,
        "M": M,
        "filename": filename,
        "E_total": state["E_total"],
        "E_initial": x._initial_total_energy,
    }


if __name__ == "__main__":
    # Parse command-line arguments
    parser = create_argument_parser("irreversible")
    args = parser.parse_args()

    # Simulation parameters
    n = args.n  # lattice size
    s = args.s  # sweeps
    r = args.r  # max bond-demon couple radius
    m = args.m  # number of sims
    validate_mode = args.validate
    use_jit = not args.no_jit

    # Determine number of cores
    if args.cores is None:
        num_cores = mp.cpu_count()
    else:
        num_cores = min(args.cores, mp.cpu_count())

    total_sims = r * m

    # Print simulation info
    print_simulation_info("irreversible", n, s, r, m, num_cores, validate_mode, use_jit)

    # Set up logging
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    log_dir = setup_logging(project_root, "irreversible")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"parallel_sim_irreversible_{timestamp}.log")

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
    )

    logging.info(
        f"Starting parallel irreversible simulation: n={n}, s={s}, r={r}, "
        f"m={m}, cores={num_cores}"
    )

    # Create data directories and build parameter list
    create_data_directories(project_root, r, "irreversible")
    sim_params = build_simulation_parameters(
        r, m, n, s, validate_mode, project_root, use_jit
    )

    # Run simulations in parallel with progress monitoring
    start_time = datetime.now()

    # Create manager and progress queue
    manager = Manager()
    progress_queue = manager.Queue()

    # Update sim_params with the actual queue
    sim_params = [
        (R, M, n, s, val, root, jit, num, progress_queue)
        for R, M, n, s, val, root, jit, num, _ in sim_params
    ]

    print(f"\nStarting {total_sims} simulations on {num_cores} cores...\n")

    # Track active simulations and completion
    active_sims = {}
    completed = 0
    sim_times = []

    with mp.Pool(processes=num_cores) as pool:
        # Start async processing
        result = pool.map_async(run_single_simulation, sim_params)

        # Monitor progress (manually update desc to avoid tqdm's postfix formatting)
        pbar = tqdm(
            total=total_sims,
            unit="sim",
            bar_format="{desc}: {percentage:3.0f}%|{bar}| {n}/{total} [{elapsed}]",
        )
        pbar.set_description("Overall Progress")
        last_refresh = time.time()

        # Process messages and handle interrupts
        completed, last_refresh, results = process_message_queue(
            progress_queue,
            active_sims,
            s,
            pbar,
            last_refresh,
            sim_times,
            total_sims,
            num_cores,
            completed,
            result,
            pool,
        )

    end_time = datetime.now()
    print_final_results(results, start_time, end_time, "irreversible")
