# Force single-threaded BLAS/LAPACK on Linux HPC to prevent thread over-subscription
# macOS/Apple Silicon has good thread management, so skip there
import os
import platform
if platform.system() == 'Linux':
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
    os.environ['NUMEXPR_NUM_THREADS'] = '1'

from inferno import Inferno
import numpy as np
import csv
from scipy.special import loggamma as logg
import math
from datetime import datetime
from pathlib import Path
import argparse
from tqdm import tqdm
import time
import sys
import multiprocessing as mp
from functools import partial
import errno
from inferno import Inferno, Sk, Su, Su0

def add_row(filename, row_data):    # appends a new row to csv file
    try:
        with open(filename, 'a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(row_data)
    except Exception as e:
         print(f"An error occurred: {e}")


def get_params():
    """Get simulation parameters from command line args or interactively."""
    parser = argparse.ArgumentParser(description='Monte Carlo simulation for statistical mechanics')
    parser.add_argument('-n', '--lattice-size', type=int, help='Lattice size (default: 10)')
    parser.add_argument('-s', '--sweeps', type=int, help='Number of sweeps (default: 100)')
    parser.add_argument('-f', '--flag', type=str, choices=['c', 'r', 'i'], help='Dynamics flag: c=combined, r=reversible, i=irreversible (default: c)')
    parser.add_argument('-r', '--radius', type=int, help='Max bond-demon couple radius (default: 11)')
    parser.add_argument('-m', '--runs', type=int, help='Number of simulation runs (default: 5)')
    parser.add_argument('-d', '--data-dir', type=str, default='data', help='Data output directory (default: data)')
    parser.add_argument('-i', '--interactive', action='store_true', help='Interactive mode')
    parser.add_argument('-np', '--no-pbar', action='store_true', help='Disable progress bar display')
    
    args = parser.parse_args()
    
    # Default values
    defaults = {
        'n': 10,        # lattice size
        's': 100,       # sweeps
        'flag': 'c',    # c=combined, r=reversible, i=irreversible
        'r': 10,        # max radius, inclusive
        'm': 6
    }
    
    # Interactive mode (only if -i flag is set)
    if args.interactive:
        print("=== Monte Carlo Simulation Configuration ===")
        print()
        
        try:
            n = input(f"Lattice size (n) [{defaults['n']}]: ").strip()
            defaults['n'] = int(n) if n else defaults['n']
            
            s = input(f"Number of sweeps (s) [{defaults['s']}]: ").strip()
            defaults['s'] = int(s) if s else defaults['s']
            
            f = input(f"Dynamics (c=combined, r=reversible, i=irreversible) [{defaults['flag']}]: ").strip()
            defaults['flag'] = f if f in ['c', 'r', 'i'] else defaults['flag']
            
            r = input(f"Max radius (r) [{defaults['r']}]: ").strip()
            defaults['r'] = int(r) if r else defaults['r']
            
            m = input(f"Number of runs (m) [{defaults['m']}]: ").strip()
            defaults['m'] = int(m) if m else defaults['m']
            
            print()
        except KeyboardInterrupt:
            print("\n\nSimulation cancelled by user.")
            sys.exit(0)
        except ValueError:
            print("\nInvalid input. Using default values")
    else:
        # Use command line args if provided
        if args.lattice_size is not None:
            defaults['n'] = args.lattice_size
        if args.sweeps is not None:
            defaults['s'] = args.sweeps
        if args.flag is not None:
            defaults['flag'] = args.flag
        if args.radius is not None:
            defaults['r'] = args.radius
        if args.runs is not None:
            defaults['m'] = args.runs
    
    print(f"Running simulation: n={defaults['n']:,}, sweeps={defaults['s']:,}, flag={defaults['flag']}, radius={defaults['r']:,}, runs={defaults['m']:,}")
    print()
    
    return (
        defaults['n'], defaults['s'], defaults['flag'], defaults['r'], defaults['m'],
        args.data_dir, not args.no_pbar
    )


def format_time(seconds):
    """Format time in days/hours/minutes"""
    if seconds >= 172800:  # More than 48 hours (2 days)
        return f"{seconds / 86400:.1f}d"
    elif seconds >= 3600:  # More than 1 hour
        return f"{seconds / 3600:.1f}h"
    elif seconds >= 60:  # More than 1 minute
        return f"{seconds / 60:.0f}m"
    else:
        return f"{seconds:.1f}s"


######################################################################################
##################################### begin sim ######################################
######################################################################################

def run_one_phase(x, dynamics_flag, active_file, s, n, t_counter, pbar_queue, PBAR_BATCH):
    """
    Run one full simulation phase (s//2 forward + s//2 reverse sweeps) for a given
    dynamics type, writing results to active_file. State carries over in-place on x.

    Every (s * 0.01) sweeps are accumulated and written as a single averaged row,
    reducing output size by 100x while preserving the overall trajectory.

    Args:
        x: Inferno instance (mutated in-place)
        dynamics_flag: 0 = reversible, 1 = irreversible
        active_file: Path to output CSV file (opened in append mode)
        s: Total sweeps (s//2 forward + s//2 reverse)
        n: Lattice size
        t_counter: Current sweep counter (continues across phases)
        pbar_queue: Shared queue for progress updates
        PBAR_BATCH: Batch size for progress IPC

    Returns:
        t_counter: Updated sweep counter after this phase
    """
    # Number of sweeps to average into one CSV row (1% of total sweeps, min 1)
    avg_window = max(1, int(s * 0.01))

    def _run_half(sweep_range, move_fn):
        nonlocal t_counter
        pbar_accum = 0
        acc_N0      = 0.0
        acc_Nx      = 0.0
        acc_S       = 0.0
        acc_E_lat   = 0.0
        acc_E_demon = 0.0
        acc_count   = 0

        with open(active_file, 'a', newline='') as f:
            writer = csv.writer(f)

            for i in sweep_range:
                move_fn(i)

                N0 = int(x.bond_count[1])
                Nx = int(x.bond_count[2])

                if N0 == 0:
                    S_conf = Su0(n, N0, Nx)
                else:
                    S_conf = Su(n, N0, Nx)

                # Store only configurational entropy (per site)
                acc_N0      += x.bond_count[1] / n * 100
                acc_Nx      += x.bond_count[2] / n * 100
                acc_S       += S_conf / n
                acc_E_lat   += x.E_lattice
                acc_E_demon += x.d_energy
                acc_count   += 1
                t_counter += 1

                if pbar_queue:
                    pbar_accum += 1
                    if pbar_accum >= PBAR_BATCH:
                        pbar_queue.put(pbar_accum)
                        pbar_accum = 0

                # Write averaged row every avg_window sweeps
                if acc_count >= avg_window:
                    h = x.d_energy_hist
                    ratios = [h[k] / h[k + 1] if h[k + 1] > 0 else 0 for k in range(4)]
                    writer.writerow([t_counter,
                                     acc_N0    / acc_count,
                                     acc_Nx    / acc_count,
                                     acc_S     / acc_count,
                                     acc_E_lat   / acc_count,
                                     acc_E_demon / acc_count,
                                     n,
                                     *ratios])
                    acc_N0 = acc_Nx = acc_S = acc_E_lat = acc_E_demon = 0.0
                    acc_count = 0
                    x.d_energy_hist[:] = 0

            # Flush any remaining accumulated sweeps
            if acc_count > 0:
                h = x.d_energy_hist
                ratios = [h[k] / h[k + 1] if h[k + 1] > 0 else 0 for k in range(4)]
                writer.writerow([t_counter,
                                 acc_N0    / acc_count,
                                 acc_Nx    / acc_count,
                                 acc_S     / acc_count,
                                 acc_E_lat   / acc_count,
                                 acc_E_demon / acc_count,
                                 n,
                                 *ratios])
                x.d_energy_hist[:] = 0

            if pbar_queue and pbar_accum:
                pbar_queue.put(pbar_accum)

    ### Forward sweeps
    def fwd_move(i):
        for _ in range(n):
            x.demon_move(dynamics_flag)

    _run_half(range(s // 2), fwd_move)
    x.d_energy_hist[:] = 0

    ### Reverse sweeps
    half = s // 2

    def rev_move(i):
        for _ in range(n):
            x.demon_reverse(dynamics_flag)

    _run_half(range(half), rev_move)

    return t_counter



def run_radius_simulations(R, n, s, flag, m, file_names, irr_files, init_files, pbar_queue):
    """
    Worker function to run all simulations for a given radius R.

    For flag='c' (combined): initializes one Inferno per radius, then for each of m
    runs alternates a full reversible phase followed by a full irreversible phase,
    with state carrying over continuously throughout. Rev and irr data are written
    to separate files indexed by M.

    For flag='r' or 'i': runs m independent simulations of the single dynamics type,
    resetting state between runs (original behaviour).

    Args:
        R: Radius value (0 to r) for bond-demon coupling
        n: Lattice size
        s: Number of sweeps (s//2 forward + s//2 reverse per phase)
        flag: Dynamics type ('c'=combined, 'r'=reversible, 'i'=irreversible)
        m: Number of rev/irr alternation cycles (combined) or independent runs
        file_names: List of file paths for reversible dynamics output
        irr_files: List of file paths for irreversible dynamics output
        init_files: List of file paths for initial/final state snapshots
        pbar_queue: Shared queue for sending progress updates to main process

    Returns:
        int: Total number of completed sweeps (for validation)
    """
    import signal
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    PBAR_BATCH = max(50, 500_000 // max(n, 1))
    data_types = ['t', 'N0 (%)', 'Nx (%)', 'S_conf/n', 'E_lattice', 'E_demon', 'n',
                  'p0/p1', 'p1/p2', 'p2/p3', 'p3/p4']

    completed_count = 0

    if flag == 'c':
        # One Inferno for the entire radius; state carries over across all phases.
        x = Inferno(n, R)

        for M in range(m):           # generate fresh order once per M
            t_rev = 0   # reset sweep counters at the start of each cycle
            t_irr = 0
            rev_filename  = file_names[R].parent / f"{file_names[R].name}_{M}.csv"
            irr_filename  = irr_files[R].parent  / f"{irr_files[R].name}_{M}.csv"

            # Create directories
            for fname in [rev_filename, irr_filename]:
                fname.parent.mkdir(parents=True, exist_ok=True)

            # Write headers (fresh file each M)
            for fname in [rev_filename, irr_filename]:
                with open(fname, 'w', newline='') as f:
                    csv.writer(f).writerow(data_types)

            # --- Reversible phase ---
            t_rev = run_one_phase(x, 0, rev_filename, s, n, t_rev, pbar_queue, PBAR_BATCH)
            x.reset()
            
            # --- Irreversible phase ---
            t_irr = run_one_phase(x, 1, irr_filename, s, n, t_irr, pbar_queue, PBAR_BATCH)
            x.reset()

    else:
        # Single-dynamics mode: m independent runs, reset between each (original behaviour).
        dynamics_flag = 0 if flag == 'r' else 1

        x = Inferno(n, R)
        for M in range(m):
            active_file = (file_names[R] if dynamics_flag == 0 else irr_files[R])
            filename = active_file.parent / f"{active_file.name}_{M}.csv"
            filename.parent.mkdir(parents=True, exist_ok=True)

            with open(filename, 'w', newline='') as f:
                csv.writer(f).writerow(data_types)

            run_one_phase(x, dynamics_flag, filename, s, n, 0, pbar_queue, PBAR_BATCH)
            x.reset()

    return completed_count


if __name__ == '__main__':
    # CRITICAL: All initialization must be inside this block!
    # On macOS/Windows, multiprocessing uses 'spawn' method which imports this module
    # in each worker process. Code outside this block would execute in every worker,
    # causing duplicate output, incorrect state, and resource conflicts.
    
    # Get parameters
    n, s, flag, r, m, data_dir, show_pbar = get_params()

    # Calculate total iterations for progress tracking
    # Note: r is the max radius, but loop goes from 0 to r inclusive, so (r+1) radii total
    total_sims = (r + 1) * m
    total_sweeps_per_sim = s  # s//2 forward + s//2 reverse

    # For combined runs, we do both reversible and irreversible, so double the sweeps
    if flag == 'c':
        total_sweeps = total_sims * s * 2  # Both rev and irr
    else:
        total_sweeps = total_sims * s  # Just one dynamics type

    # Parallelization strategy: Distribute R values (radii) across CPU cores
    # Each worker handles all M runs for one R value, writing to independent files.
    # Pool automatically queues remaining jobs when cores < (r+1).
    # Example: 8 cores, 11 radii → 8 run immediately, 3 queued, grabbed as workers finish.
    num_cores = min(mp.cpu_count(), r + 1)  # Don't create more workers than jobs
    
    print(f"Starting {total_sims:,} simulations ({r+1:,} radii × {m:,} runs)")
    print(f"Using {num_cores} CPU cores for parallel processing")
    print()

    # Use relative path from repo root
    repo_root = Path(__file__).parent.parent
    data_root = repo_root / data_dir
    init_folder = data_root / 'init_fin'

    # Ensure data root exists before writing run outputs.
    data_root.mkdir(parents=True, exist_ok=True)

    # Create timestamped folder for this run directly in data directory
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    data_folder = data_root / timestamp
    data_folder.mkdir(parents=True, exist_ok=True)

    # Determine status file location based on flag
    if flag == 'r':  # reversible only
        status_folder = data_folder / 'rev'
    elif flag == 'i':  # irreversible only
        status_folder = data_folder / 'irr'
    else:  # combined
        status_folder = data_folder
    status_folder.mkdir(parents=True, exist_ok=True)

    # Status tracking files
    status_file = status_folder / 'sim_status.txt'
    start_marker = status_folder / 'sim_started.txt'
    completion_marker = status_folder / 'sim_completed.txt'

    # Write simulation start marker with parameters
    with open(start_marker, 'w') as f:
        f.write(f"Simulation started: {datetime.now().isoformat()}\n")
        f.write(f"Parameters: n={n}, sweeps={s}, flag={flag}, radius={r}, runs={m}\n")
        f.write(f"Total simulations: {total_sims}\n")

    file_names = [data_folder / 'rev' / f'r{i}' / f'sim_data_r{i}' for i in range(r+1)]
    irr_files = [data_folder / 'irr' / f'r{i}' / f'irr_sim_data_r{i}' for i in range(r+1)]
    init_files = [init_folder / f'r{i}' / f'sim_data_r{i}' for i in range(r+1)]

    # Progress tracking
    start_time = time.time()

    def write_error_status(error_message):
        with open(status_file, 'w') as f:
            f.write("Status: ERROR\n")
            f.write(f"Time: {datetime.now().isoformat()}\n")
            f.write(f"Error: {error_message}\n")
            f.write(f"Completed sweeps: {pbar.n if pbar is not None else 0}/{total_sweeps}\n")

    pbar = None
    manager = None
    pbar_queue = None

    if show_pbar:
        def update_progress_time():
            """Update progress bar with elapsed/remaining time"""
            elapsed = time.time() - pbar_start_time
            rate = pbar.n / elapsed if elapsed > 0 else 0
            remaining = (total_sweeps - pbar.n) / rate if rate > 0 else 0
            pbar.set_postfix_str(f"[elapsed {format_time(elapsed)}|remaining {format_time(remaining)}|{rate:,.0f}sweep/s]", refresh=True)

        # Create progress bar with custom format and dynamic width
        pbar = tqdm(total=total_sweeps, desc="Progress", unit="sweep",
                    bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} {postfix}',
                    dynamic_ncols=True,  # Dynamically adjust to terminal width
                    mininterval=0.1,  # Update every 0.1 seconds
                    maxinterval=1.0,  # Force update at least every second
                    leave=True)  # Keep the bar visible after completion

        # Track for custom time display
        pbar_start_time = time.time()

        # Create a manager and queue for inter-process communication
        manager = mp.Manager()
        pbar_queue = manager.Queue()
    
    # Create worker function with fixed parameters using partial application
    worker_func = partial(run_radius_simulations, 
                         n=n, s=s, flag=flag, m=m,
                         file_names=file_names, 
                         irr_files=irr_files,
                         init_files=init_files,
                         pbar_queue=pbar_queue)
    
    # Create process pool and submit all R values
    pool = mp.Pool(processes=num_cores)
    
    try:
        results = pool.map_async(worker_func, range(r+1))
        
        if show_pbar:
            completed_sweeps = 0
            while not results.ready():
                try:
                    n_done = pbar_queue.get(timeout=0.1)
                    completed_sweeps += n_done
                    pbar.update(n_done)
                    update_progress_time()
                except Exception:
                    pass

            while not pbar_queue.empty():
                try:
                    n_done = pbar_queue.get_nowait()
                    completed_sweeps += n_done
                    pbar.update(n_done)
                    update_progress_time()
                except Exception:
                    break

            if pbar.n < total_sweeps:
                pbar.update(total_sweeps - pbar.n)
        else:
            results.wait()
        
        sweep_counts = results.get()
        pool.close()
        pool.join()
        
    except KeyboardInterrupt:
        print("\n\nInterrupted! Terminating...")
        if pbar is not None:
            pbar.close()
        pool.terminate()
        pool.join()
        
        with open(status_file, 'w') as f:
            f.write(f"Status: INTERRUPTED\n")
            f.write(f"Time: {datetime.now().isoformat()}\n")
        print(f"Status written to {status_file}")
        sys.exit(1)
    except OSError as e:
        if pbar is not None:
            pbar.close()
        pool.terminate()
        pool.join()

        if e.errno == errno.ENOSPC:
            msg = "No space left on device while writing simulation output"
            print(f"\n\nERROR: {msg}")
            print("Free disk space and re-run the simulation.")
        else:
            msg = f"OS error during simulation: {e}"
            print(f"\n\nERROR: {msg}")

        write_error_status(msg)
        print(f"Status written to {status_file}")
        sys.exit(1)
    except Exception as e:
        if pbar is not None:
            pbar.close()
        pool.terminate()
        pool.join()

        msg = f"Simulation failed: {e}"
        print(f"\n\nERROR: {msg}")
        write_error_status(msg)
        print(f"Status written to {status_file}")
        sys.exit(1)
    finally:
        if manager is not None:
            manager.shutdown()

    if pbar is not None:
        pbar.close()

    # Print summary
    total_time = time.time() - start_time
    if total_time < 60:
        time_str = f"{total_time:.1f}s"
    elif total_time < 3600:
        time_str = f"{total_time/60:.1f}m ({total_time:.0f}s)"
    else:
        time_str = f"{total_time/3600:.2f}h ({total_time/60:.1f}m)"

    avg_time_per_sim = total_time / total_sims
    print(f"\nSimulation complete!")
    print(f"  Total time: {time_str}")
    print(f"  Average: {avg_time_per_sim:.2f}s per simulation")
    print(f"  Throughput: {total_sweeps/total_time:,.0f} sweeps/sec")

    # Write completion status
    with open(completion_marker, 'w') as f:
        f.write(f"Simulation completed: {datetime.now().isoformat()}\n")
        f.write(f"Total time: {time_str}\n")
        f.write(f"Average: {avg_time_per_sim:.2f}s per simulation\n")
        f.write(f"Throughput: {total_sweeps/total_time:,.0f} sweeps/sec\n")

    with open(status_file, 'w') as f:
        f.write(f"Status: COMPLETED\n")
        f.write(f"Time: {datetime.now().isoformat()}\n")
        f.write(f"Completed: {total_sims}/{total_sims} simulations\n")

    print(f"Status written to {status_file}")