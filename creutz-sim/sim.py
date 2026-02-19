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

def add_row(filename, row_data):    # appends a new row to csv file
    try:
        with open(filename, 'a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(row_data)
    except Exception as e:
         print(f"An error occurred: {e}")

Sk = lambda N, K: logg(K + N) - logg(K+1) - logg(N) # N == lattice size, K == kinetic energy
Su = lambda N, N0, Nx, N0_exp: logg(N+1) + math.log(2**N0_exp) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1)) # N == lattice size, N0 == broken bonds, Nx == bonds between anti-aligned spins


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
    
    args = parser.parse_args()
    
    # Default values
    defaults = {
        'n': 10,        # lattice size
        's': 100,       # sweeps
        'flag': 'c',    # c=combined, r=reversible, i=irreversible
        'r': 10,        # max radius, inclusive
        'm': 5
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
    
    print(f"Running simulation: n={defaults['n']}, sweeps={defaults['s']}, flag={defaults['flag']}, radius={defaults['r']}, runs={defaults['m']}")
    print()
    
    return defaults['n'], defaults['s'], defaults['flag'], defaults['r'], defaults['m'], args.data_dir


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

def run_radius_simulations(R, n, s, flag, m, file_names, irr_files, init_files, pbar_queue):
    """
    Worker function to run all simulations for a given radius R.
    
    This function is executed in parallel by multiprocessing.Pool workers. Each worker
    handles all M runs for a single R value (radius). The pool automatically distributes
    R values across available CPU cores and queues remaining jobs.
    
    Architecture:
    - Each worker is independent (no shared state except progress queue)
    - Workers ignore SIGINT - main process handles Ctrl-C interrupts
    - Progress updates sent via shared queue for real-time progress bar
    - Each (R, M) combination writes to independent CSV files
    
    Args:
        R: Radius value (0 to r) for bond-demon coupling
        n: Lattice size
        s: Number of sweeps (s//2 forward + s//2 reverse)
        flag: Dynamics type ('c'=combined, 'r'=reversible, 'i'=irreversible)
        m: Number of independent runs per radius
        file_names: List of file paths for reversible dynamics output
        irr_files: List of file paths for irreversible dynamics output
        init_files: List of file paths for initial/final state snapshots
        pbar_queue: Shared queue for sending progress updates to main process
    
    Returns:
        int: Total number of completed sweeps (for validation)
    """
    # Ignore keyboard interrupts in worker processes - main process handles them
    import signal
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    
    completed_count = 0
    x = Inferno(n, R)
    for M in range(m):
        # Determine which dynamics to run based on flag
        if flag == 'c':  # combined
            dynamics_flags = [0, 1]
        elif flag == 'r':  # reversible only
            dynamics_flags = [0]
        else:  # 'i' - irreversible only
            dynamics_flags = [1]

        for dynamics_flag in dynamics_flags:
            filename = file_names[R].parent / f"{file_names[R].name}_{M}.csv"
            irr_filename = irr_files[R].parent / f"{irr_files[R].name}_{M}.csv"
            init_filename = init_files[R].parent / f"{init_files[R].name}_{M}.csv"
            filenames = [filename, irr_filename, init_filename]

            data_types = ['t', 'K', 'U', 'N0', 'Nx', 'S/nk', 'n'] # step counter, lattice energy, demon energy, total energy, broken bonds, anti-aligned spins, lattice size

            if dynamics_flag == 0:  # rev only
                filenames = [filename]
            elif dynamics_flag:  # irr only
                filenames = [irr_filename]
            
            # Create directories if they don't exist
            for fname in filenames:
                fname.parent.mkdir(parents=True, exist_ok=True)
            
            for fname in filenames:
                with open(fname, 'w+', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow(data_types)

            # Separate counters for each file type
            t_counter = 0
            t_irr_counter = 0

            ### Forward simulation
            for i in range(s//2):
                # flag = (i // (n//k)) % 2    # if 0, perform reversible dynamics
                # Attempt to flip each spin in lattice (full sweep)
                for j in range(n):
                    x.demon_move(dynamics_flag, i)
                
                # Calculate entropy and data ONCE per sweep (after all N moves)
                N0e = int(x.bond_count[1])
                if N0e == 0:
                    N0e = 1
                total_entropy = (Sk(n, x.d_energy) + Su(n, x.bond_count[1], x.bond_count[2], N0e))/n
                
                # Increment appropriate counter
                if dynamics_flag == 0:
                    t_counter += 1
                    t_value = t_counter
                else:
                    t_irr_counter += 1
                    t_value = t_irr_counter

                # write avg sweep results to csv
                new_row = np.array([t_value, x.d_energy/n, x.E_lattice/n, x.bond_count[1]/n, x.bond_count[2]/n, total_entropy, n])

                # Save initial and final states to init_filename, all states to appropriate file
                # if (i == 0):
                #     add_row(init_filename, new_row)

                if (dynamics_flag == 0):
                    add_row(filename, new_row)
                else:
                    add_row(irr_filename, new_row)
                
                # Signal progress to main process
                if pbar_queue:
                    pbar_queue.put(1)
                completed_count += 1

            ### Reverse simulation
            for i in range(s//2):
                # flag = (i // (n//k)) % 2    # if 0, perform reversible dynamics
                # Attempt to flip each spin in lattice (full sweep)
                for j in range(n):
                    x.demon_reverse(dynamics_flag, (s//2) - 1 - i)
                
                # Calculate entropy and data ONCE per sweep (after all N moves)
                N0_exp = int(x.bond_count[1])
                if N0_exp == 0:
                    N0_exp = 1
                total_entropy = (Sk(n, x.d_energy) + Su(n, x.bond_count[1], x.bond_count[2], N0_exp))/n
                
                # Increment appropriate counter
                if dynamics_flag == 0:
                    t_counter += 1
                    t_value = t_counter
                else:
                    t_irr_counter += 1
                    t_value = t_irr_counter

                # write avg sweep results to csv
                new_row = np.array([t_value, x.d_energy/n, x.E_lattice/n, x.bond_count[1]/n, x.bond_count[2]/n, total_entropy, n])

                # Save initial and final states to init_filename, all states to appropriate file
                # if (i == (s//2-1)):
                #     add_row(init_filename, new_row)

                if (dynamics_flag == 0):
                    add_row(filename, new_row)
                else:
                    add_row(irr_filename, new_row)
                
                # Signal progress to main process
                if pbar_queue:
                    pbar_queue.put(1)
                completed_count += 1

            # Reset "random" order for next simulation
            x.reset()
    
    return completed_count


if __name__ == '__main__':
    # CRITICAL: All initialization must be inside this block!
    # On macOS/Windows, multiprocessing uses 'spawn' method which imports this module
    # in each worker process. Code outside this block would execute in every worker,
    # causing duplicate output, incorrect state, and resource conflicts.
    
    # Get parameters
    n, s, flag, r, m, data_dir = get_params()

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
    
    print(f"Starting {total_sims} simulations ({r+1} radii × {m} runs)")
    print(f"Using {num_cores} CPU cores for parallel processing")
    print()

    # Use relative path from repo root
    repo_root = Path(__file__).parent.parent
    data_root = repo_root / data_dir
    init_folder = data_root / 'init_fin'

    # Create timestamped folder for this run directly in data directory
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    data_root.mkdir(parents=True, exist_ok=True)
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

    def update_progress_time():
        """Update progress bar with elapsed/remaining time"""
        elapsed = time.time() - pbar_start_time
        rate = pbar.n / elapsed if elapsed > 0 else 0
        remaining = (total_sweeps - pbar.n) / rate if rate > 0 else 0
        pbar.set_postfix_str(f"[elapsed {format_time(elapsed)}|remaining {format_time(remaining)}|{rate:.2f}sweep/s]", refresh=True)

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
    # Workers send progress updates (1 per sweep) through this queue.
    # Main process consumes updates and increments progress bar in real-time.
    # Queue is thread-safe and works across process boundaries.
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
        # Submit all (r+1) jobs to pool. Pool distributes across workers automatically.
        # map_async returns immediately - workers run in background.
        results = pool.map_async(worker_func, range(r+1))
        
        # Monitor progress while workers run
        # Each worker puts '1' in queue after completing a sweep.
        # We consume the queue and update the progress bar in real-time.
        completed_sweeps = 0
        while not results.ready():  # Loop until all workers finish
            try:
                # Wait 0.1s for progress update from any worker
                pbar_queue.get(timeout=0.1)
                completed_sweeps += 1
                pbar.update(1)
                update_progress_time()
            except Exception:  # Must be Exception, not bare except:!
                # Queue empty or timeout - KeyboardInterrupt must propagate up!
                pass
        
        # Get any remaining progress updates
        while not pbar_queue.empty():
            try:
                pbar_queue.get_nowait()
                completed_sweeps += 1
                pbar.update(1)
                update_progress_time()
            except Exception:
                break
        
        # Ensure progress bar is at 100%
        if pbar.n < total_sweeps:
            pbar.update(total_sweeps - pbar.n)
        
        # Get results (blocks until all complete)
        sweep_counts = results.get()
        
        # Close pool normally
        pool.close()
        pool.join()
        
    except KeyboardInterrupt:
        # Ctrl-C pressed - clean shutdown
        # Workers ignore SIGINT (set in worker function), so only main process catches it.
        # terminate() sends SIGTERM to all workers for immediate forced shutdown.
        # join() waits for workers to clean up before exiting.
        print("\n\nInterrupted! Terminating...")
        pbar.close()       # Close progress bar first to clear terminal
        pool.terminate()   # Force kill all worker processes
        pool.join()        # Wait for cleanup
        
        # Write interrupted status
        with open(status_file, 'w') as f:
            f.write(f"Status: INTERRUPTED\n")
            f.write(f"Time: {datetime.now().isoformat()}\n")
        print(f"Status written to {status_file}")
        sys.exit(1)

    # Close progress bar
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
    print(f"  Throughput: {total_sweeps/total_time:.0f} sweeps/sec")

    # Write completion status
    with open(completion_marker, 'w') as f:
        f.write(f"Simulation completed: {datetime.now().isoformat()}\n")
        f.write(f"Total time: {time_str}\n")
        f.write(f"Average: {avg_time_per_sim:.2f}s per simulation\n")
        f.write(f"Throughput: {total_sweeps/total_time:.0f} sweeps/sec\n")

    with open(status_file, 'w') as f:
        f.write(f"Status: COMPLETED\n")
        f.write(f"Time: {datetime.now().isoformat()}\n")
        f.write(f"Completed: {total_sims}/{total_sims} simulations\n")

    print(f"Status written to {status_file}")