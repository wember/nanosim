"""Parallel implementation of irreversible Creutz demon simulation.

Uses multiprocessing to run multiple independent simulations simultaneously,
significantly reducing total execution time on multi-core systems.

Supports JIT-compiled version via --jit flag for 106x additional speedup.
"""

# Dynamic import based on --jit flag
import numpy as np
import csv
import sys
from scipy.special import loggamma as logg
import os
import argparse
import logging
import json
from datetime import datetime
from tqdm import tqdm
import multiprocessing as mp
from multiprocessing import Manager
from typing import Dict, Tuple
import time
from sim_utils import format_time


def Sk_stable(N: int, K: int) -> float:
    """Calculate kinetic entropy using stable logarithm.
    
    Args:
        N: Number of sites in lattice
        K: Total demon energy (sum of all oscillator energies)
        
    Returns:
        Kinetic entropy in units of k_B
    """
    return logg(K + N) - logg(K + 1) - logg(N)


def Su_stable(N: int, N0: int, Nx: int, N0_exp: int) -> float:
    """Calculate configurational entropy using stable logarithm.
    
    Args:
        N: Number of sites in lattice
        N0: Number of broken bonds
        Nx: Number of anti-aligned neighbor pairs
        N0_exp: N0 value for exponential term (max(N0, 1) to avoid log(0))
        
    Returns:
        Configurational entropy in units of k_B
    """
    log_2_term = N0_exp * np.log(2)
    return logg(N + 1) + log_2_term - (logg(N - N0 - Nx + 1) + logg(N0 + 1) + logg(Nx + 1))


def run_single_simulation(args: Tuple[int, int, int, int, str, str, bool, int, 'mp.Queue']) -> Dict:
    """Run a single simulation for given parameters.
    
    Args:
        args: Tuple of (R, M, n, s, validate_mode, project_root, use_jit, sim_num, progress_queue)
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
    progress_queue.put({
        'type': 'start',
        'sim_num': sim_num,
        'R': R,
        'M': M
    })
    
    # Import appropriate irrInferno class
    if use_jit:
        from jit_irr_inferno import JITirrInferno as irrInferno
    else:
        from irr_inferno import irrInferno
    
    # Create irrInferno instance
    x = irrInferno(n, R+1, validate_mode=validate_mode)
    
    # Setup output paths
    folder = os.path.join(project_root, 'data', 'irr')
    if R == 0:
        output_dir = os.path.join(folder, 'r0')
        filename = os.path.join(output_dir, f'irr_sim_data_{M}.csv')
    else:
        output_dir = os.path.join(folder, f'r{R}')
        filename = os.path.join(output_dir, f'irr_sim_data_r{R}_{M}.csv')
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Pre-allocate array for all results
    all_results = np.zeros((2*s, 7), dtype=np.float64)
    
    # Forward simulation
    for i in range(s):
        # Note: JIT version does full sweep per call, original needs N calls
        if use_jit:
            x.demon_move()  # JIT: full sweep
        else:
            for j in range(n):  # Original: N calls per sweep
                x.demon_move()
        
        # Send progress update (display refresh is throttled to 1/sec in main loop)
        progress_queue.put({
            'type': 'progress',
            'sim_num': sim_num,
            'phase': 'forward',
            'sweep': i + 1,
            'total_sweeps': s
        })
        
        # Get validated state after each sweep
        state = x.get_validated_state()
        
        # Calculate entropy with stable functions
        N0e = max(1, int(state['bond_count'][1]))
        
        try:
            Sk_val = Sk_stable(n, state['E_demon_sum'])
            Su_val = Su_stable(n, state['bond_count'][1], state['bond_count'][2], N0e)
            total_entropy = (Sk_val + Su_val) / n
        except (ValueError, OverflowError):
            total_entropy = np.nan
        
        # Store results
        all_results[i] = [
            np.float64(i + 1),
            np.float64(state['E_demon_sum']) / n,
            np.float64(state['E_lattice']) / n,
            np.float64(state['bond_count'][1]) / n,
            np.float64(state['bond_count'][2]) / n,
            total_entropy,
            np.float64(n)
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
        progress_queue.put({
            'type': 'progress',
            'sim_num': sim_num,
            'phase': 'reverse',
            'sweep': i + 1,
            'total_sweeps': s
        })
        
        # Get validated state
        state = x.get_validated_state()
        
        # Calculate entropy
        N0_exp = max(1, int(state['bond_count'][1]))
        
        try:
            Sk_val = Sk_stable(n, state['E_demon_sum'])
            Su_val = Su_stable(n, state['bond_count'][1], state['bond_count'][2], N0_exp)
            total_entropy = (Sk_val + Su_val) / n
        except (ValueError, OverflowError):
            total_entropy = np.nan
        
        # Store results
        all_results[s + i] = [
            np.float64(s + i + 1),
            np.float64(state['E_demon_sum']) / n,
            np.float64(state['E_lattice']) / n,
            np.float64(state['bond_count'][1]) / n,
            np.float64(state['bond_count'][2]) / n,
            total_entropy,
            np.float64(n)
        ]
    
    # Write results
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['t', 'K', 'U', 'N0', 'Nx', 'S/nk', 'n'])
        writer.writerows(all_results)
    
    # Save metadata
    metadata = {
        'lattice_size': n,
        'sweeps': s,
        'radius': R,
        'run': M,
        'timestamp': datetime.now().isoformat(),
        'simulation_type': 'irreversible_parallel'
    }
    metadata_file = filename.replace('.csv', '_metadata.json')
    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    # Send completion status
    progress_queue.put({
        'type': 'complete',
        'sim_num': sim_num,
        'R': R,
        'M': M
    })
    
    return {
        'R': R,
        'M': M,
        'filename': filename,
        'E_total': state['E_total'],
        'E_initial': x._initial_total_energy
    }


if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Run parallel irreversible Creutz demon simulation')
    parser.add_argument('--n', type=int, default=1000000, help='Lattice size (default: 1000000)')
    parser.add_argument('--s', type=int, default=5000, help='Number of sweeps per phase (default: 5000)')
    parser.add_argument('--r', type=int, default=11, help='Max demon-coupling radius, tests 1 to r-1 (default: 11)')
    parser.add_argument('--m', type=int, default=5, help='Number of independent runs (default: 5)')
    parser.add_argument('--cores', type=int, default=None, 
                       help='Number of CPU cores to use (default: auto-detect)')
    parser.add_argument('--validate', type=str, default='off', choices=['off', 'periodic', 'frequent'],
                       help='Validation mode: off (fastest), periodic (every 100 sweeps), frequent (every sweep)')
    parser.add_argument('--jit', action='store_true',
                       help='Use JIT-compiled version for 106x speedup (requires numba)')
    args = parser.parse_args()
    
    # Simulation parameters
    n = args.n  # lattice size
    s = args.s  # sweeps
    r = args.r  # max bond-demon couple radius
    m = args.m  # number of sims
    validate_mode = args.validate
    use_jit = args.jit
    
    # Determine number of cores
    if args.cores is None:
        num_cores = mp.cpu_count()
    else:
        num_cores = min(args.cores, mp.cpu_count())
    
    total_sims = r * m
    
    print(f"Parallel irreversible simulation:")
    print(f"  Lattice size: n={n}")
    print(f"  Sweeps: s={s}")
    print(f"  Radii: R=0 to R={r-1}")
    print(f"  Runs per radius: m={m}")
    print(f"  Total simulations: {total_sims}")
    print(f"  CPU cores: {num_cores}")
    print(f"  Validation mode: {validate_mode}")
    print(f"  JIT compilation: {'ENABLED' if use_jit else 'disabled'}")
    print(f"  Parallelization: {num_cores} cores")
    
    # Set up logging
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    log_dir = os.path.join(project_root, 'logs')
    os.makedirs(log_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'parallel_sim_irreversible_{timestamp}.log')
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    
    logging.info(f"Starting parallel irreversible simulation: n={n}, s={s}, r={r}, m={m}, cores={num_cores}")
    
    # Create data directories 
    for R in range(r):
        os.makedirs(os.path.join(project_root, 'data', 'irr', f'r{R}'), exist_ok=True)
    
    # Build list of all simulation parameters with simulation numbers
    sim_params = []
    sim_num = 1
    for M in range(m):
        for R in range(r):
            sim_params.append((R, M, n, s, validate_mode, project_root, use_jit, sim_num, None))
            sim_num += 1
    
    # Run simulations in parallel with progress monitoring
    start_time = datetime.now()
    
    # Create manager and progress queue
    manager = Manager()
    progress_queue = manager.Queue()
    
    # Update sim_params with the actual queue
    sim_params = [(R, M, n, s, val, root, jit, num, progress_queue) 
                  for R, M, n, s, val, root, jit, num, _ in sim_params]
    
    print(f"\nStarting {total_sims} simulations on {num_cores} cores...\n")
    
    # Track active simulations and completion
    active_sims = {}
    completed = 0
    sim_times = []
    early_eta_shown = False  # Track if we've shown early ETA estimate
    
    with mp.Pool(processes=num_cores) as pool:
        # Start async processing
        result = pool.map_async(run_single_simulation, sim_params)
        
        # Monitor progress (manually update desc to avoid tqdm's postfix formatting)
        pbar = tqdm(total=total_sims, unit="sim",
                   bar_format='{desc}: {percentage:3.0f}%|{bar}| {n}/{total} [{elapsed}]')
        pbar.set_description("Overall Progress")
        last_refresh = time.time()
        
        try:
            while not result.ready() or not progress_queue.empty():
                try:
                    msg = progress_queue.get(timeout=0.1)
                        
                    if msg['type'] == 'start':
                        active_sims[msg['sim_num']] = {
                            'R': msg['R'],
                            'M': msg['M'],
                            'start_time': time.time(),
                            'phase': 'forward',
                            'progress': 0,
                            'phase_start_time': time.time()  # Initialize for forward phase
                        }
                        
                    elif msg['type'] == 'progress':
                        if msg['sim_num'] in active_sims:
                            # Track phase transitions
                            old_phase = active_sims[msg['sim_num']]['phase']
                            new_phase = msg['phase']
                            if old_phase != new_phase:
                                active_sims[msg['sim_num']]['phase_start_time'] = time.time()
                            
                            active_sims[msg['sim_num']]['phase'] = new_phase
                            active_sims[msg['sim_num']]['progress'] = msg['sweep']
                            
                            # Update status message with active simulations
                            active_count = len(active_sims)
                            if active_count > 0:
                                # Show first active simulation as example with percentage
                                first_sim = list(active_sims.values())[0]
                                pct = (first_sim['progress'] / s) * 100
                                
                                # Only update display if throttle allows (1 second intervals)
                                now = time.time()
                                if now - last_refresh >= 1.0:
                                    # Determine what to display based on available data
                                    if sim_times:
                                        # Use actual completion times once available (most accurate)
                                        avg_time = np.mean(sim_times)
                                        remaining = total_sims - completed
                                        eta_minutes = (remaining * avg_time / num_cores) / 60
                                        desc = f"Overall Progress (Running {active_count} | Ex: R={first_sim['R']}, M={first_sim['M']}, {first_sim['phase']} {first_sim['progress']}/{s} {pct:.0f}% | ETA: ~{format_time(eta_minutes)} remaining)"
                                    elif first_sim['progress'] >= s * 0.01 and 'phase_start_time' in first_sim:
                                        # Calculate ETA if we have 1%+ progress
                                        phase_elapsed = time.time() - first_sim['phase_start_time']
                                        phase_time_estimate = (phase_elapsed / first_sim['progress']) * s
                                        sim_time_estimate = phase_time_estimate * 2  # forward + reverse
                                        total_eta_minutes = (total_sims * sim_time_estimate / num_cores) / 60
                                        desc = f"Overall Progress (Running {active_count} | Ex: R={first_sim['R']}, M={first_sim['M']}, {first_sim['phase']} {first_sim['progress']}/{s} {pct:.0f}% | Est: ~{format_time(total_eta_minutes)} remaining)"
                                    else:
                                        # Show running status with example (before 1% threshold)
                                        desc = f"Overall Progress (Running {active_count} | Ex: R={first_sim['R']}, M={first_sim['M']}, {first_sim['phase']} {first_sim['progress']}/{s} {pct:.0f}%)"
                                    
                                    pbar.set_description(desc)
                                    last_refresh = now
                    
                    elif msg['type'] == 'complete':
                        if msg['sim_num'] in active_sims:
                            elapsed = time.time() - active_sims[msg['sim_num']]['start_time']
                            sim_times.append(elapsed)
                            del active_sims[msg['sim_num']]
                            
                            completed += 1
                            pbar.update(1)
                            # ETA now handled in progress section for continuous updates
                                
                except Exception:
                    pass
                
                # Force refresh display every second to show current status
                if time.time() - last_refresh > 1.0:
                    pbar.refresh()
                    last_refresh = time.time()
        
        except KeyboardInterrupt:
            pbar.close()
            pool.terminate()
            pool.join()
            print("\n\nSimulation interrupted by user (Ctrl-C)")
            print(f"Completed {completed}/{total_sims} simulations before interruption")
            sys.exit(130)
        
        pbar.close()
        results = result.get()
    
    end_time = datetime.now()
    elapsed = (end_time - start_time).total_seconds()
    
    # Log results
    logging.info(f"Parallel simulation complete!")
    logging.info(f"Total time: {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
    logging.info(f"Time per simulation: {elapsed/len(results):.1f} seconds")
    
    print(f"\nCompleted {len(results)} simulations in {elapsed:.1f}s ({elapsed/60:.1f} min)")
    print(f"Average time per simulation: {elapsed/len(results):.1f}s")
    
    # Verify energy conservation
    energy_errors = []
    for result in results:
        if abs(result['E_total'] - result['E_initial']) > 1e-10:
            energy_errors.append(result)
            logging.warning(f"Energy conservation error: R={result['R']}, M={result['M']}")
    
    if energy_errors:
        print(f"\nWARNING: {len(energy_errors)} simulations had energy conservation errors")
    else:
        print(f"\n✓ All simulations passed energy conservation check")
