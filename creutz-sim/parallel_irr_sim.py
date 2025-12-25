"""Parallel implementation of irreversible Creutz demon simulation.

Uses multiprocessing to run multiple independent simulations simultaneously,
significantly reducing total execution time on multi-core systems.
"""

from irr_inferno import irrInferno
import numpy as np
import csv
from scipy.special import loggamma as logg
import os
import argparse
import logging
import json
from datetime import datetime
from tqdm import tqdm
import multiprocessing as mp
from typing import Dict, Tuple


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


def run_single_simulation(args: Tuple[int, int, int, int, str, str]) -> Dict:
    """Run a single simulation for given parameters.
    
    Args:
        args: Tuple of (R, M, n, s, validate_mode, project_root)
            R: Demon-coupling radius
            M: Run number
            n: Lattice size
            s: Number of sweeps
            validate_mode: Validation mode ('off', 'periodic', 'frequent')
            project_root: Project root directory
            
    Returns:
        Dictionary with simulation results and metadata
    """
    R, M, n, s, validate_mode, project_root = args
    
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
        # Attempt to flip each spin in lattice
        for j in range(n):
            x.demon_move()
        
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
        # Attempt to flip each spin in lattice
        for j in range(n):
            x.demon_reverse()
        
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
    args = parser.parse_args()
    
    # Simulation parameters
    n = args.n  # lattice size
    s = args.s  # sweeps
    r = args.r  # max bond-demon couple radius
    m = args.m  # number of sims
    validate_mode = args.validate
    
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
    print(f"  Expected speedup: ~{min(num_cores, total_sims)}x")
    
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
    
    # Build list of all simulation parameters
    sim_params = []
    for M in range(m):
        for R in range(r):
            sim_params.append((R, M, n, s, validate_mode, project_root))
    
    # Run simulations in parallel
    start_time = datetime.now()
    
    with mp.Pool(processes=num_cores) as pool:
        results = list(tqdm(
            pool.imap(run_single_simulation, sim_params),
            total=len(sim_params),
            desc="Simulations"
        ))
    
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
