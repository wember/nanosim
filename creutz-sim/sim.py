from inferno import Inferno
import numpy as np
import csv
from scipy.special import loggamma as logg
import math
import os
import argparse
import logging
import json
from datetime import datetime
from tqdm import tqdm
from typing import Dict, Tuple

# Use stable logarithm calculations for large factorials
def Sk_stable(N: int, K: int) -> float:
    """Calculate kinetic entropy using stable logarithm.
    
    Args:
        N: Number of sites in lattice
        K: Total demon energy (sum of all oscillator energies)
        
    Returns:
        Kinetic entropy in units of k_B
        
    Note:
        Uses loggamma for numerical stability with large factorials.
        Calculates log(Ω_k) where Ω_k = (K+N-1)! / (K! * (N-1)!)
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
        
    Note:
        Uses N0_exp * log(2) instead of log(2^N0_exp) to avoid overflow.
        Calculates log(Ω_u) where Ω_u = N! * 2^N0 / ((N-N0-Nx)! * N0! * Nx!)
    """
    log_2_term = N0_exp * np.log(2)
    return logg(N + 1) + log_2_term - (logg(N - N0 - Nx + 1) + logg(N0 + 1) + logg(Nx + 1))

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Run reversible Creutz demon simulation')
parser.add_argument('--n', type=int, default=1000000, help='Lattice size (default: 1000000)')
parser.add_argument('--s', type=int, default=5000, help='Number of sweeps per phase (default: 5000)')
parser.add_argument('--r', type=int, default=11, help='Max demon-coupling radius, tests 1 to r-1 (default: 11)')
parser.add_argument('--m', type=int, default=5, help='Number of independent runs (default: 5)')
parser.add_argument('--validate', type=str, default='off', choices=['off', 'periodic', 'frequent'],
                   help='Validation mode: off (fastest), periodic (every 100 sweeps), frequent (every sweep)')
args = parser.parse_args()

# Simulation parameters
n = args.n  # lattice size
s = args.s  # sweeps
r = args.r  # max bond-demon couple radius
m = args.m  # number of sims
validate_mode = args.validate  # validation mode

print(f"Reversible simulation with: n={n}, s={s}, r={r}, m={m}, validate={validate_mode}")

# Set up logging
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
log_dir = os.path.join(project_root, 'logs')
os.makedirs(log_dir, exist_ok=True)

timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
log_file = os.path.join(log_dir, f'sim_reversible_{timestamp}.log')

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)

logging.info(f"Starting reversible simulation: n={n}, s={s}, r={r}, m={m}")
logging.info(f"Log file: {log_file}")

# Use relative path from project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
folder = os.path.join(project_root, 'data') + os.sep

# Create data directories if they don't exist
for R in range(r):
    os.makedirs(os.path.join(project_root, 'data', f'r{R}'), exist_ok=True)

status_file = f"{folder}sim_status.csv"
file_names = [f'{folder}r0/sim_data',
              f'{folder}r1/sim_data_r1',
              f'{folder}r2/sim_data_r2',
              f'{folder}r3/sim_data_r3',
              f'{folder}r4/sim_data_r4',
              f'{folder}r5/sim_data_r5',
              f'{folder}r6/sim_data_r6',
              f'{folder}r7/sim_data_r7',
              f'{folder}r8/sim_data_r8',
              f'{folder}r9/sim_data_r9',
              f'{folder}r10/sim_data_r10']

for M in tqdm(range(m), desc="Runs", position=0):
    for R in tqdm(range(r), desc="Radii", position=1, leave=False):
        logging.info(f"Starting run M={M}, radius R={R}")
        x = Inferno(n, R+1, validate_mode=validate_mode)
        filename = f"{file_names[R]}_{M}.csv"

        # Save metadata
        metadata = {
            'lattice_size': n,
            'sweeps': s,
            'radius': R,
            'run': M,
            'timestamp': datetime.now().isoformat(),
            'simulation_type': 'reversible'
        }
        metadata_file = filename.replace('.csv', '_metadata.json')
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)

        # Pre-allocate array for all results - use float64 for final calculations
        all_results = np.zeros((2*s, 7), dtype=np.float64)

        # Forward simulation
        for i in tqdm(range(s), desc="Forward", position=2, leave=False):
            # Attempt to flip each spin in lattice
            for j in range(n):
                x.demon_move()

            # Get validated state after each sweep
            state = x.get_validated_state()

            # Calculate entropy with stable functions
            N0e = max(1, int(state['bond_count'][1]))

            # Use stable entropy calculations
            try:
                Sk_val = Sk_stable(n, state['E_demon_sum'])
                Su_val = Su_stable(n, state['bond_count'][1], state['bond_count'][2], N0e)
                total_entropy = (Sk_val + Su_val) / n
            except (ValueError, OverflowError) as e:
                logging.warning(f"Entropy calculation error at sweep {i}: {e}")
                total_entropy = np.nan

            # Store results - division by n done in float64
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
        for i in tqdm(range(s), desc="Reverse", position=2, leave=False):
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
            except (ValueError, OverflowError) as e:
                logging.warning(f"Entropy calculation error at reverse sweep {i}: {e}")
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

        # Write all results at once
        with open(filename, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['t', 'K', 'U', 'N0', 'Nx', 'S/nk', 'n'])
            writer.writerows(all_results)

        logging.info(f"Completed run M={M}, radius R={R}, output: {filename}")

logging.info("Reversible simulation complete!")

        print(f"Completed M={M}, R={R} - saved to {filename}")
        print(f"  Final energy check: {state['E_total']} == {x._initial_total_energy}")
