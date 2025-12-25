#!/usr/bin/env python3
"""
Wrapper script to run irreversible simulation with small parameters for quick testing.
Parameters: n=1000, s=100, r=3, m=2
"""
from irr_inferno import irrInferno
import numpy as np
import csv
from scipy.special import loggamma as logg
import math
import os

# Use stable logarithm calculations for large factorials
def Sk_stable(N, K):
    """Stable entropy calculation for kinetic energy"""
    return logg(K + N) - logg(K + 1) - logg(N)

def Su_stable(N, N0, Nx, N0_exp):
    """Stable entropy calculation for potential energy"""
    # Use N0_exp * log(2) instead of log(2^N0_exp) to avoid overflow
    log_2_term = N0_exp * np.log(2)
    return logg(N + 1) + log_2_term - (logg(N - N0 - Nx + 1) + logg(N0 + 1) + logg(Nx + 1))

# Small parameters for quick testing
n = 1000     # lattice size
s = 100      # sweeps
r = 3        # max bond-demon couple radius (tests R=1,2)
m = 2        # number of sims

# Use relative path from project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
folder = os.path.join(project_root, 'data') + os.sep

# Create data directories if they don't exist
for R in range(r):
    os.makedirs(os.path.join(project_root, 'data', 'irr', f'r{R}'), exist_ok=True)

file_names = [f'{folder}irr/r0/irr_sim_data',
              f'{folder}irr/r1/irr_sim_data_r1',
              f'{folder}irr/r2/irr_sim_data_r2',
              f'{folder}irr/r3/irr_sim_data_r3',
              f'{folder}irr/r4/irr_sim_data_r4',
              f'{folder}irr/r5/irr_sim_data_r5',
              f'{folder}irr/r6/irr_sim_data_r6',
              f'{folder}irr/r7/irr_sim_data_r7',
              f'{folder}irr/r8/irr_sim_data_r8',
              f'{folder}irr/r9/irr_sim_data_r9',
              f'{folder}irr/r10/irr_sim_data_r10']

for M in range(m):
    for R in range(r):
        x = irrInferno(n, R+1)
        filename = f"{file_names[R]}_{M}.csv"

        # Pre-allocate array for all results - use float64 for final calculations
        all_results = np.zeros((2*s, 7), dtype=np.float64)

        print(f"Starting M={M}, R={R}")

        # Forward simulation
        for i in range(s):
            for j in range(n):
                x.demon_move()

            state = x.get_validated_state()
            N0e = max(1, int(state['bond_count'][1]))

            try:
                Sk_val = Sk_stable(n, state['E_demon_sum'])
                Su_val = Su_stable(n, state['bond_count'][1], state['bond_count'][2], N0e)
                total_entropy = (Sk_val + Su_val) / n
            except (ValueError, OverflowError) as e:
                print(f"Warning: Entropy calculation error at sweep {i}: {e}")
                total_entropy = 0.0

            all_results[i] = [i, state['K'], state['U'], state['N0'], state['Nx'], total_entropy, n]

        # Reverse simulation
        for i in range(s):
            for j in range(n):
                x.demon_reverse()

            state = x.get_validated_state()
            N0e = max(1, int(state['bond_count'][1]))

            try:
                Sk_val = Sk_stable(n, state['E_demon_sum'])
                Su_val = Su_stable(n, state['bond_count'][1], state['bond_count'][2], N0e)
                total_entropy = (Sk_val + Su_val) / n
            except (ValueError, OverflowError) as e:
                print(f"Warning: Entropy calculation error at sweep {s+i}: {e}")
                total_entropy = 0.0

            all_results[s + i] = [s + i, state['K'], state['U'], state['N0'], state['Nx'], total_entropy, n]

        # Write all results at once
        with open(filename, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['t', 'K', 'U', 'N0', 'Nx', 'S/nk', 'n'])
            writer.writerows(all_results)

        print(f"Completed M={M}, R={R}")

print("Small simulation complete!")
