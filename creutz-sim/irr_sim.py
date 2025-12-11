from irr_inferno import irrInferno
import numpy as np
import csv
from scipy.special import loggamma as logg
import math
import socket

# Use stable logarithm calculations for large factorials
def Sk_stable(N, K):
    """Stable entropy calculation for kinetic energy"""
    return logg(K + N) - logg(K + 1) - logg(N)

def Su_stable(N, N0, Nx, N0_exp):
    """Stable entropy calculation for potential energy"""
    # Use N0_exp * log(2) instead of log(2^N0_exp) to avoid overflow
    log_2_term = N0_exp * np.log(2)
    return logg(N + 1) + log_2_term - (logg(N - N0 - Nx + 1) + logg(N0 + 1) + logg(Nx + 1))

# lattice size
n = 1000000
# sweeps
s = 5000
# max bond-demon couple radius
r = 11
# number of sims
m = 5

folder = "/Users/winry/Documents/ASU/thesis/dev/data/"

host = socket.gethostname()
if host != 'Luli.local':
    folder = '/home/wember/2025thesis/nanosim/data/'

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

        print(f"Starting irreversible sim M={M}, R={R}")

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
            except (ValueError, OverflowError) as e:
                print(f"Warning: Entropy calculation error at sweep {i}: {e}")
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

            # Progress indicator every 100 sweeps
            if (i + 1) % 100 == 0:
                print(f"  Forward sweep {i+1}/{s} complete")
                # Verify energy conservation
                assert state['E_total'] == x._initial_total_energy, \
                    f"Energy conservation violated: {state['E_total']} != {x._initial_total_energy}"

        # Reverse simulation
        print(f"  Starting reverse simulation")
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
            except (ValueError, OverflowError) as e:
                print(f"Warning: Entropy calculation error at reverse sweep {i}: {e}")
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

            # Progress indicator
            if (i + 1) % 100 == 0:
                print(f"  Reverse sweep {i+1}/{s} complete")
                # Verify energy conservation
                assert state['E_total'] == x._initial_total_energy, \
                    f"Energy conservation violated: {state['E_total']} != {x._initial_total_energy}"

        # Write all results at once (much faster than row-by-row)
        with open(filename, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['t', 'K', 'U', 'N0', 'Nx', 'S/nk', 'n'])
            writer.writerows(all_results)

        print(f"Completed M={M}, R={R} - saved to {filename}")
        print(f"  Final energy check: {state['E_total']} == {x._initial_total_energy}")
