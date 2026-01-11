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
    parser.add_argument('-f', '--flag', type=int, choices=[0, 1], help='Dynamics flag: 0=reversible, 1=irreversible (default: 0)')
    parser.add_argument('-r', '--radius', type=int, help='Max bond-demon couple radius (default: 11)')
    parser.add_argument('-m', '--runs', type=int, help='Number of simulation runs (default: 5)')
    parser.add_argument('-i', '--interactive', action='store_true', help='Interactive mode')
    
    args = parser.parse_args()
    
    # Default values
    defaults = {
        'n': 10,
        's': 100,
        'flag': 0,
        'r': 11,
        'm': 5
    }
    
    # Interactive mode
    if args.interactive or not any([args.lattice_size, args.sweeps, args.flag is not None, args.radius, args.runs]):
        print("=== Monte Carlo Simulation Configuration ===")
        print()
        
        try:
            n = input(f"Lattice size [{defaults['n']}]: ").strip()
            defaults['n'] = int(n) if n else defaults['n']
            
            s = input(f"Number of sweeps [{defaults['s']}]: ").strip()
            defaults['s'] = int(s) if s else defaults['s']
            
            f = input(f"Dynamics (0=reversible, 1=irreversible) [{defaults['flag']}]: ").strip()
            defaults['flag'] = int(f) if f else defaults['flag']
            
            r = input(f"Max radius [{defaults['r']}]: ").strip()
            defaults['r'] = int(r) if r else defaults['r']
            
            m = input(f"Number of runs [{defaults['m']}]: ").strip()
            defaults['m'] = int(m) if m else defaults['m']
            
            print()
        except (ValueError, KeyboardInterrupt):
            print("\nUsing default values")
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
    
    return defaults['n'], defaults['s'], defaults['flag'], defaults['r'], defaults['m']


n, s, flag, r, m = get_params()
k = 100  # number of sweeps before switching dynamics

# Calculate total iterations for progress tracking
total_sims = r * m
total_sweeps_per_sim = s  # s//2 forward + s//2 reverse

print(f"Starting {total_sims} simulations ({r} radii × {m} runs)")
print()

# Use relative path from repo root
repo_root = Path(__file__).parent.parent
data_folder = repo_root / 'data'
init_folder = repo_root / 'init-fin'

status_file = data_folder / 'sim_status.csv'

file_names = [data_folder / f'r{i}' / f'sim_data_r{i}' for i in range(11)]
irr_files = [data_folder / 'irr' / f'r{i}' / f'irr_sim_data_r{i}' for i in range(11)]
init_files = [init_folder / f'r{i}' / f'sim_data_r{i}' for i in range(11)]

# Progress tracking
sim_counter = 0
start_time = time.time()
total_sweeps = total_sims * s
completed_sweeps = 0

pbar = tqdm(total=total_sweeps, desc="Progress", unit="sweep", 
            bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [elapsed {elapsed}|remaining {remaining}|{rate_fmt}]')

for M in range(m):
    for R in range(r):
        x = Inferno(n, R+1)

        filename = file_names[R].parent / f"{file_names[R].name}_{M}.csv"
        irr_filename = irr_files[R].parent / f"{irr_files[R].name}_{M}.csv"
        init_filename = init_files[R].parent / f"{init_files[R].name}_{M}.csv"
        filenames = [filename, irr_filename, init_filename]

        data_types = ['t', 'K', 'U', 'N0', 'Nx', 'S/nk', 'n'] # step counter, lattice energy, demon energy, total energy, broken bonds, anti-aligned spins, lattice size

        # rev only
        filenames = [filename]
        # # irr only
        # filenames = [irr_filename]
        
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
            # Update progress description periodically
            if i % max(1, (s//2)//10) == 0:  # Update ~10 times during phase
                pbar.set_description(f"R={R}/{r-1}, M={M}/{m-1} [fwd {i}/{s//2}]")
            # flag = (i // (n//k)) % 2    # if 0, perform reversible dynamics
            data = np.zeros(5)
            # Attempt to flip each spin in lattice
            for j in range(n):
                x.demon_move(flag)
                # Calculate total entropy
                N0e = int(x.bond_count[1])
                if N0e == 0:
                    N0e = 1
                total_entropy = (Sk(n, sum(x.E_demon)) + Su(n, x.bond_count[1], x.bond_count[2], N0e))/n
                # Add results to totals
                data += [sum(x.E_demon), x.E_lattice, x.bond_count[1]/n, x.bond_count[2]/n, total_entropy]

            # Increment appropriate counter
            if flag == 0:
                t_counter += 1
                t_value = t_counter
            else:
                t_irr_counter += 1
                t_value = t_irr_counter

            # write avg sweep results to csv
            new_row = np.array([t_value, data[0]/n, data[1]/n, data[2]/n, data[3]/n, data[4]/n, n])

            # Save initial and final states to init_filename, all states to appropriate file
            # if (i == 0):
            #     add_row(init_filename, new_row)

            if (flag == 0):
                add_row(filename, new_row)
            else:
                add_row(irr_filename, new_row)
            
            # Update sweep progress
            completed_sweeps += 1
            pbar.update(1)

        ### Reverse simulation
        for i in range(s//2):
            # Update progress description periodically
            if i % max(1, (s//2)//10) == 0:  # Update ~10 times during phase
                pbar.set_description(f"R={R}/{r-1}, M={M}/{m-1} [rev {i}/{s//2}]")
            # flag = (i // (n//k)) % 2    # if 0, perform reversible dynamics
            data = np.zeros(5)
            # Attempt to flip each spin in lattice (full sweep)
            for j in range(n):
                x.demon_reverse(flag)
                # Calculate total entropy
                N0_exp = int(x.bond_count[1])
                if N0_exp == 0:
                    N0_exp = 1
                total_entropy = (Sk(n, sum(x.E_demon)) + Su(n, x.bond_count[1], x.bond_count[2], N0_exp))/n
                # Add results to totals
                data += [sum(x.E_demon), x.E_lattice, x.bond_count[1]/n, x.bond_count[2]/n, total_entropy]

            # Increment appropriate counter
            if flag == 0:
                t_counter += 1
                t_value = t_counter
            else:
                t_irr_counter += 1
                t_value = t_irr_counter

            # write avg sweep results to csv
            new_row = np.array([t_value, data[0]/n, data[1]/n, data[2]/n, data[3]/n, data[4]/n, n])

            # Save initial and final states to init_filename, all states to appropriate file
            # if (i == (s//2-1)):
            #     add_row(init_filename, new_row)

            if (flag == 0):
                add_row(filename, new_row)
            else:
                add_row(irr_filename, new_row)
            
            # Update sweep progress
            completed_sweeps += 1
            pbar.update(1)
        
        # Update simulation counter
        sim_counter += 1

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
