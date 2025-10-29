from irr_inferno import irrInferno
import numpy as np
import csv
from scipy.special import loggamma as logg
import math

def add_row(filename, row_data):    # appends a new row to csv file
    try:
        with open(filename, 'a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(row_data)
    except Exception as e:
         print(f"An error occurred: {e}")

Sk = lambda N, K: logg(K + N) - logg(K+1) - logg(N) # N == lattice size, K == kinetic energy
Su = lambda N, N0, Nx, N0_exp: logg(N+1) + math.log(2**N0_exp) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1)) # N == lattice size, N0 == broken bonds, Nx == bonds between anti-aligned spins

# lattice size
n=1000000
# sweeps
s = 10000
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

        init = np.array([x.lattice])

        temp_file = 'temp_sim.csv'
        temp_rev = 'temp_rev.csv'
        file_path = f"{file_names[R]}_{M}.csv"

        data_types = ['t', 'K', 'U', 'E', 'N0', 'Nx', 'S/nk', 'n'] # step counter, lattice energy, demon energy, total energy, broken bonds, anti-aligned spins, lattice size
        for files in [file_path]:
            with open(files, 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(data_types)

        for i in range(s):
            with open(temp_file, 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(['K', 'U', 'E', 'N0', 'Nx', 'S/nk', 'Smax'])
            # Attempt to flip each spin in lattice
            for j in range(n):
                x.demon_move()
                # Calculate total entropy
                N0e = int(x.bond_count[0])
                if N0e == 0:
                    N0e = 1
                total_entropy = (Sk(n, sum(x.E_demon)) + Su(n, x.bond_count[0], x.bond_count[1], N0e))/n
                # Write results to temp file
                new_row = [sum(x.E_demon), x.E_lattice, sum(x.E_demon) + x.E_lattice, x.bond_count[0]/n, x.bond_count[1]/n, int(total_entropy*1000)]
                add_row(temp_file, new_row)
            # write avg sweep results to csv
            with open(temp_file, 'r') as f:
                reader = csv.reader(f)
                next(reader)  # Skip the header row
                data = np.array(list(reader), dtype=float)
            new_row = [i+1, np.mean(data[:, 0]), np.mean(data[:, 1]), np.mean(data[:, 2]), np.mean(data[:, 3]), np.mean(data[:, 4]), np.mean(data[:, 5])/1000, n]
            add_row(file_path, new_row)

        ### Reverse simulation
        for i in range(s):
            with open(temp_rev, 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(['K', 'U', 'E', 'N0', 'Nx', 'S/nk', 'Smax'])
            # Attempt to flip each spin in lattice
            for j in range(n):
                x.demon_reverse()
                # Calculate total entropy
                N0_exp = int(x.bond_count[0])
                if N0_exp == 0:
                    N0_exp = 1
                total_entropy = (Sk(n, sum(x.E_demon)) + Su(n, x.bond_count[0], x.bond_count[1], N0_exp))/n
                # Write results to temp file
                new_row = [sum(x.E_demon), x.E_lattice, sum(x.E_demon) + x.E_lattice, x.bond_count[0]/n, x.bond_count[1]/n, int(total_entropy*1000)]
                add_row(temp_rev, new_row)
            # write avg sweep results to csv
            with open(temp_rev, 'r') as f:
                reader = csv.reader(f)
                next(reader)  # Skip the header row
                data = np.array(list(reader), dtype=float)
            new_row = [s+i, np.mean(data[:, 0]), np.mean(data[:, 1]), np.mean(data[:, 2]), np.mean(data[:, 3]), np.mean(data[:, 4]), np.mean(data[:, 5])/1000, n]
            add_row(file_path, new_row)
        print(f"R{R} complete")
    print(f'############################## Sim #{M+1} complete #################################### ')
