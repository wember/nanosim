from inferno import Inferno
import numpy as np
import csv
from scipy.special import loggamma as logg
import math
import socket
from datetime import datetime

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
n=10000
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


with open(status_file, 'w+', newline='') as file:
    writer = csv.writer(file)

add_row(status_file, f"{datetime.now()}\t Begin sim")

for M in range(m):
    for R in range(r):
        x = Inferno(n, R+1)

        filename = f"{file_names[R]}_{M}.csv"

        data_types = ['t', 'K', 'U', 'N0', 'Nx', 'S/nk', 'n'] # step counter, lattice energy, demon energy, total energy, broken bonds, anti-aligned spins, lattice size

        with open(filename, 'w+', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(data_types)

        for i in range(s):
            data = np.array([0.0,0.0,0.0,0.0,0.0])
            # Attempt to flip each spin in lattice
            for j in range(n):
                x.demon_move()
                # Calculate total entropy
                N0e = int(x.bond_count[1])
                if N0e == 0:
                    N0e = 1
                total_entropy = (Sk(n, sum(x.E_demon)) + Su(n, x.bond_count[1], x.bond_count[2], N0e))/n
                # Add results to totals
                data += [float(sum(x.E_demon)), float(x.E_lattice), x.bond_count[1]/n, x.bond_count[2]/n, total_entropy]
            # write avg sweep results to csv
            new_row = np.array([i+1, data[0]/n, data[1]/n, data[2]/n, data[3]/n, data[4]/n, n])
            add_row(filename, new_row)
            add_row(status_file, f"{datetime.now()}\t sweep {i} complete.")

        ### Reverse simulation
        for i in range(s):
            data = np.array([0.0,0.0,0.0,0.0,0.0])
            # Attempt to flip each spin in lattice
            for j in range(n):
                x.demon_reverse()
                # Calculate total entropy
                N0_exp = int(x.bond_count[1])
                if N0_exp == 0:
                    N0_exp = 1
                total_entropy = (Sk(n, sum(x.E_demon)) + Su(n, x.bond_count[1], x.bond_count[2], N0_exp))/n
                # Add results to totals
                data += [float(sum(x.E_demon)), float(x.E_lattice), x.bond_count[1]/n, x.bond_count[2]/n, total_entropy]
            # write avg sweep results to csv
            new_row = np.array([s+i+1, data[0]/n, data[1]/n, data[2]/n, data[3]/n, data[4]/n, n])
            add_row(filename, new_row)
