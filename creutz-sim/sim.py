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


n=10            # lattice size
s = 100         # sweeps
flag = 0        # 0 for rev, 1 for irr dynamics
k = 100         # number of sweeps before switching dynamics
r = 11          # max bond-demon couple radius
m = 5           # number of sims




folder = "/Users/winry/Documents/ASU/thesis/dev/data/"
init_folder = "/Users/winry/Documents/ASU/thesis/dev/init-fin/"

host = socket.gethostname()
if host != 'Luli':
  folder = '/home/wember/2025thesis/nanosim/data/'

status_file = f"{folder}sim_status.csv"

file_names = [f'{folder}r{i}/sim_data_r{i}' for i in range(11)]
irr_files = [f'{folder}irr/r{i}/irr_sim_data_r{i}' for i in range(11)]
init_files = [f'{init_folder}r{i}/sim_data_r{i}' for i in range(11)]

for M in range(m):
    for R in range(r):
        x = Inferno(n, R+1)

        filename = f"{file_names[R]}_{M}.csv"
        irr_filename = f"{irr_files[R]}_{M}.csv"
        init_filename = f"{init_files[R]}_{M}.csv"
        filenames = [filename, irr_filename, init_filename]

        data_types = ['t', 'K', 'U', 'N0', 'Nx', 'S/nk', 'n'] # step counter, lattice energy, demon energy, total energy, broken bonds, anti-aligned spins, lattice size

        # rev only
        filenames = [filename]
        # # irr only
        # filenames = [irr_filename]
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

        ### Reverse simulation
        for i in range(s//2):
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
