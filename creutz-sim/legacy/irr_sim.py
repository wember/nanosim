"""
Copyright (c) 2026 Winry Ember
Licensed under the MIT License.
See LICENSE file in the project root for full license information.
"""

import argparse
import csv
import math
import os
import sys
from datetime import datetime

import numpy as np
from scipy.special import loggamma as logg

# Ensure parent directory is in sys.path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from inferno_irr import irrInferno


def add_row(filename, row_data):  # appends a new row to csv file
    try:
        with open(filename, "a", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(row_data)
    except Exception as e:
        print(f"An error occurred: {e}")


Sk = (
    lambda N, K: logg(K + N) - logg(K + 1) - logg(N)
)  # N == lattice size, K == kinetic energy
Su = (
    lambda N, N0, Nx, N0_exp: logg(N + 1)
    + math.log(2**N0_exp)
    - (logg(N - N0 - Nx + 1) + logg(N0 + 1) + logg(Nx + 1))
)  # N == lattice size, N0 == broken bonds, Nx == bonds between anti-aligned spins

# Parse command line arguments
parser = argparse.ArgumentParser(description="Run legacy irreversible simulation")
parser.add_argument("--n", type=int, default=10000, help="Lattice size")
parser.add_argument("--s", type=int, default=10000, help="Number of sweeps")
parser.add_argument("--r", type=int, default=11, help="Max radius")
parser.add_argument("--m", type=int, default=5, help="Number of runs")
args = parser.parse_args()

n = args.n
s = args.s
r = args.r
m = args.m

# Use relative path from project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
folder = os.path.join(project_root, "data") + os.sep

# Create data directories if they don't exist
for R in range(r):
    os.makedirs(os.path.join(project_root, "data", "irr", f"r{R}"), exist_ok=True)

status_file = f"{folder}irr_sim_status.csv"
file_names = [f"{folder}irr/r{R}/irr_sim_data_r{R}" for R in range(r)]


with open(status_file, "w+", newline="") as file:
    writer = csv.writer(file)

add_row(status_file, f"{datetime.now()}\t Begin irr_sim")


for M in range(m):
    for R in range(r):
        x = irrInferno(n, R + 1)

        filename = f"{file_names[R]}_{M}.csv"

        data_types = [
            "t",
            "K",
            "U",
            "N0",
            "Nx",
            "S/nk",
            "n",
        ]
        # step counter, lattice energy, demon energy, total energy, broken bonds,
        # anti-aligned spins, lattice size

        with open(filename, "w+", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(data_types)

        for i in range(s):
            data = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
            # Attempt to flip each spin in lattice
            for j in range(n):
                x.demon_move()
                # Calculate total entropy
                N0e = int(x.bond_count[1])
                if N0e == 0:
                    N0e = 1
                total_entropy = (
                    Sk(n, sum(x.E_demon)) + Su(n, x.bond_count[1], x.bond_count[2], N0e)
                ) / n
                # Add results to totals
                data += [
                    float(sum(x.E_demon)),
                    float(x.E_lattice),
                    x.bond_count[1] / n,
                    x.bond_count[2] / n,
                    total_entropy,
                ]
            # write avg sweep results to csv
            new_row = np.array(
                [i, data[0] / n, data[1] / n, data[2] / n, data[3] / n, data[4] / n, n]
            )
            add_row(filename, new_row)

        # Reverse simulation
        for i in range(s):
            data = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
            # Attempt to flip each spin in lattice
            for j in range(n):
                x.demon_reverse()
                # Calculate total entropy
                N0_exp = int(x.bond_count[1])
                if N0_exp == 0:
                    N0_exp = 1
                total_entropy = (
                    Sk(n, sum(x.E_demon))
                    + Su(n, x.bond_count[1], x.bond_count[2], N0_exp)
                ) / n
                # Add results to totals
                data += [
                    float(sum(x.E_demon)),
                    float(x.E_lattice),
                    x.bond_count[1] / n,
                    x.bond_count[2] / n,
                    total_entropy,
                ]
            # write avg sweep results to csv
            new_row = np.array(
                [
                    s + i,
                    data[0] / n,
                    data[1] / n,
                    data[2] / n,
                    data[3] / n,
                    data[4] / n,
                    n,
                ]
            )
            add_row(filename, new_row)
