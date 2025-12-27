"""
Minimal working example of the Creutz demon simulation.

This script demonstrates the basic usage of the Inferno class
for a reversible microcanonical Monte Carlo simulation.
"""

import os
import sys

# Add creutz-sim to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))

import numpy as np
from inferno import Inferno

# Create a small lattice
print("Creating 1D Ising lattice with N=100 sites...")
N = 100
R = 3  # Demon coupling radius
x = Inferno(N, R)

print(f"Initial state:")
print(f"  Total energy: {x.E_total}")
print(f"  Lattice energy: {x.E_lattice}")
print(f"  Demon energy: {np.sum(x.E_demon)}")
print(f"  Energy conservation: {x.E_total == x.E_lattice + np.sum(x.E_demon)}")

# Run forward simulation
print(f"\nRunning forward simulation (10 sweeps)...")
for sweep in range(10):
    # One sweep = N attempted moves
    for _ in range(N):
        x.demon_move()

    if sweep % 5 == 0:
        state = x.get_validated_state()
        K = state["E_demon_sum"] / N
        U = state["E_lattice"] / N
        print(
            f"  Sweep {sweep}: E_total={state['E_total']:.1f}, "
            + f"<K>={K:.3f}, <U>={U:.3f}"
        )

# Run reverse simulation
print(f"\nRunning reverse simulation (10 sweeps)...")
for sweep in range(10):
    for _ in range(N):
        x.demon_reverse()

    if sweep % 5 == 0:
        state = x.get_validated_state()
        K = state["E_demon_sum"] / N
        U = state["E_lattice"] / N
        print(
            f"  Sweep {sweep}: E_total={state['E_total']:.1f}, "
            + f"<K>={K:.3f}, <U>={U:.3f}"
        )

print(f"\nFinal state:")
print(f"  Total energy: {x.E_total}")
print(f"  Energy conserved: {np.abs(x.E_total - 2*N) < 1e-10}")
print(f"\n✓ Simulation complete!")
