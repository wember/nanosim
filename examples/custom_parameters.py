"""
Example showing how to customize simulation parameters.

Demonstrates running simulations with different lattice sizes,
radii, and analyzing the results.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))

import numpy as np
from inferno import Inferno
from inferno_irr import irrInferno
from scipy.special import loggamma as logg


def calculate_entropy(N, K, N0, Nx):
    """Calculate total entropy per site."""
    Sk = logg(K + N) - logg(K + 1) - logg(N)
    N0_exp = max(N0, 1)
    Su = (
        logg(N + 1)
        + N0_exp * np.log(2)
        - (logg(N - N0 - Nx + 1) + logg(N0 + 1) + logg(Nx + 1))
    )
    return (Sk + Su) / N


print("Comparing different lattice sizes and radii...\n")

# Test different configurations
configs = [
    (50, 2, "Small lattice, small radius"),
    (50, 10, "Small lattice, large radius"),
    (200, 2, "Large lattice, small radius"),
    (200, 10, "Large lattice, large radius"),
]

for N, R, description in configs:
    print(f"{description} (N={N}, R={R}):")

    # Reversible
    x_rev = Inferno(N, R)
    for _ in range(100):  # 100 moves
        x_rev.demon_move()
    state_rev = x_rev.get_validated_state()
    K_rev = state_rev["E_demon_sum"] / N
    U_rev = state_rev["E_lattice"] / N

    # Irreversible
    x_irr = irrInferno(N, R)
    for _ in range(100):
        x_irr.demon_move()
    state_irr = x_irr.get_validated_state()
    K_irr = state_irr["E_demon_sum"] / N
    U_irr = state_irr["E_lattice"] / N

    print(f"  Reversible:   <K>={K_rev:.4f}, <U>={U_rev:.4f}")
    print(f"  Irreversible: <K>={K_irr:.4f}, <U>={U_irr:.4f}")
    print()

print("✓ Comparison complete!")
