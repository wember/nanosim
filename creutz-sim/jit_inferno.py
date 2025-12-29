"""
Copyright (c) 2026 Winry Ember
Licensed under the MIT License.
See LICENSE file in the project root for full license information.

JIT-optimized Inferno class using Numba-compiled functions.

This module provides a drop-in replacement for the Inferno class with
3-10x performance improvement via Numba JIT compilation of hot functions.

Performance comparison:
- Original Inferno: ~8.15s for N=10,000, sweeps=100
- JIT Inferno: ~0.8-2.7s (expected)

Usage:
    from jit_inferno import JITInferno as Inferno
"""

import numpy as np
from jit_functions import demon_move_jit
from sim_utils import SimulationBase, Sk, Su, Su0  # noqa: F401


class JITInferno(SimulationBase):
    """
    JIT-optimized Inferno simulation with reversible dynamics.

    Uses Numba-compiled functions for spin flips and bond changes.
    Maintains identical physics and API to original Inferno class.
    """

    def __init__(self, N, R, validate_mode="off"):
        """
        Initialize JIT-optimized Inferno simulation.

        Args:
            N: Number of lattice sites
            R: Demon coupling radius
            validate_mode: Validation frequency
                'off' - No validation (fastest, production mode)
                'periodic' - Validate every 100 sweeps (default for testing)
                'frequent' - Validate every sweep (debug mode, slowest)
        """
        total_energy = 2 * N

        self.N = N
        self.R = R
        self.order, self.rev_order = SimulationBase.initialize_order_arrays(N)
        self.radius_spin = np.random.randint(0, R, size=N) * np.random.choice(
            [-1, 1], size=N
        )
        self.rev_radius_spin = np.flip(self.radius_spin)
        self.radius_bond = np.flip(self.radius_spin)
        self.rev_radius_bond = self.radius_spin.copy()

        # Initialize lattice, bonds, and bond counts using base class helper
        self.lattice, self.bonds, self.bond_count = SimulationBase.initialize_lattice(N)

        self.E_lattice = np.sum(self.bonds, dtype=np.int64)

        # Initialize demon energy using base class helper
        self.E_demon, self.E_demon_sum, self.d_energy = (
            SimulationBase.initialize_demon_energy(N, total_energy, self.E_lattice)
        )
        self.E_total = self.E_lattice + self.E_demon_sum

        # Setup validation using base class helper
        self.setup_validation_mode(N, validate_mode, total_energy)

        # Pre-compute neighbor indices using base class helper
        self.right_neighbor, self.left_neighbor = SimulationBase.setup_neighbor_arrays(
            N
        )

    def demon_move(self):
        """
        Perform one complete sweep using JIT-compiled functions.

        This replaces the N calls to spin_flip/bond_change with a single
        JIT-compiled loop, providing massive speedup (3-10x).
        """
        # Call JIT-compiled sweep
        self.E_demon_sum, self.E_lattice, self.bond_count = demon_move_jit(
            self.lattice,
            self.bonds,
            self.bond_count,
            self.E_demon,
            self.E_demon_sum,
            self.E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            self.order,
            self.radius_spin,
            self.radius_bond,
            self.N,
            self.R,
        )

        # Periodic validation (only if enabled)
        self.perform_periodic_validation()

    def demon_reverse(self):
        """
        Perform one reverse sweep using JIT-compiled functions.

        Uses reversed order and radii for time-reversible dynamics.
        """
        # Call JIT-compiled sweep with reversed parameters
        self.E_demon_sum, self.E_lattice, self.bond_count = demon_move_jit(
            self.lattice,
            self.bonds,
            self.bond_count,
            self.E_demon,
            self.E_demon_sum,
            self.E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            self.rev_order,
            self.rev_radius_bond,  # Reversed!
            self.rev_radius_spin,  # Reversed!
            self.N,
            self.R,
        )

        # Periodic validation (only if enabled)
        self.perform_periodic_validation()


# Convenience alias for drop-in replacement
Inferno = JITInferno
