"""JIT-optimized irrInferno class using Numba-compiled functions.

This module provides a drop-in replacement for the irrInferno class with
3-10x performance improvement via Numba JIT compilation.

Performance comparison:
- Original irrInferno: ~15.78s for N=10,000, sweeps=100
- JIT irrInferno: ~1.6-5.3s (expected)

Usage:
    from jit_inferno_irr import JITirrInferno as irrInferno
"""

import numpy as np
from jit_functions import demon_move_irr_jit
from sim_utils import SimulationBase, Sk, Su, Su0  # noqa: F401


class JITirrInferno(SimulationBase):
    """
    JIT-optimized irreversible Inferno simulation.

    Uses Numba-compiled functions with RNG for truly random dynamics.
    Maintains identical physics and API to original irrInferno class.
    """

    def __init__(self, N, R, validate_mode="off"):
        """
        Initialize JIT-optimized irreversible Inferno simulation.

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

        # Random seed counter for reproducibility
        self._seed_counter = 0

    def demon_move(self):
        """
        Perform one complete sweep with random radii (JIT-compiled).

        Generates random radii on-the-fly for true irreversibility.
        """
        # Generate unique seed for this sweep
        seed = np.random.randint(0, 2**31) + self._seed_counter
        self._seed_counter += 1

        # Call JIT-compiled sweep with RNG
        self.E_demon_sum, self.E_lattice, self.bond_count = demon_move_irr_jit(
            self.lattice,
            self.bonds,
            self.bond_count,
            self.E_demon,
            self.E_demon_sum,
            self.E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            self.order,
            self.N,
            self.R,
            seed,
        )

        # Periodic validation (only if enabled)
        self.perform_periodic_validation()

    def demon_reverse(self):
        """
        Perform one reverse sweep with random radii (JIT-compiled).

        Uses reversed order but still generates random radii.
        """
        # Generate unique seed for this sweep
        seed = np.random.randint(0, 2**31) + self._seed_counter
        self._seed_counter += 1

        # Call JIT-compiled sweep with reversed order
        self.E_demon_sum, self.E_lattice, self.bond_count = demon_move_irr_jit(
            self.lattice,
            self.bonds,
            self.bond_count,
            self.E_demon,
            self.E_demon_sum,
            self.E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            self.rev_order,  # Reversed order
            self.N,
            self.R,
            seed,
        )

        # Periodic validation (only if enabled)
        self.perform_periodic_validation()


# Convenience alias for drop-in replacement
irrInferno = JITirrInferno
