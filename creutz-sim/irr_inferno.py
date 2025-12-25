import numpy as np
import random
import math

from sim_utils import Sk, Su, Su0, SimulationBase


class irrInferno(SimulationBase):
    """
        Optimized irrInferno class with roundoff error prevention
        Uses truly random radius selection (no pre-generated arrays)
    """

    def __init__(self, N, R, validate_mode='off'):
        """
        Initialize irrInferno simulation.
        
        Args:
            N: Number of lattice sites
            R: Demon coupling radius
            validate_mode: Validation frequency
                'off' - No validation (fastest, production mode)
                'periodic' - Validate every 100 sweeps (default for testing)
                'frequent' - Validate every sweep (debug mode, slowest)
        """
        total_energy = 2*N

        self.N = N
        self.order, self.rev_order = SimulationBase.initialize_order_arrays(N)
        self.R = R

        # Initialize lattice, bonds, and bond counts using base class helper
        self.lattice, self.bonds, self.bond_count = SimulationBase.initialize_lattice(N)

        self.E_lattice = np.sum(self.bonds, dtype=np.int64)
        
        # Initialize demon energy using base class helper
        self.E_demon, self.E_demon_sum, self.d_energy = SimulationBase.initialize_demon_energy(
            N, total_energy, self.E_lattice
        )
        self.E_total = self.E_lattice + self.E_demon_sum

        # Setup validation using base class helper
        self.setup_validation_mode(N, validate_mode, total_energy)

        # Indices for rolling
        self.order_idx = 0
        self.rev_order_idx = 0

        # Pre-compute neighbor indices using base class helper
        self.right_neighbor, self.left_neighbor = SimulationBase.setup_neighbor_arrays(N)

    def demon_move(self):
        """
            Move the demon with truly random radius selection
            This version generates random radius on-the-fly for irreversibility
        """
        a = self.order[self.order_idx]

        # Generate random radius and direction for spin flip
        radius_spin = np.random.randint(0, self.R) * (2 * np.random.randint(0, 2) - 1)
        bonds_changed = self.spin_flip(a, (a + radius_spin) % self.N)
        if bonds_changed:
            self.update_bonds_incremental(a)

        # Generate random radius and direction for bond change
        radius_bond = np.random.randint(0, self.R) * (2 * np.random.randint(0, 2) - 1)
        self.bond_change(a, (a + radius_bond) % self.N)

        # Move to next in order
        self.order_idx += 1
        if self.order_idx >= self.N:
            self.order_idx = 0

        # Periodic validation (only if enabled)
        self.perform_periodic_validation()

    def demon_reverse(self):
        """
            Reverse order with random radius selection
        """
        a = self.rev_order[self.rev_order_idx]

        # Generate random radius for bond change
        radius_bond = np.random.randint(0, self.R) * (2 * np.random.randint(0, 2) - 1)
        self.bond_change(a, (a + radius_bond) % self.N)

        # Generate random radius for spin flip
        radius_spin = np.random.randint(0, self.R) * (2 * np.random.randint(0, 2) - 1)
        bonds_changed = self.spin_flip(a, (a + radius_spin) % self.N)
        if bonds_changed:
            self.update_bonds_incremental(a)

        # Move to next in order
        self.rev_order_idx += 1
        if self.rev_order_idx >= self.N:
            self.rev_order_idx = 0

        # Periodic validation (only if enabled)
        self.perform_periodic_validation()
