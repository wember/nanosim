import numpy as np
import random
import math

from sim_utils import Sk, Su, Su0, SimulationBase


class Inferno(SimulationBase):
    """
        Optimized Inferno class with roundoff error prevention
    """

    def __init__(self, N, R, validate_mode='off'):
        """
        Initialize Inferno simulation.
        
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
        self.radius_spin = np.random.randint(0, R, size=N)*np.random.choice([-1, 1], size=N)
        self.rev_radius_spin = np.flip(self.radius_spin)
        self.radius_bond = np.flip(self.radius_spin)
        self.rev_radius_bond = self.radius_spin.copy()

        # Initialize lattice, bonds, and bond counts using base class helper
        self.lattice, self.bonds, self.bond_count = SimulationBase.initialize_lattice(N)

        self.E_lattice = np.sum(self.bonds, dtype=np.int64)  # Use int64 for sum
        
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
        self.radius_spin_idx = 0
        self.radius_bond_idx = 0
        self.rev_radius_spin_idx = 0
        self.rev_radius_bond_idx = 0

        # Pre-compute neighbor indices using base class helper
        self.right_neighbor, self.left_neighbor = SimulationBase.setup_neighbor_arrays(N)

    def demon_move(self):
        """
            Move the demon with periodic validation
        """
        a = self.order[self.order_idx]

        # Attempt to flip spin
        bonds_changed = self.spin_flip(a, (a + self.radius_spin[self.radius_spin_idx]) % self.N)
        if bonds_changed:
            self.update_bonds_incremental(a)

        self.radius_spin_idx += 1
        if self.radius_spin_idx >= self.N:
            self.radius_spin_idx = 0

        # Attempt to change bond
        self.bond_change(a, (a + self.radius_bond[self.radius_bond_idx]) % self.N)
        self.radius_bond_idx += 1
        if self.radius_bond_idx >= self.N:
            self.radius_bond_idx = 0

        # Move to next in order
        self.order_idx += 1
        if self.order_idx >= self.N:
            self.order_idx = 0

        # Periodic validation (only if enabled)
        self.perform_periodic_validation()

    def demon_reverse(self):
        """
            Reverse order with periodic validation
        """
        a = self.rev_order[self.rev_order_idx]

        # Attempt to change bond
        self.bond_change(a, (a + self.rev_radius_bond[self.rev_radius_bond_idx]) % self.N)
        self.rev_radius_bond_idx += 1
        if self.rev_radius_bond_idx >= self.N:
            self.rev_radius_bond_idx = 0

        # Attempt to flip spin
        bonds_changed = self.spin_flip(a, (a + self.rev_radius_spin[self.rev_radius_spin_idx]) % self.N)
        if bonds_changed:
            self.update_bonds_incremental(a)

        self.rev_radius_spin_idx += 1
        if self.rev_radius_spin_idx >= self.N:
            self.rev_radius_spin_idx = 0

        # Move to next in order
        self.rev_order_idx += 1
        if self.rev_order_idx >= self.N:
            self.rev_order_idx = 0

        # Periodic validation (only if enabled)
        self.perform_periodic_validation()
