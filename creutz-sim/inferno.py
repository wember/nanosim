"""Reversible Creutz Demon Simulation (Inferno)

Implements the reversible version of the Creutz demon algorithm for 1D Ising lattice.
This class uses pre-computed radius arrays to enable exact time-reversibility testing.

Key features for reversibility:
- Pre-computed random walk patterns (radius_spin, radius_bond)
- Reversed arrays (rev_radius_spin, rev_radius_bond) for backward traversal
- Forward (demon_move) and reverse (demon_reverse) phases that mirror each other
- Can return to initial state after forward then reverse traversal

The reversibility property allows testing the second law of thermodynamics by
observing whether entropy returns to its initial value during reverse phase.

Inherits validation, energy conservation, and bond operations from SimulationBase.
"""

import numpy as np
from sim_utils import SimulationBase, Sk, Su, Su0  # noqa: F401


class Inferno(SimulationBase):
    """
    Reversible Creutz demon simulation with pre-computed radius arrays.

    This class implements the reversible version of the Creutz demon algorithm.
    The key to reversibility is using fixed random walk patterns (radius arrays)
    and their exact reversals, allowing the system to be run forward then backward
    to return to the initial state.

    Reversibility mechanism:
    - radius_spin: Pre-computed random radii for spin flips (forward)
    - rev_radius_spin: Flipped radius_spin for reverse traversal
    - radius_bond: Flipped radius_spin (for bond changes, forward)
    - rev_radius_bond: Copy of radius_spin (for bond changes, reverse)

    The specific pattern ensures demon_reverse() exactly reverses demon_move().

    Attributes:
        N: Number of lattice sites
        R: Maximum demon-coupling radius
        lattice: Spin configuration array (±1)
        bonds: Bond state array (-1=aligned, 0=broken, 1=anti-aligned)
        bond_count: Array [N0, N1, Nx] tracking bond counts
        E_demon: Array of demon energies (one per site)
        E_lattice: Total lattice energy
        E_total: Conserved total energy (lattice + demon)
        order: Randomized forward traversal order
        rev_order: Reversed traversal order (flipped order array)
        radius_spin, radius_bond: Forward phase radius arrays
        rev_radius_spin, rev_radius_bond: Reverse phase radius arrays
    """

    def __init__(self, N, R, validate_mode="off"):
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
        total_energy = 2 * N

        self.N = N
        self.order, self.rev_order = SimulationBase.initialize_order_arrays(N)

        # Pre-compute radius arrays for reversibility
        # radius_spin: random radii in range [-R+1, R-1] for forward spin flips
        self.radius_spin = np.random.randint(0, R, size=N) * np.random.choice(
            [-1, 1], size=N
        )
        # rev_radius_spin: flipped radius_spin for reverse spin flips
        self.rev_radius_spin = np.flip(self.radius_spin)
        # radius_bond: flipped radius_spin for forward bond changes
        self.radius_bond = np.flip(self.radius_spin)
        # rev_radius_bond: copy of radius_spin for reverse bond changes
        # This specific pattern ensures exact reversibility
        self.rev_radius_bond = self.radius_spin.copy()

        # Initialize lattice, bonds, and bond counts using base class helper
        self.lattice, self.bonds, self.bond_count = SimulationBase.initialize_lattice(N)

        self.E_lattice = np.sum(self.bonds, dtype=np.int64)  # Use int64 for sum

        # Initialize demon energy using base class helper
        self.E_demon, self.E_demon_sum, self.d_energy = (
            SimulationBase.initialize_demon_energy(N, total_energy, self.E_lattice)
        )
        self.E_total = self.E_lattice + self.E_demon_sum

        # Setup validation using base class helper
        self.setup_validation_mode(N, validate_mode, total_energy)

        # Index tracking for cycling through arrays
        # These wrap around to 0 after reaching N
        self.order_idx = 0
        self.rev_order_idx = 0
        self.radius_spin_idx = 0
        self.radius_bond_idx = 0
        self.rev_radius_spin_idx = 0
        self.rev_radius_bond_idx = 0

        # Pre-compute neighbor indices using base class helper
        self.right_neighbor, self.left_neighbor = SimulationBase.setup_neighbor_arrays(
            N
        )

    def demon_move(self):
        """
        Perform one reversible Monte Carlo move (forward phase).

        For each lattice site (in order):
        1. Attempt spin flip using pre-computed radius from radius_spin
        2. Update bond states if spin was flipped
        3. Attempt bond change using pre-computed radius from radius_bond

        Uses fixed radius arrays to enable exact reversal by demon_reverse().
        The sequence (spin flip → bond change) is mirrored in reverse phase.

        Updates all index counters and performs periodic validation if enabled.
        """
        a = self.order[self.order_idx]

        # Attempt to flip spin
        bonds_changed = self.spin_flip(
            a, (a + self.radius_spin[self.radius_spin_idx]) % self.N
        )
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
        Perform one reversible Monte Carlo move (reverse phase).

        Exactly reverses demon_move() by:
        1. Using reversed traversal order (rev_order)
        2. Performing bond change BEFORE spin flip (opposite sequence)
        3. Using reversed radius arrays (rev_radius_bond, rev_radius_spin)

        The specific combination of reversed order, reversed sequence, and
        reversed radii ensures this exactly undoes demon_move(), returning
        the system to its previous state.

        Updates all reverse index counters and performs periodic validation if enabled.
        """
        a = self.rev_order[self.rev_order_idx]

        # Attempt to change bond
        self.bond_change(
            a, (a + self.rev_radius_bond[self.rev_radius_bond_idx]) % self.N
        )
        self.rev_radius_bond_idx += 1
        if self.rev_radius_bond_idx >= self.N:
            self.rev_radius_bond_idx = 0

        # Attempt to flip spin
        bonds_changed = self.spin_flip(
            a, (a + self.rev_radius_spin[self.rev_radius_spin_idx]) % self.N
        )
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
