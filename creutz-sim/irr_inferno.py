"""Irreversible Creutz Demon Simulation (irrInferno)

Implements an irreversible version of the Creutz demon algorithm for 1D Ising lattice.
Unlike the reversible Inferno class, this generates truly random demon-coupling radii
on-the-fly during each move, making the dynamics irreversible.

Key differences from reversible version:
- Generates new random radii for each spin flip and bond change
- No pre-computed radius arrays (radius_spin, radius_bond)
- Cannot reverse back to initial state even with reversed order
- Used to study thermodynamic irreversibility

Inherits validation, energy conservation, and bond operations from SimulationBase.
"""

import numpy as np
import random
import math

from sim_utils import Sk, Su, Su0, SimulationBase


class irrInferno(SimulationBase):
    """
    Irreversible Creutz demon simulation with on-the-fly random radius generation.
    
    This class implements the irreversible version of the Creutz demon algorithm.
    Unlike the reversible Inferno class which uses pre-computed radius arrays,
    this generates new random radii for each move, ensuring true irreversibility.
    
    The random radii are generated as: radius = randint(0, R) * (±1)
    where ±1 is randomly chosen, giving radii in range [-R+1, R-1].
    
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
        Perform one irreversible Monte Carlo move with random demon coupling.
        
        This is the forward phase move. For each lattice site (in order):
        1. Generate random radius and attempt spin flip
        2. Update bond states if spin was flipped
        3. Generate new random radius and attempt bond change
        
        The random radius generation ensures this move cannot be exactly
        reversed by demon_reverse(), creating irreversible dynamics.
        
        Updates order_idx to cycle through lattice sites.
        Performs periodic validation if enabled.
        """
        a = self.order[self.order_idx]

        # Generate random radius in range [-R+1, R-1] with random sign
        # This is regenerated each call, ensuring irreversibility
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
        Perform one irreversible Monte Carlo move in reverse order.
        
        This is the reverse phase move. Uses reversed traversal order but:
        1. Generates NEW random radii (not the same as demon_move)
        2. Performs bond change BEFORE spin flip (opposite order)
        
        Despite using reversed order, the newly generated random radii mean
        this does NOT reverse the forward phase - the dynamics remain irreversible.
        
        Updates rev_order_idx to cycle through lattice sites.
        Performs periodic validation if enabled.
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
