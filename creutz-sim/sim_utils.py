"""
Shared utilities for Monte Carlo simulations.

This module contains common functionality shared between Inferno and irrInferno
classes, including entropy calculations, validation methods, and initialization
helpers.
"""

import numpy as np
from scipy.special import loggamma as logg


# Entropy calculation functions (high-precision log calculations)
Sk = lambda N, K: logg(K + N) - logg(K+1) - logg(N)
Su = lambda N, N0, Nx: logg(N+1) + np.log(2**(N0)) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1))
Su0 = lambda N, N0, Nx: logg(N+1) + np.log(2**(N0+1)) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1))


class SimulationBase:
    """
    Base class for Monte Carlo simulations with common validation and
    initialization methods.
    
    This class provides shared functionality for energy conservation,
    bond count validation, and state retrieval that is common to both
    reversible (Inferno) and irreversible (irrInferno) simulations.
    """
    
    @staticmethod
    def initialize_order_arrays(N):
        """
        Initialize randomized order arrays for lattice traversal.
        
        Creates a shuffled array for forward traversal and its reverse
        for backward traversal (used in reversibility tests).
        
        Args:
            N: Number of lattice sites
            
        Returns:
            tuple: (order, rev_order) - forward and reverse traversal arrays
        """
        order = np.arange(N)
        np.random.shuffle(order)
        rev_order = np.flip(order)
        
        return order, rev_order
    
    @staticmethod
    def initialize_lattice(N):
        """
        Initialize a lattice with half spins up (+1) and half down (-1).
        
        Args:
            N: Number of lattice sites
            
        Returns:
            tuple: (lattice, bonds, bond_count) arrays
        """
        lattice = np.concatenate((np.ones(N//2, dtype=np.int8),
                                 (-1)*np.ones(N//2, dtype=np.int8)))
        bonds = np.ones(N, dtype=np.int8)*(-1)
        bonds[[N//2-1, -1]] = 1
        bond_count = np.array([N-2, 0, 2], dtype=np.int64)
        
        return lattice, bonds, bond_count
    
    @staticmethod
    def initialize_demon_energy(N, total_energy, E_lattice):
        """
        Initialize demon energy array by distributing remaining energy randomly.
        
        Args:
            N: Number of lattice sites (demons)
            total_energy: Total energy to conserve
            E_lattice: Current lattice energy
            
        Returns:
            tuple: (E_demon array, E_demon_sum, d_energy)
        """
        d_energy = total_energy - E_lattice
        result = np.zeros(N, dtype=np.int64)
        indices = np.random.randint(0, N, size=d_energy)
        np.add.at(result, indices, 1)
        
        return result, np.int64(d_energy), d_energy
    
    @staticmethod
    def setup_neighbor_arrays(N):
        """
        Pre-compute neighbor index arrays for periodic boundaries.
        
        This optimization replaces modulo operations with direct array lookups.
        
        Args:
            N: Number of lattice sites
            
        Returns:
            tuple: (right_neighbor, left_neighbor) arrays
        """
        right_neighbor = np.arange(1, N+1, dtype=np.int32) % N
        left_neighbor = np.arange(-1, N-1, dtype=np.int32) % N
        
        return right_neighbor, left_neighbor
    
    def setup_validation_mode(self, N, validate_mode, initial_total_energy):
        """
        Configure validation parameters based on mode.
        
        Args:
            N: Number of lattice sites (used to calculate check intervals)
            validate_mode: 'off', 'periodic', or 'frequent'
            initial_total_energy: Total energy for conservation checks
        """
        self._initial_total_energy = np.int64(initial_total_energy)
        self._validate_mode = validate_mode
        self._check_counter = 0
        
        # Configure validation interval based on mode
        if validate_mode == 'off':
            self._check_interval = float('inf')  # Never validate automatically
        elif validate_mode == 'frequent':
            self._check_interval = N  # Every sweep (debug mode)
        else:  # 'periodic' or default
            self._check_interval = 100 * N  # Every 100 sweeps
    
    def perform_periodic_validation(self):
        """
        Perform periodic validation check if enabled.
        
        Should be called at the end of demon_move() and demon_reverse().
        Increments counter and validates when interval is reached.
        """
        if self._validate_mode != 'off':
            self._check_counter += 1
            if self._check_counter >= self._check_interval:
                self._check_counter = 0
                self.validate_energy_conservation()
                self.validate_bond_counts()
    
    def validate_energy_conservation(self):
        """
        Periodic energy conservation check to catch drift.
        
        Recalculates total energy from scratch and compares with tracked value.
        If drift is detected, corrects the cached values and prints warning.
        
        Returns:
            bool: True if energy is conserved, False if drift was detected and corrected
        """
        current_total = self.E_lattice + self.E_demon_sum
        if current_total != self._initial_total_energy:
            # Recalculate from scratch
            actual_lattice = np.sum(self.bonds, dtype=np.int64)
            actual_demon = np.sum(self.E_demon, dtype=np.int64)

            print(f"WARNING: Energy drift detected!")
            print(f"  Expected total: {self._initial_total_energy}")
            print(f"  Tracked total: {current_total}")
            print(f"  Actual total: {actual_lattice + actual_demon}")
            print(f"  Drift: {current_total - self._initial_total_energy}")

            # Correct the cached values
            self.E_lattice = actual_lattice
            self.E_demon_sum = actual_demon
            self.d_energy = self.E_demon_sum

            return False
        return True
    
    def validate_bond_counts(self):
        """
        Periodic bond count validation to catch drift.
        
        Recalculates bond counts from scratch and compares with tracked values.
        If drift is detected, corrects the cached values and prints warning.
        
        Returns:
            bool: True if bond counts match, False if drift was detected and corrected
        """
        actual_counts = np.bincount(self.bonds + 1, minlength=3).astype(np.int64)
        if not np.array_equal(actual_counts, self.bond_count):
            print(f"WARNING: Bond count drift detected!")
            print(f"  Tracked: {self.bond_count}")
            print(f"  Actual: {actual_counts}")

            # Correct the cached values
            self.bond_count = actual_counts
            return False
        return True
    
    def get_validated_state(self):
        """
        Return current state with validation.
        
        Always recalculates values from scratch rather than using cached values.
        Use this method when you need guaranteed accurate values.
        
        Returns:
            dict: Dictionary with keys:
                - E_lattice: Lattice energy
                - E_demon_sum: Total demon energy
                - E_total: Total energy (lattice + demon)
                - bond_count: Array of bond counts [broken, neutral, anti-aligned]
                - d_energy: Demon energy (same as E_demon_sum)
        """
        # Force validation - recalculate from scratch
        actual_lattice = np.sum(self.bonds, dtype=np.int64)
        actual_demon = np.sum(self.E_demon, dtype=np.int64)
        actual_bond_counts = np.bincount(self.bonds + 1, minlength=3).astype(np.int64)

        # Return validated values
        return {
            'E_lattice': actual_lattice,
            'E_demon_sum': actual_demon,
            'E_total': actual_lattice + actual_demon,
            'bond_count': actual_bond_counts,
            'd_energy': actual_demon
        }
    
    def spin_flip(self, a, i):
        """
        Attempt to flip the spin at site a using demon i.
        
        Uses all integer arithmetic to prevent floating point drift.
        
        Args:
            a: Lattice site index
            i: Demon index
            
        Returns:
            bool: True if spin was flipped (bonds may have changed), False otherwise
        """
        s = self.lattice[a]
        d = self.E_demon[i]

        # Calculate energy cost (integer arithmetic)
        nb = (self.lattice[self.right_neighbor[a]] * abs(self.bonds[a]) +
              self.lattice[self.left_neighbor[a]] * abs(self.bonds[self.left_neighbor[a]]))

        cost = 2 * s * nb  # Always an integer
        bonds_changed = False

        if cost < 0 or cost <= d:
            s *= -1
            # Update energies using integer arithmetic only
            self.E_demon[i] -= cost
            self.E_demon_sum -= cost
            self.d_energy -= cost
            self.E_lattice += cost
            self.lattice[a] = s
            bonds_changed = True

        return bonds_changed
    
    def update_bonds_incremental(self, a):
        """
        Update bond states after spin flip at site a.
        
        Updates both the right bond (at site a) and left bond (at left_neighbor[a])
        using careful integer counting to maintain bond_count accuracy.
        
        Args:
            a: Lattice site index where spin was flipped
        """
        # Update right bond
        if self.bonds[a] != 0:
            old_bond = self.bonds[a]
            new_bond = np.int8(-1 if self.lattice[a] == self.lattice[self.right_neighbor[a]] else 1)

            if old_bond != new_bond:
                self.bonds[a] = new_bond
                # Update counts with explicit indexing
                old_idx = 0 if old_bond == -1 else 2
                new_idx = 0 if new_bond == -1 else 2
                self.bond_count[old_idx] -= 1
                self.bond_count[new_idx] += 1

        # Update left bond
        left_idx = self.left_neighbor[a]
        if self.bonds[left_idx] != 0:
            old_bond = self.bonds[left_idx]
            new_bond = np.int8(-1 if self.lattice[a] == self.lattice[left_idx] else 1)

            if old_bond != new_bond:
                self.bonds[left_idx] = new_bond
                old_idx = 0 if old_bond == -1 else 2
                new_idx = 0 if new_bond == -1 else 2
                self.bond_count[old_idx] -= 1
                self.bond_count[new_idx] += 1
    
    def bond_change(self, a, i):
        """
        Attempt to change the bond state at site a using demon i.
        
        Can either break an existing bond (0 → ±1) or create a bond (±1 → 0).
        Uses integer arithmetic only to prevent drift.
        
        Args:
            a: Lattice site index
            i: Demon index
        """
        s = self.lattice[a]
        b = self.bonds[a]
        d = self.E_demon[i]
        n = self.lattice[self.right_neighbor[a]]

        cost = -1 if s == n else 1  # Always ±1

        if b == 0 and d - cost >= 0:
            # Update energies
            self.E_lattice += cost
            self.E_demon[i] -= cost
            self.E_demon_sum -= cost
            self.d_energy -= cost
            self.bonds[a] = np.int8(cost)

            # Update bond_count: broken → aligned/misaligned
            self.bond_count[1] -= 1
            self.bond_count[0 if cost == -1 else 2] += 1

        elif b != 0 and d + cost >= 0:
            # Update energies
            self.E_lattice -= cost
            self.E_demon[i] += cost
            self.E_demon_sum += cost
            self.d_energy += cost

            # Update bond_count
            old_idx = 0 if self.bonds[a] == -1 else 2
            self.bond_count[old_idx] -= 1
            self.bond_count[1] += 1

            self.bonds[a] = 0

        # Update left neighbor bond
        left_idx = self.left_neighbor[a]
        if self.bonds[left_idx] != 0:
            old_bond = self.bonds[left_idx]
            new_bond = np.int8(-1 if self.lattice[a] == self.lattice[left_idx] else 1)

            if old_bond != new_bond:
                self.bonds[left_idx] = new_bond
                old_idx = 0 if old_bond == -1 else 2
                new_idx = 0 if new_bond == -1 else 2
                self.bond_count[old_idx] -= 1
                self.bond_count[new_idx] += 1
