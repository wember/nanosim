import numpy as np
import random
import math

from scipy.special import loggamma as logg

# Use high-precision log calculations
Sk = lambda N, K: logg(K + N) - logg(K+1) - logg(N)
Su = lambda N, N0, Nx: logg(N+1) + np.log(2**(N0)) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1))
Su0 = lambda N, N0, Nx: logg(N+1) + np.log(2**(N0+1)) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1))

class Inferno:
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
        a = np.arange(N)
        np.random.shuffle(a)

        total_energy = 2*N

        self.N = N
        self.order = a
        self.rev_order = np.flip(a)
        self.radius_spin = np.random.randint(0, R, size=N)*np.random.choice([-1, 1], size=N)
        self.rev_radius_spin = np.flip(self.radius_spin)
        self.radius_bond = np.flip(self.radius_spin)
        self.rev_radius_bond = self.radius_spin.copy()

        self.lattice = np.concatenate((np.ones(N//2, dtype=np.int8),
                                       (-1)*np.ones(N//2, dtype=np.int8)))
        self.bonds = np.ones(N, dtype=np.int8)*(-1)
        self.bonds[[N//2-1, -1]] = 1

        # Initialize bond counts incrementally
        self.bond_count = np.array([N-2, 0, 2], dtype=np.int64)

        self.E_lattice = np.sum(self.bonds, dtype=np.int64)  # Use int64 for sum
        self.d_energy = total_energy - self.E_lattice

        # Energy distribution
        result = np.zeros(N, dtype=np.int64)  # Use int64 for demon energies
        indices = np.random.randint(0, N, size=self.d_energy)
        np.add.at(result, indices, 1)

        self.E_demon = result
        self.E_demon_sum = np.int64(self.d_energy)  # Explicit int64
        self.E_total = self.E_lattice + self.E_demon_sum

        # Store initial values for validation
        self._initial_total_energy = np.int64(total_energy)
        self._validate_mode = validate_mode
        self._check_counter = 0
        
        # Configure validation interval based on mode
        if validate_mode == 'off':
            self._check_interval = float('inf')  # Never validate automatically
        elif validate_mode == 'frequent':
            self._check_interval = N  # Every sweep (debug mode)
        else:  # 'periodic' or default
            self._check_interval = 100 * N  # Every 100 sweeps

        # Indices for rolling
        self.order_idx = 0
        self.rev_order_idx = 0
        self.radius_spin_idx = 0
        self.radius_bond_idx = 0
        self.rev_radius_spin_idx = 0
        self.rev_radius_bond_idx = 0

    def validate_energy_conservation(self):
        """
            Periodic energy conservation check to catch drift
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
            Periodic bond count validation to catch drift
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

    def spin_flip(self, a, i):
        """
            Attempt to flip the spin - all integer arithmetic
        """
        s = self.lattice[a]
        d = self.E_demon[i]

        # All operations in integer arithmetic
        nb = (self.lattice[(a+1) % self.N] * abs(self.bonds[a]) +
              self.lattice[(a-1) % self.N] * abs(self.bonds[(a-1) % self.N]))

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
            Update bonds with careful integer counting
        """
        # Update right bond
        if self.bonds[a] != 0:
            old_bond = self.bonds[a]
            new_bond = np.int8(-1 if self.lattice[a] == self.lattice[(a+1) % self.N] else 1)

            if old_bond != new_bond:
                self.bonds[a] = new_bond
                # Update counts with explicit indexing
                old_idx = 0 if old_bond == -1 else 2
                new_idx = 0 if new_bond == -1 else 2
                self.bond_count[old_idx] -= 1
                self.bond_count[new_idx] += 1

        # Update left bond
        left_idx = (a-1) % self.N
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
            Attempt to change the bond - integer arithmetic only
        """
        s = self.lattice[a]
        b = self.bonds[a]
        d = self.E_demon[i]
        n = self.lattice[(a+1) % self.N]

        cost = -1 if s == n else 1  # Always ±1

        if b == 0 and d - cost >= 0:
            # Update energies
            self.E_lattice += cost
            self.E_demon[i] -= cost
            self.E_demon_sum -= cost
            self.d_energy -= cost
            self.bonds[a] = np.int8(cost)

            # Update bond_count: broken -> aligned/misaligned
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
        left_idx = (a-1) % self.N
        if self.bonds[left_idx] != 0:
            old_bond = self.bonds[left_idx]
            new_bond = np.int8(-1 if self.lattice[a] == self.lattice[left_idx] else 1)

            if old_bond != new_bond:
                self.bonds[left_idx] = new_bond
                old_idx = 0 if old_bond == -1 else 2
                new_idx = 0 if new_bond == -1 else 2
                self.bond_count[old_idx] -= 1
                self.bond_count[new_idx] += 1

    def demon_move(self):
        """
            Move the demon with periodic validation
        """
        a = self.order[self.order_idx]

        # Attempt to flip spin
        bonds_changed = self.spin_flip(a, (a + self.radius_spin[self.radius_spin_idx]) % self.N)
        if bonds_changed:
            self.update_bonds_incremental(a)

        self.radius_spin_idx = (self.radius_spin_idx + 1) % self.N

        # Attempt to change bond
        self.bond_change(a, (a + self.radius_bond[self.radius_bond_idx]) % self.N)
        self.radius_bond_idx = (self.radius_bond_idx + 1) % self.N

        # Move to next in order
        self.order_idx = (self.order_idx + 1) % self.N

        # Periodic validation (only if enabled)
        if self._validate_mode != 'off':
            self._check_counter += 1
            if self._check_counter >= self._check_interval:
                self._check_counter = 0
                self.validate_energy_conservation()
                self.validate_bond_counts()

    def demon_reverse(self):
        """
            Reverse order with periodic validation
        """
        a = self.rev_order[self.rev_order_idx]

        # Attempt to change bond
        self.bond_change(a, (a + self.rev_radius_bond[self.rev_radius_bond_idx]) % self.N)
        self.rev_radius_bond_idx = (self.rev_radius_bond_idx + 1) % self.N

        # Attempt to flip spin
        bonds_changed = self.spin_flip(a, (a + self.rev_radius_spin[self.rev_radius_spin_idx]) % self.N)
        if bonds_changed:
            self.update_bonds_incremental(a)

        self.rev_radius_spin_idx = (self.rev_radius_spin_idx + 1) % self.N

        # Move to next in order
        self.rev_order_idx = (self.rev_order_idx + 1) % self.N

        # Periodic validation (only if enabled)
        if self._validate_mode != 'off':
            self._check_counter += 1
            if self._check_counter >= self._check_interval:
                self._check_counter = 0
                self.validate_energy_conservation()
                self.validate_bond_counts()

    def get_validated_state(self):
        """
            Return current state with validation
        """
        # Force validation
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
