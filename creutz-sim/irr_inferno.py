import numpy as np
import random
import math

from scipy.special import loggamma as logg

# Use high-precision log calculations
Sk = lambda N, K: logg(K + N) - logg(K+1) - logg(N)
Su = lambda N, N0, Nx: logg(N+1) + np.log(2**(N0)) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1))
Su0 = lambda N, N0, Nx: logg(N+1) + np.log(2**(N0+1)) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1))

class irrInferno:
    """
        Optimized irrInferno class with roundoff error prevention
        Uses truly random radius selection (no pre-generated arrays)
    """

    def __init__(self, N, R):
        a = np.arange(N)
        np.random.shuffle(a)

        total_energy = 2*N

        self.N = N
        self.order = a
        self.rev_order = np.flip(a)
        self.R = R

        self.lattice = np.concatenate((np.ones(N//2, dtype=np.int8),
                                       (-1)*np.ones(N//2, dtype=np.int8)))
        self.bonds = np.ones(N, dtype=np.int8)*(-1)
        self.bonds[[N//2-1, -1]] = 1

        # Initialize bond counts incrementally
        self.bond_count = np.array([N-2, 0, 2], dtype=np.int64)

        self.E_lattice = np.sum(self.bonds, dtype=np.int64)
        self.d_energy = total_energy - self.E_lattice

        # More efficient energy distribution
        result = np.zeros(N, dtype=np.int64)
        indices = np.random.randint(0, N, size=self.d_energy)
        np.add.at(result, indices, 1)

        self.E_demon = result
        self.E_demon_sum = np.int64(self.d_energy)
        self.E_total = self.E_lattice + self.E_demon_sum

        # Store initial values for validation
        self._initial_total_energy = np.int64(total_energy)
        self._check_counter = 0
        self._check_interval = N  # Check every sweep

        # Indices for rolling
        self.order_idx = 0
        self.rev_order_idx = 0

    def validate_energy_conservation(self):
        """Periodic energy conservation check"""
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
        """Periodic bond count validation"""
        actual_counts = np.bincount(self.bonds + 1, minlength=3).astype(np.int64)
        if not np.array_equal(actual_counts, self.bond_count):
            print(f"WARNING: Bond count drift detected!")
            print(f"  Tracked: {self.bond_count}")
            print(f"  Actual: {actual_counts}")

            self.bond_count = actual_counts
            return False
        return True

    def spin_flip(self, a, i):
        """Attempt to flip the spin - all integer arithmetic"""
        s = self.lattice[a]
        d = self.E_demon[i]

        nb = (self.lattice[(a+1) % self.N] * abs(self.bonds[a]) +
              self.lattice[(a-1) % self.N] * abs(self.bonds[(a-1) % self.N]))

        cost = 2 * s * nb
        bonds_changed = False

        if cost < 0 or cost <= d:
            s *= -1
            self.E_demon[i] -= cost
            self.E_demon_sum -= cost
            self.d_energy -= cost
            self.E_lattice += cost
            self.lattice[a] = s
            bonds_changed = True

        return bonds_changed

    def update_bonds_incremental(self, a):
        """Update bonds with careful integer counting"""
        # Update right bond
        if self.bonds[a] != 0:
            old_bond = self.bonds[a]
            new_bond = np.int8(-1 if self.lattice[a] == self.lattice[(a+1) % self.N] else 1)

            if old_bond != new_bond:
                self.bonds[a] = new_bond
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
        """Attempt to change the bond - integer arithmetic only"""
        s = self.lattice[a]
        b = self.bonds[a]
        d = self.E_demon[i]
        n = self.lattice[(a+1) % self.N]

        cost = -1 if s == n else 1

        if b == 0 and d - cost >= 0:
            self.E_lattice += cost
            self.E_demon[i] -= cost
            self.E_demon_sum -= cost
            self.d_energy -= cost
            self.bonds[a] = np.int8(cost)

            # Update bond_count: broken -> aligned/misaligned
            self.bond_count[1] -= 1
            self.bond_count[0 if cost == -1 else 2] += 1

        elif b != 0 and d + cost >= 0:
            self.E_lattice -= cost
            self.E_demon[i] += cost
            self.E_demon_sum += cost
            self.d_energy += cost

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
            Move the demon with truly random radius selection
            This version generates random radius on-the-fly for irreversibility
        """
        a = self.order[self.order_idx]

        # Generate random radius and direction for spin flip
        radius_spin = np.random.randint(0, self.R) * np.random.choice([-1, 1])
        bonds_changed = self.spin_flip(a, (a + radius_spin) % self.N)
        if bonds_changed:
            self.update_bonds_incremental(a)

        # Generate random radius and direction for bond change
        radius_bond = np.random.randint(0, self.R) * np.random.choice([-1, 1])
        self.bond_change(a, (a + radius_bond) % self.N)

        # Move to next in order
        self.order_idx = (self.order_idx + 1) % self.N

        # Periodic validation
        self._check_counter += 1
        if self._check_counter >= self._check_interval:
            self._check_counter = 0
            self.validate_energy_conservation()
            self.validate_bond_counts()

    def demon_reverse(self):
        """
            Reverse order with random radius selection
        """
        a = self.rev_order[self.rev_order_idx]

        # Generate random radius for bond change
        radius_bond = np.random.randint(0, self.R) * np.random.choice([-1, 1])
        self.bond_change(a, (a + radius_bond) % self.N)

        # Generate random radius for spin flip
        radius_spin = np.random.randint(0, self.R) * np.random.choice([-1, 1])
        bonds_changed = self.spin_flip(a, (a + radius_spin) % self.N)
        if bonds_changed:
            self.update_bonds_incremental(a)

        # Move to next in order
        self.rev_order_idx = (self.rev_order_idx + 1) % self.N

        # Periodic validation
        self._check_counter += 1
        if self._check_counter >= self._check_interval:
            self._check_counter = 0
            self.validate_energy_conservation()
            self.validate_bond_counts()

    def get_validated_state(self):
        """Return current state with validation"""
        actual_lattice = np.sum(self.bonds, dtype=np.int64)
        actual_demon = np.sum(self.E_demon, dtype=np.int64)
        actual_bond_counts = np.bincount(self.bonds + 1, minlength=3).astype(np.int64)

        return {
            'E_lattice': actual_lattice,
            'E_demon_sum': actual_demon,
            'E_total': actual_lattice + actual_demon,
            'bond_count': actual_bond_counts,
            'd_energy': actual_demon
        }
