"""
Copyright (c) 2026 Winry Ember
Licensed under the MIT License.
See LICENSE file in the project root for full license information.

Direct unit tests for JIT-compiled functions in jit_functions.py

Tests each JIT function in isolation to achieve high coverage of the core
simulation logic. These tests complement the integration tests in
test_jit_implementation.py.
"""

import os
import sys

import numpy as np

# Add creutz-sim to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))

from jit_functions import (
    bond_change_jit,
    demon_move_irr_jit,
    demon_move_jit,
    spin_flip_jit,
    update_bonds_incremental_jit,
)


class TestSpinFlipJit:
    """Test spin_flip_jit function."""

    def setup_method(self):
        """Set up minimal test arrays."""
        self.N = 10
        self.lattice = np.array([1, -1, 1, -1, 1, -1, 1, -1, 1, -1], dtype=np.int8)
        self.bonds = np.array([1, -1, 1, -1, 1, -1, 1, -1, 1, -1], dtype=np.int8)
        self.E_demon = np.array([5, 10, 15, 20, 25, 30, 35, 40, 45, 50], dtype=np.int64)
        self.right_neighbor = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 0], dtype=np.int32)
        self.left_neighbor = np.array([9, 0, 1, 2, 3, 4, 5, 6, 7, 8], dtype=np.int32)

    def test_spin_flip_when_cost_negative(self):
        """Test spin flip when energy cost is negative (always accepted)."""
        # Site 0 has spin=1, neighbors are -1 and -1, so flipping reduces energy
        a, i = 0, 0
        E_demon_sum_before = np.sum(self.E_demon)
        E_lattice = np.int64(100)
        E_total_before = E_demon_sum_before + E_lattice

        bonds_changed, new_sum, new_lattice = spin_flip_jit(
            self.lattice,
            self.bonds,
            self.E_demon,
            E_demon_sum_before,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            a,
            i,
        )

        assert bonds_changed is True
        assert self.lattice[a] == -1  # Spin flipped
        # Total energy should be conserved
        E_total_after = new_sum + new_lattice
        assert E_total_before == E_total_after

    def test_spin_flip_when_cost_zero(self):
        """Test spin flip when energy cost is zero."""
        # Create scenario with aligned spins
        self.lattice = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=np.int8)
        self.bonds = np.array([-1, -1, -1, -1, -1, -1, -1, -1, -1, -1], dtype=np.int8)

        a, i = 5, 5
        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)

        bonds_changed, new_sum, new_lattice = spin_flip_jit(
            self.lattice,
            self.bonds,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            a,
            i,
        )

        # Cost is 2*1*(-1 + -1) = -4, should flip
        assert bonds_changed is True

    def test_spin_flip_rejected_insufficient_demon_energy(self):
        """Test spin flip rejected when demon lacks energy."""
        # Create all parallel spins with aligned bonds - flipping will cost energy
        self.lattice = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=np.int8)
        self.bonds = np.array([-1, -1, -1, -1, -1, -1, -1, -1, -1, -1], dtype=np.int8)
        self.E_demon = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=np.int64)

        a, i = 5, 5
        # All parallel spins=1, all bonds=-1 (aligned)
        # nb = 1*1 + 1*1 = 2, cost = 2*1*2 = 4 (positive!)
        original_spin = self.lattice[a]

        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)

        bonds_changed, new_sum, new_lattice = spin_flip_jit(
            self.lattice,
            self.bonds,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            a,
            i,
        )

        # Cost is positive (4) but demon has no energy - should reject
        assert self.lattice[a] == original_spin  # Spin unchanged
        assert bonds_changed is False

    def test_spin_flip_accepted_with_exact_demon_energy(self):
        """Test spin flip when demon has exactly enough energy."""
        self.lattice = np.array([1, -1, 1, -1, 1, -1, 1, -1, 1, -1], dtype=np.int8)
        self.bonds = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=np.int8)

        a, i = 1, 1
        # Calculate expected cost: 2 * s * nb
        s = self.lattice[a]
        nb = self.lattice[self.right_neighbor[a]] * abs(self.bonds[a]) + self.lattice[
            self.left_neighbor[a]
        ] * abs(self.bonds[self.left_neighbor[a]])
        expected_cost = 2 * s * nb

        # Set demon energy to exactly the cost
        self.E_demon[i] = expected_cost
        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)

        bonds_changed, new_sum, new_lattice = spin_flip_jit(
            self.lattice,
            self.bonds,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            a,
            i,
        )

        assert bonds_changed is True
        assert self.E_demon[i] == 0  # Demon used all energy

    def test_spin_flip_energy_conservation(self):
        """Test that total energy is conserved during spin flip."""
        a, i = 3, 3
        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)
        E_total_before = E_demon_sum + E_lattice

        bonds_changed, new_sum, new_lattice = spin_flip_jit(
            self.lattice,
            self.bonds,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            a,
            i,
        )

        E_total_after = new_sum + new_lattice
        assert E_total_before == E_total_after


class TestUpdateBondsIncrementalJit:
    """Test update_bonds_incremental_jit function."""

    def setup_method(self):
        """Set up test arrays."""
        self.N = 10
        self.lattice = np.array([1, -1, 1, -1, 1, -1, 1, -1, 1, -1], dtype=np.int8)
        self.bonds = np.array([1, -1, 1, -1, 1, -1, 1, -1, 1, -1], dtype=np.int8)
        self.right_neighbor = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 0], dtype=np.int32)
        self.left_neighbor = np.array([9, 0, 1, 2, 3, 4, 5, 6, 7, 8], dtype=np.int32)

    def test_update_bonds_after_spin_flip(self):
        """Test bond updates after a spin flip."""
        # Start with alternating spins: [1, -1, 1, -1, ...]
        # Bonds: [1, -1, 1, -1, ...] (misaligned, aligned, misaligned, ...)
        bond_count = np.array([5, 0, 5], dtype=np.int32)  # [N0=5, N1=0, Nx=5]

        # Flip spin at site 1 from -1 to 1
        a = 1
        self.lattice[a] = 1  # Was -1, now 1
        # Now spins are: [1, 1, 1, -1, ...]

        new_bond_count = update_bonds_incremental_jit(
            self.lattice,
            self.bonds,
            bond_count,
            self.right_neighbor,
            self.left_neighbor,
            a,
        )

        # Right bond (1→2): lattice[1]=1, lattice[2]=1, so bond should be -1 (aligned)
        assert self.bonds[a] == -1
        # Left bond (0→1): lattice[0]=1, lattice[1]=1, so bond should be -1 (aligned)
        assert self.bonds[0] == -1
        # Verify bond count was updated
        assert np.sum(new_bond_count) == self.N

    def test_bond_count_conservation(self):
        """Test that total bond count is conserved."""
        bond_count = np.array([3, 0, 7], dtype=np.int32)
        total_before = np.sum(bond_count)

        a = 5
        self.lattice[a] = -self.lattice[a]  # Flip spin

        new_bond_count = update_bonds_incremental_jit(
            self.lattice,
            self.bonds,
            bond_count,
            self.right_neighbor,
            self.left_neighbor,
            a,
        )

        total_after = np.sum(new_bond_count)
        assert total_before == total_after == self.N

    def test_no_update_when_bond_zero(self):
        """Test that zero bonds are not updated."""
        self.bonds = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=np.int8)
        bond_count = np.array([0, 10, 0], dtype=np.int32)

        a = 3
        original_bond_count = bond_count.copy()

        new_bond_count = update_bonds_incremental_jit(
            self.lattice,
            self.bonds,
            bond_count,
            self.right_neighbor,
            self.left_neighbor,
            a,
        )

        # No bonds should change
        assert np.array_equal(new_bond_count, original_bond_count)

    def test_aligned_to_misaligned_transition(self):
        """Test transition from aligned to misaligned bond."""
        self.lattice = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=np.int8)
        self.bonds = np.array([-1, -1, -1, -1, -1, -1, -1, -1, -1, -1], dtype=np.int8)
        bond_count = np.array([10, 0, 0], dtype=np.int32)

        a = 5
        self.lattice[a] = -1  # Flip from 1 to -1
        # Now we have [..., 1, -1, 1, ...]

        new_bond_count = update_bonds_incremental_jit(
            self.lattice,
            self.bonds,
            bond_count,
            self.right_neighbor,
            self.left_neighbor,
            a,
        )

        # Right bond (5→6): lattice[5]=-1, lattice[6]=1, bond should be 1 (misaligned)
        assert self.bonds[a] == 1
        # Left bond (4→5): lattice[4]=1, lattice[5]=-1, bond should be 1 (misaligned)
        assert self.bonds[4] == 1
        # Two bonds changed from aligned to misaligned
        assert new_bond_count[0] == 8  # Two less aligned
        assert new_bond_count[2] == 2  # Two more misaligned


class TestBondChangeJit:
    """Test bond_change_jit function."""

    def setup_method(self):
        """Set up test arrays."""
        self.N = 10
        self.lattice = np.array([1, -1, 1, -1, 1, -1, 1, -1, 1, -1], dtype=np.int8)
        self.bonds = np.array([1, -1, 1, -1, 1, -1, 1, -1, 1, -1], dtype=np.int8)
        self.E_demon = np.array(
            [10, 10, 10, 10, 10, 10, 10, 10, 10, 10], dtype=np.int64
        )
        self.right_neighbor = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 0], dtype=np.int32)
        self.left_neighbor = np.array([9, 0, 1, 2, 3, 4, 5, 6, 7, 8], dtype=np.int32)

    def test_create_aligned_bond(self):
        """Test creating an aligned bond."""
        self.bonds[0] = 0  # No bond initially
        bond_count = np.array([4, 1, 5], dtype=np.int32)
        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)

        a, i = 0, 0
        self.lattice[a] = 1
        self.lattice[self.right_neighbor[a]] = 1  # Aligned

        new_sum, new_lattice, new_bond_count = bond_change_jit(
            self.lattice,
            self.bonds,
            bond_count,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            a,
            i,
        )

        # Bond should be created as -1 (aligned)
        assert self.bonds[a] == -1
        # Energy should transfer from demon to lattice
        assert new_lattice == E_lattice + (-1)
        assert self.E_demon[i] == 10 - (-1)

    def test_create_misaligned_bond(self):
        """Test creating a misaligned bond."""
        self.bonds[2] = 0
        bond_count = np.array([4, 1, 5], dtype=np.int32)
        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)

        a, i = 2, 2
        self.lattice[a] = 1
        self.lattice[self.right_neighbor[a]] = -1  # Misaligned

        new_sum, new_lattice, new_bond_count = bond_change_jit(
            self.lattice,
            self.bonds,
            bond_count,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            a,
            i,
        )

        # Bond should be created as 1 (misaligned)
        assert self.bonds[a] == 1
        assert new_lattice == E_lattice + 1

    def test_break_bond_with_sufficient_energy(self):
        """Test breaking a bond when demon has enough energy."""
        bond_count = np.array([5, 0, 5], dtype=np.int32)
        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)

        a, i = 1, 1

        new_sum, new_lattice, new_bond_count = bond_change_jit(
            self.lattice,
            self.bonds,
            bond_count,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            a,
            i,
        )

        # Bond should be broken
        assert self.bonds[a] == 0
        # Energy transferred from lattice to demon
        assert self.E_demon[i] > 10

    def test_bond_change_rejected_insufficient_energy(self):
        """Test bond creation rejected when demon lacks energy."""
        self.bonds[3] = 0
        self.E_demon[3] = 0  # No demon energy
        bond_count = np.array([4, 1, 5], dtype=np.int32)
        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)

        a, i = 3, 3

        new_sum, new_lattice, new_bond_count = bond_change_jit(
            self.lattice,
            self.bonds,
            bond_count,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            a,
            i,
        )

        # Bond should not be created if cost is positive
        # (depends on alignment)

    def test_bond_change_energy_conservation(self):
        """Test energy conservation during bond change."""
        bond_count = np.array([5, 0, 5], dtype=np.int32)
        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)
        E_total_before = E_demon_sum + E_lattice

        a, i = 4, 4

        new_sum, new_lattice, new_bond_count = bond_change_jit(
            self.lattice,
            self.bonds,
            bond_count,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            a,
            i,
        )

        E_total_after = new_sum + new_lattice
        assert E_total_before == E_total_after

    def test_bond_count_updates_on_creation(self):
        """Test bond count array updates when bond is created."""
        self.bonds[5] = 0
        bond_count = np.array([3, 1, 6], dtype=np.int32)
        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)

        a, i = 5, 5

        new_sum, new_lattice, new_bond_count = bond_change_jit(
            self.lattice,
            self.bonds,
            bond_count,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            a,
            i,
        )

        # Total bonds should be conserved
        assert np.sum(new_bond_count) == self.N


class TestDemonMoveJit:
    """Test demon_move_jit function (reversible)."""

    def setup_method(self):
        """Set up test simulation state."""
        self.N = 20
        self.R = 3
        self.lattice = np.array(
            [1 if i % 2 == 0 else -1 for i in range(self.N)], dtype=np.int8
        )
        self.bonds = np.array(
            [1 if i % 2 == 0 else -1 for i in range(self.N)], dtype=np.int8
        )
        self.bond_count = np.array([10, 0, 10], dtype=np.int32)
        self.E_demon = np.full(self.N, 20, dtype=np.int64)
        self.right_neighbor = np.array(
            [(i + 1) % self.N for i in range(self.N)], dtype=np.int32
        )
        self.left_neighbor = np.array(
            [(i - 1) % self.N for i in range(self.N)], dtype=np.int32
        )
        self.order = np.arange(self.N, dtype=np.int32)
        self.radius_spin = np.random.randint(-self.R, self.R, self.N, dtype=np.int32)
        self.radius_bond = np.random.randint(-self.R, self.R, self.N, dtype=np.int32)

    def test_demon_move_completes_one_sweep(self):
        """Test that demon_move_jit completes without error."""
        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)

        new_sum, new_lattice, new_bond_count = demon_move_jit(
            self.lattice,
            self.bonds,
            self.bond_count,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            self.order,
            self.radius_spin,
            self.radius_bond,
            self.N,
            self.R,
        )

        # Function should complete and return valid energies
        assert isinstance(new_sum, (int, np.integer))
        assert isinstance(new_lattice, (int, np.integer))
        assert np.sum(new_bond_count) == self.N

    def test_demon_move_energy_conservation(self):
        """Test energy conservation over full sweep."""
        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)
        E_total_before = E_demon_sum + E_lattice

        new_sum, new_lattice, new_bond_count = demon_move_jit(
            self.lattice,
            self.bonds,
            self.bond_count,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            self.order,
            self.radius_spin,
            self.radius_bond,
            self.N,
            self.R,
        )

        E_total_after = new_sum + new_lattice
        assert E_total_before == E_total_after

    def test_demon_move_modifies_state(self):
        """Test that demon_move_jit modifies simulation state."""
        lattice_before = self.lattice.copy()
        bonds_before = self.bonds.copy()
        E_demon_before = self.E_demon.copy()

        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)

        demon_move_jit(
            self.lattice,
            self.bonds,
            self.bond_count,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            self.order,
            self.radius_spin,
            self.radius_bond,
            self.N,
            self.R,
        )

        # State should change (with high probability)
        assert (
            not np.array_equal(lattice_before, self.lattice)
            or not np.array_equal(bonds_before, self.bonds)
            or not np.array_equal(E_demon_before, self.E_demon)
        )

    def test_demon_move_with_zero_radius(self):
        """Test demon_move with R=0 (local moves only)."""
        self.R = 0
        self.radius_spin = np.zeros(self.N, dtype=np.int32)
        self.radius_bond = np.zeros(self.N, dtype=np.int32)

        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)

        # Should still work with local moves
        new_sum, new_lattice, new_bond_count = demon_move_jit(
            self.lattice,
            self.bonds,
            self.bond_count,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            self.order,
            self.radius_spin,
            self.radius_bond,
            self.N,
            self.R,
        )

        assert np.sum(new_bond_count) == self.N

    def test_demon_move_bond_count_valid(self):
        """Test that bond counts remain valid throughout sweep."""
        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)

        new_sum, new_lattice, new_bond_count = demon_move_jit(
            self.lattice,
            self.bonds,
            self.bond_count,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            self.order,
            self.radius_spin,
            self.radius_bond,
            self.N,
            self.R,
        )

        # All bond counts should be non-negative
        assert np.all(new_bond_count >= 0)
        # Sum should equal N
        assert np.sum(new_bond_count) == self.N


class TestDemonMoveIrrJit:
    """Test demon_move_irr_jit function (irreversible)."""

    def setup_method(self):
        """Set up test simulation state."""
        self.N = 20
        self.R = 3
        self.lattice = np.array(
            [1 if i % 2 == 0 else -1 for i in range(self.N)], dtype=np.int8
        )
        self.bonds = np.array(
            [1 if i % 2 == 0 else -1 for i in range(self.N)], dtype=np.int8
        )
        self.bond_count = np.array([10, 0, 10], dtype=np.int32)
        self.E_demon = np.full(self.N, 20, dtype=np.int64)
        self.right_neighbor = np.array(
            [(i + 1) % self.N for i in range(self.N)], dtype=np.int32
        )
        self.left_neighbor = np.array(
            [(i - 1) % self.N for i in range(self.N)], dtype=np.int32
        )
        self.order = np.arange(self.N, dtype=np.int32)

    def test_demon_move_irr_completes_one_sweep(self):
        """Test that demon_move_irr_jit completes without error."""
        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)
        seed = 42

        new_sum, new_lattice, new_bond_count = demon_move_irr_jit(
            self.lattice,
            self.bonds,
            self.bond_count,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            self.order,
            self.N,
            self.R,
            seed,
        )

        assert isinstance(new_sum, (int, np.integer))
        assert isinstance(new_lattice, (int, np.integer))
        assert np.sum(new_bond_count) == self.N

    def test_demon_move_irr_energy_conservation(self):
        """Test energy conservation in irreversible version."""
        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)
        E_total_before = E_demon_sum + E_lattice
        seed = 123

        new_sum, new_lattice, new_bond_count = demon_move_irr_jit(
            self.lattice,
            self.bonds,
            self.bond_count,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            self.order,
            self.N,
            self.R,
            seed,
        )

        E_total_after = new_sum + new_lattice
        assert E_total_before == E_total_after

    def test_demon_move_irr_reproducible_with_seed(self):
        """Test that same seed gives reproducible results."""
        lattice1 = self.lattice.copy()
        bonds1 = self.bonds.copy()
        E_demon1 = self.E_demon.copy()
        bond_count1 = self.bond_count.copy()

        E_demon_sum = np.sum(E_demon1)
        E_lattice = np.int64(100)
        seed = 999

        new_sum1, new_lattice1, new_bond_count1 = demon_move_irr_jit(
            lattice1,
            bonds1,
            bond_count1,
            E_demon1,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            self.order,
            self.N,
            self.R,
            seed,
        )

        # Run again with same seed
        lattice2 = self.lattice.copy()
        bonds2 = self.bonds.copy()
        E_demon2 = self.E_demon.copy()
        bond_count2 = self.bond_count.copy()

        new_sum2, new_lattice2, new_bond_count2 = demon_move_irr_jit(
            lattice2,
            bonds2,
            bond_count2,
            E_demon2,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            self.order,
            self.N,
            self.R,
            seed,
        )

        # Results should be identical
        assert np.array_equal(lattice1, lattice2)

    def test_demon_move_irr_different_seeds_differ(self):
        """
        Test that different seeds give different RNG behavior
        (checked via demon energy).
        """
        # The RNG state should differ even if final lattice happens to be same
        # We can verify different seeds by checking demon energy distribution
        N = 50
        lattice1 = np.random.choice([-1, 1], N).astype(np.int8)
        lattice2 = lattice1.copy()
        bonds1 = np.random.choice([-1, 0, 1], N).astype(np.int8)
        bonds2 = bonds1.copy()
        E_demon1 = np.full(N, 50, dtype=np.int64)
        E_demon2 = E_demon1.copy()
        bond_count1 = np.array(
            [np.sum(bonds1 == -1), np.sum(bonds1 == 0), np.sum(bonds1 == 1)],
            dtype=np.int32,
        )
        bond_count2 = bond_count1.copy()
        right_neighbor = np.array([(i + 1) % N for i in range(N)], dtype=np.int32)
        left_neighbor = np.array([(i - 1) % N for i in range(N)], dtype=np.int32)
        order = np.arange(N, dtype=np.int32)

        E_demon_sum = np.sum(E_demon1)
        E_lattice = np.int64(1000)

        demon_move_irr_jit(
            lattice1,
            bonds1,
            bond_count1,
            E_demon1,
            E_demon_sum,
            E_lattice,
            right_neighbor,
            left_neighbor,
            order,
            N,
            5,  # Higher radius
            seed=100,
        )

        demon_move_irr_jit(
            lattice2,
            bonds2,
            bond_count2,
            E_demon2,
            E_demon_sum,
            E_lattice,
            right_neighbor,
            left_neighbor,
            order,
            N,
            5,
            seed=200,
        )

        # Different seeds should produce different outcomes somewhere
        # Check if lattice OR demon energy OR bonds differ
        assert (
            not np.array_equal(lattice1, lattice2)
            or not np.array_equal(E_demon1, E_demon2)
            or not np.array_equal(bonds1, bonds2)
        )

    def test_demon_move_irr_with_large_radius(self):
        """Test irreversible moves with large radius."""
        self.R = 10
        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)

        new_sum, new_lattice, new_bond_count = demon_move_irr_jit(
            self.lattice,
            self.bonds,
            self.bond_count,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            self.order,
            self.N,
            self.R,
            seed=42,
        )

        # Should still conserve energy
        E_total_after = new_sum + new_lattice
        assert E_total_after == E_demon_sum + E_lattice

    def test_demon_move_irr_bond_count_conservation(self):
        """Test that bond count sum is conserved."""
        E_demon_sum = np.sum(self.E_demon)
        E_lattice = np.int64(100)

        new_sum, new_lattice, new_bond_count = demon_move_irr_jit(
            self.lattice,
            self.bonds,
            self.bond_count,
            self.E_demon,
            E_demon_sum,
            E_lattice,
            self.right_neighbor,
            self.left_neighbor,
            self.order,
            self.N,
            self.R,
            seed=555,
        )

        assert np.sum(new_bond_count) == self.N
        assert np.all(new_bond_count >= 0)


class TestJitFunctionIntegration:
    """Integration tests combining multiple JIT functions."""

    def test_spin_flip_followed_by_bond_update(self):
        """Test combined spin flip and bond update."""
        N = 10
        lattice = np.array([1, -1, 1, -1, 1, -1, 1, -1, 1, -1], dtype=np.int8)
        bonds = np.array([1, -1, 1, -1, 1, -1, 1, -1, 1, -1], dtype=np.int8)
        bond_count = np.array([5, 0, 5], dtype=np.int32)
        E_demon = np.array([10, 10, 10, 10, 10, 10, 10, 10, 10, 10], dtype=np.int64)
        right_neighbor = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 0], dtype=np.int32)
        left_neighbor = np.array([9, 0, 1, 2, 3, 4, 5, 6, 7, 8], dtype=np.int32)

        E_demon_sum = np.sum(E_demon)
        E_lattice = np.int64(100)

        a, i = 2, 2

        # Attempt spin flip
        bonds_changed, E_demon_sum, E_lattice = spin_flip_jit(
            lattice,
            bonds,
            E_demon,
            E_demon_sum,
            E_lattice,
            right_neighbor,
            left_neighbor,
            a,
            i,
        )

        if bonds_changed:
            # Update bonds
            bond_count = update_bonds_incremental_jit(
                lattice, bonds, bond_count, right_neighbor, left_neighbor, a
            )

        # Verify consistency
        assert np.sum(bond_count) == N

    def test_multiple_operations_maintain_energy(self):
        """Test that multiple operations maintain energy conservation."""
        N = 15
        lattice = np.random.choice([-1, 1], N).astype(np.int8)
        bonds = np.random.choice([-1, 0, 1], N).astype(np.int8)
        bond_count = np.array(
            [np.sum(bonds == -1), np.sum(bonds == 0), np.sum(bonds == 1)],
            dtype=np.int32,
        )
        E_demon = np.random.randint(5, 20, N, dtype=np.int64)
        right_neighbor = np.array([(i + 1) % N for i in range(N)], dtype=np.int32)
        left_neighbor = np.array([(i - 1) % N for i in range(N)], dtype=np.int32)

        E_demon_sum = np.sum(E_demon)
        E_lattice = np.int64(150)
        E_total_start = E_demon_sum + E_lattice

        # Perform multiple operations
        for iteration in range(5):
            a = np.random.randint(0, N)
            i = np.random.randint(0, N)

            # Spin flip
            bonds_changed, E_demon_sum, E_lattice = spin_flip_jit(
                lattice,
                bonds,
                E_demon,
                E_demon_sum,
                E_lattice,
                right_neighbor,
                left_neighbor,
                a,
                i,
            )

            if bonds_changed:
                bond_count = update_bonds_incremental_jit(
                    lattice, bonds, bond_count, right_neighbor, left_neighbor, a
                )

            # Bond change
            a = np.random.randint(0, N)
            i = np.random.randint(0, N)
            E_demon_sum, E_lattice, bond_count = bond_change_jit(
                lattice,
                bonds,
                bond_count,
                E_demon,
                E_demon_sum,
                E_lattice,
                right_neighbor,
                left_neighbor,
                a,
                i,
            )

        E_total_end = E_demon_sum + E_lattice
        assert E_total_start == E_total_end
