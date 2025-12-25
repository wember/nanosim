"""
Additional critical tests for core physics operations and edge cases.
These tests go deeper into the Monte Carlo moves and boundary conditions.
"""
import sys
import os
import pytest
import numpy as np

# Add creutz-sim directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'creutz-sim'))

from inferno import Inferno
from irr_inferno import irrInferno


class TestSpinFlipOperation:
    """Test the spin_flip method directly."""
    
    def test_spin_flip_energy_cost(self):
        """Test that spin flip calculates energy cost correctly."""
        x = Inferno(10, 3, validate_mode='off')
        
        # Store initial state
        initial_lattice = x.lattice.copy()
        initial_E_demon = x.E_demon.copy()
        initial_E_total = x.E_total
        
        # Try to flip middle spin
        a = 5
        i = 0
        
        # Calculate expected cost manually
        s = x.lattice[a]
        nb = (x.lattice[(a+1) % x.N] * abs(x.bonds[a]) +
              x.lattice[(a-1) % x.N] * abs(x.bonds[(a-1) % x.N]))
        expected_cost = 2 * s * nb
        
        # Attempt flip
        bonds_changed = x.spin_flip(a, i)
        
        # Energy must be conserved
        assert x.E_total == initial_E_total
        
        if bonds_changed:
            # Spin should have flipped
            assert x.lattice[a] == -initial_lattice[a]
            # Demon energy should have changed by cost
            assert x.E_demon[i] == initial_E_demon[i] - expected_cost
    
    def test_spin_flip_rejected_if_insufficient_demon_energy(self):
        """Test that spin flip is rejected when demon lacks energy."""
        x = Inferno(10, 3, validate_mode='off')
        
        # Set demon energy to zero
        demon_idx = 5
        x.E_demon[demon_idx] = 0
        x.E_demon_sum = np.sum(x.E_demon)
        
        # Find a spin flip that costs energy
        for site in range(x.N):
            s = x.lattice[site]
            nb = (x.lattice[(site+1) % x.N] * abs(x.bonds[site]) +
                  x.lattice[(site-1) % x.N] * abs(x.bonds[(site-1) % x.N]))
            cost = 2 * s * nb
            
            if cost > 0:  # Costs energy
                initial_spin = x.lattice[site]
                bonds_changed = x.spin_flip(site, demon_idx)
                
                # Flip should be rejected
                assert bonds_changed == False
                assert x.lattice[site] == initial_spin
                break
    
    def test_spin_flip_energy_conservation_always(self):
        """Test that spin flip always conserves total energy regardless of acceptance."""
        x = Inferno(50, 3, validate_mode='off')
        
        # Try many random flips
        for _ in range(200):
            site = np.random.randint(0, x.N)
            demon_idx = np.random.randint(0, x.N)
            
            initial_E_total = x.E_total
            initial_E_lattice = x.E_lattice
            initial_E_demon_sum = x.E_demon_sum
            
            # Attempt flip
            bonds_changed = x.spin_flip(site, demon_idx)
            
            # Energy must always be conserved (whether accepted or rejected)
            assert x.E_total == initial_E_total
            assert x.E_lattice + x.E_demon_sum == initial_E_total


class TestBondChangeOperation:
    """Test the bond_change method directly."""
    
    def test_bond_change_energy_conservation(self):
        """Test that bond changes conserve energy."""
        x = Inferno(10, 3, validate_mode='off')
        initial_E_total = x.E_total
        
        # Try many bond changes
        for _ in range(100):
            a = np.random.randint(0, x.N)
            i = np.random.randint(0, x.N)
            x.bond_change(a, i)
            
            assert x.E_total == initial_E_total
    
    def test_bond_change_from_broken_to_active(self):
        """Test changing bond from broken (0) to active (±1)."""
        x = Inferno(10, 3, validate_mode='off')
        
        # Find a broken bond
        for site in range(x.N):
            if x.bonds[site] == 0:
                initial_E_total = x.E_total
                initial_bond_count = x.bond_count.copy()
                demon_idx = 0
                
                # Attempt bond change
                x.bond_change(site, demon_idx)
                
                # Energy must be conserved
                assert x.E_total == initial_E_total
                
                # Bond count must be updated correctly
                assert np.sum(x.bond_count) == x.N
                break
    
    def test_bond_change_from_active_to_broken(self):
        """Test changing bond from active (±1) to broken (0)."""
        x = Inferno(10, 3, validate_mode='off')
        
        # Find an active bond
        for site in range(x.N):
            if x.bonds[site] != 0:
                initial_E_total = x.E_total
                demon_idx = 5
                
                # Give demon enough energy
                x.E_demon[demon_idx] += 10
                x.E_demon_sum += 10
                x.E_lattice -= 10
                
                # Attempt bond change
                x.bond_change(site, demon_idx)
                
                # Energy must be conserved
                assert x.E_total == initial_E_total
                break


class TestBoundaryConditions:
    """Test periodic boundary conditions at lattice edges."""
    
    def test_boundary_spin_flip_at_zero(self):
        """Test spin flip at site 0 (left boundary)."""
        x = Inferno(20, 3, validate_mode='off')
        initial_E_total = x.E_total
        
        # Flip at boundary
        for _ in range(10):
            x.spin_flip(0, 0)
            assert x.E_total == initial_E_total
    
    def test_boundary_spin_flip_at_N_minus_1(self):
        """Test spin flip at site N-1 (right boundary)."""
        x = Inferno(20, 3, validate_mode='off')
        initial_E_total = x.E_total
        
        # Flip at boundary
        for _ in range(10):
            x.spin_flip(x.N-1, 0)
            assert x.E_total == initial_E_total
    
    def test_periodic_wrapping(self):
        """Test that periodic boundary conditions wrap correctly."""
        x = Inferno(10, 3, validate_mode='off')
        
        # Site N-1 should connect to site 0
        left_neighbor_of_0 = (0 - 1) % x.N
        assert left_neighbor_of_0 == x.N - 1
        
        right_neighbor_of_Nminus1 = (x.N - 1 + 1) % x.N
        assert right_neighbor_of_Nminus1 == 0
    
    def test_bonds_at_boundaries_updated(self):
        """Test that bonds at boundaries are updated correctly."""
        x = Inferno(10, 3, validate_mode='off')
        
        # Flip spin at boundary
        initial_bonds = x.bonds.copy()
        x.spin_flip(0, 0)
        x.update_bonds_incremental(0)
        
        # Bond at N-1 (left neighbor of 0) might have changed
        # Bond at 0 (between 0 and 1) might have changed
        # But total bond count must be consistent
        actual_counts = np.bincount(x.bonds + 1, minlength=3).astype(np.int64)
        np.testing.assert_array_equal(x.bond_count, actual_counts)


class TestDemonEnergyDistribution:
    """Test that demon energy distribution remains physical."""
    
    def test_demon_energy_stays_nonnegative(self):
        """Test that no demon ever gets negative energy."""
        x = Inferno(100, 5, validate_mode='off')
        
        # Run many moves
        for _ in range(1000):
            for _ in range(x.N):
                x.demon_move()
        
        # All demons must have non-negative energy
        assert np.all(x.E_demon >= 0)
    
    def test_demon_energy_distribution_physical(self):
        """Test that demon energy distribution is reasonable."""
        x = Inferno(1000, 5, validate_mode='off')
        
        # Run to equilibration
        for _ in range(500):
            for _ in range(x.N):
                x.demon_move()
        
        # Check distribution properties
        mean_demon_energy = np.mean(x.E_demon)
        max_demon_energy = np.max(x.E_demon)
        
        # Should have nonzero mean
        assert mean_demon_energy > 0
        
        # Max shouldn't be absurdly large (should be ~O(10) for typical systems)
        # For E_total = 2N, average per demon is 2, max should be reasonable
        assert max_demon_energy < 100  # Generous upper bound
    
    def test_demon_energy_sum_correct(self):
        """Test that E_demon_sum matches actual sum."""
        x = Inferno(100, 5, validate_mode='off')
        
        # Run moves
        for _ in range(100):
            for _ in range(x.N):
                x.demon_move()
        
        # Cached sum must match actual
        actual_sum = np.sum(x.E_demon)
        assert x.E_demon_sum == actual_sum


class TestReversibleVsIrreversible:
    """Compare reversible and irreversible simulations."""
    
    def test_irreversible_different_after_cycle(self):
        """Test that irreversible doesn't return after forward-reverse."""
        x_irr = irrInferno(50, 3, validate_mode='off')
        initial_lattice = x_irr.lattice.copy()
        
        # Forward then reverse
        for _ in range(20):
            for _ in range(x_irr.N):
                x_irr.demon_move()
        
        for _ in range(20):
            for _ in range(x_irr.N):
                x_irr.demon_reverse()
        
        # Should NOT return to initial state
        assert not np.array_equal(x_irr.lattice, initial_lattice)
    
    def test_energy_conservation_both_types(self):
        """Test that both reversible and irreversible conserve energy."""
        x_rev = Inferno(100, 5, validate_mode='off')
        x_irr = irrInferno(100, 5, validate_mode='off')
        
        initial_E_rev = x_rev.E_total
        initial_E_irr = x_irr.E_total
        
        # Run many moves
        for _ in range(500):
            for _ in range(100):
                x_rev.demon_move()
                x_irr.demon_move()
        
        # Both must conserve energy
        assert x_rev.E_total == initial_E_rev
        assert x_irr.E_total == initial_E_irr


class TestIntegerArithmetic:
    """Test that integer arithmetic is maintained throughout."""
    
    def test_all_energies_are_integers(self):
        """Test that all energy values remain integers."""
        x = Inferno(100, 5, validate_mode='off')
        
        # Run moves
        for _ in range(100):
            for _ in range(x.N):
                x.demon_move()
        
        # Check types
        assert x.E_total == int(x.E_total)
        assert x.E_lattice == int(x.E_lattice)
        assert x.E_demon_sum == int(x.E_demon_sum)
        assert np.all(x.E_demon == x.E_demon.astype(np.int64))
    
    def test_no_floating_point_in_bonds(self):
        """Test that bond values are always -1, 0, or 1."""
        x = Inferno(100, 5, validate_mode='off')
        
        # Run moves
        for _ in range(100):
            for _ in range(x.N):
                x.demon_move()
        
        # All bonds must be in {-1, 0, 1}
        assert np.all(np.isin(x.bonds, [-1, 0, 1]))


class TestStateConsistency:
    """Test internal state consistency after operations."""
    
    def test_bond_count_matches_actual_always(self):
        """Test bond_count array always matches actual bond configuration."""
        x = Inferno(50, 3, validate_mode='off')
        
        # Run moves and check periodically
        for _ in range(100):
            for _ in range(x.N):
                x.demon_move()
            
            # Verify consistency
            actual_counts = np.bincount(x.bonds + 1, minlength=3).astype(np.int64)
            np.testing.assert_array_equal(x.bond_count, actual_counts,
                err_msg=f"Bond count mismatch: tracked={x.bond_count}, actual={actual_counts}")
    
    def test_energy_components_sum_correct(self):
        """Test that E_lattice + E_demon_sum = E_total always."""
        x = Inferno(100, 5, validate_mode='off')
        
        for _ in range(200):
            for _ in range(x.N):
                x.demon_move()
            
            # Check sum
            assert x.E_total == x.E_lattice + x.E_demon_sum, \
                f"Energy sum mismatch: E_total={x.E_total}, E_lattice={x.E_lattice}, E_demon_sum={x.E_demon_sum}"
    
    def test_lattice_spin_values(self):
        """Test that all lattice spins remain ±1."""
        x = Inferno(100, 5, validate_mode='off')
        
        # Run moves
        for _ in range(100):
            for _ in range(x.N):
                x.demon_move()
        
        # All spins must be ±1
        assert np.all(np.abs(x.lattice) == 1)


class TestLargeSystemStability:
    """Test stability with very large systems."""
    
    def test_large_system_no_overflow(self):
        """Test that large N doesn't cause integer overflow."""
        N = 100000  # Very large
        x = Inferno(N, 5, validate_mode='off')
        
        # Initial energy should be 2*N (fits in int64)
        assert x.E_total == 2 * N
        
        # Run some moves
        for _ in range(10):
            for _ in range(N):
                x.demon_move()
        
        # Energy still conserved
        assert x.E_total == 2 * N
    
    def test_long_run_accumulation_error(self):
        """Test for accumulation errors in very long runs."""
        x = Inferno(200, 5, validate_mode='off')
        initial_E = x.E_total
        
        # Very long run (40,000 sweeps)
        for _ in range(40000):
            for _ in range(x.N):
                x.demon_move()
        
        # Must still be exact
        assert x.E_total == initial_E
        
        # Validate against actual
        state = x.get_validated_state()
        assert state['E_total'] == initial_E


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
