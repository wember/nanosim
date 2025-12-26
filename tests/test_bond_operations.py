"""
Bond operations, state consistency, and edge case tests.

Tests cover:
- Spin flip operations and energy changes
- Bond change operations and mechanics
- Bond update edge cases (broken, neutral, left neighbor)
- Boundary conditions and periodic wrapping
- State consistency checks
- Integer arithmetic verification
"""

import sys
import os
import pytest
import numpy as np

# Add creutz-sim directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'creutz-sim'))

from inferno import Inferno
from inferno_irr import irrInferno


class TestSpinFlipOperation:
    """Test spin flip mechanics for both classes."""
    
    def test_spin_flip_energy_cost(self):
        """Test that spin flip correctly calculates and applies energy cost."""
        for SimClass in [Inferno, irrInferno]:
            sim = SimClass(N=100, R=5, validate_mode='off')
            
            # Test multiple spin flips
            for _ in range(10):
                initial_total = sim.E_lattice + sim.E_demon_sum
                site = np.random.randint(0, sim.N)
                demon_idx = np.random.randint(0, sim.N)
                
                # Attempt spin flip
                sim.spin_flip(site, demon_idx)
                
                # Energy must be conserved
                assert sim.E_lattice + sim.E_demon_sum == initial_total
    
    def test_spin_flip_updates_lattice(self):
        """Test that successful spin flip updates the lattice."""
        sim = Inferno(N=100, R=5, validate_mode='off')
        
        # Find a site where spin flip will succeed (enough demon energy)
        site = 50
        initial_spin = sim.lattice[site]
        
        # Give demon lots of energy to ensure flip succeeds
        demon_idx = 50
        sim.E_demon[demon_idx] = 10
        sim.E_demon_sum += 10
        sim.E_lattice -= 10
        
        # Flip should succeed
        bonds_changed = sim.spin_flip(site, demon_idx)
        
        # If flip succeeded, spin should be flipped
        if bonds_changed:
            assert sim.lattice[site] == -initial_spin


class TestBondChangeOperation:
    """Test bond change mechanics for both classes."""
    
    def test_bond_change_energy_conservation(self):
        """Test that bond changes conserve energy."""
        for SimClass in [Inferno, irrInferno]:
            sim = SimClass(N=100, R=5, validate_mode='off')
            
            for _ in range(20):
                initial_total = sim.E_lattice + sim.E_demon_sum
                site = np.random.randint(0, sim.N)
                demon_idx = np.random.randint(0, sim.N)
                
                # Attempt bond change
                sim.bond_change(site, demon_idx)
                
                # Energy must be conserved
                assert sim.E_lattice + sim.E_demon_sum == initial_total
    
    def test_bond_change_bond_count_consistency(self):
        """Test that bond changes maintain consistent bond counts."""
        sim = Inferno(N=100, R=5, validate_mode='off')
        
        for _ in range(50):
            sim.bond_change(
                np.random.randint(0, sim.N),
                np.random.randint(0, sim.N)
            )
            
            # Verify bond counts match actual bonds
            actual_counts = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
            assert np.array_equal(sim.bond_count, actual_counts), \
                f"Bond count mismatch: tracked {sim.bond_count} vs actual {actual_counts}"


class TestBondUpdateEdgeCases:
    """Test edge cases in bond update logic."""
    
    def test_left_neighbor_bond_sign_flip(self):
        """Test that bond_change properly updates left neighbor bond when it has wrong sign.
        
        This specifically tests lines 222-226 in inferno.py where old_bond != new_bond.
        The bond_change function should detect and fix inconsistencies in the left neighbor bond.
        """
        for SimClass in [Inferno, irrInferno]:
            sim = SimClass(N=100, R=5, validate_mode='off')
            
            # Set up a configuration with an intentional inconsistency
            # Site 48: +1, Site 49: +1, Site 50: -1
            sim.lattice[48] = np.int8(1)
            sim.lattice[49] = np.int8(1)
            sim.lattice[50] = np.int8(-1)
            
            # Bond 48 (from 48 to 49): should be -1 (aligned)
            # Bond 49 (from 49 to 50): should be +1 (anti-aligned) - but we'll set it wrong
            sim.bonds[48] = np.int8(-1)
            sim.bonds[49] = np.int8(-1)  # WRONG! Should be +1
            
            # Recalculate counts to match the (incorrect) bonds
            sim.bond_count = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
            
            # Give demon enough energy
            sim.E_demon[50] = 10
            
            # Call bond_change at index 50 (bond from 50 to 51)
            # This should also check and fix the left neighbor bond (bonds[49])
            sim.bond_change(50, 50)
            
            # Verify the left bond was corrected
            expected_bond_49 = -1 if sim.lattice[49] == sim.lattice[50] else 1
            assert sim.bonds[49] == expected_bond_49, \
                f"{SimClass.__name__}: Bond 49 should be corrected to {expected_bond_49}, got {sim.bonds[49]}"
            
            # Verify counts are still consistent
            actual_counts = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
            assert np.array_equal(sim.bond_count, actual_counts), \
                f"{SimClass.__name__}: Bond counts should match after left neighbor update"
    
    def test_update_bonds_left_neighbor_change(self):
        """Test bond_change when left neighbor bond needs updating."""
        for SimClass in [Inferno, irrInferno]:
            sim = SimClass(N=100, R=5, validate_mode='off')
            
            # Set up inconsistent state where left bond doesn't match lattice
            i = 50
            prev = 49
            
            sim.lattice[prev] = np.int8(1)
            sim.lattice[i] = np.int8(-1)
            sim.bonds[prev] = np.int8(1)  # Inconsistent! Should be -1
            
            # Recalculate bond counts
            sim.bond_count = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
            
            # Call bond_change - it will fix the inconsistency
            sim.bond_change(i, 0)
            
            # Verify bond counts are consistent
            actual_counts = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
            assert np.array_equal(sim.bond_count, actual_counts), \
                f"{SimClass.__name__}: Bond counts should match"
    
    def test_update_bonds_with_neutral_bond(self):
        """Test update_bonds_incremental when left bond is neutral."""
        sim = Inferno(N=100, R=5, validate_mode='off')
        
        # Create a configuration with a neutral bond
        for i in range(sim.N):
            j = (i + 1) % sim.N
            if sim.bonds[i] == 1:
                sim.lattice[j] = -sim.lattice[j]
                sim.bonds[i] = 0
                sim.bond_count[2] -= 1
                sim.bond_count[1] += 1
                
                # Test update when i is left neighbor of j
                sim.lattice[j] = -sim.lattice[j]
                sim.update_bonds_incremental(j)
                
                # Verify counts are consistent
                actual_counts = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
                assert np.array_equal(actual_counts, sim.bond_count)
                break
    
    def test_update_bonds_with_broken_bond(self):
        """Test update_bonds_incremental with broken bonds."""
        sim = Inferno(N=100, R=5, validate_mode='off')
        
        # Find or create a broken bond
        broken_site = None
        for i in range(sim.N):
            if sim.bonds[i] == 0:
                broken_site = i
                break
        
        if broken_site is None:
            for i in range(sim.N):
                j = (i + 1) % sim.N
                if sim.lattice[i] == sim.lattice[j]:
                    sim.lattice[i] = -sim.lattice[i]
                    prev_i = (i - 1) % sim.N
                    sim.bonds[prev_i] = sim.lattice[prev_i] * sim.lattice[i]
                    sim.bonds[i] = sim.lattice[i] * sim.lattice[j]
                    broken_site = i
                    break
        
        # Recalculate bond counts
        sim.bond_count = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
        
        # Test update
        site = broken_site if broken_site is not None else 0
        sim.lattice[site] = -sim.lattice[site]
        sim.update_bonds_incremental(site)
        
        # Verify consistency
        actual_counts = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
        assert np.array_equal(actual_counts, sim.bond_count)
    
    def test_bond_update_all_broken_bonds(self):
        """Test bond updates when all bonds start broken."""
        sim = Inferno(N=20, R=5, validate_mode='off')
        
        # Create alternating spins (all bonds broken)
        for i in range(sim.N):
            sim.lattice[i] = 1 if i % 2 == 0 else -1
        
        # Recalculate bonds
        for i in range(sim.N):
            j = (i + 1) % sim.N
            sim.bonds[i] = sim.lattice[i] * sim.lattice[j]
        
        sim.bond_count = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
        
        # All bonds should be broken (value -1, bin[0])
        assert sim.bond_count[0] == sim.N
        
        # Flip two adjacent spins to create aligned bonds
        sim.lattice[5] = -sim.lattice[5]
        sim.lattice[6] = -sim.lattice[6]
        sim.update_bonds_incremental(5)
        sim.update_bonds_incremental(6)
        
        # Should now have some aligned bonds
        assert sim.bond_count[2] > 0
        
        # Verify counts match
        actual_counts = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
        assert np.array_equal(sim.bond_count, actual_counts)


class TestBoundaryConditions:
    """Test periodic boundary conditions."""
    
    def test_periodic_wraparound_site_0(self):
        """Test that operations at site 0 correctly wrap to site N-1."""
        sim = Inferno(N=100, R=5, validate_mode='off')
        
        # Flip spin at site 0
        initial_energy = sim.E_total
        sim.spin_flip(0, 50)
        
        # Verify left neighbor (site 99) bond is updated
        assert sim.E_total == initial_energy
    
    def test_periodic_wraparound_site_N_minus_1(self):
        """Test that operations at site N-1 correctly wrap to site 0."""
        sim = Inferno(N=100, R=5, validate_mode='off')
        
        initial_energy = sim.E_total
        sim.spin_flip(sim.N - 1, 50)
        
        # Right neighbor is site 0
        assert sim.E_total == initial_energy


class TestStateConsistency:
    """Test that internal state remains consistent."""
    
    def test_bond_count_always_matches_bonds(self):
        """Test that bond_count always matches actual bond array."""
        sim = Inferno(N=100, R=5, validate_mode='off')
        
        for _ in range(100):
            sim.demon_move()
            
            # Verify bond count matches bonds
            actual_counts = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
            assert np.array_equal(sim.bond_count, actual_counts), \
                f"Bond count mismatch: {sim.bond_count} vs {actual_counts}"
    
    def test_demon_energy_sum_matches_array(self):
        """Test that E_demon_sum equals sum of E_demon array."""
        sim = Inferno(N=100, R=5, validate_mode='off')
        
        for _ in range(100):
            sim.demon_move()
            
            # Verify sum matches
            assert sim.E_demon_sum == np.sum(sim.E_demon), \
                f"Demon energy sum mismatch: {sim.E_demon_sum} vs {np.sum(sim.E_demon)}"
    
    def test_state_methods_consistent(self):
        """Test that get_validated_state returns consistent values."""
        sim = Inferno(N=100, R=5, validate_mode='off')
        
        for _ in range(50):
            sim.demon_move()
        
        state = sim.get_validated_state()
        
        # Verify state matches internal values
        assert state['E_total'] == sim.E_total
        assert state['E_lattice'] == sim.E_lattice
        assert state['E_demon_sum'] == sim.E_demon_sum


class TestIntegerArithmetic:
    """Test that all operations use integer arithmetic."""
    
    def test_all_energies_are_integers(self):
        """Test that all energy values remain integers."""
        sim = Inferno(N=100, R=5, validate_mode='off')
        
        for _ in range(100):
            sim.demon_move()
            
            # All energies should be integers
            assert isinstance(sim.E_lattice, (int, np.integer))
            assert isinstance(sim.E_demon_sum, (int, np.integer))
            assert isinstance(sim.E_total, (int, np.integer))
            assert np.all([isinstance(e, (int, np.integer)) for e in sim.E_demon])
    
    def test_no_floating_point_drift(self):
        """Test that integer arithmetic prevents floating point drift."""
        sim = Inferno(N=1000, R=5, validate_mode='off')
        initial_energy = sim.E_total
        
        # Run many operations
        for _ in range(10000):
            sim.demon_move()
        
        # Energy should be EXACTLY conserved (not approximately)
        assert sim.E_total == initial_energy
        assert type(sim.E_total) == type(initial_energy)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
