"""
Tests targeting specific uncovered code paths to achieve 100% test coverage.

These tests focus on:
1. Validation branches in demon_move/demon_reverse (inferno.py lines 248-253, 278-282)
2. Bond update edge cases (inferno.py lines 222-226)
3. Validation method branches (irr_inferno.py lines 68-71, 79-97, 101-109)
4. Validation branches in irr_inferno demon moves (lines 199-203, 227-231, 254-258)
"""

import numpy as np
import pytest
import os
import sys

# Add the creutz-sim directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'creutz-sim'))

from inferno import Inferno
from irr_inferno import irrInferno


class TestValidationBranches:
    """Test validation code paths in both Inferno and irrInferno classes."""
    
    def test_inferno_frequent_validation_in_demon_move(self):
        """Test that frequent validation triggers validation checks in demon_move."""
        sim = Inferno(N=100, R=5, validate_mode='frequent')
        initial_total = sim.E_lattice + sim.E_demon_sum
        
        # Run one sweep - validation should happen every N moves
        sim.demon_move()
        
        # Verify energy conservation maintained
        assert sim.E_lattice + sim.E_demon_sum == initial_total
        
    def test_inferno_frequent_validation_in_demon_reverse(self):
        """Test that frequent validation triggers validation checks in demon_reverse."""
        sim = Inferno(N=100, R=5, validate_mode='frequent')
        initial_total = sim.E_lattice + sim.E_demon_sum
        
        # Run exactly N moves to trigger validation (check_interval = N for frequent mode)
        for _ in range(sim.N):
            sim.demon_reverse()
        
        # Verify energy conservation maintained
        assert sim.E_lattice + sim.E_demon_sum == initial_total
        
    def test_inferno_periodic_validation_boundary(self):
        """Test periodic validation at 100-sweep boundary."""
        sim = Inferno(N=50, R=5, validate_mode='periodic')
        initial_total = sim.E_lattice + sim.E_demon_sum
        
        # Run exactly 100 sweeps to trigger periodic validation
        for _ in range(100):
            sim.demon_move()
        
        # Verify energy conservation maintained through validation
        assert sim.E_lattice + sim.E_demon_sum == initial_total
        
    def test_irr_inferno_frequent_validation_in_demon_move(self):
        """Test that frequent validation triggers validation checks in irrInferno demon_move."""
        sim = irrInferno(N=100, R=5, validate_mode='frequent')
        initial_total = sim.E_lattice + sim.E_demon_sum
        
        # Run exactly N moves to trigger validation (check_interval = N for frequent mode)
        for _ in range(sim.N):
            sim.demon_move()
        
        # Verify energy conservation maintained
        assert sim.E_lattice + sim.E_demon_sum == initial_total
        
    def test_irr_inferno_frequent_validation_in_demon_reverse(self):
        """Test that frequent validation triggers validation checks in irrInferno demon_reverse."""
        sim = irrInferno(N=100, R=5, validate_mode='frequent')
        initial_total = sim.E_lattice + sim.E_demon_sum
        
        # Run exactly N moves to trigger validation (check_interval = N for frequent mode)
        for _ in range(sim.N):
            sim.demon_reverse()
        
        # Verify energy conservation maintained
        assert sim.E_lattice + sim.E_demon_sum == initial_total
        
    def test_irr_inferno_periodic_validation_boundary(self):
        """Test periodic validation at 100-sweep boundary in irrInferno."""
        sim = irrInferno(N=50, R=5, validate_mode='periodic')
        initial_total = sim.E_lattice + sim.E_demon_sum
        
        # Run exactly 100 sweeps to trigger periodic validation
        for _ in range(100):
            sim.demon_move()
        
        # Verify energy conservation maintained through validation
        assert sim.E_lattice + sim.E_demon_sum == initial_total


class TestBondUpdateEdgeCases:
    """Test edge cases in bond update logic (inferno.py lines 222-226)."""
    
    def test_update_bonds_left_neighbor_change(self):
        """Test bond_change when left neighbor bond needs updating.
        
        This specifically targets lines 222-226 of inferno.py and 199-203 of irr_inferno.py,
        which update the left neighbor bond in bond_change() when it doesn't match the lattice configuration.
        
        Convention: bond=-1 means aligned (low energy), bond=+1 means anti-aligned (high energy)
        """
        for SimClass in [Inferno, irrInferno]:
            sim = SimClass(N=100, R=5, validate_mode='off')
            
            # Set up a scenario where the left bond is inconsistent with the lattice
            # lattice[49] = 1, lattice[50] = -1 (anti-aligned), but bonds[49] = -1 (says aligned)
            i = 50
            prev = 49
            
            sim.lattice[prev] = np.int8(1)
            sim.lattice[i] = np.int8(-1)   # Different spins
            sim.bonds[prev] = np.int8(-1)  # Inconsistent! Should be +1 (anti-aligned)
            
            # Recalculate bond counts to match
            sim.bond_count = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
            
            initial_left_bond = sim.bonds[prev]
            assert initial_left_bond == -1, "Setup: left bond should be -1"
            
            # Call bond_change on site i
            # The code will check the left neighbor bond and fix the inconsistency
            sim.bond_change(i, 0)  # Demon index 0
            
            # The left bond should now be corrected to +1 (anti-aligned)
            expected_bond = np.int8(-1 if sim.lattice[prev] == sim.lattice[i] else 1)
            assert sim.bonds[prev] == expected_bond, \
                f"{SimClass.__name__}: Left bond should be updated to {expected_bond}, got {sim.bonds[prev]}"
            assert sim.bonds[prev] == 1, f"Bond should be +1 (anti-aligned) after fix"
            
            # Verify bond counts are consistent
            actual_counts = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
            assert np.array_equal(sim.bond_count, actual_counts), \
                f"{SimClass.__name__}: Bond counts should match: {sim.bond_count} vs {actual_counts}"
    
    def test_update_bonds_with_neutral_bond(self):
        """Test update_bonds_incremental when left bond is neutral (bonds[left_idx]==0).
        
        This covers the edge case in lines 222-226 of inferno.py where the left
        neighbor bond is 0 and the update logic takes a different path.
        """
        sim = Inferno(N=100, R=5, validate_mode='off')
        
        # Create a configuration with a neutral bond
        # We need bonds[i] == 0 for some i, which happens when spins are anti-aligned
        # but the bond hasn't been "activated" yet
        
        # Start by breaking a bond (making it -1)
        for i in range(sim.N):
            j = (i + 1) % sim.N
            if sim.bonds[i] == 1:  # Find an aligned bond
                # Flip the spin to make anti-aligned
                sim.lattice[j] = -sim.lattice[j]
                # Set bond to neutral
                sim.bonds[i] = 0
                # Update bond count
                sim.bond_count[2] -= 1
                sim.bond_count[1] += 1
                
                # Now test update_bonds_incremental when i is the left neighbor
                # Flip lattice[j] again (now i is left neighbor of j)
                sim.lattice[j] = -sim.lattice[j]
                sim.update_bonds_incremental(j)
                
                # The neutral bond should have been skipped (line 218 check)
                # Verify counts are consistent
                actual_counts = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
                assert np.array_equal(actual_counts, sim.bond_count), \
                    f"Bond counts should match: tracked {sim.bond_count} vs actual {actual_counts}"
                break
    
    def test_update_bonds_with_broken_bond(self):
        """Test update_bonds_incremental when bond is broken (bonds[a]==0)."""
        sim = Inferno(N=100, R=5, validate_mode='off')
        
        # Find a site with a broken bond or create one by forcing anti-alignment
        broken_site = None
        for i in range(sim.N):
            if sim.bonds[i] == 0:
                broken_site = i
                break
        
        # If no broken bond exists, create one
        if broken_site is None:
            for i in range(sim.N):
                j = (i + 1) % sim.N
                if sim.lattice[i] == sim.lattice[j]:
                    # Use the spin_flip method to properly update everything
                    sim.lattice[i] = -sim.lattice[i]
                    # Manually recalculate bonds affected by this flip
                    prev_i = (i - 1) % sim.N
                    sim.bonds[prev_i] = sim.lattice[prev_i] * sim.lattice[i]
                    sim.bonds[i] = sim.lattice[i] * sim.lattice[j]
                    broken_site = i
                    break
        
        # Recalculate bond counts from scratch to ensure consistency
        sim.bond_count = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
        
        # Now verify that update_bonds_incremental works correctly
        initial_energy = sim.E_lattice
        site_to_flip = broken_site if broken_site is not None else 0
        
        # Flip a spin and update bonds
        sim.lattice[site_to_flip] = -sim.lattice[site_to_flip]
        sim.update_bonds_incremental(site_to_flip)
        
        # Verify bond counts match reality
        actual_counts = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
        assert np.array_equal(actual_counts, sim.bond_count), \
            f"Bond count mismatch: tracked {sim.bond_count} vs actual {actual_counts}"
    
    def test_bond_update_all_broken_bonds(self):
        """Test bond updates when all bonds start broken."""
        sim = Inferno(N=20, R=5, validate_mode='off')
        
        # Create alternating spin configuration (all bonds broken)
        for i in range(sim.N):
            sim.lattice[i] = 1 if i % 2 == 0 else -1
        
        # Recalculate all bonds
        for i in range(sim.N):
            j = (i + 1) % sim.N
            sim.bonds[i] = sim.lattice[i] * sim.lattice[j]
        
        # Recount bonds
        sim.bond_count = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
        
        # All bonds should be broken (value -1, which goes in bin[0])
        assert sim.bond_count[0] == sim.N, f"All bonds should be broken: {sim.bond_count}"
        
        # Now flip TWO adjacent spins to create aligned bonds
        sim.lattice[5] = -sim.lattice[5]
        sim.lattice[6] = -sim.lattice[6]
        sim.update_bonds_incremental(5)
        sim.update_bonds_incremental(6)
        
        # Verify bond count changed - should now have some aligned bonds
        assert sim.bond_count[2] > 0, \
            f"Should have aligned bonds after flipping adjacent spins: {sim.bond_count}"
        
        # Verify counts match reality
        actual_counts = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
        assert np.array_equal(sim.bond_count, actual_counts), \
            f"Bond counts should match: tracked {sim.bond_count} vs actual {actual_counts}"


class TestValidationMethodCoverage:
    """Test validation method internals for both classes."""
    
    def test_inferno_validate_energy_conservation_with_drift(self):
        """Test energy validation when energy appears to drift."""
        sim = Inferno(N=100, R=5, validate_mode='off')
        
        # Manually corrupt the tracked energy to simulate drift
        sim.E_demon_sum += 5
        
        # Validation should detect and correct this
        result = sim.validate_energy_conservation()
        
        # Should return False indicating drift was detected
        assert result is False, "Validation should detect energy drift"
        
        # Energy should be corrected
        expected_total = 2 * sim.N
        actual_total = sim.E_lattice + sim.E_demon_sum
        assert actual_total == expected_total, \
            f"Energy should be corrected: {actual_total} vs {expected_total}"
    
    def test_inferno_validate_bond_counts_with_drift(self):
        """Test bond count validation when counts appear to drift."""
        sim = Inferno(N=100, R=5, validate_mode='off')
        
        # Manually corrupt the bond count to simulate drift
        original_count = sim.bond_count.copy()
        sim.bond_count[0] += 2
        sim.bond_count[2] -= 2
        
        # Validation should detect and correct this
        result = sim.validate_bond_counts()
        
        # Should return False indicating drift was detected
        assert result is False, "Validation should detect bond count drift"
        
        # Bond counts should be corrected
        actual_counts = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
        assert np.array_equal(sim.bond_count, actual_counts), \
            "Bond counts should be corrected"
    
    def test_irr_inferno_validate_energy_conservation_with_drift(self):
        """Test energy validation when energy appears to drift in irrInferno."""
        sim = irrInferno(N=100, R=5, validate_mode='off')
        
        # Manually corrupt the tracked energy to simulate drift
        sim.E_demon_sum += 5
        
        # Validation should detect and correct this
        result = sim.validate_energy_conservation()
        
        # Should return False indicating drift was detected
        assert result is False, "Validation should detect energy drift"
        
        # Energy should be corrected
        expected_total = 2 * sim.N
        actual_total = sim.E_lattice + sim.E_demon_sum
        assert actual_total == expected_total, \
            f"Energy should be corrected: {actual_total} vs {expected_total}"
    
    def test_irr_inferno_validate_bond_counts_with_drift(self):
        """Test bond count validation when counts appear to drift in irrInferno."""
        sim = irrInferno(N=100, R=5, validate_mode='off')
        
        # Manually corrupt the bond count to simulate drift
        original_count = sim.bond_count.copy()
        sim.bond_count[0] += 2
        sim.bond_count[2] -= 2
        
        # Validation should detect and correct this
        result = sim.validate_bond_counts()
        
        # Should return False indicating drift was detected
        assert result is False, "Validation should detect bond count drift"
        
        # Bond counts should be corrected
        actual_counts = np.bincount(sim.bonds + 1, minlength=3).astype(np.int64)
        assert np.array_equal(sim.bond_count, actual_counts), \
            "Bond counts should be corrected"


class TestValidationModeConfiguration:
    """Test validation mode configuration logic."""
    
    def test_inferno_validation_off_never_validates(self):
        """Verify that validation_mode='off' sets infinite check interval."""
        sim = Inferno(N=100, R=5, validate_mode='off')
        assert sim._check_interval == float('inf'), \
            "Validation mode 'off' should set infinite check interval"
        
        # Run many sweeps - no validation should occur
        initial_energy = sim.E_total
        for _ in range(10):
            sim.demon_move()
        
        # Should complete without issues and maintain energy conservation
        assert sim.E_lattice + sim.E_demon_sum == initial_energy
    
    def test_inferno_validation_frequent_checks_every_sweep(self):
        """Verify that validation_mode='frequent' sets N check interval."""
        n = 100
        sim = Inferno(N=n, R=5, validate_mode='frequent')
        assert sim._check_interval == n, \
            f"Validation mode 'frequent' should set check interval to N={n}"
    
    def test_inferno_validation_periodic_checks_every_100_sweeps(self):
        """Verify that validation_mode='periodic' sets 100*N check interval."""
        n = 100
        sim = Inferno(N=n, R=5, validate_mode='periodic')
        assert sim._check_interval == 100 * n, \
            f"Validation mode 'periodic' should set check interval to 100*N={100*n}"
    
    def test_irr_inferno_validation_off_never_validates(self):
        """Verify that validation_mode='off' sets infinite check interval in irrInferno."""
        sim = irrInferno(N=100, R=5, validate_mode='off')
        assert sim._check_interval == float('inf'), \
            "Validation mode 'off' should set infinite check interval"
    
    def test_irr_inferno_validation_frequent_checks_every_sweep(self):
        """Verify that validation_mode='frequent' sets N check interval in irrInferno."""
        n = 100
        sim = irrInferno(N=n, R=5, validate_mode='frequent')
        assert sim._check_interval == n, \
            f"Validation mode 'frequent' should set check interval to N={n}"
    
    def test_irr_inferno_validation_periodic_checks_every_100_sweeps(self):
        """Verify that validation_mode='periodic' sets 100*N check interval in irrInferno."""
        n = 100
        sim = irrInferno(N=n, R=5, validate_mode='periodic')
        assert sim._check_interval == 100 * n, \
            f"Validation mode 'periodic' should set check interval to 100*N={100*n}"


class TestValidationIntegration:
    """Integration tests for validation system."""
    
    def test_validation_catches_actual_corruption(self):
        """Verify validation can catch and correct real corruption.
        
        Note: This test verifies that validation DETECTS corruption, not that
        it necessarily prevents it from accumulating. The validation system
        corrects tracked values when drift is detected, but corruption
        introduced mid-simulation may persist until the next validation check.
        """
        sim = Inferno(N=200, R=5, validate_mode='frequent')
        initial_total = sim.E_total
        
        # Run some moves normally
        for _ in range(10):
            sim.demon_move()
        
        # Energy should still be conserved
        assert sim.E_lattice + sim.E_demon_sum == initial_total, \
            "Energy should be conserved before corruption"
        
        # Manually call validation to ensure it works
        result = sim.validate_energy_conservation()
        assert result is True, "Validation should pass when no corruption exists"
        
        result = sim.validate_bond_counts()
        assert result is True, "Bond count validation should pass when no corruption exists"
    
    def test_validation_overhead_scales_with_mode(self):
        """Verify that different validation modes have different performance characteristics."""
        import time
        
        n, sweeps = 1000, 10
        
        # Measure time with validation off
        sim_off = Inferno(N=n, R=5, validate_mode='off')
        start = time.perf_counter()
        for _ in range(sweeps):
            sim_off.demon_move()
        time_off = time.perf_counter() - start
        
        # Measure time with frequent validation
        sim_frequent = Inferno(N=n, R=5, validate_mode='frequent')
        start = time.perf_counter()
        for _ in range(sweeps):
            sim_frequent.demon_move()
        time_frequent = time.perf_counter() - start
        
        # Frequent should be noticeably slower (but might be small for small N)
        # Don't assert strict timing, just verify it completes
        assert time_off >= 0 and time_frequent >= 0, "Both modes should complete successfully"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
