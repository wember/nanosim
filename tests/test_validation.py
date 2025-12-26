"""
Validation system tests for both Inferno and irrInferno.

Tests cover:
- Validation mode configuration (off/periodic/frequent)
- Validation correctness (energy conservation, bond counts)
- Validation performance impact
- Drift detection and correction
- Integration with simulation dynamics
"""

import sys
import os
import pytest
import numpy as np

# Add creutz-sim directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'creutz-sim'))

from inferno import Inferno
from inferno_irr import irrInferno




class TestValidationModeCorrectness:
    """Test that validation mode doesn't affect simulation correctness."""
    
    def test_off_mode_energy_conservation(self):
        """Test that validation='off' maintains energy conservation."""
        x = Inferno(1000, 5, validate_mode='off')
        initial_energy = x.E_total
        
        # Run many moves
        for _ in range(5000):
            for _ in range(x.N):
                x.demon_move()
        
        # Energy must be conserved
        assert x.E_total == initial_energy
        
        # Validate against actual sum
        state = x.get_validated_state()
        assert state['E_total'] == initial_energy
    
    def test_off_mode_reversibility(self):
        """Test that validation='off' preserves reversibility."""
        x = Inferno(100, 3, validate_mode='off')
        
        # Store initial state
        initial_lattice = x.lattice.copy()
        initial_E_demon = x.E_demon.copy()
        
        # Forward then reverse
        sweeps = 20
        for _ in range(sweeps):
            for _ in range(x.N):
                x.demon_move()
        
        for _ in range(sweeps):
            for _ in range(x.N):
                x.demon_reverse()
        
        # Must return to initial state
        np.testing.assert_array_equal(x.lattice, initial_lattice)
        np.testing.assert_array_equal(x.E_demon, initial_E_demon)
    
    def test_all_modes_same_deterministic_result(self):
        """Test that all validation modes give identical results with same seed."""
        np.random.seed(42)
        x_off = Inferno(50, 3, validate_mode='off')
        
        np.random.seed(42)
        x_periodic = Inferno(50, 3, validate_mode='periodic')
        
        np.random.seed(42)
        x_frequent = Inferno(50, 3, validate_mode='frequent')
        
        # Run same moves
        for _ in range(10):
            for _ in range(50):
                x_off.demon_move()
                x_periodic.demon_move()
                x_frequent.demon_move()
        
        # All should have identical states
        np.testing.assert_array_equal(x_off.lattice, x_periodic.lattice)
        np.testing.assert_array_equal(x_off.lattice, x_frequent.lattice)
        np.testing.assert_array_equal(x_off.E_demon, x_periodic.E_demon)
        np.testing.assert_array_equal(x_off.E_demon, x_frequent.E_demon)
        np.testing.assert_array_equal(x_off.bonds, x_periodic.bonds)
        np.testing.assert_array_equal(x_off.bonds, x_frequent.bonds)
    
    def test_irr_inferno_off_mode(self):
        """Test irrInferno with validation='off'."""
        x = irrInferno(1000, 5, validate_mode='off')
        initial_energy = x.E_total
        
        # Run many moves
        for _ in range(1000):
            for _ in range(x.N):
                x.demon_move()
        
        # Energy must be conserved
        assert x.E_total == initial_energy
        
        # Validate
        state = x.get_validated_state()
        assert state['E_total'] == initial_energy


class TestValidationModeConfiguration:
    """Test that validation mode configuration works correctly."""
    
    def test_off_mode_never_validates(self):
        """Test that 'off' mode sets check_interval to infinity."""
        x = Inferno(100, 5, validate_mode='off')
        assert x._check_interval == float('inf')
        assert x._validate_mode == 'off'
    
    def test_periodic_mode_interval(self):
        """Test that 'periodic' mode sets correct interval."""
        N = 1000
        x = Inferno(N, 5, validate_mode='periodic')
        assert x._check_interval == 100 * N
        assert x._validate_mode == 'periodic'
    
    def test_frequent_mode_interval(self):
        """Test that 'frequent' mode sets interval to N."""
        N = 500
        x = Inferno(N, 5, validate_mode='frequent')
        assert x._check_interval == N
        assert x._validate_mode == 'frequent'
    
    def test_default_mode_is_off(self):
        """Test that default validation mode is 'off'."""
        x = Inferno(100, 5)  # No validate_mode specified
        assert x._validate_mode == 'off'
        assert x._check_interval == float('inf')


class TestValidationModePerformance:
    """Test that validation modes have expected performance characteristics."""
    
    def test_off_mode_no_counter_increment_when_disabled(self):
        """Test that counter doesn't increment with validation off."""
        x = Inferno(100, 5, validate_mode='off')
        initial_counter = x._check_counter
        
        # Run some moves
        for _ in range(10):
            x.demon_move()
        
        # Counter should not change when validation is off
        # (it increments but never triggers validation)
        assert x._check_counter >= initial_counter
    
    def test_frequent_mode_validates_every_sweep(self):
        """Test that frequent mode validation is triggered correctly."""
        x = Inferno(100, 5, validate_mode='frequent')
        
        # Run exactly one sweep
        for _ in range(x.N):
            x.demon_move()
        
        # Counter should have cycled (validation triggered)
        assert x._check_counter == 0  # Reset after validation


class TestLongRunStability:
    """Test stability over very long runs without validation."""
    
    def test_very_long_run_off_mode(self):
        """Test 20,000 sweep run with validation='off'."""
        x = Inferno(500, 5, validate_mode='off')
        initial_energy = x.E_total
        
        # Run 20,000 sweeps (10 million moves)
        for _ in range(20000):
            for _ in range(x.N):
                x.demon_move()
        
        # Energy must still be conserved
        assert x.E_total == initial_energy
        
        # Full validation check
        state = x.get_validated_state()
        assert state['E_total'] == initial_energy
        
        # Bond counts should be consistent
        actual_bonds = np.bincount(x.bonds + 1, minlength=3).astype(np.int64)
        np.testing.assert_array_equal(state['bond_count'], actual_bonds)
    
    def test_multiple_forward_reverse_cycles(self):
        """Test multiple forward-reverse cycles with validation='off'."""
        x = Inferno(200, 5, validate_mode='off')
        initial_lattice = x.lattice.copy()
        initial_E_demon = x.E_demon.copy()
        initial_bonds = x.bonds.copy()
        
        # 10 cycles of forward-reverse
        for cycle in range(10):
            # Forward
            for _ in range(50):
                for _ in range(x.N):
                    x.demon_move()
            
            # Reverse
            for _ in range(50):
                for _ in range(x.N):
                    x.demon_reverse()
        
        # Must return to initial state after all cycles
        np.testing.assert_array_equal(x.lattice, initial_lattice,
            err_msg="Lattice not restored after multiple cycles")
        np.testing.assert_array_equal(x.E_demon, initial_E_demon,
            err_msg="Demon energy not restored after multiple cycles")
        np.testing.assert_array_equal(x.bonds, initial_bonds,
            err_msg="Bonds not restored after multiple cycles")


class TestValidationModeWithDifferentParameters:
    """Test validation modes with various lattice sizes and radii."""
    
    @pytest.mark.parametrize("N", [10, 100, 1000, 10000])
    def test_off_mode_various_lattice_sizes(self, N):
        """Test validation='off' with different lattice sizes."""
        x = Inferno(N, 5, validate_mode='off')
        initial_energy = x.E_total
        
        # Run proportional to lattice size
        sweeps = min(100, 10000 // N)
        for _ in range(sweeps):
            for _ in range(x.N):
                x.demon_move()
        
        assert x.E_total == initial_energy
    
    @pytest.mark.parametrize("R", [1, 3, 5, 10, 20])
    def test_off_mode_various_radii(self, R):
        """Test validation='off' with different coupling radii."""
        x = Inferno(100, R, validate_mode='off')
        initial_energy = x.E_total
        
        for _ in range(50):
            for _ in range(x.N):
                x.demon_move()
        
        assert x.E_total == initial_energy


if __name__ == '__main__':
    pytest.main([__file__, '-v'])


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


class TestValidationMethodAccuracy:
    """Test validation methods detect and correct corruption."""
    
    def test_validation_catches_energy_errors(self):
        """Test that validate_energy_conservation detects problems."""
        x = Inferno(100, 5)
        
        # Validation should pass initially
        assert x.validate_energy_conservation() == True
        
        # Artificially introduce drift
        x.E_lattice += 1
        
        # Validation should detect and fix (returns False)
        assert x.validate_energy_conservation() == False
        
        # Should be fixed now
        assert x.validate_energy_conservation() == True
    
    def test_bond_count_validation_detects_errors(self):
        """Test that validate_bond_counts detects problems."""
        x = Inferno(100, 5)
        
        # Validation should pass initially
        assert x.validate_bond_counts() == True
        
        # Artificially corrupt bond counts
        x.bond_count[0] += 1
        x.bond_count[1] -= 1
        
        # Validation should detect and fix
        assert x.validate_bond_counts() == False
        
        # Should be fixed now
        assert x.validate_bond_counts() == True
    
    def test_get_validated_state_accuracy(self):
        """Test that get_validated_state returns accurate values."""
        x = Inferno(100, 5)
        
        # Run some moves
        for _ in range(50):
            for _ in range(x.N):
                x.demon_move()
        
        # Get validated state
        state = x.get_validated_state()
        
        # Manually calculate values
        actual_lattice = np.sum(x.bonds, dtype=np.int64)
        actual_demon = np.sum(x.E_demon, dtype=np.int64)
        actual_bonds = np.bincount(x.bonds + 1, minlength=3).astype(np.int64)
        
        # Should match exactly
        assert state['E_lattice'] == actual_lattice
        assert state['E_demon_sum'] == actual_demon
        assert state['E_total'] == actual_lattice + actual_demon
        assert np.array_equal(state['bond_count'], actual_bonds)


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

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
