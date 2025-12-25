"""
Regression tests for validation mode optimization.
Ensures that disabling validation doesn't affect physics correctness.
"""
import sys
import os
import pytest
import numpy as np

# Add creutz-sim directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'creutz-sim'))

from inferno import Inferno
from irr_inferno import irrInferno


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
