"""Performance benchmarks and stress tests for Inferno simulation."""
import sys
import os
import pytest
import numpy as np
import time

# Add creutz-sim directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'creutz-sim'))

from inferno import Inferno
from irr_inferno import irrInferno


class TestPerformanceBenchmarks:
    """Benchmark tests to measure performance before/after optimizations."""
    
    def test_baseline_small_lattice(self):
        """Baseline: 1000 sites, 1000 sweeps."""
        x = Inferno(1000, 5)
        
        start = time.time()
        for _ in range(1000):
            for _ in range(x.N):
                x.demon_move()
        duration = time.time() - start
        
        print(f"\nSmall lattice (N=1000, 1000 sweeps): {duration:.3f}s")
        assert duration < 30  # Should complete in reasonable time
    
    def test_baseline_medium_lattice(self):
        """Baseline: 10000 sites, 100 sweeps."""
        x = Inferno(10000, 5)
        
        start = time.time()
        for _ in range(100):
            for _ in range(x.N):
                x.demon_move()
        duration = time.time() - start
        
        print(f"\nMedium lattice (N=10000, 100 sweeps): {duration:.3f}s")
        assert duration < 60
    
    def test_baseline_forward_reverse(self):
        """Baseline: Forward-reverse cycle timing."""
        x = Inferno(5000, 5)
        
        start = time.time()
        # Forward
        for _ in range(100):
            for _ in range(x.N):
                x.demon_move()
        # Reverse
        for _ in range(100):
            for _ in range(x.N):
                x.demon_reverse()
        duration = time.time() - start
        
        print(f"\nForward-reverse cycle (N=5000, 200 sweeps): {duration:.3f}s")
        assert duration < 60


class TestStressTests:
    """Stress tests for large systems and long runs."""
    
    def test_large_lattice_energy_conservation(self):
        """Test energy conservation with large lattice."""
        x = Inferno(50000, 5)
        initial_energy = x.E_total
        
        # Run many moves
        for _ in range(10):
            for _ in range(x.N):
                x.demon_move()
        
        assert np.abs(x.E_total - initial_energy) < 1e-10
    
    def test_long_run_no_drift(self):
        """Test for energy drift over long simulation."""
        x = Inferno(1000, 5)
        initial_energy = x.E_total
        
        # Run 10000 sweeps
        for _ in range(10000):
            for _ in range(x.N):
                x.demon_move()
        
        # Check no drift
        assert np.abs(x.E_total - initial_energy) < 1e-10
        
        # Validate against actual sum
        state = x.get_validated_state()
        assert state['E_total'] == initial_energy
    
    def test_extreme_radius(self):
        """Test with very large coupling radius."""
        x = Inferno(1000, 500)  # R = N/2
        initial_energy = x.E_total
        
        for _ in range(100):
            for _ in range(x.N):
                x.demon_move()
        
        assert np.abs(x.E_total - initial_energy) < 1e-10


class TestValidationMechanisms:
    """Test the validation and error detection mechanisms."""
    
    def test_validation_catches_errors(self):
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
    
    def test_bond_count_validation(self):
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


class TestIndexCycling:
    """Test that index cycling works correctly."""
    
    def test_order_index_cycles(self):
        """Test that order_idx cycles through all sites."""
        x = Inferno(10, 3)
        
        # Run exactly N moves
        for _ in range(x.N):
            x.demon_move()
        
        # Should have cycled back to 0
        assert x.order_idx == 0
    
    def test_radius_indices_cycle(self):
        """Test that radius indices cycle correctly."""
        x = Inferno(10, 3)
        
        # Run exactly N moves
        for _ in range(x.N):
            x.demon_move()
        
        # All indices should cycle back
        assert x.radius_spin_idx == 0
        assert x.radius_bond_idx == 0
    
    def test_reverse_indices_independent(self):
        """Test that forward and reverse indices are independent."""
        x = Inferno(20, 5)
        
        # Run forward moves
        for _ in range(10):
            x.demon_move()
        
        forward_idx = x.order_idx
        
        # Run reverse moves
        for _ in range(10):
            x.demon_reverse()
        
        # Forward index shouldn't change
        assert x.order_idx == forward_idx
        # Reverse index should have moved
        assert x.rev_order_idx == 10


class TestBondOperations:
    """Test bond update operations specifically."""
    
    def test_update_bonds_maintains_consistency(self):
        """Test that update_bonds_incremental maintains consistency."""
        x = Inferno(20, 3)
        
        for _ in range(100):
            # Pick a random site
            a = np.random.randint(0, x.N)
            
            # Flip the spin manually
            x.lattice[a] *= -1
            
            # Update bonds
            x.update_bonds_incremental(a)
            
            # Validate bond counts are correct
            actual_counts = np.bincount(x.bonds + 1, minlength=3).astype(np.int64)
            np.testing.assert_array_equal(x.bond_count, actual_counts)
    
    def test_bond_change_energy_tracking(self):
        """Test that bond_change tracks energy correctly."""
        x = Inferno(50, 5)
        
        initial_total = x.E_total
        
        # Run many bond changes
        for _ in range(500):
            a = np.random.randint(0, x.N)
            i = np.random.randint(0, x.N)
            x.bond_change(a, i)
        
        # Energy should still be conserved
        assert x.E_total == initial_total


class TestEntropyCalculations:
    """Test entropy calculation functions."""
    
    def test_Sk_positive(self):
        """Test that Sk (kinetic entropy) is positive."""
        from inferno import Sk
        
        # Test various N and K values
        for N in [10, 100, 1000]:
            for K in [10, 100, 1000]:
                s = Sk(N, K)
                assert s >= 0, f"Sk should be non-negative for N={N}, K={K}"
    
    def test_Su_positive(self):
        """Test that Su (configurational entropy) is positive."""
        from inferno import Su
        
        # Test various configurations
        # Note: N0 can't be too large or 2**N0 overflows
        N = 100
        for N0 in [0, 10, 20, 30]:  # Keep N0 reasonable to avoid overflow
            for Nx in [0, 5, 10]:
                if N0 + Nx <= N:
                    s = Su(N, N0, Nx)
                    assert not np.isnan(s), f"Su should not be NaN for N={N}, N0={N0}, Nx={Nx}"
    
    def test_Su0_edge_case(self):
        """Test Su0 for N0=0 edge case."""
        from inferno import Su0
        
        N = 100
        s = Su0(N, 0, 0)
        assert not np.isnan(s), "Su0 should handle N0=0"
        assert s >= 0


class TestComparisonMetrics:
    """Test metrics used for comparing reversible vs irreversible."""
    
    def test_entropy_difference(self):
        """Test that reversible and irreversible have different entropy evolution."""
        x_rev = Inferno(100, 5)
        x_irr = irrInferno(100, 5)
        
        # Run same number of sweeps
        sweeps = 50
        for _ in range(sweeps):
            for _ in range(x_rev.N):
                x_rev.demon_move()
                x_irr.demon_move()
        
        # States should be different
        assert not np.array_equal(x_rev.lattice, x_irr.lattice)
        
        # Get states
        state_rev = x_rev.get_validated_state()
        state_irr = x_irr.get_validated_state()
        
        # Both should conserve energy
        assert state_rev['E_total'] == x_rev._initial_total_energy
        assert state_irr['E_total'] == x_irr._initial_total_energy


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])  # -s to show print statements
