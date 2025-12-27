"""Performance benchmarks and stress tests for Inferno simulation."""

import os
import sys
import time

import numpy as np
import pytest

# Add creutz-sim directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))

from inferno import Inferno


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
        """Test for energy drift over a substantial (but faster) simulation."""
        x = Inferno(1000, 5)
        initial_energy = x.E_total

        # Run 2000 sweeps (reduced from 10000 for speed)
        for _ in range(2000):
            for _ in range(x.N):
                x.demon_move()

        # Check no drift
        assert np.abs(x.E_total - initial_energy) < 1e-10

        # Validate against actual sum
        state = x.get_validated_state()
        assert state["E_total"] == initial_energy

    def test_extreme_radius(self):
        """Test with very large coupling radius."""
        x = Inferno(1000, 500)  # R = N/2
        initial_energy = x.E_total

        for _ in range(100):
            for _ in range(x.N):
                x.demon_move()

        assert np.abs(x.E_total - initial_energy) < 1e-10


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


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])  # -s to show print statements
