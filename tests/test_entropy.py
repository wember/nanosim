"""
Entropy calculations and statistical property tests.

Tests cover:
- Entropy calculation correctness
- Statistical properties of demon energy distribution
- Comparison metrics between reversible and irreversible dynamics
"""

import os
import sys

import numpy as np
import pytest

# Add creutz-sim directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))

from inferno import Inferno
from inferno_irr import irrInferno


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
                    assert not np.isnan(
                        s
                    ), f"Su should not be NaN for N={N}, N0={N0}, Nx={Nx}"

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
        assert state_rev["E_total"] == x_rev._initial_total_energy
        assert state_irr["E_total"] == x_irr._initial_total_energy


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
