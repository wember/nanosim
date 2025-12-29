"""
Copyright (c) 2026 Winry Ember
Licensed under the MIT License.
See LICENSE file in the project root for full license information.

Unit tests for reversible Inferno simulation class."""

import os
import sys

import numpy as np
import pytest

# Add creutz-sim directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))

from inferno import Inferno


class TestInfernoInitialization:
    """Test Inferno class initialization."""

    def test_basic_initialization(self):
        """Test that Inferno initializes without errors."""
        x = Inferno(100, 5)
        assert x.N == 100
        assert x._initial_total_energy == 200

    def test_energy_conservation_on_init(self):
        """Test that total energy equals 2*N on initialization."""
        for N in [10, 50, 100, 500]:
            x = Inferno(N, 5)
            assert x.E_total == 2 * N
            assert x.E_total == x.E_lattice + np.sum(x.E_demon)

    def test_lattice_structure(self):
        """Test that lattice has correct size and spin values."""
        x = Inferno(100, 5)
        assert len(x.lattice) == 100
        assert np.all(np.abs(x.lattice) == 1)  # All spins are ±1

    def test_demon_energy_nonnegative(self):
        """Test that all demon energies are non-negative."""
        x = Inferno(100, 5)
        assert np.all(x.E_demon >= 0)


class TestInfernoStateMethods:
    """Test state inspection methods."""

    def test_get_validated_state(self):
        """Test that get_validated_state returns correct structure."""
        x = Inferno(100, 5)
        state = x.get_validated_state()

        assert "E_total" in state
        assert "E_lattice" in state
        assert "E_demon_sum" in state
        assert "bond_count" in state

    def test_state_consistency(self):
        """Test that state values are consistent."""
        x = Inferno(100, 5)
        state = x.get_validated_state()

        # Energy should sum correctly
        assert (
            np.abs(state["E_total"] - (state["E_lattice"] + state["E_demon_sum"]))
            < 1e-10
        )

        # Per-site values can be calculated
        K = state["E_demon_sum"] / x.N
        U = state["E_lattice"] / x.N
        assert K >= 0  # Demon energy should be non-negative
        assert -2 <= U <= 2  # Lattice energy per site bounded by ±2


class TestInfernoEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_small_lattice(self):
        """Test with very small lattice."""
        x = Inferno(4, 1)
        assert x.E_total == 8

        # Should run without errors
        for _ in range(10):
            x.demon_move()

    def test_large_radius(self):
        """Test with radius equal to lattice size."""
        x = Inferno(20, 20)

        # Should run without errors
        for _ in range(50):
            x.demon_move()

        assert x.E_total == 40

    def test_zero_radius(self):
        """Test with radius of 1 (minimum coupling)."""
        x = Inferno(50, 1)

        for _ in range(100):
            x.demon_move()

        assert x.E_total == 100


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
