"""
Copyright (c) 2026 Winry Ember
Licensed under the MIT License.
See LICENSE file in the project root for full license information.

Unit tests for irreversible irrInferno simulation class.
"""

import os
import sys

import pytest

# Add creutz-sim directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))

from inferno_irr import irrInferno


class TestIrrInfernoInitialization:
    """Test irrInferno class initialization."""

    def test_basic_initialization(self):
        """Test that irrInferno initializes without errors."""
        x = irrInferno(100, 5)
        assert x.N == 100
        assert x._initial_total_energy == 200

    def test_energy_conservation_on_init(self):
        """Test that total energy equals 2*N on initialization."""
        for N in [10, 50, 100, 500]:
            x = irrInferno(N, 5)
            assert x.E_total == 2 * N


class TestIrrInfernoComparison:
    """Compare irrInferno with Inferno behavior."""

    def test_same_initialization_structure(self):
        """Test that both classes initialize with same structure."""
        from inferno import Inferno

        x_rev = Inferno(100, 5)
        x_irr = irrInferno(100, 5)

        assert x_rev.N == x_irr.N
        assert x_rev.E_total == x_irr.E_total
        assert len(x_rev.lattice) == len(x_irr.lattice)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
