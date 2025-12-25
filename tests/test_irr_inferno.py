"""Unit tests for irreversible irrInferno simulation class."""
import sys
import os
import pytest
import numpy as np

# Add creutz-sim directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'creutz-sim'))

from irr_inferno import irrInferno


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


class TestIrrInfernoEnergyConservation:
    """Test energy conservation during simulation."""
    
    def test_energy_conservation_forward(self):
        """Test energy conservation during forward moves."""
        x = irrInferno(100, 5)
        initial_energy = x.E_total
        
        # Run 1000 forward moves
        for _ in range(1000):
            x.demon_move()
        
        assert np.abs(x.E_total - initial_energy) < 1e-10
    
    def test_energy_conservation_reverse(self):
        """Test energy conservation during reverse moves."""
        x = irrInferno(100, 5)
        initial_energy = x.E_total
        
        # Run 1000 reverse moves
        for _ in range(1000):
            x.demon_reverse()
        
        assert np.abs(x.E_total - initial_energy) < 1e-10


class TestIrrInfernoIrreversibility:
    """Test that dynamics are indeed irreversible."""
    
    def test_irreversibility(self):
        """Test that forward-reverse does NOT return to initial state."""
        x = irrInferno(50, 3)
        
        # Store initial state
        initial_lattice = x.lattice.copy()
        
        # Run forward then reverse
        sweeps = 10
        for _ in range(sweeps):
            for _ in range(x.N):
                x.demon_move()
        
        for _ in range(sweeps):
            for _ in range(x.N):
                x.demon_reverse()
        
        # State should be different (irreversible)
        # Note: There's a very small chance they could be equal by coincidence
        # But for 50 sites and 10 sweeps, this is extremely unlikely
        assert not np.array_equal(x.lattice, initial_lattice), \
            "Irreversible simulation should not return to initial state"


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


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
