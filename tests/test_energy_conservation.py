"""
Energy conservation and reversibility tests for both Inferno and irrInferno.

Tests cover:
- Energy conservation during forward/reverse moves
- Reversibility in the Inferno class (time-reversal symmetry)
- Irreversibility in the irrInferno class (broken symmetry)
- Comparison between reversible and irreversible dynamics
"""

import sys
import os
import pytest
import numpy as np

# Add creutz-sim directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'creutz-sim'))

from inferno import Inferno
from irr_inferno import irrInferno


class TestInfernoEnergyConservation:
    """Test energy conservation in reversible Inferno class."""
    
    def test_energy_conservation_forward(self):
        """Test energy conservation during forward moves."""
        x = Inferno(100, 5)
        initial_energy = x.E_total
        
        # Run 1000 forward moves
        for _ in range(1000):
            x.demon_move()
        
        assert x.E_total == initial_energy
        assert x.E_lattice + np.sum(x.E_demon) == initial_energy
    
    def test_energy_conservation_reverse(self):
        """Test energy conservation during reverse moves."""
        x = Inferno(100, 5)
        initial_energy = x.E_total
        
        # Run 1000 reverse moves
        for _ in range(1000):
            x.demon_reverse()
        
        assert x.E_total == initial_energy
        assert x.E_lattice + np.sum(x.E_demon) == initial_energy
    
    def test_energy_conservation_forward_reverse(self):
        """Test energy conservation over combined forward and reverse moves."""
        x = Inferno(100, 5)
        initial_energy = x.E_total
        
        # Interleave forward and reverse moves
        for _ in range(100):
            for _ in range(10):
                x.demon_move()
            for _ in range(10):
                x.demon_reverse()
        
        assert x.E_total == initial_energy


class TestInfernoReversibility:
    """Test reversibility (time-reversal symmetry) in Inferno class."""
    
    def test_reversibility_small(self):
        """Test that forward then reverse returns to initial state."""
        x = Inferno(50, 5)
        
        # Save initial state
        initial_lattice = x.lattice.copy()
        initial_E_demon = x.E_demon.copy()
        initial_bonds = x.bonds.copy()
        
        # Run forward then reverse
        for _ in range(100):
            x.demon_move()
        
        for _ in range(100):
            x.demon_reverse()
        
        # Should return to initial state
        np.testing.assert_array_equal(x.lattice, initial_lattice)
        np.testing.assert_array_equal(x.E_demon, initial_E_demon)
        np.testing.assert_array_equal(x.bonds, initial_bonds)
    
    def test_reversibility_different_radii(self):
        """Test reversibility with different demon coupling radii."""
        for R in [1, 3, 7]:
            x = Inferno(100, R)
            initial_lattice = x.lattice.copy()
            initial_E_demon = x.E_demon.copy()
            
            # Forward and reverse
            for _ in range(100):
                x.demon_move()
            for _ in range(100):
                x.demon_reverse()
            
            np.testing.assert_array_equal(x.lattice, initial_lattice,
                                        err_msg=f"Reversibility failed for R={R}")
            np.testing.assert_array_equal(x.E_demon, initial_E_demon,
                                        err_msg=f"Demon energy reversibility failed for R={R}")


class TestIrrInfernoEnergyConservation:
    """Test energy conservation in irreversible irrInferno class."""
    
    def test_energy_conservation_forward(self):
        """Test energy conservation during forward moves."""
        x = irrInferno(100, 5)
        initial_energy = x.E_total
        
        # Run 1000 forward moves
        for _ in range(1000):
            x.demon_move()
        
        assert x.E_total == initial_energy
        assert x.E_lattice + np.sum(x.E_demon) == initial_energy
    
    def test_energy_conservation_reverse(self):
        """Test energy conservation during reverse moves."""
        x = irrInferno(100, 5)
        initial_energy = x.E_total
        
        # Run 1000 reverse moves
        for _ in range(1000):
            x.demon_reverse()
        
        assert x.E_total == initial_energy
        assert x.E_lattice + np.sum(x.E_demon) == initial_energy


class TestIrrInfernoIrreversibility:
    """Test that irrInferno breaks time-reversal symmetry."""
    
    def test_irreversibility(self):
        """Test that forward then reverse does NOT return to initial state.
        
        This is a statistical test - we expect the state to be different
        because the random radii are regenerated on each move.
        """
        x = irrInferno(50, 5)
        
        # Save initial state
        initial_lattice = x.lattice.copy()
        
        # Run forward then reverse (same number of moves)
        for _ in range(100):
            x.demon_move()
        
        for _ in range(100):
            x.demon_reverse()
        
        # Should NOT return to initial state (statistical test)
        # Allow small probability of returning by chance
        differences = np.sum(x.lattice != initial_lattice)
        assert differences > 0, \
            "State should differ after forward-reverse cycle (irreversible dynamics)"


class TestReversibleVsIrreversible:
    """Compare reversible and irreversible dynamics."""
    
    def test_same_initialization_structure(self):
        """Test that both classes initialize with same structure."""
        x_rev = Inferno(100, 5)
        x_irr = irrInferno(100, 5)
        
        # Same lattice size
        assert x_rev.N == x_irr.N
        
        # Same total energy
        assert x_rev.E_total == x_irr.E_total
        
        # Same energy components structure
        assert len(x_rev.E_demon) == len(x_irr.E_demon)
        assert len(x_rev.bonds) == len(x_irr.bonds)
    
    def test_energy_distribution_similar(self):
        """Test that initial energy distributions are statistically similar."""
        np.random.seed(42)
        x_rev = Inferno(1000, 5)
        
        np.random.seed(42)
        x_irr = irrInferno(1000, 5)
        
        # Same total energy
        assert x_rev.E_total == x_irr.E_total
        
        # Similar mean demon energy
        assert abs(np.mean(x_rev.E_demon) - np.mean(x_irr.E_demon)) < 0.5
    
    def test_dynamics_preserve_energy_both_classes(self):
        """Test that both classes preserve energy identically."""
        for SimClass in [Inferno, irrInferno]:
            sim = SimClass(100, 5)
            initial_energy = sim.E_total
            
            # Run mixed forward/reverse
            for _ in range(50):
                sim.demon_move()
                sim.demon_reverse()
            
            assert sim.E_total == initial_energy, \
                f"{SimClass.__name__} failed to conserve energy"


class TestLargeSystemStability:
    """Test energy conservation in large systems over long runs."""
    
    @pytest.mark.slow
    def test_inferno_large_system(self):
        """Test energy conservation in large Inferno system."""
        x = Inferno(N=100000, R=5, validate_mode='off')
        initial_energy = x.E_total
        
        # Run 40,000 sweeps (4 billion moves)
        for _ in range(40000):
            x.demon_move()
        
        # Energy should be exactly conserved (integer arithmetic)
        assert x.E_total == initial_energy
        assert x.E_lattice + x.E_demon_sum == initial_energy
    
    @pytest.mark.slow
    def test_irr_inferno_large_system(self):
        """Test energy conservation in large irrInferno system."""
        x = irrInferno(N=100000, R=5, validate_mode='off')
        initial_energy = x.E_total
        
        # Run 40,000 sweeps
        for _ in range(40000):
            x.demon_move()
        
        # Energy should be exactly conserved
        assert x.E_total == initial_energy
        assert x.E_lattice + x.E_demon_sum == initial_energy


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
