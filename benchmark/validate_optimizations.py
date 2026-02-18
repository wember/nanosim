#!/usr/bin/env python3
"""
Validation script to verify optimizations maintain correctness.
Tests:
1. count_bonds() produces same results as before
2. Energy conservation holds
3. Small simulation runs without errors
"""

import numpy as np
import sys
sys.path.insert(0, '../creutz-sim')
from inferno import Inferno

def test_count_bonds():
    """Test that count_bonds() correctly counts bond types."""
    print("Testing count_bonds()...")
    
    # Create a small system
    x = Inferno(100, 5)
    
    # Count manually to verify
    manual_count = [0, 0, 0]
    for b in x.bonds:
        if b == -1:
            manual_count[0] += 1
        elif b == 0:
            manual_count[1] += 1
        elif b == 1:
            manual_count[2] += 1
    
    # Compare with method result
    assert np.array_equal(x.bond_count, manual_count), \
        f"count_bonds mismatch: got {x.bond_count}, expected {manual_count}"
    
    print(f"  ✓ count_bonds() correct: {x.bond_count}")

def test_energy_conservation():
    """Test that total energy is conserved during simulation."""
    print("\nTesting energy conservation...")
    
    N = 100
    x = Inferno(N, 5)
    E_total_init = x.E_total
    
    # Run a few demon moves
    for i in range(100):
        x.demon_move(flag=0, sweep_count=i)
        E_total_current = x.E_lattice + sum(x.E_demon)
        assert E_total_current == E_total_init, \
            f"Energy not conserved at step {i}: {E_total_current} != {E_total_init}"
    
    print(f"  ✓ Energy conserved over 100 steps: E_total = {E_total_init}")

def test_reversibility():
    """Test that forward + reverse returns to initial state."""
    print("\nTesting reversibility...")
    
    N = 100
    x = Inferno(N, 5)
    
    # Save initial state
    lattice_init = x.lattice.copy()
    bonds_init = x.bonds.copy()
    E_demon_init = x.E_demon.copy()
    E_lattice_init = x.E_lattice
    
    # Forward sweep
    for i in range(N):
        x.demon_move(flag=0, sweep_count=i)
    
    # Reverse sweep
    for i in range(N-1, -1, -1):
        x.demon_reverse(flag=0, sweep_count=i)
    
    # Check if returned to initial state
    assert np.array_equal(x.lattice, lattice_init), "Lattice not reversed!"
    assert np.array_equal(x.bonds, bonds_init), "Bonds not reversed!"
    assert np.array_equal(x.E_demon, E_demon_init), "E_demon not reversed!"
    assert x.E_lattice == E_lattice_init, "E_lattice not reversed!"
    
    print(f"  ✓ Reversibility validated: forward + reverse returns to initial state")

def test_spin_flip():
    """Test spin flip logic with various scenarios."""
    print("\nTesting spin_flip()...")
    
    x = Inferno(100, 5)
    E_total_before = x.E_lattice + sum(x.E_demon)
    
    # Test a few spin flips
    for _ in range(50):
        a = np.random.randint(0, x.N)
        i = np.random.randint(0, x.N)
        x.spin_flip(a, i)
        
        E_total_after = x.E_lattice + sum(x.E_demon)
        assert E_total_after == E_total_before, \
            f"spin_flip broke energy conservation: {E_total_after} != {E_total_before}"
    
    print(f"  ✓ spin_flip() maintains energy conservation")

def test_bond_change():
    """Test bond change logic."""
    print("\nTesting bond_change()...")
    
    x = Inferno(100, 5)
    E_total_before = x.E_lattice + sum(x.E_demon)
    
    # Test a few bond changes
    for _ in range(50):
        a = np.random.randint(0, x.N)
        i = np.random.randint(0, x.N)
        x.bond_change(a, i)
        
        E_total_after = x.E_lattice + sum(x.E_demon)
        assert E_total_after == E_total_before, \
            f"bond_change broke energy conservation: {E_total_after} != {E_total_before}"
    
    print(f"  ✓ bond_change() maintains energy conservation")

def main():
    """Run all validation tests."""
    print("="*60)
    print("Validating Inferno Optimizations")
    print("="*60)
    
    try:
        test_count_bonds()
        test_energy_conservation()
        test_reversibility()
        test_spin_flip()
        test_bond_change()
        
        print("\n" + "="*60)
        print("✅ All validation tests passed!")
        print("="*60)
        return 0
        
    except AssertionError as e:
        print(f"\n❌ Validation failed: {e}")
        return 1
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
