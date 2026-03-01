# Force single-threaded BLAS/LAPACK on Linux HPC to prevent thread over-subscription
# macOS/Apple Silicon has good thread management, so skip there
import os
import platform
if platform.system() == 'Linux':
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
    os.environ['NUMEXPR_NUM_THREADS'] = '1'

import numpy as np
import random
import math
from numba import njit

# test #
from scipy.special import factorial as f
from scipy.special import loggamma as logg
Sk = lambda N, K: logg(K + N) - logg(K+1) - logg(N) # N == lattice size, K == kinetic energy
Su = lambda N, N0, Nx: logg(N+1) + np.log(2**(N0)) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1)) # N == lattice size, N0 == broken bonds, Nx == bonds between anti-aligned spins
Su0 = lambda N, N0, Nx: logg(N+1) + np.log(2**(N0+1)) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1)) # N == lattice size, N0 == broken bonds, Nx == bonds between anti-aligned spins


################################################################################
# JIT-Compiled Functions (Hot Path Optimization)
################################################################################

@njit
def spin_flip_jit(lattice, bonds, E_demon, E_lattice, d_energy, a, i, N):
    """
    JIT-compiled spin flip with Metropolis criterion.
    
    Args:
        lattice: array of spin values
        bonds: array of bond values
        E_demon: array of demon energies
        E_lattice: lattice energy
        d_energy: total demon energy
        a: lattice site index
        i: demon index
        N: lattice size
    
    Returns:
        (E_lattice, d_energy): updated energies
    """
    s = lattice[a]
    
    # Calculate energy based on nearest neighbors
    nb = lattice[(a+1)%N]*abs(bonds[a%N]) + lattice[(a-1)%N]*abs(bonds[(a-1)%N])
    cost = 2*s*nb
    
    # Metropolis acceptance: flip if favorable or if demon has energy
    if cost < 0 or cost <= E_demon[i]:
        s *= -1
        E_demon[i] -= cost
        d_energy -= cost
        E_lattice += cost
        lattice[a] = s
        
        # Update bonds
        if bonds[a] != 0:
            if lattice[a] == lattice[(a+1)%N]:
                bonds[a] = -1
            else:
                bonds[a] = 1
        
        if bonds[(a-1)%N] != 0:
            if lattice[a] == lattice[(a-1)%N]:
                bonds[(a-1)%N] = -1
            else:
                bonds[(a-1)%N] = 1
    
    return E_lattice, d_energy


@njit
def bond_change_jit(lattice, bonds, E_demon, E_lattice, d_energy, a, i, N):
    """
    JIT-compiled bond creation/breaking.
    
    Args:
        lattice: array of spin values
        bonds: array of bond values
        E_demon: array of demon energies
        E_lattice: lattice energy
        d_energy: total demon energy
        a: lattice site index
        i: demon index
        N: lattice size
    
    Returns:
        (E_lattice, d_energy): updated energies
    """
    s = lattice[a]
    b = bonds[a]
    n = lattice[(a+1)%N]
    cost = -1 if s == n else 1
    
    # If bond is broken, attempt to remake
    if (b == 0) and (E_demon[i] >= cost):
        E_lattice += cost
        E_demon[i] -= cost
        d_energy -= cost
        bonds[a] = cost
    
    # If bond is made, attempt to break
    elif (b != 0) and (E_demon[i] + cost >= 0):
        E_lattice -= cost
        E_demon[i] += cost
        d_energy += cost
        bonds[a] = 0
    
    # Update neighbor bond if it exists
    if bonds[(a-1)%N] != 0:
        if lattice[a] == lattice[(a-1)%N]:
            bonds[(a-1)%N] = -1
        else:
            bonds[(a-1)%N] = 1
    
    return E_lattice, d_energy


@njit
def count_bonds_jit(bonds):
    """
    JIT-compiled bond counting.
    
    Args:
        bonds: array of bond values
    
    Returns:
        bond_count: array [count_of_-1, count_of_0, count_of_1]
    """
    bond_count = np.zeros(3, dtype=np.int64)
    for b in bonds:
        if b == -1:
            bond_count[0] += 1
        elif b == 0:
            bond_count[1] += 1
        elif b == 1:
            bond_count[2] += 1
    return bond_count


################################################################################
# Inferno Class
################################################################################

class Inferno:
    """
        Inferno:
            - Main class for implementing microcanonical Monte Carlo simulation

        :instance methods:
            calc_E_lat - calculates the energy of a given latice configuration
            demon_move - updates the lattice by moving the demon around
    """

    def __init__(self,N, R):
        """
            :params:
                N - size of lattice
                lattice - state of lattice
                bonds - state of the bonds
                E_lattice - energy of lattice
                E_demon - energy of the demon
        """
        # every integer from 0-N placed in a random order
        a = np.arange(N)
        np.random.shuffle(a)

        total_energy = N//2 # total energy of the system

        self.N = N
        self.order = a
        self.rev_order = np.flip(a)
        self.order_idx = 0  # Index pointer for order array
        self.rev_order_idx = 0  # Index pointer for rev_order array
        self.radius = R
        self.radius_cycle = 2 * R + 1
        self.R_counter = 0
        # self.radius_spin = self.rev_radius_bond = np.random.randint(0, R, size=N)*np.random.choice([-1, 1], size=N)
        # self.rev_radius_spin = self.radius_bond = np.flip(self.radius_spin)
        self.lattice = np.concatenate((np.ones(N//2, dtype=int), (-1)*np.ones(N//2, dtype=int)))
        self.bonds = np.ones(N, dtype=int)*(-1)
        self.bonds[[N//2-1, -1]] = 1
        self.bond_count = np.ones(3, dtype=int)
        self.count_bonds()
        self.E_lattice = sum(self.bonds)

        self.d_energy = total_energy - self.E_lattice
        # randomly assign energy to einstein oscillators
        result = np.zeros(N, dtype=int)
        for i in range(self.d_energy):
            result[random.randint(0,N-1)] += 1

        self.E_demon = np.array(result)
        self.E_total = self.E_lattice + sum(self.E_demon)

################################################################################
################################################################################
################################################################################
    def reset(self):
        """
            Resets the "random" order for reversible simulations
        """     
        a = np.arange(self.N)
        np.random.shuffle(a)
        self.order = a
        self.rev_order = np.flip(a)

    def spin_flip(self, a, i):
        """
            Attempt to flip the spin of a given lattice site (JIT-optimized)
        """
        self.E_lattice, self.d_energy = spin_flip_jit(
            self.lattice, self.bonds, self.E_demon,
            self.E_lattice, self.d_energy, a, i, self.N
        )

    def bond_change(self, a, i):
        """
            Attempt to change the bond given lattice site (JIT-optimized)
        """
        self.E_lattice, self.d_energy = bond_change_jit(
            self.lattice, self.bonds, self.E_demon,
            self.E_lattice, self.d_energy, a, i, self.N
        )

    def count_bonds(self):
        """
            Updates the bond-count array of number of aligned (-1), broken (0), and misaligned (1) bonds (JIT-optimized)
        """
        self.bond_count = count_bonds_jit(self.bonds)

    def demon_move(self, flag, sweep_count):
        """
            "Randomly" move the demon around and flip spins & change bonds
        """
        lattice = self.lattice
        bonds = self.bonds
        E_demon = self.E_demon
        E_lattice = self.E_lattice
        d_energy = self.d_energy
        N = self.N

        a = self.order[self.order_idx]
        R = (self.R_counter % self.radius_cycle) - self.radius
        self.R_counter += 1

        # If irr flag is on, generate a random number instead
        if (flag != 0):
            a = np.random.randint(0, N)
            if self.radius != 0:
                R = np.random.randint(0, self.radius)

        # Attempt to flip spin
        E_lattice, d_energy = spin_flip_jit(
            lattice, bonds, E_demon, E_lattice, d_energy, a, (a + R) % N, N
        )

        R = (self.R_counter % self.radius_cycle) - self.radius
        self.R_counter += 1

        # If irr flag is on, generate a random number instead
        if (flag != 0):
            if self.radius != 0:
                R = np.random.randint(0, self.radius)

        # Attempt to change bond
        E_lattice, d_energy = bond_change_jit(
            lattice, bonds, E_demon, E_lattice, d_energy, a, (a + R) % N, N
        )

        self.E_lattice = E_lattice
        self.d_energy = d_energy

        # Update bond count
        self.bond_count = count_bonds_jit(self.bonds)

        # Advance index pointer (replaces np.roll)
        if (flag == 0):
            self.order_idx = (self.order_idx + 1) % self.N

    def demon_reverse(self, flag, sweep_count):
        """
            In reverse order, flip spins & change bonds
        """
        lattice = self.lattice
        bonds = self.bonds
        E_demon = self.E_demon
        E_lattice = self.E_lattice
        d_energy = self.d_energy
        N = self.N

        a = self.rev_order[self.rev_order_idx]
        self.R_counter -= 1
        R = (self.R_counter % self.radius_cycle) - self.radius

        # If irr flag is on, generate a random number instead
        if (flag != 0):
            a = np.random.randint(0, N)
            if self.radius != 0:
                R = np.random.randint(0, self.radius)

        # Attempt to change bond
        E_lattice, d_energy = bond_change_jit(
            lattice, bonds, E_demon, E_lattice, d_energy, a, (a + R) % N, N
        )

        self.R_counter -= 1
        R = (self.R_counter % self.radius_cycle) - self.radius

        # If irr flag is on, generate a random number instead
        if (flag != 0):
            if self.radius != 0:
                R = np.random.randint(0, self.radius)

        # Attempt to flip spin
        E_lattice, d_energy = spin_flip_jit(
            lattice, bonds, E_demon, E_lattice, d_energy, a, (a + R) % N, N
        )

        self.E_lattice = E_lattice
        self.d_energy = d_energy

        # Update bond count
        self.bond_count = count_bonds_jit(self.bonds)

        # Advance index pointer (replaces np.roll)
        if (flag == 0):
            self.rev_order_idx = (self.rev_order_idx + 1) % self.N
