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
        self.r_idx = 0
        self.radius = R
        offsets = np.array([0] + list(range(1, R+1)) + list(range(-R, 0))).reshape(-1, 1)
        self.d_order = (a + offsets)  % N
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
        a = np.arange(self.N)
        np.random.shuffle(a)
        self.order = a
        self.rev_order = np.flip(a)
        # Rebuild d_order to match new order
        offsets = np.array([0] + list(range(1, self.radius+1)) + list(range(-self.radius, 0))).reshape(-1, 1)
        self.d_order = (a + offsets) % self.N
        self.order_idx = 0  # Also reset index pointers
        self.r_idx = 0

        # reset demon energy distribution
        self.d_energy = self.N//2 - self.E_lattice

        result = np.zeros(self.N, dtype=int)
        for i in range(self.d_energy):
            result[random.randint(0, self.N-1)] += 1

        self.E_demon = np.array(result)
        self.E_total = self.E_lattice + np.sum(self.E_demon)

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

    def _choose_rev_pair(self):
        """
        Deterministic reversible site and two demon indices for the CURRENT step.
        Returns (a, b1, b2) where b1 is used for spin flip and b2 for bond change.
        Do not advance indices here.
        """
        a  = self.order[self.order_idx]
        b1 = self.d_order[self.r_idx % (self.radius * 2 + 1)][self.order_idx]
        b2 = self.d_order[(self.r_idx + 1) % (self.radius * 2 + 1)][self.order_idx]
        return a, b1, b2


    def _choose_irr_pair(self):
        """
        Random irreversible site and two independently sampled demon indices.
        Returns (a, b1, b2) where b1 is used for spin flip and b2 for bond change.
        """
        rand_idx = np.random.randint(0, self.N)
        a = self.order[rand_idx]

        if self.radius == 0:
            b1 = a
            b2 = a
        else:
            local_idx1 = np.random.randint(0, self.radius * 2 + 1)
            local_idx2 = np.random.randint(0, self.radius * 2 + 1)
            b1 = self.d_order[local_idx1][rand_idx]
            b2 = self.d_order[local_idx2][rand_idx]

        return a, b1, b2


    def _advance_rev_forward(self):
        self.order_idx = (self.order_idx + 1) % self.N
        self.r_idx = (self.r_idx + 2) % (self.radius * 2 + 1)


    def _advance_rev_backward(self):
        self.order_idx = (self.order_idx - 1) % self.N
        self.r_idx = (self.r_idx - 2) % (self.radius * 2 + 1)


    def demon_move(self, flag, sweep_count):
        """
        One forward step.

        Reversible:
        - choose one deterministic (site, demon) pair
        - apply spin/bond to that SAME pair
        - then advance indices once

        Irreversible:
        - choose one random (site, demon) pair
        - randomly choose spin-first or bond-first
        - use SAME pair for both sub-attempts
        """
        if flag == 0:
            # Reversible: b1 for spin flip, b2 for bond change
            a, b1, b2 = self._choose_rev_pair()

            # Keep forward ordering fixed for reversibility
            self.E_lattice, self.d_energy = spin_flip_jit(
                self.lattice, self.bonds, self.E_demon,
                self.E_lattice, self.d_energy, a, b1, self.N
            )
            self.E_lattice, self.d_energy = bond_change_jit(
                self.lattice, self.bonds, self.E_demon,
                self.E_lattice, self.d_energy, a, b2, self.N
            )

            self._advance_rev_forward()

        else:
            # Irreversible: independently sampled b1 for spin, b2 for bond
            a, b1, b2 = self._choose_irr_pair()
            spin_first = np.random.randint(0, 2) == 0

            if spin_first:
                self.E_lattice, self.d_energy = spin_flip_jit(
                    self.lattice, self.bonds, self.E_demon,
                    self.E_lattice, self.d_energy, a, b1, self.N
                )
                self.E_lattice, self.d_energy = bond_change_jit(
                    self.lattice, self.bonds, self.E_demon,
                    self.E_lattice, self.d_energy, a, b2, self.N
                )
            else:
                self.E_lattice, self.d_energy = bond_change_jit(
                    self.lattice, self.bonds, self.E_demon,
                    self.E_lattice, self.d_energy, a, b2, self.N
                )
                self.E_lattice, self.d_energy = spin_flip_jit(
                    self.lattice, self.bonds, self.E_demon,
                    self.E_lattice, self.d_energy, a, b1, self.N
                )

        self.bond_count = count_bonds_jit(self.bonds)
        self.E_total = self.E_lattice + np.sum(self.E_demon)


    def demon_reverse(self, flag, sweep_count):
        """
        One reverse step.

        Reversible:
        - step indices backward once to recover the PREVIOUS forward pair
        - undo in reverse sub-order: bond then spin
        - use SAME pair for both sub-attempts

        Irreversible:
        - there is no microscopic reverse path, so treat as another
            stochastic step but with reversed sub-order bias if desired.
        """
        if flag == 0:
            # Move back to the pair used in the previous forward step
            self._advance_rev_backward()
            a, b1, b2 = self._choose_rev_pair()

            # Reverse sub-order of the forward move: bond (b2) then spin (b1)
            self.E_lattice, self.d_energy = bond_change_jit(
                self.lattice, self.bonds, self.E_demon,
                self.E_lattice, self.d_energy, a, b2, self.N
            )
            self.E_lattice, self.d_energy = spin_flip_jit(
                self.lattice, self.bonds, self.E_demon,
                self.E_lattice, self.d_energy, a, b1, self.N
            )

        else:
            # Irreversible "reverse half" is not truly invertible;
            # use the same one-pair-per-step rule with reversed sub-order bias.
            a, b1, b2 = self._choose_irr_pair()
            spin_first = np.random.randint(0, 2) == 0

            # Optional: bias toward opposite order from demon_move()
            if spin_first:
                self.E_lattice, self.d_energy = bond_change_jit(
                    self.lattice, self.bonds, self.E_demon,
                    self.E_lattice, self.d_energy, a, b2, self.N
                )
                self.E_lattice, self.d_energy = spin_flip_jit(
                    self.lattice, self.bonds, self.E_demon,
                    self.E_lattice, self.d_energy, a, b1, self.N
                )
            else:
                self.E_lattice, self.d_energy = spin_flip_jit(
                    self.lattice, self.bonds, self.E_demon,
                    self.E_lattice, self.d_energy, a, b1, self.N
                )
                self.E_lattice, self.d_energy = bond_change_jit(
                    self.lattice, self.bonds, self.E_demon,
                    self.E_lattice, self.d_energy, a, b2, self.N
                )

        self.bond_count = count_bonds_jit(self.bonds)
        self.E_total = self.E_lattice + np.sum(self.E_demon)