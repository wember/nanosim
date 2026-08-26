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
from numba import njit

# test #
from scipy.special import factorial as f
from scipy.special import loggamma as logg

Sk = lambda N, K: logg(K + N) - logg(K+1) - logg(N) # N == lattice size, K == kinetic energy
Su  = lambda N, N0, Nx: logg(N+1) + N0       * np.log(2) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1))
Su0 = lambda N, N0, Nx: logg(N+1) + (N0 + 1) * np.log(2) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1))

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

    # # Update neighbor bond if it exists
    # if bonds[(a-1)%N] != 0:
    #     if lattice[a] == lattice[(a-1)%N]:
    #         bonds[(a-1)%N] = -1
    #     else:
    #         bonds[(a-1)%N] = 1

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


@njit
def run_sweep_fwd_jit(lattice, bonds, E_demon, d_order, order, order_type,
                      E_lattice, d_energy, sweep_row, order_idx,
                      N, n_demon_rows, flag, d_energy_hist):
    """
    Run one complete forward sweep (N steps) inside a single JIT frame.

    Eliminates the per-step Python→Numba boundary crossing and moves
    count_bonds to once per sweep instead of once per step.
    """
    for _ in range(N):
        if flag == 0:                               # reversible
            a          = order[order_idx]
            row1       = sweep_row
            b1         = d_order[row1][order_idx]
            b2         = b1                          # same demon for both moves
            spin_first = (order_type[order_idx] == 0)
        else:                                       # irreversible
            rand_idx   = np.random.randint(0, N)
            a          = order[rand_idx]
            loc1       = np.random.randint(0, n_demon_rows)
            b1         = d_order[loc1][rand_idx]
            b2         = b1                          # same demon for both moves
            spin_first = np.random.randint(0, 2) == 0

        if spin_first:
            E_lattice, d_energy = spin_flip_jit(
                lattice, bonds, E_demon, E_lattice, d_energy, a, b1, N)
            E_lattice, d_energy = bond_change_jit(
                lattice, bonds, E_demon, E_lattice, d_energy, a, b2, N)
        else:
            E_lattice, d_energy = bond_change_jit(
                lattice, bonds, E_demon, E_lattice, d_energy, a, b2, N)
            E_lattice, d_energy = spin_flip_jit(
                lattice, bonds, E_demon, E_lattice, d_energy, a, b1, N)

        d_energy_hist[E_demon[b1]] += 1
        d_energy_hist[E_demon[b2]] += 1

        if flag == 0:
            order_idx = (order_idx + 1) % N
            if order_idx == 0:
                sweep_row = (sweep_row + 1) % n_demon_rows

    bond_count = count_bonds_jit(bonds)
    E_total    = E_lattice + np.sum(E_demon)
    return E_lattice, d_energy, sweep_row, order_idx, bond_count, E_total


@njit
def run_sweep_rev_jit(lattice, bonds, E_demon, d_order, order, order_type,
                      E_lattice, d_energy, sweep_row, order_idx,
                      N, n_demon_rows, flag, d_energy_hist):
    """
    Run one complete reverse sweep (N steps) inside a single JIT frame.
    """
    for _ in range(N):
        if flag == 0:                               # reversible: undo one forward step
            if order_idx == 0:
                sweep_row = (sweep_row - 1) % n_demon_rows
            order_idx  = (order_idx - 1) % N

            a          = order[order_idx]
            row1       = sweep_row
            b1         = d_order[row1][order_idx]
            b2         = b1                          # same demon for both moves
            spin_first = (order_type[order_idx] == 0)

            if spin_first:                          # undo in opposite sub-order
                E_lattice, d_energy = bond_change_jit(
                    lattice, bonds, E_demon, E_lattice, d_energy, a, b2, N)
                E_lattice, d_energy = spin_flip_jit(
                    lattice, bonds, E_demon, E_lattice, d_energy, a, b1, N)
            else:
                E_lattice, d_energy = spin_flip_jit(
                    lattice, bonds, E_demon, E_lattice, d_energy, a, b1, N)
                E_lattice, d_energy = bond_change_jit(
                    lattice, bonds, E_demon, E_lattice, d_energy, a, b2, N)
        else:                                       # irreversible: another random step
            rand_idx   = np.random.randint(0, N)
            a          = order[rand_idx]
            loc1       = np.random.randint(0, n_demon_rows)
            b1         = d_order[loc1][rand_idx]
            b2         = b1                          # same demon for both moves
            spin_first = np.random.randint(0, 2) == 0

            if spin_first:
                E_lattice, d_energy = spin_flip_jit(
                    lattice, bonds, E_demon, E_lattice, d_energy, a, b1, N)
                E_lattice, d_energy = bond_change_jit(
                    lattice, bonds, E_demon, E_lattice, d_energy, a, b2, N)
            else:
                E_lattice, d_energy = bond_change_jit(
                    lattice, bonds, E_demon, E_lattice, d_energy, a, b2, N)
                E_lattice, d_energy = spin_flip_jit(
                    lattice, bonds, E_demon, E_lattice, d_energy, a, b1, N)

        d_energy_hist[E_demon[b1]] += 1
        d_energy_hist[E_demon[b2]] += 1

    bond_count = count_bonds_jit(bonds)
    E_total    = E_lattice + np.sum(E_demon)
    return E_lattice, d_energy, sweep_row, order_idx, bond_count, E_total


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
        self.order_idx = 0  # Index pointer for order array
        self.r_idx = 0
        self.radius = R
        self.n_demon_rows = self.radius + 1
        offsets = np.array(list(range(0, R+1))).reshape(-1, 1)
        self.d_order = (a + offsets)  % N
        row_perm = np.random.permutation(self.n_demon_rows)
        self.d_order = self.d_order[row_perm]
        self.order_type = np.random.randint(0, 2, size=self.N) # 0 = spin-first, 1 = bond-first

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

        self.total_energy = total_energy
        self.d_energy_hist = np.zeros(self.N + self.total_energy + 1, dtype=np.int64)
        self.sweep_row = 0  # Which d_order row is active; increments each full sweep

################################################################################
################################################################################
################################################################################
    def reset(self):
        a = np.arange(self.N)
        np.random.shuffle(a)
        self.order = a
        # Rebuild d_order to match new order
        offsets = np.array(list(range(0, self.radius+1))).reshape(-1, 1)
        self.d_order = (a + offsets) % self.N
        self.order_idx = 0  # Also reset index pointers
        self.r_idx = 0
        self.sweep_row = 0
        self.order_type = np.random.randint(0, 2, size=self.N)

        # reset demon energy distribution
        result = np.zeros(self.N, dtype=int)
        for i in range(self.d_energy):
            result[random.randint(0, self.N-1)] += 1

        row_perm = np.random.permutation(self.n_demon_rows)
        self.d_order = self.d_order[row_perm]

        self.E_demon = np.array(result)
        self.E_total = self.E_lattice + np.sum(self.E_demon)

        self.d_energy_hist[:] = 0
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

    def do_sweep(self, flag):
        """
        Run one full forward sweep (N steps) in a single JIT call.
        count_bonds is called once at the end of the sweep, not per step.
        """
        (self.E_lattice, self.d_energy, self.sweep_row, self.order_idx,
         self.bond_count, self.E_total) = run_sweep_fwd_jit(
            self.lattice, self.bonds, self.E_demon,
            self.d_order, self.order, self.order_type,
            self.E_lattice, self.d_energy, self.sweep_row, self.order_idx,
            self.N, self.n_demon_rows, flag, self.d_energy_hist
        )

    def do_sweep_reverse(self, flag):
        """
        Run one full reverse sweep (N steps) in a single JIT call.
        count_bonds is called once at the end of the sweep, not per step.
        """
        (self.E_lattice, self.d_energy, self.sweep_row, self.order_idx,
         self.bond_count, self.E_total) = run_sweep_rev_jit(
            self.lattice, self.bonds, self.E_demon,
            self.d_order, self.order, self.order_type,
            self.E_lattice, self.d_energy, self.sweep_row, self.order_idx,
            self.N, self.n_demon_rows, flag, self.d_energy_hist
        )

    def _choose_rev_pair(self):
        a = self.order[self.order_idx]
        row1 = self.sweep_row
        b1 = self.d_order[row1][self.order_idx]
        b2 = b1  # same demon services both the spin flip and bond change
        return a, b1, b2


    def _choose_irr_pair(self):
        """
        Random irreversible site and a single sampled demon index that
        services both the spin flip and bond change attempts this step.
        Returns (a, b1, b2) with b1 == b2.
        """
        rand_idx = np.random.randint(0, self.N)
        a = self.order[rand_idx]

        local_idx1 = np.random.randint(0, self.radius + 1)
        b1 = self.d_order[local_idx1][rand_idx]
        b2 = b1  # same demon services both moves

        return a, b1, b2

    def _advance_rev_forward(self):
        self.order_idx = (self.order_idx + 1) % self.N
        if self.order_idx == 0:  # completed a full sweep
            self.sweep_row = (self.sweep_row + 1) % self.n_demon_rows

    def _advance_rev_backward(self):
        if self.order_idx == 0:  # about to cross a sweep boundary going back
            self.sweep_row = (self.sweep_row - 1) % self.n_demon_rows
        self.order_idx = (self.order_idx - 1) % self.N


    def demon_move(self, flag):
        """
        One forward step.

        Reversible:
        - choose one deterministic (site, demon) pair; the same demon
          services both the spin flip and bond change attempts
        - use the stored order_type at this step:
            0 -> spin then bond
            1 -> bond then spin
        - then advance indices once

        Irreversible:
        - choose one random (site, demon) pair; the same demon
          services both the spin flip and bond change attempts
        - choose spin-first or bond-first randomly each step
        """
        if flag == 0:
            a, b1, b2 = self._choose_rev_pair()

            # Fixed reversible substep order for this step
            spin_first = (self.order_type[self.order_idx] == 0)

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

            self._advance_rev_forward()

        else:
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

        self.d_energy_hist[self.E_demon[b1]] += 1
        self.d_energy_hist[self.E_demon[b2]] += 1


    def demon_reverse(self, flag):
        """
        One reverse step.

        Reversible:
        - move indices backward first to recover the previous forward step
        - read the SAME stored order_type used on that forward step
        - undo in the opposite sub-order

        Irreversible:
        - not microscopically invertible, so just use another random step
        with random substep order
        """
        if flag == 0:
            # Go back to the exact reversible step we want to undo
            self._advance_rev_backward()
            a, b1, b2 = self._choose_rev_pair()

            # This is the order that was used during the forward step
            spin_first = (self.order_type[self.order_idx] == 0)

            # Undo in the opposite order
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

        else:
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

        self.d_energy_hist[self.E_demon[b1]] += 1
        self.d_energy_hist[self.E_demon[b2]] += 1