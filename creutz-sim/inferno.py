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

        # Reversible demon scheduling parameters.
        # Keep the map deterministic/invertible, but avoid simple lock-step
        # resonance between order_idx and demon-row order.
        self.r_step = self._coprime_step(self.n_demon_rows, prefer=2)
        self.r_sweep_step = self._coprime_step(self.n_demon_rows, prefer=max(2, self.n_demon_rows // 2))
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
        self.d_energy_hist = np.zeros(self.total_energy + 1, dtype=np.int64)

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

    def _coprime_step(self, modulus, prefer=2):
        """
        Choose a deterministic step size that is coprime to `modulus`.
        Prefer a value larger than 1 when possible so the reversible demon-row
        schedule does not simply march 0, 1, 2, ... in lock-step with order_idx.
        """
        if modulus <= 1:
            return 0

        start = prefer % modulus
        if start == 0:
            start = 1

        for delta in range(modulus):
            cand = (start + delta) % modulus
            if cand == 0:
                continue
            if math.gcd(cand, modulus) == 1:
                return cand

        return 1


    def _choose_rev_pair(self, sweep_count):
        a = self.order[self.order_idx]

        if self.n_demon_rows == 1:
            row1 = 0
        else:
            sweep_phase = (sweep_count * self.r_sweep_step) % self.n_demon_rows
            row1 = (sweep_phase + self.r_step * self.order_idx) % self.n_demon_rows

        row2 = (row1 + 1) % self.n_demon_rows
        b1 = self.d_order[row1][self.order_idx]
        b2 = self.d_order[row2][self.order_idx]
        return a, b1, b2


    def _choose_irr_pair(self):
        """
        Random irreversible site and two independently sampled demon indices.
        Returns (a, b1, b2) where b1 is used for spin flip and b2 for bond change.
        """
        rand_idx = np.random.randint(0, self.N)
        a = self.order[rand_idx]

        local_idx1 = np.random.randint(0, self.radius + 1)
        local_idx2 = np.random.randint(0, self.radius + 1)
        b1 = self.d_order[local_idx1][rand_idx]
        b2 = self.d_order[local_idx2][rand_idx]

        return a, b1, b2

    def _advance_rev_forward(self):
        self.order_idx = (self.order_idx + 1) % self.N

    def _advance_rev_backward(self):
        self.order_idx = (self.order_idx - 1) % self.N


    def demon_move(self, flag, sweep_count):
        """
        One forward step.

        Reversible:
        - choose one deterministic (site, spin_demon, bond_demon) triple
        - use the stored order_type at this step:
            0 -> spin then bond
            1 -> bond then spin
        - then advance indices once

        Irreversible:
        - choose one random (site, spin_demon, bond_demon) triple
        - choose spin-first or bond-first randomly each step
        """
        if flag == 0:
            a, b1, b2 = self._choose_rev_pair(sweep_count)

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
        self.bond_count = count_bonds_jit(self.bonds)
        self.E_total = self.E_lattice + np.sum(self.E_demon)


    def demon_reverse(self, flag, sweep_count):
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
            a, b1, b2 = self._choose_rev_pair(sweep_count)

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
        self.bond_count = count_bonds_jit(self.bonds)
        self.E_total = self.E_lattice + np.sum(self.E_demon)