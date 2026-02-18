import numpy as np
import random
import math


# test #
from scipy.special import factorial as f
from scipy.special import loggamma as logg
Sk = lambda N, K: logg(K + N) - logg(K+1) - logg(N) # N == lattice size, K == kinetic energy
Su = lambda N, N0, Nx: logg(N+1) + np.log(2**(N0)) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1)) # N == lattice size, N0 == broken bonds, Nx == bonds between anti-aligned spins
Su0 = lambda N, N0, Nx: logg(N+1) + np.log(2**(N0+1)) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1)) # N == lattice size, N0 == broken bonds, Nx == bonds between anti-aligned spins

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
        self.radius = R
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
            Attempt to flip the spin of a given lattice site
        """
        # Grab the lattice site spin value
        s = self.lattice[a]
        # Calculate the energy of the configuration based on
        # nearest neighbors
        nb = self.lattice[(a+1)%self.N]*abs(self.bonds[(a)%self.N]) + self.lattice[(a-1)%self.N]*abs(self.bonds[(a-1)%self.N])
        # Check the cost of flipping the spin
        cost = 2*s*nb
        # If energetically favorable, flip and add energy to demon
        if cost < 0:
            s *= -1
            # Notice we substract the cost to maintain net0 energy
            self.E_demon[i] -= cost
            self.d_energy -= cost
            self.E_lattice += cost
        # If it costs energy, only flip if demon has enough energy
        elif cost <= self.E_demon[i]:
            s *= -1
            self.E_demon[i] -= cost
            self.d_energy -= cost
            self.E_lattice += cost

        # Update spin
        self.lattice[a] = s

        # Update bond of lattice site and of leftmost neighbor
        if (self.bonds[a] != 0):
            if (self.lattice[a] == self.lattice[(a+1)%self.N]):
                self.bonds[a] = -1
            else:
                self.bonds[a] = 1


        if (self.bonds[(a-1)%self.N] != 0):
            if (self.lattice[a] == self.lattice[(a-1)%self.N]):
                self.bonds[(a-1)%self.N] = -1
            else:
                self.bonds[(a-1)%self.N] = 1

    def bond_change(self, a, i):
        """
            Attempt to change the bond given lattice site
        """
        # Grab the lattice site spin, bond value
        s = self.lattice[a]
        b = self.bonds[a]
        # Grab value of bonded neighbor
        n = self.lattice[(a+1)%self.N]
        # Check the cost of breaking the bond
        cost = -1 if s == n else 1

        # if bond is broken, attempt to remake
        if (b == 0) and (self.E_demon[i] >= cost):
            self.E_lattice += cost
            self.E_demon[i] -= cost
            self.d_energy -= cost
            self.bonds[a] = cost

        # if bond is made, attempt to break
        elif (b != 0) and (self.E_demon[i] + cost >= 0):
            self.E_lattice -= cost
            self.E_demon[i] += cost
            self.d_energy += cost
            self.bonds[a] = 0

        # Update neighbor bond if it exists
        if (self.bonds[(a-1)%self.N] != 0):
            if (self.lattice[a] == self.lattice[(a-1)%self.N]):
                self.bonds[(a-1)%self.N] = -1
            else:
                self.bonds[(a-1)%self.N] = 1

    def count_bonds(self):
        """
            Updates the bond-count array of number of aligned (-1), broken (0), and misaligned (1) bonds
        """
        # Simple loop is faster than np.unique and JIT-compatible
        self.bond_count[0] = 0  # count of -1 (aligned)
        self.bond_count[1] = 0  # count of 0 (broken)
        self.bond_count[2] = 0  # count of 1 (misaligned)
        
        for b in self.bonds:
            if b == -1:
                self.bond_count[0] += 1
            elif b == 0:
                self.bond_count[1] += 1
            elif b == 1:
                self.bond_count[2] += 1

    def demon_move(self, flag, sweep_count):
        """
            "Randomly" move the demon around and flip spins & change bonds
        """
        a = self.order[0]
        radius_cycle = 2 * self.radius + 1
        R = (self.R_counter % radius_cycle) - self.radius
        self.R_counter += 1
        # If irr flag is on, generate a random number instead
        if (flag != 0):
            a = np.random.randint(0, self.N)
            if self.radius != 0:
                R = np.random.randint(0, self.radius)

        # Attempt to flip spin
        self.spin_flip(a, (a + R)%self.N)

        R = (self.R_counter % radius_cycle) - self.radius
        self.R_counter += 1
        # If irr flag is on, generate a random number instead
        if (flag != 0):
            if self.radius != 0:
                R = np.random.randint(0, self.radius)
        # Attempt to change bond
        self.bond_change(a, (a + R)%self.N)

        # Update bond count
        self.count_bonds()

        # Move first element in order to back
        if (flag == 0):
            self.order = np.roll(self.order, -1)

    def demon_reverse(self, flag, sweep_count):
        """
            In reverse order, flip spins & change bonds
        """
        a = self.rev_order[0]
        radius_cycle = 2 * self.radius + 1 
        self.R_counter -= 1
        R = (self.R_counter % radius_cycle) - self.radius
        # If irr flag is on, generate a random number instead
        if (flag != 0):
            a = np.random.randint(0, self.N)
            if self.radius != 0:
                R = np.random.randint(0, self.radius)
        # Attempt to change bond
        self.bond_change(a, (a + R)%self.N)

        self.R_counter -= 1
        R = (self.R_counter % radius_cycle) - self.radius
        # If irr flag is on, generate a random number instead
        if (flag != 0):
            if self.radius != 0:
                R = np.random.randint(0, self.radius)
        # Attempt to flip spin
        self.spin_flip(a, (a + R)%self.N)

        # Update bond count
        self.count_bonds()

        # Move first element in order to back
        if (flag == 0):
            self.rev_order = np.roll(self.rev_order, -1)
