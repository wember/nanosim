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

        total_energy = 2*N # total energy of the system

        self.N = N
        self.order = a
        self.rev_order = np.flip(a)
        self.radius_spin = self.rev_radius_bond = np.random.randint(0, R, size=N)*np.random.choice([-1, 1], size=N)
        self.rev_radius_spin = self.radius_bond = np.flip(self.radius_spin)
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
                ###            TEST SETUP           ###
################################################################################
        # a = np.arange(N)
        # total_energy = N # total energy of the system
        # self.N = N
        # self.order = a
        # self.rev_order = np.flip(a)
        # self.lattice = np.array([1,1])
        # self.bonds = [-1,-1]
        # self.bond_count = np.ones(3, dtype=int)
        # self.count_bonds()
        # self.E_lattice = sum(self.bonds)
        # self.E_demon = np.array([4,0])
        # self.d_energy = total_energy - self.E_lattice
        # self.E_total = self.E_lattice + sum(self.E_demon)

################################################################################
################################################################################
################################################################################

    # def calc_E_lat(self,lattice,N):
    #     """
    #         Calculate energy of the lattice configuration.
    #     """
    #     ETOT = 0
    #
    #     # Loop over the entire lattice calculating nearest neighbor interactions
    #     for a in range(N):
    #         # Grab the lattice site spin value
    #         s =  lattice[a]
    #         # Calculate the energy of the configuration based on
    #         # nearest neighbors
    #         nb = lattice[(a+1)%N] + lattice[(a-1)%N]
    #         # running sum of energy of Ising latus
    #         ETOT += 2 * s * nb
    #
    #     # Update the value of the lattice energy
    #     return ETOT

    def spin_flip(self, a, i):
        # print(self.lattice, "lattice: ", self.E_lattice,  "demon: ", self.E_demon,  self.d_energy, "total: ", self.E_lattice+sum(self.E_demon))
        # print(self.bonds, self.bond_count)
        # ### entropy calc for testing
        # self.count_bonds()
        #
        # print(" N0 | Nx | U/J | Su/k  | K/J | Sk/k  | Se/k | p(u,k)")
        # if (self.bond_count[1] == 0):
        #     uk = (f(self.N)*2**(self.bond_count[1]+1))/(f(self.N-self.bond_count[1]-self.bond_count[2])*f(self.bond_count[1])*f(self.bond_count[2]))
        #     SUtest = Su0(self.N, self.bond_count[1], self.bond_count[2])
        # else:
        #     uk = (f(self.N)*2**(self.bond_count[1]))/(f(self.N-self.bond_count[1]-self.bond_count[2])*f(self.bond_count[1])*f(self.bond_count[2]))
        #     SUtest = Su(self.N, self.bond_count[1], self.bond_count[2])
        # sk = f(self.d_energy+self.N-1)/(f(self.d_energy)*f(self.N-1))
        # # print(f" {self.bond_count[1]}  |  {self.bond_count[2]} | {self.E_lattice}  |ln({uk})|  {self.d_energy}  |ln({sk})|ln({sk*uk})")
        #
        # print("Su : ", SUtest, math.log(uk))
        # print("Sk : ", Sk(self.N, self.d_energy), math.log(sk))
        # print("----------------------")
        """
            Attempt to flip the spin of a given lattice site
        """
        # Grab the lattice site spin value
        s =  self.lattice[a]
        d = self.E_demon[i]
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
        # Otherwise, pass
        else:
            pass

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
        # print(self.lattice, "lattice: ", self.E_lattice, "demon: ", self.E_demon,  self.d_energy, "total: ", self.E_lattice+sum(self.E_demon))
        # print(self.bonds, self.bond_count)
        # ### entropy calc for testing
        # self.count_bonds()
        #
        # print(" N0 | Nx | U/J | Su/k  | K/J | Sk/k  | Se/k | p(u,k)")
        # if (self.bond_count[1] == 0):
        #     uk = (f(self.N)*2**(self.bond_count[1]+1))/(f(self.N-self.bond_count[1]-self.bond_count[2])*f(self.bond_count[1])*f(self.bond_count[2]))
        #     SUtest = Su0(self.N, self.bond_count[1], self.bond_count[2])
        # else:
        #     uk = (f(self.N)*2**(self.bond_count[1]))/(f(self.N-self.bond_count[1]-self.bond_count[2])*f(self.bond_count[1])*f(self.bond_count[2]))
        #     SUtest = Su(self.N, self.bond_count[1], self.bond_count[2])
        # sk = f(self.d_energy+self.N-1)/(f(self.d_energy)*f(self.N-1))
        # # print(f" {self.bond_count[1]}  |  {self.bond_count[2]} | {self.E_lattice}  |ln({uk})|  {self.d_energy}  |ln({sk})|ln({sk*uk})")
        #
        # print("Su : ", SUtest, math.log(uk))
        # print("Sk : ", Sk(self.N, self.d_energy), math.log(sk))
        # print("----------------------")
        """
            Attempt to change the bond given lattice site
        """
        # Grab the lattice site spin, bond value, and demon energy
        s =  self.lattice[a]
        b =  self.bonds[a]
        d = self.E_demon[i]
        # Grab value of bonded neighbor
        n = self.lattice[(a+1)%self.N]
        # Check the cost of breaking the bond
        if (s == n):
            cost = -1
        else:
            cost = 1

        # if bond is broken, attempt to remake
        if (b == 0) and (d - cost >= 0):
            b = cost

            if (self.bonds[a] == 0):
                self.E_lattice += cost
                self.E_demon[i] -= cost
                self.d_energy -= cost
                self.bonds[a] = b

        # if bond is made, attempt to break
        elif (d + cost >= 0):
            b = 0

            if (self.bonds[a] != 0):
                self.E_lattice -= cost
                self.E_demon[i] += cost
                self.d_energy += cost
                self.bonds[a] = b

        else:
            pass

        if (self.bonds[(a-1)%self.N] != 0):
            if (self.lattice[a] == self.lattice[(a-1)%self.N]):
                self.bonds[(a-1)%self.N] = -1
            else:
                self.bonds[(a-1)%self.N] = 1

    def count_bonds(self):
        """
            Updates the bond-count array of number of aligned (-1), broken (0), and misaligned (-1) bonds
        """
        unique, counts = np.unique(self.bonds, return_counts=True)
        list = dict(zip(unique.astype(str), counts.astype(str)))
        index = 0
        for i in [-1,0,1]:
            bond_type = str(i)
            if i in unique:
                self.bond_count[index] = list[bond_type]
            else:
                self.bond_count[index] = 0
            index += 1


    def demon_move(self):
        """
            "Randomly" move the demon around and flip spins & change bonds
        """
        a = self.order[0]

        # Attempt to flip spin
        self.spin_flip(a, (a + self.radius_spin[0])%self.N)
        self.radius_spin = np.roll(self.radius_spin, -1)

        # Attempt to change bond
        self.bond_change(a, (a + self.radius_bond[0])%self.N)
        self.radius_bond = np.roll(self.radius_bond, -1)

        # Update bond count
        self.count_bonds()

        # Move first element in order to back
        self.order = np.roll(self.order, -1)

    def demon_reverse(self):
        """
            In reverse order, flip spins & change bonds
        """
        a = self.rev_order[0]

        # Attempt to change bond
        self.bond_change(a, (a + self.rev_radius_bond[0])%self.N)
        self.rev_radius_bond = np.roll(self.rev_radius_bond, -1)

        # Attempt to flip spin
        self.spin_flip(a, (a + self.rev_radius_spin[0])%self.N)
        self.rev_radius_spin = np.roll(self.rev_radius_spin, -1)

        # Update bond count
        self.count_bonds()

        # Move first element in order to back
        self.rev_order = np.roll(self.rev_order, -1)
