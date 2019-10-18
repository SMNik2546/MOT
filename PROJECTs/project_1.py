# Project: Molecular Orbital Theory
# SMNik2546

import numpy as np
from numpy import linalg as lig
import sympy as smp
import matplotlib.pyplot as plt
from operator import itemgetter
import collections

"""
Preface:
Part A - Huckel(H) theory for pi-molecular orbitals.
For a given matrix H, with Alpha diagonals and Beta off-diagonals, we will determine the Eigenvalues and Eigenvectors
for the Huckel Effective Hamiltonian and use them to create an energy level diagram for the electronic configuration of the molecule.
"""

a, b = smp.symbols('a, b')  # Generate symbols so we can use them in our fancy matrix


class Molecule:
    """This class represents a molecule with a Huckel Matrix and associated methods"""

    def __init__(self, name, H, num_pi_electrons, num_carbons, num_double_bonds):
        self.H = H
        self.name = name
        self.num_carbons = num_carbons
        self.eigenvalues = []
        self.eigenvectors = []
        self.normalized_eigenvectors = []
        self.eigval_eigvect = []         # Associate Eigenvalue with its Eigenvect
        self.eigval_multiplicity = []
        self.num_pi_electrons = num_pi_electrons
        self.deloc_energy = 0.0
        self.alpha = None
        self.beta = None
        self.charge_density = []
        self.bond_order = []
        self.num_additional_connections = 0
        self.con = []
        self.num_double_bonds = num_double_bonds

        # This is our internal data structure that contains [(eig_value, # of Electrons)]
        self.e_per_energy_lvl = []
        self.e_per_eigen_vect = []

    def __str__(self):
        """Create the string representation for the molecule"""
        return "---" + self.name + "--- \n" + str(self.H) + "\n" \
               + "Charge Density :" + str(self.charge_density) + "\n" \
               + "Delocalization Energy :" + str(self.deloc_energy) + "\n" \
               + "Bond Order :" + str(self.bond_order) + "\n"

    def set_constants(self, al, be):
        """This function substitutes numeric value in place for Alpha and Beta in the Huckel Matrix"""

        # Replace our H matrix with numerical values
        self.H = np.matrix([[al if x == a else be if x == b else x for x in i] for i in self.H.tolist()])

        # Store the resonance integral values
        self.alpha = al
        self.beta = be

    def generate_H(self):
        """generates H matrix for linear carbon chain"""
        N = self.num_carbons
        H = np.zeros((N, N)).tolist()

        for i in range(N):
            H[i][i] = a
            if i == 0:
                H[i][i + 1] = b
            elif i == N - 1:
                H[i][i - 1] = b
            else:
                H[i][i + 1] = b
                H[i][i - 1] = b

        self.H = np.matrix(H)

    def add_connections(self, con):
        """adds con to the Huckel matrix, takes a list of lists and inserts beta into H at the
        specified coordinates
        con: [[coord1, coord2]...]
        """
        self.con = con
        for i in range(len(con)):
            self.H[con[i][0] - 1, con[i][1] - 1] = b
            self.H[con[i][1] - 1, con[i][0] - 1] = b
        self.num_additional_connections += len(con)

    def delete_connections(self, con):
        """deletes con to the Huckel matrix, takes a list of tuples and inserts zero into H at the
        specified coordinates
        con: [[coord1, coord2]...]
        """

        for i in range(len(con)):
            self.H[con[i][0] - 1, con[i][1] - 1] = 0
            self.H[con[i][1] - 1, con[i][0] - 1] = 0
        self.num_additional_connections -= len(con)

    def generate_eigen(self):
        """Finds the eigenvalue and eigenvector for the H matrix"""
        assert(self.alpha is not None and self.beta is not None)    # We are no longer using symbolic evals

        # Reset our arrays
        self.eigenvalues = []
        self.eigenvectors = []
        self.eigval_multiplicity = []
        self.eigval_eigvect = []

        e_vals, e_vects = lig.eig(self.H)                     # Generate using Numpy's eigenvals
        e_vals = np.around(e_vals, decimals=5)                    # Round these eigen values
        freq_dict = collections.Counter(e_vals)                   # Counts up # of Eigenvalue and their multiplicity

        # Store the eigenvalues along with their multiplicities. Note: Eigenvalues in this array is unique
        self.eigval_multiplicity = sorted(freq_dict.items(), key=lambda x: x[0])

        # Associate each eigenvalues with their eigenvectors. Note: Eigenvalues may repeat in this array
        for i in range(len(e_vals)):
            eigenvalue = e_vals[i]
            mega_tuple = tuple([eigenvalue, e_vects[:, i]])
            self.eigval_eigvect.append(mega_tuple)

        # Sort our eigval_eigvect list, in ascending eigenvalue order
        self.eigval_eigvect.sort(key=lambda x: x[0])

        # Store just the eigenvalues and the eigenvectors in the appropriate arrays
        for eig_set in self.eigval_eigvect:
            self.eigenvalues = eig_set[0]
            self.eigenvectors.append(eig_set[1].tolist())

        # Reset the array which will hold the number of electrons per eigenvalue
        self.e_per_energy_lvl = []

        # Compile our electron per energy level array
        electrons_available = self.num_pi_electrons
        for eig_val_multiplicity in self.eigval_multiplicity:
            eig_val = eig_val_multiplicity[0]
            multiplicity = eig_val_multiplicity[1]
            if electrons_available >= (multiplicity * 2):
                self.e_per_energy_lvl.append(tuple([eig_val, 2 * multiplicity]))
                electrons_available -= 2 * multiplicity
            else:
                self.e_per_energy_lvl.append(tuple([eig_val, electrons_available]))
                electrons_available -= electrons_available

        # Count up the number of electrons associated with each eigenvector by associating it with the eigenvalue
        self.e_per_eigen_vect = []
        e_per_energy_lvl_copy = [x for x in self.e_per_energy_lvl]      # Get a temporary copy to work with

        for eig_set in self.eigval_eigvect:
            eig_val = eig_set[0]            # Eigenvalue
            eig_vect = eig_set[1]           # Eigenvectors
            electrons = 0

            for electron_index in range(len(e_per_energy_lvl_copy)):
                if e_per_energy_lvl_copy[electron_index][0] == eig_val:
                    # We have a match of eigenvalues
                    if e_per_energy_lvl_copy[electron_index][1] > 2:
                        # pull out 2
                        electrons = 2
                        e_per_energy_lvl_copy[electron_index] = tuple([e_per_energy_lvl_copy[electron_index][0],
                                                                       e_per_energy_lvl_copy[electron_index][1] - 2])
                    else:
                        # Pull out as many as it has
                        electrons = e_per_energy_lvl_copy[electron_index][1]
                        e_per_energy_lvl_copy[electron_index] = tuple([e_per_energy_lvl_copy[electron_index][0], 0])

            self.e_per_eigen_vect.append(tuple([eig_vect, electrons]))

    def find_nodes(self):
        """Finds the nodes for the Eigenvector"""

        def find_nodes_helper(eigenvector):
            """ A helper function that count the number of nodes for a specific eigenvector"""

            nodes = 0  # Number of Nodes
            for i in range(len(eigenvector) - 1):
                if eigenvector[i] * eigenvector[i + 1] < 0:
                    nodes += 1  # We have a node
            return nodes

        result = []
        for eigvector in self.eigenvectors:
            result.append(find_nodes_helper(eigvector))

        # Result = [num_nodes, ...]
        return result

    def energy_level_plot(self):
        """Plots the energy levels and denote spin for the electrons"""
        assert(self.alpha is not None and self.beta is not None)  # Make sure alpha and beta are defined

        # Get the unicode arrow
        down_arrow = u'$\u2193$'
        up_arrow = u'$\u2191$'

        electrons_used = self.num_pi_electrons
        max_multiplicity = max(self.eigval_multiplicity, key=itemgetter(1))[1]      # Count the maximum multiplicity

        # Draw the graphs by iterating through each eigenvalue and its associated multiplicity
        for eig in self.eigval_multiplicity:
            eig_val = eig[0]
            if eig[1] == 1:
                plt.axhline(eig[0])  # Draw the eigenvalues as lines on the graph
            else:
                for i in range(eig[1]):
                    plt.plot([i, i + 0.95], 2 * [eig_val],
                             color=np.random.rand(3, ))

            # Fill up to two electrons per level per line, from bottom up
            if eig[1] == 1:
                if electrons_used > 1:
                    plt.plot(0.2 * max_multiplicity, eig_val, linestyle='none', marker=up_arrow, markersize=15)
                    plt.plot(0.8 * max_multiplicity, eig_val, linestyle='none', marker=down_arrow, markersize=15)
                    electrons_used -= 2
                elif electrons_used == 1:
                    plt.plot(0.2 * max_multiplicity, eig_val, linestyle='none', marker=up_arrow, markersize=15)
                    electrons_used -= 1

                else:
                    pass
            else:
                for i in range(eig[1]):
                    # Add all the up arrows
                    if electrons_used >= 1:
                        plt.plot(0.2 + i, eig_val, linestyle='none', marker=up_arrow, markersize=15)
                        electrons_used -= 1
                for i in range(eig[1]):
                    # Add all the up arrows
                    if electrons_used >= 1:
                        plt.plot(0.8 + i, eig_val, linestyle='none', marker=down_arrow, markersize=15)
                        electrons_used -= 1

        # Draw the graphs
        plt.title('Energy Level Plot for ' + str(self.name))
        plt.xlim(0, max_multiplicity)  # Format the Graph
        plt.xticks([])  # Hide the x-axes
        plt.ylabel('Energy')
        plt.show()

    def find_deloc_energy(self):
        """Calculates the delocalization energy by going through our augmented array and adding up all the energies"""

        deloc_energy = 0.0
        for e, num_electrons in self.e_per_energy_lvl:
            deloc_energy += e * num_electrons
        # We calculated the total Pi Electron Energy

        # Subtract Pi Electrons in Isolated Bonds
        deloc_energy -= (2 * a + 2 * b) * self.num_double_bonds

        if self.alpha is not None and self.beta is not None:
            # If we have the alpha and beta constants
            deloc_energy = deloc_energy.subs(a, self.alpha).subs(b, self.beta)

        self.deloc_energy = smp.N(deloc_energy)     # Set it nicely

        return self.deloc_energy

    def find_charge_density(self):
        """finds the charge density of Pi electrons for each carbon atom in the molecule"""
        charge_density = []

        for c in range(self.num_carbons):                           # For each carbon atom
            charge_sum = 0.0
            for eig_e in self.e_per_eigen_vect:
                """Iterate through each eigenvector and its number of electrons and add up the charge sum
                according to that formula """
                eig_vect_val = (eig_e[0][c].tolist())[0][0]
                charge_sum += eig_e[1] * (abs(eig_vect_val) ** 2)

            charge_density.append(charge_sum)

        self.charge_density = charge_density        # list of charge densities for carbon atoms 1-n

    def find_bond_order(self):
        """finds the bond order of molecule, stores as list of values beginning at C1"""

        bond_order = []

        for c in range(self.num_carbons - 1 + self.num_additional_connections):
            # For each carbon atom - 1 plus extra bonds between carbons

            bond_sum = 1.0

            # Iterate through our eigenvector and their associated number of electrons
            for eigv_e in self.e_per_eigen_vect:
                num_elec = eigv_e[1]

                if c < self.num_carbons - 1:
                    eigvector_at_c = (eigv_e[0][c].tolist())[0][0]
                    eigvector_at_c1= (eigv_e[0][c + 1].tolist())[0][0]
                    bond_sum += num_elec * eigvector_at_c * eigvector_at_c1
                else:
                    first_eigv_index = self.con[c % (self.num_carbons-1)][0]-1
                    eigv_1 = (eigv_e[0][first_eigv_index].tolist())[0][0]
                    second_eigv_index = self.con[c % (self.num_carbons-1)][1] - 1
                    eigv_2 = (eigv_e[0][second_eigv_index].tolist())[0][0]
                    bond_sum += num_elec * eigv_1 * eigv_2

            bond_order.append(bond_sum)

        self.bond_order = bond_order            # list of charge densities for carbon atoms 1-n
