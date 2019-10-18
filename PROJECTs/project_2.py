"""PART 2
Preface: We extend our Molecule class to be able to calculate the various properties of types of Graphene
"""

from project_1 import Molecule, a, b
import numpy as np
import matplotlib.pyplot as plt
import copy


class Carbon:
    def __init__(self, x, y, magnitudes, id):
        self.pos = (x, y)
        self.psi_magnitudes = magnitudes
        self.id = id        # Which carbon atom it is?


class Graphene(Molecule):
    def __init__(self, *args):
        super(Graphene, self).__init__(*args)
        self.carbons = [] # List of Carbon Atoms
        self.prev_id = 0 # Stores the id of the most recently generated carbon
        self.current_row = 1 # Stores the row of hexagons that the most recently generated carbon belongs to

    def generate_carbon(self, id, m, n, position_set):  # m is number of hexagon rows, n is number of hexagon columns
        """Mint a new Carbon atom"""
        mags = [(abs(eig[id - 1][0]) ** 2) for eig in self.eigenvectors]
        position = position_set
        if position_set == 'find':
            position = self.find_position(id, m, n)
        return Carbon(position[0], position[1], mags, id)

    def zig_zipper(self, m, n):
        """creates row of linking carbons to wrap zigzag"""
        a = 0.5 # Short leg of 30-60-90 triangle
        b = float((3**(1/2.0))/2.0) # long leg of 30-60-90 triangle
        carb_iter = self.num_carbons - (2*n + 1)
        pull_tab = copy.copy(self.carbons[carb_iter])
        self.carbons[carb_iter + 1].pos = (pull_tab.pos[0], pull_tab.pos[1]-2*b)
        t = copy.copy(self.carbons[carb_iter + 1])
        if m % 2 == 0:
            dist = 2
            for i in range(2*n-1):
                self.carbons[(carb_iter+2+i)].pos = (t.pos[0]-dist, t.pos[1])
                t = copy.copy(self.carbons[(carb_iter+2+i)])
                if dist == 2:
                    dist = 1
                else:
                    dist = 2

        else:
            dist = 2
            for i in range(2*n-1):
                self.carbons[(carb_iter+2+i)].pos = (t.pos[0]+dist, t.pos[1])
                t = copy.copy(self.carbons[(carb_iter+2+i)])
                if dist == 2:
                    dist = 1
                else:
                    dist = 2

    def find_position(self, id, m, n):
        """Finds the x, y coord of a specific carbon atom"""
        prev = self.prev_id
        xprev = 0
        yprev = 0
        if self.prev_id != 0:
            xprev = self.carbons[prev-1].pos[0]
            yprev = self.carbons[prev-1].pos[1]
        row = self.current_row
        a = 0.5 # Short leg of 30-60-90 triangle
        b = float((3**(1/2.0))/2.0) # long leg of 30-60-90 triangle

        if row == 1:
            if prev == 0:
                xpos = 0.0
                ypos = 0.0
            elif prev == 4 * n:
                xpos = xprev - a
                ypos = yprev - b  #(sqrt(3))/2
                self.current_row += 1
            else:
                if prev % 2 == 0:
                    xpos = xprev + 1
                    ypos = yprev
                else:
                    xpos = xprev + a
                    if id % 4 == 0:
                        ypos = yprev - b
                    else:
                        ypos = yprev + b
        elif row == 2:

            if prev == 4*n + (row-1)*(4*n - 1):
                xpos = xprev + a
                ypos = yprev - b
                self.current_row += 1
            else:
                if id % 2 == 0:
                    xpos = xprev - 1
                    ypos = yprev
                else:
                    xpos = xprev - a
                    if (id+row-1) % 4 == 0:
                        ypos = yprev - b
                    else:
                        ypos = yprev + b
        elif row == 3:

            if prev == 4*n + (row-1)*(4*n - 1):
                xpos = xprev
                ypos = yprev - 2*b
                self.current_row += 1
            elif prev == 4*n + (row-1)*(4*n -1) - 1:
                xpos = xprev + a
                ypos = yprev + b
            else:
                if prev % 2 == 0:
                    xpos = xprev + 1
                    ypos = yprev
                else:
                    xpos = xprev + a
                    if (id+row-1) % 4 == 0:
                        ypos = yprev - b
                    else:
                        ypos = yprev + b
        elif row == m + 1:

            if id == 4*n + 2*(4*n -1) + (row-3)*(2*n +2):
                if row % 2 == 0:
                    xpos = xprev - a
                    ypos = yprev + b
                else:
                    xpos = xprev + a
                    ypos = yprev + b
            elif prev == 4*n + 2*(4*n -1) + (row-4)*(4*n) + 1:
                ypos = yprev - b
                if row % 2 == 0:
                    xpos = xprev - a
                else:
                    xpos = xprev + a
            else:
                if row % 2 == 0:
                    if prev % 2 == 0:
                        xpos = xprev - 1
                        ypos = yprev
                    else:
                        xpos = xprev - 2
                        ypos = yprev
                else:
                    if prev % 2 == 0:
                        xpos = xprev + 1
                        ypos = yprev
                    else:
                        xpos = xprev + 2
                        ypos = yprev
        else:
            if prev == 4*n + 2*(4*n -1) + (row-3)*(4*n):
                xpos = xprev
                ypos = yprev - 2*b
                self.current_row += 1
            elif prev == 4*n + 2*(4*n -1) + (row-3)*(4*n) - 1:
                ypos = yprev + b
                if row % 2 == 0:
                    xpos = xprev - a
                else:
                    xpos = xprev + a
            elif prev == 4*n + 2*(4*n -1) + (row-4)*(4*n) + 1:
                ypos = yprev - b
                if row % 2 == 0:
                    xpos = xprev - a
                else:
                    xpos = xprev + a
            else:
                if prev % 2 == 0:
                    if row % 2 == 0:
                        xpos = xprev - 1
                        ypos = yprev
                    else:
                        xpos = xprev + 1
                        ypos = yprev
                else:
                    if row % 2 == 0:
                        xpos = xprev - a
                        if (prev - 1) % 4 == 0:
                            ypos = yprev - b
                        else:
                            ypos = yprev + b

        self.prev_id = id
        return (xpos, ypos)

    def generate_carbons(self, m, n):
        """generates the carbon atoms for a graphene molecule of type 'self.name' and of block size mxn"""
        for i in range(self.num_carbons):
            self.carbons.append(self.generate_carbon(i+1, m, n, 'find'))

    def plot_lattice(self, index):
        """plots the graphene lattice"""
        x_list = []
        y_list = []
        for c in self.carbons:
            x_list.append(c.pos[0])
            y_list.append(c.pos[1])
        plt.scatter(x_list, y_list)
        for i in range(self.num_carbons):
            for j in range(self.num_carbons):
                if self.H[i, j] == 1:
                    plt.plot([self.carbons[i].pos[0], self.carbons[j].pos[0]], [self.carbons[i].pos[1], self.carbons[j].pos[1]])

        for c in self.carbons:
            circ = plt.Circle((c.pos[0], c.pos[1]), c.psi_magnitudes[index] * 5)
            plt.gcf().gca().add_artist(circ)

        # Hide the axes in our plot
        plt.xticks([])  # Hide the x-axes
        plt.yticks([])
        plt.title("Molecule Structure for " + self.name + " with energy level " + str(index) )
        plt.show()
