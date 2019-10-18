import numpy as np
from project_1 import Molecule, a, b
# BuckyBall
bucky = Molecule("Buckminsterfullerene", np.matrix([]), 60, 60, 30)
bucky.generate_H()
bucky.add_connections([[1, 5], [1, 9], [2, 12], [3, 15], [4, 18], [6, 20], [7, 22], [8, 25], [10, 26], [11, 29],
                       [13, 30], [14, 33], [16, 34], [17, 37], [19, 38], [21, 40], [23, 42], [24, 44], [27, 45],
                       [28, 47], [31, 48], [32, 50], [35, 51], [36, 53], [39, 54], [41, 55], [43, 57], [46, 58],
                       [49, 59], [52, 60], [56, 60]])
bucky.set_constants(0, -1)
bucky.generate_eigen()
bucky.find_deloc_energy()
bucky.energy_level_plot()
bucky.find_charge_density()
bucky.find_bond_order()
print(bucky)
