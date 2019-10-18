import numpy as np
from project_1 import Molecule, a, b
#Butadiene

butadiene_H = np.matrix([[a, b, 0, 0], [b, a, b, 0],
                               [0, b, a, b], [0, 0, b, a]])
butadiene = Molecule("Butadine", butadiene_H, 4, 4, 2)

butadiene.set_constants(0, -1)
butadiene.generate_eigen()
butadiene.find_deloc_energy()
butadiene.energy_level_plot()
butadiene.find_charge_density()
butadiene.find_bond_order()
print(butadiene)
