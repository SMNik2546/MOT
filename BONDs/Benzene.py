import numpy as np
from project_1 import Molecule, a, b
# Benzene
benzene = Molecule("Benzene", np.matrix([]), 6, 6, 3)
benzene.generate_H()
benzene.add_connections([[1, 6]])
benzene.set_constants(0, -1)
benzene.generate_eigen()
benzene.find_deloc_energy()
benzene.energy_level_plot()
benzene.find_charge_density()
benzene.find_bond_order()
print(benzene)
