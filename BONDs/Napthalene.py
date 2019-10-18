import numpy as np
from project_1 import Molecule, a, b

# Napthalene
napthalene = Molecule("Napthalene", np.matrix([]), 10, 10, 5)
napthalene.generate_H()
napthalene.add_connections([[5, 10], [1, 6]])
napthalene.set_constants(0, -1)
napthalene.generate_eigen()
napthalene.find_deloc_energy()
napthalene.energy_level_plot()
napthalene.find_charge_density()
napthalene.find_bond_order()
print(napthalene)
