from project_1 import Molecule, a, b
import numpy as np

# Toluene
toluene = Molecule("Toluene", np.matrix([]), 7, 7, 3)
toluene.generate_H()
toluene.add_connections([[1, 6]])
toluene.set_constants(0, -1)
toluene.generate_eigen()
toluene.find_deloc_energy()
toluene.energy_level_plot()
toluene.find_charge_density()
toluene.find_bond_order()
print(toluene)
