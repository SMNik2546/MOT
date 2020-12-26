import numpy as np
import unittest
import project1
#import Benzene

class Test(unittest.TestCase):
    def test_1(self):
        assert project1.Molecule
        assert project1.Molecule.find_bond_order
        assert project1.Molecule.find_charge_density
        assert project1.Molecule.set_constants
        assert project1.Molecule.generate_eigen
        assert project1.Molecule.energy_level_plot
        assert project1.Molecule.find_deloc_energy
        assert project1.Molecule.energy_level_plot
        assert project1.Molecule.add_connections
        assert project1.Molecule.generate_H
        assert project1.Molecule.find_nodes
        assert project1.Molecule.delete_connections
        assert project1.Molecule.__str__
        assert project1.Molecule.__init__
        
