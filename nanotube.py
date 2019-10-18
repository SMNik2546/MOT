from Armchair import armchair
# Wrap this up in nanotube
armchair.add_connections([[1, 12], [23, 34], [35, 42]])

armchair.delete_connections([[34, 35], [37, 38], [39, 40]])
armchair.set_constants(0, 1)
armchair.generate_eigen()
armchair.find_charge_density()
armchair.find_deloc_energy()
armchair.generate_carbons(3, 3)
print(armchair)
