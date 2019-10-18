from project_2 import Graphene
import numpy as np
# Zigzag
zigzag = Graphene('ZigZag', np.matrix([]), 48, 48, 18)
zigzag.generate_H()
zigzag.add_connections([[1, 22], [4, 21], [5, 18], [8, 17], [9, 14], [15, 32], [16, 29], [19, 28],
                        [20, 25], [33, 35], [31, 37], [30, 38], [27, 39], [26, 40], [24, 42], [13, 34]])
zigzag.delete_connections([[34, 35], [37, 38], [39, 40], [42, 43], [43, 44], [44, 45], [45, 46], [46, 47], [47, 48]])
zigzag.add_connections([[41, 43], [40, 44], [39, 45], [38, 46], [37, 47], [36, 48], [43, 2], [44, 3], [45, 6], [46, 7],
                        [47, 10], [48, 11]])
zigzag.set_constants(0, 1)
zigzag.generate_eigen()
zigzag.generate_carbons(3, 3)
zigzag.zig_zipper(3, 3)
zigzag.find_charge_density()
zigzag.find_deloc_energy()
print(zigzag)
