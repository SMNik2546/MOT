from project_2 import Graphene
import numpy as np
# Armchair
armchair = Graphene('Armchair', np.matrix([]), 42, 42, 13)
armchair.generate_H()
armchair.add_connections([[1, 22], [4, 21], [5, 18], [8, 17], [9, 14], [15, 32], [16, 29], [19, 28],
                          [20, 25], [33, 35], [31, 37], [30, 38], [27, 39], [26, 40], [24, 42], [13, 34]])
