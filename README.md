# softwar-project
# Libraries:(numpy , sympy , matplotlib.pyplot)
# Project subject: Molecular Orbital Theory (SMNik2546)
# In each Project's file (#) lines explain the category of each step & explain each level of project.(BLUE=codes structure, PINK=category)

" Preface: Part A - Huckel(H) theory for pi-molecular orbitals.
For a given matrix H, with Alpha diagonals and Beta off-diagonals, we will determine the Eigenvalues and Eigenvectors
for the Huckel Effective Hamiltonian and use them to create an energy level diagram for the electronic configuration of the molecule."
" In this case, we want to creat the class that represent a molecule with a Huckel Matrix and associated methods"
" During the code, after defining Eigenvalue and Eigenvector, the string representation for the molecule is described.(Charge Density, deloc Energy & Bond Order)"
" Then, definition numerical matrix elements is required by codes for linear carbon chain."
" So, denoting spin for electrons from the energy levels clarifies the plot of each bonds. Besids, we find the bond order of molecule."


" Now, our class is ready, afterward, in Project two, we extend our molecule class to be able to calculate the various properties of Graphene."
" Hexagonal structure of carbon is used for generating the chain."
" creates row of linking carbons to wrap zigzag"
" Finds the x, y coord of a specific carbon atom"
" generates the carbon atoms for a graphene molecule of type 'self.name' and of block size mxn"
" Now, the structure of each project is ready for calculating different carbon band structure in many bulk. But, befor that we need consider Armchair as a Graphene molecule matrix calls from Project2. Zigzag matrix as a numerical calculation of neighbor connection in Graphen molecule for charge density and deloc energy and plot them. Nanotube works as a carbon tube for molecule to prepaire the graphs."

" Different molecules are Benzene, BuckBall, Butadiene, Napthalene and Toluene. Code of each one involves generate_eigen, find_deloc_energy, energy_level_plot, find_charge_density and find_bond_order that call from Project1." 




