
from ket import Ket2D
from geometric_operators import GeometricOperators

import numpy as np
import sympy as sp
import math

geo = GeometricOperators()
ket = Ket2D()

#first, compute the spin network standard basis, given by

v1_ket = ket.compute_ket((0, 1), (2, 3))
v2_ket = ket.compute_ket((0, 3), (1, 2))

v1_vec = ket.compute_vector([(0, 1), (2, 3)])
v2_vec = ket.compute_vector([(0, 3), (1, 2)])

print(f"v1 ket: {v1_ket}", f"   v1 vect: {v1_vec.real}")
print(f"v2 ket: {v2_ket}", f"   v2 vect: {v2_vec.real}")

# Let's compute now all the geometric quantities of this quantum tetrahedron in this standard basis

# Ciclo per i valori di i e j
for i in range(1, 4):
    for j in range(i, 4):  # j parte da i per evitare duplicati simmetrici
        D_ij = geo.diedral_matrix((i, j), [v1_vec, v2_vec])

        if i == j:
            print(f"\nArea of the {i}-upper face, D_{i}{j}:")
        else:
            print(f"\nMatrice Diedrale D_{i}{j}:")

        geo.print_matrix(D_ij)

print(f"\nArea of the down face, D_00:")
D_00 = geo.diedral_matrix((0, 0), [v1_vec, v2_vec])
geo.print_matrix(D_00)


print(f"\nVolume of the tetrahedron, Vol:")
Vol_std = geo.commutatore([v1_vec, v2_vec])
geo.print_matrix(np.abs(Vol_std))

#now we change towards the orth basis

orth_basis = [(1/2, 0), (math.sqrt(3)/6, math.sqrt(3)/3)]

for i in range(1, 4):
    for j in range(i, 4):  # j parte da i per evitare duplicati simmetrici
        D_ij = geo.diedral_matrix((i, j), [v1_vec, v2_vec])
        D_ij_orth = geo.change_basis(D_ij, orth_basis)

        if i == j:
            print(f"\nArea of the {i}-upper face, D_{i}{j}_orth:")
        else:
            print(f"\nMatrice Diedrale D_{i}{j}_orth:")

        geo.print_matrix(D_ij_orth)

print(f"\nArea of the down face, D_00_orth:")
D_00 = geo.diedral_matrix((0, 0), [v1_vec, v2_vec])
geo.print_matrix(geo.change_basis(D_00, orth_basis))


print(f"\nVolume of the tetrahedron, Vol_orth:")
Vol_std = geo.commutatore([v1_vec, v2_vec])
Vol_orth = geo.change_basis(Vol_std, orth_basis)
geo.print_matrix(np.abs(Vol_orth))
