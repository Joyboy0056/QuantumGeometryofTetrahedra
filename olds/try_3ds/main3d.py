import math
import numpy as np

from ket import Ket3D
from geometric_operators import GeometricOperators1

ket_3d = Ket3D()
GeoOps1 = GeometricOperators1()

# Vettori di base per Inv(1, 1, 1, 1)
print("\nSpin networks per (1, 1, 1, 1)")

# Qui lo spazio delle spin networks è 3-dimensionale
#V1 = ket_3d.compute_ket(([(0, 3), (1, 2), (4, 7), (5, 6)]))
#V2 = ket_3d.compute_ket(([(0, 7), (1, 6), (2, 5), (3, 4)]))
#V3 = ket_3d.compute_ket(([(0, 7), (1, 2), (3, 4), (5, 6)]))

#aggiorniamoli
V1 = ket_3d.compute_ket(([(0, 2), (1, 3), (4, 6), (5, 7)]))
V2 = ket_3d.compute_ket(([(0, 6), (1, 7), (2, 4), (3, 5)]))
V3 = ket_3d.compute_ket(([(0, 6), (1, 3), (2, 4), (5, 7)]))

print("v1 = ", V1)
print("v2 = ", V2)
print("v3 = ", V3)

#v1 = ket_3d.compute_vector([(0, 3), (1, 2), (4, 7), (5, 6)])
#v2 = ket_3d.compute_vector([(0, 7), (1, 6), (2, 5), (3, 4)])
#v3 = ket_3d.compute_vector([(0, 7), (1, 2), (3, 4), (5, 6)])
#aggiorniamoli
v1 = ket_3d.compute_vector(([(0, 2), (1, 3), (4, 6), (5, 7)]))
v2 = ket_3d.compute_vector(([(0, 6), (1, 7), (2, 4), (3, 5)]))
v3 = ket_3d.compute_vector(([(0, 6), (1, 3), (2, 4), (5, 7)]))

print(np.shape(v1))

D_12 = GeoOps1.diedral_matrix((1, 2), [v1, v2, v3])
GeoOps1.print_matrix(D_12)
print(np.shape(D_12))

D_12_op = GeoOps1.diedral_operator(1, 2)
print(np.shape(D_12_op))
D_12_op_v1 = GeoOps1.apply_operator_to_vector(D_12_op, v1)
D_12_op_v2 = GeoOps1.apply_operator_to_vector(D_12_op, v2)
D_12_op_v3 = GeoOps1.apply_operator_to_vector(D_12_op, v3)

print(ket_3d.vector_to_ket(D_12_op_v1))

#print(GeoOps1.parse_ket_expression(ket_3d.vector_to_ket(D_12_op_v1))) #funziona bene

print(GeoOps1.find_combination_coefficients(ket_3d.vector_to_ket(D_12_op_v1), V1, V2, V3))
print(GeoOps1.find_combination_coefficients(ket_3d.vector_to_ket(D_12_op_v2), V1, V2, V3))
print(GeoOps1.find_combination_coefficients(ket_3d.vector_to_ket(D_12_op_v3), V1, V2, V3))