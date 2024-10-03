# lqg_framework/main.py
# Esempio di utilizzo della classe ket
from ket import Ket2D, Ket3D

# Esempio di utilizzo ket2D per lo spazio Inv(1/2, 1/2, 1/2, 1/2)
ket_2d = Ket2D()
v1 = ket_2d.contract_with_epsilon((0, 1), (2, 3))  # A-B e C-D
print("\nRisultato della contrazione (2D A-B e C-D):", v1)

v2 = ket_2d.contract_with_epsilon((0, 3), (1, 2))
print("\nRisultato della contrazione (2D A-D e B-C):", v2)


# Vettori di base per Inv(1, 1, 1, 1)

ket = Ket3D()

# Contrai i vettori con epsilon
v1_terms = ket.contract_with_epsilon([(0, 3), (1, 2), (4, 7), (5, 6)])
v1_mapped = ket.map_to_C3(v1_terms)

print("\nv1:", v1_mapped)

v2_terms = ket.contract_with_epsilon([(0, 7), (1, 6), (2, 5), (3, 4)])
v2_mapped = ket.map_to_C3(v2_terms)

print("\nv2:", v2_mapped)

v3_terms = ket.contract_with_epsilon([(0, 7), (1, 2), (3, 4), (5, 6)])
v3_mapped = ket.map_to_C3(v3_terms)

print("\nv3:", v3_mapped)
