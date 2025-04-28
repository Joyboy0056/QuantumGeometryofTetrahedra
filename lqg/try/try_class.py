from TensorVspaces import Vspace, TensorProductVspace, pretty_ket
import sympy as sp

from sympy import symbols, pprint
from sympy.physics.quantum import Ket

# Esempio di utilizzo
R2 = Vspace([symbols('e1'), symbols('e2')])

pprint(R2.basis)
pprint(R2.dimension)

kets = [Ket('+'), Ket('-')]
C2 = Vspace(kets)
pprint(C2.basis)
pprint(C2.dimension)
pprint(C2.vec_basis)

# Creazione di nuovi elementi nello spazio
v1 = C2.get_element(2, -1)
v2 = C2.get_element(0, 4)

# Stampa delle rappresentazioni
pprint(v1)
pprint(v2)

# res = Matrix.zeros(C2.dimension, 1)
# coeffs = (2, -1)
# for c in coeffs:
#     for vec in C2.vec_basis:
#         res += c * vec
# pprint(res)

pprint(C2.elements_dict)

pprint(TensorProductVspace(C2, C2, C2, C2).bases)

pprint(TensorProductVspace(C2, C2, C2, C2).basis)

kets = [Ket('++'), sp.sqrt(2) * Ket('+-'), Ket('--')]
C3 = Vspace(kets)

pprint(TensorProductVspace(C3, C3, C3, C3).total_dimension)
pprint(TensorProductVspace(C3, C3, C3, C3).basis)

pprint(TensorProductVspace(C3, C2, C3, C2).total_dimension)
pprint(TensorProductVspace(C3, C2, C3, C2).basis)

ground_state = Vspace(TensorProductVspace(C2, C2, C2, C2).basis)

v = ground_state.get_element(1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0)
pprint(v)


w = ground_state.get_element(1,0,0,1,0,1,0,0,1,0,1,0,0,0,0,1)
#w = ground_state.get_element(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
pprint(pretty_ket(w))
print(type(w))