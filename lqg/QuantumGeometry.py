from TensorVspaces import Vspace, TensorProductVspace
import sympy as sp
from sympy import symbols, pprint
from sympy.physics.quantum import Ket, TensorProduct
import itertools

class SpinTetrahedron:
    def __init__(self, j1, j2, j3, j4):

        self.j1 = j1
        self.j2 = j3
        self.j3 = j3
        self.j4 = j4
        self.dimension = (2 * j1 + 1) * (2 * j2 + 1) * (2 * j3 + 1) * (2 * j4 + 1)

        self.support1 = Vspace([Ket(A) for A in range(int(2 * j1 +1))])
        self.support2 = Vspace([Ket(B) for B in range(int(2 * j2 +1))])
        self.support3 = Vspace([Ket(C) for C in range(int(2 * j3 +1))])
        self.support4 = Vspace([Ket(D) for D in range(int(2 * j4 +1))])

        self.total_support = TensorProductVspace(self.support1, self.support2, self.support3, self.support4)

GroundState = SpinTetrahedron(1/2, 1/2, 1/2, 1/2).total_support

pprint(GroundState.basis)

# pprint(itertools.product(*GroundState.basis))

# for basis in itertools.product(*GroundState.bases):
#    pprint(TensorProduct(*basis))

pprint(SpinTetrahedron(1,1,1,1).total_support.basis)







