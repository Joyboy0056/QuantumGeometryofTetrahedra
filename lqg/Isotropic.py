from TensorVspaces import Vspace, TensorProductVspace, pretty_ket

from sympy import pprint
from sympy.physics.quantum import Ket, TensorProduct

C2 = Vspace([Ket('+'), Ket('-')])

GroundState = TensorProductVspace(C2, C2, C2, C2)

print('\nGround state (1/2, 1/2, 1/2, 1/2) tetrahedron basis:\n')
pprint([pretty_ket(GroundState.basis[j]) for j in range(len(GroundState.basis))])

epsilon = {
                ('+', '+'): 0,
                ('-', '-'): 0,
                ('+', '-'): 1,
                ('-', '+'): -1
            }

v1 = sum(
    epsilon[(A, B)] * epsilon[(C, D)] * TensorProduct(Ket(A), Ket(B), Ket(C), Ket(D))
    for A in ['+', '-'] for B in ['+', '-'] for C in ['+', '-'] for D in ['+', '-']
)
v2 = sum(
    epsilon[(A, D)] * epsilon[(B, C)] * TensorProduct(Ket(A), Ket(B), Ket(C), Ket(D))
    for A in ['+', '-'] for B in ['+', '-'] for C in ['+', '-'] for D in ['+', '-']
)

Inv = [v1, v2]

print('\nIsotropic subspace of the spin (1/2, 1/2, 1/2, 1/2) tetrahedron:\n')
pprint([pretty_ket(Inv[j]) for j in range(len(Inv))])

