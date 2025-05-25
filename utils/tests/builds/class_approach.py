import sympy as sp
import sympy.physics.quantum as spq

class C:
    def __init__(self, n):
        self.dimension = int(n) #the complex dimension
        self.real_dimension = int(2*n)
        self.canonical_basis = [sp.symbols(f'e{j+1}') for j in range(self.dimension)]
        self.canonical_basis_val = [sp.eye(self.dimension)[:, j] for j in range(self.dimension)]
        self.basis_dict = dict(zip(self.canonical_basis, self.canonical_basis_val))


print(C(3).dimension)
print(C(3).real_dimension)
sp.pprint(C(5).canonical_basis_val)
sp.pprint(C(4).canonical_basis)
sp.pprint(C(3).basis_dict)


class SpinTetrahedron:
    def __init__(self, j1, j2, j3, j4):

        self.j1 = j1
        self.j2 = j2
        self.j3 = j3
        self.j4 = j4
        self.dimension = (2*self.j1+1)*(2*self.j2+1)*(2*self.j3+1)*(2*self.j4+1)

        self.support1 = C(2 * j1 + 1)
        self.support2 = C(2 * j2 + 1)
        self.support3 = C(2 * j3 + 1)
        self.support4 = C(2 * j4 + 1)

        if self.j1 == self.j2 and self.j2 == self.j3 and self.j3 == self.j4:
            if self.j1 == 1/2:
                print('helloo')

ground_state = SpinTetrahedron(1/2, 1/2, 1/2, 1/2)

print(ground_state.dimension)
print(ground_state.support4.dimension)




