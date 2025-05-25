from src.utilia import *

from src.spin_jjjj_tetrahedra.Spin1HalfTetrahedron import GroundState
from src.spin_jjjj_tetrahedra.Spin1Tetrahedron import Spin1Tetrahedron
from src.spin_jjjj_tetrahedra.Spin3HalfTetrahedron import Spin3_2Tetrahedron


class SpinTetrahedron:
    def __init__(self, j1, j2, j3, j4):
        """Inizialize a `wrapper` class for a general spin tetrahedron
            ρ = ρ^j1 ⊗ ρ^j2 ⊗ ρ^j3 ⊗ ρ^j4 : Spin(3) -> End(V_j1 ⊗ V_j2 ⊗ V_j3 ⊗ V_j4)
        """
        if j1 == j2 == j3 == j4 == 1/2:
            self.impl = GroundState()

        elif j1 == j4 == 1 and j2 == j3 == 1/2:
            self.impl = FirstExcitedState()

        elif j1 == j2 == j3 == j4 == 1:
            self.impl = Spin1Tetrahedron()

        elif j1 == j4 == 3/2 and j2 == j3 == 1:
            self.impl = ThirdExcitedState()

        elif j1 == j2 == j3 == j4 == 3/2:
            self.impl = Spin3_2Tetrahedron

        else:
            raise NotImplementedError(f"Spin values {j1},{j2},{j3},{j4} not yet covered.")

    def __getattr__(self, name):
        """Delegate all the attributes and methods to the actual wrapped class"""
        return getattr(self.impl, name)



class FirstExcitedState:
    def __init__(self):
        """Inizialize a class for the first excited state spin tetrahedron
            ρ = ρ^1 ⊗ ρ^1/2 ⊗ ρ^1/2 ⊗ ρ^1 : Spin(3) -> End(C^3 ⊗ C^2 ⊗ C^2 ⊗ C^3)
        """
        self.Vj1 = Vspace([Ket('++'), sqrt(2)*Ket('+-'), Ket('--')])
        self.Vj2 = Vspace([Ket('+'), Ket('-')])
        self.Vj3 = self.Vj2
        self.Vj4 = self.Vj1
        self.supp = TensorProductVspace(self.Vj1, self.Vj2, self.Vj3, self.Vj4)
        # Now it inherits all the attributes, instances and methods of a TensorProductVspace object



class ThirdExcitedState:
    def __init__(self):
        """Inizialize a class for the third excited state spin tetrahedron
            ρ = ρ^3/2 ⊗ ρ^1 ⊗ ρ^1 ⊗ ρ^3/2 : Spin(3) -> End(C^4 ⊗ C^3 ⊗ C^3 ⊗ C^4)
        """
        self.Vj1 = Vspace([Ket('++'), Ket('+-'), Ket('-+'), Ket('--')])
        self.Vj2 = Vspace([Ket('++'), sqrt(2)*Ket('+-'), Ket('--')])
        self.Vj3 = self.Vj2
        self.Vj4 = self.Vj1
        self.supp = TensorProductVspace(self.Vj1, self.Vj2, self.Vj3, self.Vj4)
        # Now it inherits all the attributes, instances and methods of a TensorProductVspace object



