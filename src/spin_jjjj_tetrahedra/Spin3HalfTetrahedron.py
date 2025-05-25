from src.utilia import *

class Spin3_2Tetrahedron:
    def __init__(self):
        """Inizialize a class for the fourth excited state spin tetrahedron
            ρ = ρ^3/2 ⊗ ρ^3/2 ⊗ ρ^3/2 ⊗ ρ^3/2 : Spin(3) -> End(C^4 ⊗ C^4 ⊗ C^4 ⊗ C^4)
        """
        self.Vj1 = Vspace([Ket('++'), Ket('+-'), Ket('-+'), Ket('--')])
        self.Vj2 = self.Vj1
        self.Vj3 = self.Vj2
        self.Vj4 = self.Vj1
        self.supp = TensorProductVspace(self.Vj1, self.Vj2, self.Vj3, self.Vj4)
        # Now it inherits all the attributes, instances and methods of a TensorProductVspace object