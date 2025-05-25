from sympy.physics.quantum import *
from sympy import *
from src.utilia import Vspace, TensorProductVspace, pket, epsilon
from src.QuantumGeometry import extract_coefficients

class Spin1Tetrahedron:
    def __init__(self):
        """Inizialize a class for the second excited state spin tetrahedron
            ρ = ρ^1 ⊗ ρ^1 ⊗ ρ^1 ⊗ ρ^1 : Spin(3) -> End(C^3 ⊗ C^3 ⊗ C^3 ⊗ C^3)
        """
        self.Vj1 = Vspace([Ket('++'), sqrt(2)*Ket('+-'), Ket('--')])
        self.Vj2 = self.Vj1
        self.Vj3 = self.Vj2
        self.Vj4 = self.Vj3
        self.supp = TensorProductVspace(self.Vj1, self.Vj2, self.Vj3, self.Vj4)
        # Now it inherits all the attributes, instances and methods of a TensorProductVspace object
        self.Inv = [
            sum(
                epsilon[(a1, b1)] * epsilon[(a2, b2)] * epsilon[(c1, d1)] * epsilon[(c2, d2)] * TensorProduct(
                    Ket(a1 + a2),
                    Ket(b1 + b2),
                    Ket(c1 + c2),
                    Ket(d1 + d2))
                for a1 in ['+', '-'] for b1 in ['+', '-'] for c1 in ['+', '-'] for d1 in ['+', '-']
                for a2 in ['+', '-'] for b2 in ['+', '-'] for c2 in ['+', '-'] for d2 in ['+', '-']
            ),
            sum(
                epsilon[(a1, d1)] * epsilon[(a2, d2)] * epsilon[(b1, c1)] * epsilon[(b2, c2)] * TensorProduct(
                    Ket(a1 + a2),
                    Ket(b1 + b2),
                    Ket(c1 + c2),
                    Ket(d1 + d2))
                for a1 in ['+', '-'] for b1 in ['+', '-'] for c1 in ['+', '-'] for d1 in ['+', '-']
                for a2 in ['+', '-'] for b2 in ['+', '-'] for c2 in ['+', '-'] for d2 in ['+', '-']
            ),
            sum(
                epsilon[(a1, d1)] * epsilon[(a2, b2)] * epsilon[(b1, c1)] * epsilon[(c2, d2)] * TensorProduct(
                    Ket(a1 + a2),
                    Ket(b1 + b2),
                    Ket(c1 + c2),
                    Ket(d1 + d2))
                for a1 in ['+', '-'] for b1 in ['+', '-'] for c1 in ['+', '-'] for d1 in ['+', '-']
                for a2 in ['+', '-'] for b2 in ['+', '-'] for c2 in ['+', '-'] for d2 in ['+', '-']
            )
        ] # This is Inv(𝜌) subspace of C^3 ⊗ C^3 ⊗ C^3 ⊗ C^3 made of isotropic (ket) vectors

    def Pauli(self, ket):
        """A function modeling the action of the standard Pauli matrices σ = [σ₁, σ₂, σ₃] on C2
        :param ket: a Ket object from sympy.physics.quantum.Ket
        """
        if not isinstance(ket, Ket):
            raise TypeError("Input must be a Ket")
        sigma1, sigma2, sigma3 = Ket(''), Ket(''), Ket('')
        if ket == Ket('++'):
        #if ket.args[0] == '+':
            sigma1, sigma2, sigma3 = -I*Ket('+-'), Ket('+-'), -I*Ket('++')
        elif ket == Ket('+-') or ket == Ket('-+'):
        #elif ket.args[0] == '-':
            sigma1, sigma2, sigma3 = -I*Ket('++') -I*Ket('--'), -Ket('++') -Ket('--'), 0*Ket('++') + 0*Ket('+-') +0*Ket('--')
        elif ket == Ket('--'):
        #elif ket.args[0] == '-':
            sigma1, sigma2, sigma3 = -I*Ket('+-'), -Ket('+-'), I*Ket('--')

        return [sigma1, sigma2, sigma3]

    def tau(self, ket):
        if not isinstance(ket, Ket):
            raise TypeError("Input must be a Ket")
        return [simplify(-I/2*res) for res in self.Pauli(ket)]

    def lie_operator(self, indexes, ket):
        """Function for the action of the Lie operator L_(i,a) on a spinnet state
            i,a vary in {0,1,2,3} and {1,2,3} respectively.
        :param indexes: a tuple (i,a) for the a-th Pauli-tau matrix applied to the i-th component of the ket
        :param ket: a sympy.Ket object
        """
        i, a = indexes
        if a == 0:
            raise ValueError("The index 'a' ranges in {1,2,3}")
        elif i == 4:
            raise ValueError("The index 'i' ranges in {0,1,2,3}")

        # Se ket è Add estendo linearmente
        if isinstance(ket, Add):
            return Add(*[self.lie_operator(indexes, term) for term in ket.args])

        # Se ket è una moltiplicazione coefficiente * ket
        coeff = S.One
        true_ket = ket
        if isinstance(ket, Mul):
            for factor in ket.args:
                if isinstance(factor, Ket):
                    true_ket = factor
                else:
                    coeff *= factor
        elif isinstance(ket, Ket):
            true_ket = ket
        else:
            raise TypeError(f"Expected Ket or Mul of Ket, got {type(ket)}")
        # Ora true_ket è sicuramente un Ket

        label = str(true_ket.label[0])  # <-- label?
        # if len(label) % 2 != 0:
        #     raise ValueError("Label must have even number of symbols to form pairs.")

        pairs = [(label[j], label[j + 1]) for j in range(0, len(label), 2)]
        pauli_action = self.Pauli(Ket(f'{pairs[i][0]}{pairs[i][1]}'))[a-1]

        # Ricostruiamo il nuovo Ket usando il risultato della Pauli action
        result = S.Zero
        for term in pauli_action.args if isinstance(pauli_action, Add) else [pauli_action]:
            if isinstance(term, Mul):
                factors = term.args
                local_coeff = S.One
                local_ket = None
                for factor in factors:
                    if isinstance(factor, Ket):
                        local_ket = factor
                    else:
                        local_coeff *= factor
            else:
                local_coeff = S.One
                local_ket = term

            label_term = str(local_ket.label[0])
            new_pairs = list(pairs)
            new_pairs[i] = (label_term[0], label_term[1])
            new_label = ''.join([x + y for (x, y) in new_pairs])
            result += local_coeff * Ket(new_label)
        return simplify(expand(coeff * result))


    def dihedral(self, indixes, ket):
        """Function for dihedral operators on spinnet ground states
                D𝛼𝛽 ∶= 𝐿_(𝛼,𝑎) 𝛿𝑎𝑏 𝐿_(𝛽,𝑏) ∈ Aut(𝒦_Γ)
        """
        alpha, beta = indixes
        return sum(self.lie_operator((alpha, a), self.lie_operator((beta, a), ket)) for a in range(1, 4))