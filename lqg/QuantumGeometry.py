from soupsieve.pretty import pretty

from TensorVspaces import Vspace, TensorProductVspace
from sympy.physics.quantum import Ket, TensorProduct
from sympy import I, Add, Mul, sqrt, S, simplify, expand, symbols, Eq, solve, Matrix

from math import prod

epsilon = {
                ('+', '+'): 0,
                ('-', '-'): 0,
                ('+', '-'): 1,
                ('-', '+'): -1
            }

class SpinTetrahedron:
    def __init__(self, j1, j2, j3, j4):
        """Inizialize a class for a general spin tetrahedron
            ρ = ρ^j1 ⊗ ρ^j2 ⊗ ρ^j3 ⊗ ρ^j4 : Spin(3) -> End(V_j1 ⊗ V_j2 ⊗ V_j3 ⊗ V_j4)
        """
        self.j1 = j1
        self.j2 = j3
        self.j3 = j3
        self.j4 = j4
        self.dimension = (2 * j1 + 1) * (2 * j2 + 1) * (2 * j3 + 1) * (2 * j4 + 1)

        # V_ji support spaces of the factors ρ^ji
        self.Vj1 = Vspace([Ket(A) for A in range(int(2 * j1 +1))])
        self.Vj2 = Vspace([Ket(B) for B in range(int(2 * j2 +1))])
        self.Vj3 = Vspace([Ket(C) for C in range(int(2 * j3 +1))])
        self.Vj4 = Vspace([Ket(D) for D in range(int(2 * j4 +1))])

        self.supp = TensorProductVspace(self.Vj1, self.Vj2, self.Vj3, self.Vj4) # Support space of ρ

        # Left: Generalize computation of Inv subspace (then operators) for general Spin Tetrahedron (quite challenging)


class GroundState:
    def __init__(self):
        """Inizialize a class for the ground state spin tetrahedron
            ρ = ρ^1/2 ⊗ ρ^1/2 ⊗ ρ^1/2 ⊗ ρ^1/2 : Spin(3) -> End(C^2 ⊗ C^2 ⊗ C^2 ⊗ C^2)
        """
        self.Vj1 = Vspace([Ket('+'), Ket('-')])
        self.Vj2 = self.Vj1
        self.Vj3 = self.Vj2
        self.Vj4 = self.Vj3
        self.supp = TensorProductVspace(self.Vj1, self.Vj2, self.Vj3, self.Vj4)
        # Now it inherits all the attributes, instances and methods of a TensorProductVspace object
        self.Inv = [
            sum(
                epsilon[(A, B)] * epsilon[(C, D)] * TensorProduct(Ket(A), Ket(B), Ket(C), Ket(D))
                for A in ['+', '-'] for B in ['+', '-'] for C in ['+', '-'] for D in ['+', '-']
            ),
            sum(
            epsilon[(A, D)] * epsilon[(B, C)] * TensorProduct(Ket(A), Ket(B), Ket(C), Ket(D))
            for A in ['+', '-'] for B in ['+', '-'] for C in ['+', '-'] for D in ['+', '-']
            )
        ] # This is Inv(𝜌) subspace of C^2 ⊗ C^2 ⊗ C^2 ⊗ C^2 made of isotropic (ket) vectors

    def Pauli(self, ket):
        """A function modeling the action of the standard Pauli matrices σ = [σ₁, σ₂, σ₃] on C2
        :param ket: a Ket object from sympy.physics.quantum.Ket
        """
        if not isinstance(ket, Ket):
            raise TypeError("Input must be a Ket")
        sigma1, sigma2, sigma3 = Ket(''), Ket(''), Ket('')
        if ket == Ket('+'):
        #if ket.args[0] == '+':
            sigma1, sigma2, sigma3 = Ket('-'), I*Ket('-'), Ket('+')
        elif ket == Ket('-'):
        #elif ket.args[0] == '-':
            sigma1, sigma2, sigma3 = Ket('+'), -I * Ket('+'), Ket('-')

        return [sigma1, sigma2, sigma3]

    def tau(self, ket):
        if not isinstance(ket, Ket):
            raise TypeError("Input must be a Ket")
        return [-I/2*res for res in self.Pauli(ket)]

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

        if isinstance(ket, Ket) or isinstance(ket, Mul):
            res = Ket('')
            if len(ket.args) > 1:
                # ket.args could be (-1/2, I, |->)
                coeff = prod(ket.args[:-1])  # Prodotto dei primi n-1 elementi
                ket = ket.args[-1]
            elif len(ket.args) == 1:
                coeff, ket = (1, Ket(ket.args[0]))
            labels = [lab for lab in str(ket.label[0])]

            single_ket = self.tau(Ket(str(ket.args[0])[i]))[a-1]

            if len(single_ket.args) == 1:
                sket = single_ket.args[1]
                labels[2] = sket.args[0]
                res += coeff * Ket(f'{labels[0]}{labels[1]}{labels[2]}{labels[3]}')
            elif len(single_ket.args) > 1:
                coef = prod(single_ket.args[:-1])  # Prodotto dei primi n-1 elementi
                sket = single_ket.args[-1]
                labels[i] = sket.args[0]
                res += coef * coeff * Ket(f'{labels[0]}{labels[1]}{labels[2]}{labels[3]}')


        elif isinstance(ket, Add):
            deeplist = [ket.args[j].args for j in range(len(ket.args))]
            # deeplist is [(coeff, ket)] or [(ket,)]
            res = Ket('')
            for t in deeplist:
                if len(t) > 1:
                    coeff = prod(t[:-1])  # Prodotto dei primi n-1 elementi
                    ket = t[-1]
                    # print(coeff, ket)
                elif len(t) == 1:
                    coeff = 1
                    ket = Ket(t[0])
                    # print(coeff, ket)
                labels = [lab for lab in str(ket.label[0])]

                single_ket = self.tau(Ket(str(ket.args[0])[i]))[a-1]

                if len(single_ket.args) == 1:
                    sket = single_ket.args[1]
                    print(sket)
                    labels[2] = sket.args[0]
                    res += coeff * Ket(f'{labels[0]}{labels[1]}{labels[2]}{labels[3]}')
                elif len(single_ket.args) > 1:
                    coef = prod(single_ket.args[:-1]) # Prodotto dei primi n-1 elementi
                    sket = single_ket.args[-1]
                    labels[i] = sket.args[0]
                    res += coef * coeff * Ket(f'{labels[0]}{labels[1]}{labels[2]}{labels[3]}')

        elif isinstance(ket, TensorProduct):
            raise TypeError("If you wanna apply Lie ops to sympy.TensorProduct, please use 'pretty_ket' on it before.")

        else:
            raise ValueError("This function handles sympy.Ket, sympy.Mul or sympy.Add only.")

        return res - Ket('')

    def dihedral(self, indixes, ket):
        """Function for dihedral operators on spinnet ground states
                D𝛼𝛽 ∶= 𝐿_(𝛼,𝑎) 𝛿𝑎𝑏 𝐿_(𝛽,𝑏) ∈ Aut(𝒦_Γ)
        """
        alpha, beta = indixes
        return sum(self.lie_operator((alpha, a), self.lie_operator((beta, a), ket)) for a in range(1, 4))


    # LEFT: dihedral_matrix


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

        label = str(true_ket.label[0])  # <-- label corretto, senza fare str(ket.args[0])
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
