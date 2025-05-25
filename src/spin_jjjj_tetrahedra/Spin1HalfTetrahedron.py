from sympy.physics.quantum import *
from sympy import *
from src.utilia import Vspace, TensorProductVspace, pket, epsilon


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
        
        self.Inv_tensor = [
            sum(
                epsilon[(A, B)] * epsilon[(C, D)] * TensorProduct(Ket(A), Ket(B), Ket(C), Ket(D))
                for A in ['+', '-'] for B in ['+', '-'] for C in ['+', '-'] for D in ['+', '-']
            ),
            sum(
            epsilon[(A, D)] * epsilon[(B, C)] * TensorProduct(Ket(A), Ket(B), Ket(C), Ket(D))
            for A in ['+', '-'] for B in ['+', '-'] for C in ['+', '-'] for D in ['+', '-']
            )
        ] # This is Inv(𝜌) subspace of C^2 ⊗ C^2 ⊗ C^2 ⊗ C^2 made of isotropic (ket) vectors

        # Convert to pretty ket format for easier manipulation
        self.Inv = [pket(inv) for inv in self.Inv_tensor]
        
        # Explicit expansion for coefficient extraction
        self.Inv_expanded = [expand(inv) for inv in self.Inv]

        print('\nGround state (1/2, 1/2, 1/2, 1/2) tetrahedron invariant basis vectors, aka spin networks:\n')
        for i, inv in enumerate(self.Inv):
            print(f"|Inv_{i}⟩ =")
            pprint(inv)
            


    def Pauli(self, ket_symbol):
        """A function modeling the action of the standard Pauli matrices σ = [σ₁, σ₂, σ₃] on C2
        :param ket_symbol: '+' or '-' character
        :return: [σ₁|ket⟩, σ₂|ket⟩, σ₃|ket⟩]
        """
        if ket_symbol == '+':
            return [Ket('-'), I*Ket('-'), Ket('+')]
        elif ket_symbol == '-':
            return [Ket('+'), -I*Ket('+'), -Ket('-')]
        else:
            raise ValueError(f"Invalid ket symbol: {ket_symbol}")


    def tau(self, ket_symbol):
        """Action of τᵢ = -i/2 σᵢ matrices"""
        pauli_results = self.Pauli(ket_symbol)
        return [-I/2 * result for result in pauli_results]

    def lie_operator(self, indexes, ket):
        """Function for the action of the Lie operator 𝐿_(i,a) on a spinnet state
            i,a vary in {0,1,2,3} and {1,2,3} respectively.

        :param indexes: a tuple (i,a) for the a-th Pauli-tau matrix applied to the i-th component of the ket
        :param ket: a sympy.Ket object
        :return: 𝐿_(i,a)|ket⟩
        """
        i, a = indexes
        if not (1 <= a <= 3):
            raise ValueError("The index 'a' ranges in {1,2,3}")
        if not (0 <= i <= 3):
            raise ValueError("The index 'i' ranges in {0,1,2,3}")

        # Handling linear combinations
        if isinstance(ket, Add):
            return Add(*[self.lie_operator((i, a), term) for term in ket.args])
        
        # Extract coefficient and pure ket
        coeff = S.One
        pure_ket = ket
        if isinstance(ket, Mul):
            coeff_factors = []
            for factor in ket.args:
                if isinstance(factor, Ket):
                    pure_ket = factor
                else:
                    coeff_factors.append(factor)
            coeff = Mul(*coeff_factors) if coeff_factors else S.One
        
        if not isinstance(pure_ket, Ket):
            raise TypeError(f"Expected Ket, got {type(pure_ket)}")
        
        # Get the label string (should be 4 characters for ground state)
        label = str(pure_ket.label[0])
        if len(label) != 4:
            raise ValueError(f"Expected 4-character label, got {label}")
        
        # Apply τ matrix to the i-th component
        tau_result = self.tau(label[i])[a-1]
        
        # Build new ket(s)
        result = S.Zero
        if isinstance(tau_result, Add):
            terms = tau_result.args
        else:
            terms = [tau_result]
        
        for term in terms:
            # Extract coefficient and new ket symbol from tau result
            tau_coeff = S.One
            new_symbol = None
            
            if isinstance(term, Mul):
                for factor in term.args:
                    if isinstance(factor, Ket):
                        new_symbol = str(factor.label[0])
                    else:
                        tau_coeff *= factor
            elif isinstance(term, Ket):
                new_symbol = str(term.label[0])
            else:
                tau_coeff = term
                continue
            
            if new_symbol is not None:
                # Construct new label by replacing i-th component
                new_label = label[:i] + new_symbol + label[i+1:]
                result += tau_coeff * Ket(new_label)
        
        return simplify(coeff * result)

    def dihedral(self, indexes, ket):
        """Function for dihedral operators on spinnet ground states
                D𝛼𝛽 ∶= 𝐿_(𝛼,𝑎) 𝛿𝑎𝑏 𝐿_(𝛽,𝑏) ∈ Aut(𝒦_Γ)
        
        :param alpha: first index (0,1,2,3)
        :param beta: second index (0,1,2,3)  
        :param ket: target ket
        :return: D_{αβ}|ket⟩
        """
        alpha, beta = indexes
        result = S.Zero
        for a in range(1, 4):  # a ∈ {1,2,3}
            # Apply L_(β,a) first, then L_(α,a)
            temp = self.lie_operator((beta, a), ket)
            result += self.lie_operator((alpha, a), temp)
        return simplify(result)
    

    def __getitem__(self, idx):
        """Magic method to make dihedral callable as `dihedral[i,j](ket)`"""
        i, j = idx
        D = [
            [lambda ket, i=i, j=j: self.dihedral((i, j), ket) for j in range(4)]
            for i in range(4)
        ]
        return lambda ket: self.dihedral(i, j, ket)


    # By Claude Sonnet 4
    def dihedral_matrix(self, alpha, beta):
        """Compute the matrix representation of D_{αβ} on the invariant subspace
        :param alpha: first index (0,1,2,3)
        :param beta: second index (0,1,2,3)
        :return: 2x2 Matrix representing D_{αβ} on Inv(ρ)
        """

        from src.QuantumGeometry import extract_coefficients

        # print(f"\nComputing D_{{{alpha},{beta}}} matrix...")
        
        # Apply D_{αβ} to both basis vectors of the invariant subspace
        D_results = []
        for i, basis_vec in enumerate(self.Inv):
            result = self.dihedral((alpha, beta), basis_vec)
            result = expand(result)
            D_results.append(result)
            # print(f"D_{{{alpha},{beta}}}|Inv_{i}⟩ = {result}")
        
        # Extract coefficients: D_results[j] = Σᵢ M[i,j] * Inv[i]
        matrix_elements = Matrix.zeros(2, 2)
        for j in range(2):  # for each result
            coeffs = extract_coefficients(D_results[j], self.Inv_expanded)
            # print(f"Coefficients for D_{{{alpha},{beta}}}|Inv_{j}⟩: {coeffs}")
            for i in range(2):  # for each basis vector
                matrix_elements[i, j] = coeffs[i]
        
        # print(f"Matrix representation:\n{matrix_elements}")
        return matrix_elements
    

    def volume_squared(self, coeff=1):
        """Compute the volume squared operator c * |[D_{1,3}, D_{1,2}]|"""
        
        def matrix_abs_elementwise(M):
            """Apply |·| to each entry as |M|ij = |Mij|"""
            from sympy import Abs
            return M.applyfunc(Abs)
        
        def commutator(A, B):
            return A*B-B*A

        D13, D12 = self.dihedral_matrix(1,3), self.dihedral_matrix(1,2)     
        
        return coeff * matrix_abs_elementwise(commutator(D13, D12))