from sympy.physics.quantum import *
from sympy import *
from src.utilia import Vspace, TensorProductVspace, pket, epsilon, bcolors


class Spin1Tetrahedron:
    def __init__(self):
        """Initialize a class for the second excited state spin tetrahedron
            ρ = ρ^1 ⊗ ρ^1 ⊗ ρ^1 ⊗ ρ^1 : Spin(3) -> End(C^3 ⊗ C^3 ⊗ C^3 ⊗ C^3)
        """
        # Define the basis for each C^3 space with standard j=1 representation
        # |++⟩ ≅ |1⟩, √2|+-⟩ ≅ |0⟩, |--⟩ ≅ |-1⟩
        self.Vj1 = Vspace([Ket('++'), sqrt(2)*Ket('+-'), Ket('--')])
        self.Vj2 = self.Vj1
        self.Vj3 = self.Vj2
        self.Vj4 = self.Vj3
        self.supp = TensorProductVspace(self.Vj1, self.Vj2, self.Vj3, self.Vj4)
        
        # Define the invariant subspace basis (spin networks)
        self.Inv_tensor = [
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
        ]  # This is Inv(ρ) subspace of C^3 ⊗ C^3 ⊗ C^3 ⊗ C^3 made of isotropic (ket) vectors

        # Convert to pretty ket format for easier manipulation
        self.Inv = [pket(inv) for inv in self.Inv_tensor]
        
        # Explicit expansion for coefficient extraction
        self.Inv_expanded = [expand(inv) for inv in self.Inv]

        print(f'\nSpin-1 tetrahedron Inv(ρ)-basis {bcolors.ITALIC}aka{bcolors.ENDC} {bcolors.BOLD}spin networks{bcolors.ENDC}:\n')
        for i, inv in enumerate(self.Inv):
            print(f"|Inv_{i}⟩ =")
            pprint(inv)
            print()

    def Pauli3D(self, ket_symbol):
        """Action of the 3D Pauli matrices L^i on C^3 with basis {|++⟩, √2|+-⟩, |--⟩}
        
        For j=1 representation, the correct generators are:
        L^1 = (1/√2) * [[0,1,0], [1,0,1], [0,1,0]]
        L^2 = (i/√2) * [[0,-1,0], [1,0,-1], [0,1,0]]  
        L^3 = [[1,0,0], [0,0,0], [0,0,-1]]
        
        :param ket_symbol: '++', '+-', or '--' (corresponding to |1⟩, |0⟩, |-1⟩)
        :return: [L^1|ket⟩, L^2|ket⟩, L^3|ket⟩]
        """
        match ket_symbol:
            case '++':  # |1⟩ state
                return [
                    (1/sqrt(2)) * sqrt(2)*Ket('+-'),     # L^1|++⟩ = (1/√2) * √2|+-⟩ = |+-⟩
                    (I/sqrt(2)) * (-sqrt(2)*Ket('+-')),  # L^2|++⟩ = (i/√2) * (-√2|+-⟩) = -i|+-⟩
                    Ket('++')                             # L^3|++⟩ = |++⟩
                ]
            case '+-' | '-+':  # |0⟩ state  
                return [
                    (1/sqrt(2)) * (Ket('++') + Ket('--')),    # L^1|+-⟩ = (1/√2) * (|++⟩ + |--⟩)
                    (I/sqrt(2)) * (Ket('++') - Ket('--')),    # L^2|+-⟩ = (i/√2) * (|++⟩ - |--⟩)
                    S.Zero                                     # L^3|+-⟩ = 0
                ]
            case '--':  # |-1⟩ state
                return [
                    (1/sqrt(2)) * sqrt(2)*Ket('+-'),     # L^1|--⟩ = (1/√2) * √2|+-⟩ = |+-⟩
                    (I/sqrt(2)) * sqrt(2)*Ket('+-'),     # L^2|--⟩ = (i/√2) * √2|+-⟩ = i|+-⟩
                    -Ket('--')                            # L^3|--⟩ = -|--⟩
                ]
            case _:
                raise ValueError(f"Invalid ket symbol: {ket_symbol}")

    def tau3D(self, ket_symbol):
        """Action of τᵢ = -i/2 Lᵢ matrices for j=1 representation"""
        pauli_results = self.Pauli3D(ket_symbol)
        return [-I/2 * result for result in pauli_results]

    def lie_operator(self, indexes, ket):
        """Function for the action of the Lie operator 𝐋_(i,a) on a spin-1 spinnet state
            i,a vary in {0,1,2,3} and {1,2,3} respectively.

        :param indexes: a tuple (i,a) for the a-th Pauli-tau matrix applied to the i-th component of the ket
        :param ket: a sympy.Ket object or linear combination
        :return: 𝐋_(i,a)|ket⟩
        """
        i, a = indexes
        if not (1 <= a <= 3):
            raise ValueError("The index 'a' ranges in {1,2,3}")
        if not (0 <= i <= 3):
            raise ValueError("The index 'i' ranges in {0,1,2,3}")

        # Handling linear combinations
        if isinstance(ket, Add):
            return Add(*[self.lie_operator((i, a), term) for term in ket.args])
        
        # Handle the case where ket is a Mul object (coefficient * linear combination)
        if isinstance(ket, Mul):
            # Check if any factor is an Add (linear combination)
            add_factor = None
            coeff_factors = []
            
            for factor in ket.args:
                if isinstance(factor, Add):
                    add_factor = factor
                else:
                    coeff_factors.append(factor)
            
            if add_factor is not None:
                # Apply lie_operator to the linear combination and multiply by coefficient
                coeff = Mul(*coeff_factors) if coeff_factors else S.One
                return coeff * self.lie_operator((i, a), add_factor)
        
        # Extract coefficient and pure ket
        coeff = S.One
        pure_ket = ket
        if isinstance(ket, Mul):
            coeff_factors = []
            ket_factor = None
            for factor in ket.args:
                if isinstance(factor, Ket):
                    ket_factor = factor
                else:
                    coeff_factors.append(factor)
            
            if ket_factor is not None:
                pure_ket = ket_factor
                coeff = Mul(*coeff_factors) if coeff_factors else S.One
            else:
                # No Ket found in Mul - this shouldn't happen in normal usage
                raise TypeError(f"Expected Ket in Mul expression, got {ket}")
        
        if not isinstance(pure_ket, Ket):
            raise TypeError(f"Expected Ket, got {type(pure_ket)}")
        
        # Get the label string (should be 8 characters for spin-1: 4 components × 2 chars each)
        label = str(pure_ket.label[0])
        if len(label) != 8:
            raise ValueError(f"Expected 8-character label for spin-1, got '{label}' with length {len(label)}")
        
        # Extract the i-th component (2 characters)
        component = label[2*i:2*i+2]
        
        # Apply τ matrix to the i-th component
        tau_result = self.tau3D(component)[a-1]
        
        # Build new ket(s)
        result = S.Zero
        if tau_result == S.Zero:
            return S.Zero
            
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
            elif term == S.Zero:
                continue
            else:
                tau_coeff = term
                continue
            
            if new_symbol is not None:
                # Construct new label by replacing i-th component
                new_label = label[:2*i] + new_symbol + label[2*i+2:]
                result += tau_coeff * Ket(new_label)
        
        return simplify(coeff * result)

    def dihedral(self, indexes, ket):
        """Function for dihedral operators on spin-1 spinnet states
                D𝛼𝛽 ∶= 𝐋_(𝛼,𝑎) 𝛿𝑎𝑏 𝐋_(𝛽,𝑏) ∈ Aut(𝒦_Γ)
        
        :param indexes: tuple (alpha, beta) with indices (0,1,2,3)
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

    def dihedral_matrix(self, alpha, beta, basis=None):
        """Compute the matrix representation of D_{αβ} on the invariant subspace
        :param alpha: first index (0,1,2,3)
        :param beta: second index (0,1,2,3)
        :param basis: list of basis (pket) vectors wrt express the matrix.
                        By default basis=Inv(ρ) spinnets basis
        :return: 3x3 Matrix representing D_{αβ} on basis
        """
        from src.utilia import extract_coefficients
        if basis is None:
            basis = self.Inv
        basis_expanded = [expand(vec) for vec in basis]

        print(f"\nComputing D_{{{alpha},{beta}}} matrix...")
        
        # Apply D_{αβ} to all basis vectors of the invariant subspace
        D_results = []
        for i, basis_vec in enumerate(basis):
            result = self.dihedral((alpha, beta), basis_vec)
            result = expand(result)
            D_results.append(result)
            print(f"D_{{{alpha},{beta}}}|Inv_{i}⟩ = {result}")
        
        # Extract coefficients: D_results[j] = Σᵢ M[i,j] * Inv[i]
        n = len(basis)
        matrix_elements: Matrix = Matrix.zeros(n, n)
        for j in range(n):  # for each result
            coeffs = extract_coefficients(D_results[j], basis_expanded)
            print(f"Coefficients for D_{{{alpha},{beta}}}|Inv_{j}⟩: {coeffs}")
            for i in range(n):  # for each basis vector
                matrix_elements[i, j] = coeffs[i]

        # Force symbolic representation and simplify
        matrix_elements = matrix_elements.applyfunc(lambda x: nsimplify(simplify(x), rational=False))
        
        pprint(matrix_elements)
        return matrix_elements

    def volume_squared(self, coeff=2/9, basis=None):
        """Compute the volume squared operator c * |[D_{1,3}, D_{1,2}]|"""
        
        def matrix_abs_elementwise(M: Matrix):
            """Apply |·| to each entry as |M|ij = |Mij|"""
            from sympy import Abs
            return M.applyfunc(Abs)
        
        def commutator(A, B):
            return A*B-B*A
            
        if basis is None:
            basis = self.Inv
        
        D13, D12 = self.dihedral_matrix(0, 2, basis), self.dihedral_matrix(0, 1, basis)     
        
        print(f"\nComputing V^2 matrix...")
        vol_squared = coeff * I * matrix_abs_elementwise(commutator(D13, D12))
        pprint(vol_squared)
        return vol_squared

    def change_basis_matrix(self, matrix, old_basis, new_basis):
        """
        Convert a matrix from one basis to another using the change of basis formula:
        M_new = P^(-1) * M_old * P
        where P is the change of basis matrix from new_basis to old_basis
        
        :param matrix: Matrix in the old basis representation
        :param old_basis: List of basis vectors for the original representation
        :param new_basis: List of basis vectors for the target representation  
        :return: Matrix in the new basis representation
        """
        from src.utilia import extract_coefficients
        
        if len(old_basis) != len(new_basis):
            raise ValueError("Old and new basis must have the same dimension")
        
        n = len(old_basis)
        
        # Expand all basis vectors for coefficient extraction
        old_basis_expanded = [expand(vec) for vec in old_basis]
        new_basis_expanded = [expand(vec) for vec in new_basis]
        
        # Build change of basis matrix P: new_basis[i] = Σⱼ P[j,i] * old_basis[j]
        P: Matrix = Matrix.zeros(n, n)
        for i in range(n):  # for each new basis vector
            coeffs = extract_coefficients(new_basis_expanded[i], old_basis_expanded)
            for j in range(n):  # coefficient of each old basis vector
                P[j, i] = coeffs[j]
        
        # Apply change of basis formula: M_new = P^(-1) * M_old * P
        P_inv = P.inv()
        matrix_new: Matrix = P_inv * matrix * P
        
        # Simplify and convert to symbolic form
        matrix_new = matrix_new.applyfunc(lambda x: nsimplify(simplify(x), rational=False))
        
        return matrix_new