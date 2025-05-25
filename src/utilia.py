from sympy import Add, Mul, S, Matrix
from sympy.physics.quantum import Ket, TensorProduct
from sympy import expand, symbols, solve, simplify, sqrt, Eq
import itertools

# ============== Classes ===============

# Colors class
class bcolors:
    BLACK = "\033[30m"
    RED = "\033[91m"
    GREEN = "\033[92m"
    YELLOW = "\033[33m"
    BLUE = "\033[94m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    BORDEAUX = "\033[38;5;52m"         # Bordeaux (scuro)
    PURPLE = "\033[38;5;129m"          # Viola
    MIDNIGHT_BLUE = "\033[38;5;17m"    # Blu notte
    DARK_GREEN = "\033[38;5;22m"       # Verde scuro
    DARK_GREY = "\033[90m"             # Grigio scuro
    ORANGE = "\033[38;5;208m"          # Arancione
    DARK_ORANGE = "\033[38;5;166m"     # Arancione scuro
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
    ENDC = "\033[0m"
    ITALIC = "\033[3m"                 # Corsivo 


# Vector spaces classes
class Vspace:
    """A class for a vector space"""
    def __init__(self, basis):
        """
            Inizialize a vector space with its basis.
            :param basis: list of SymPy vectors for the basis.
        """
        self.basis = basis
        self.dimension = len(basis)

        # Genera i vettori non simbolici per la base canonica
        self.vec_basis = [Matrix.eye(self.dimension)[:, j] for j in range(self.dimension)]

        # Registra la relazione in un dizionario
        self.basis_dict = {self.basis[i]: self.vec_basis[i] for i in range(self.dimension)}

        # Dizionario per memorizzare la rappresentazione degli elementi generati
        self.elements_dict = {}

    def get_element(self, *coeffs):
        """Returns the linear combination of the basis vectors with the inpunt coefficients
            and records both its symbolic and vectorial representation in the dict "self.elements_dict"
        """
        if len(coeffs) != self.dimension:
            raise ValueError("Il numero di coefficienti deve corrispondere alla dimensione dello spazio vettoriale.")

        # Creazione dell'elemento simbolico
        res = sum(c * v for c, v in zip(coeffs, self.basis))

        # Creazione dell'elemento vettoriale
        res_vec = sum((c * vec for c, vec in zip(coeffs, self.vec_basis)),
                                start=Matrix.zeros(self.dimension, 1))

        # Registra l'elemento nel dizionario
        self.elements_dict[res] = res_vec

        return res

    # def get_elements_dict(self):
    #     """Restituisce il dizionario con tutte le rappresentazioni degli elementi generati."""
    #     return self.elements_dict
    # non serve, si aggiorna automaticamente con get_element(coeffs)


class TensorProductVspace:
    """
    Class for the tensorial product of k vector spaces.
    """

    def __init__(self, *Vspaces):
        """
            Inizialize the tensor product space of k vector spaces.
            :param Vspaces: Objects of the Vspace class above.
        """
        self.Vspaces = Vspaces
        self.k = len(Vspaces)  # Numero di fattori
        self.dimensions = [V.dimension for V in Vspaces]  # lista con le dimensioni dei fattori
        self.dimension = 1
        for dim in self.dimensions:
            self.dimension *= dim  # dimensione del prodotto

        self.bases = [V.basis for V in Vspaces]
        self.vec_bases = [V.vec_basis for V in Vspaces]

        # Costruiamo la base del prodotto tensoriale
        #self.basis = list(itertools.product(*self.bases))  # Base simbolica come tuple di vettori
        self.vec_basis = list(itertools.product(*self.vec_bases))  # Base vettoriale come tuple di matrici

        self.basis = [TensorProduct(*basis) for basis in itertools.product(*self.bases)]

        # # Genera la base vettoriale combinata
        # self.vec_basis = [Matrix.tensor_product(vecs) for vecs in product((v.vec_basis for v in Vspaces))]

        # # Dizionario di associazione simbolica-vettoriale
        # self.basis_dict = {self.basis[i]: self.vec_basis[i] for i in range(len(self.basis))}

        # Dizionario per elementi generati
        self.elements_dict = {}

    # def get_element(self, *coeffs):
    #    pass

    def get_tensor_basis(self): # Unused
        res = []
        for tuple in list(itertools.product(*self.bases)):
            res.append(TensorProduct(*tuple))
        return res


# =================== AUXs utilia ====================

epsilon = {
                ('+', '+'): 0,
                ('-', '-'): 0,
                ('+', '-'): 1,
                ('-', '+'): -1
            }

def pket(expr):
    """
    Pretty transforms an expression made of Ket, TensorProduct, Add or Mul in a single Ket
    e.g. 2*|+⟩ ⊗ |-⟩ ⊗ 3*|+⟩ ⊗ |-⟩ = 6*|+-+-⟩
    """
    if expr == S.Zero:
        return S.Zero

    if isinstance(expr, Ket):
        # If already a single Ket with the full string, just return it
        if isinstance(expr.label[0], str):
            return expr
        # Otherwise, build the full string
        aux = ''.join(single_ket.label[0] for single_ket in expr.args)
        return Ket(aux)


    elif isinstance(expr, TensorProduct):
        aux = ''
        for component in expr.args:
            if isinstance(component, Ket):
                aux += f'{component.label[0]}'
            else:
                raise ValueError(f"Unexpected component in TensorProduct: {component}")
        return Ket(aux)

    elif isinstance(expr, Add):
        terms = []
        for term in expr.args:
            terms.append(pket(term))
        return Add(*terms)

    elif isinstance(expr, Mul):
        coeff = S.One
        ket_expr = None
        for factor in expr.args:
            if isinstance(factor, (Ket, TensorProduct)):
                ket_expr = factor
            else:
                coeff *= factor
        if ket_expr is None:
            raise ValueError("Mul does not contain a Ket or TensorProduct.")
        new_ket = pket(ket_expr)
        if coeff == 1:
            return new_ket
        else:
            return coeff * new_ket

    else:
        raise TypeError("Input must be a Ket, TensorProduct, Add or Mul of them.")
    

# Claude function for coefficients extraction
def extract_coefficients(expr, basis_vectors):
    """Extract coefficients when expressing expr as linear combination of basis_vectors
    Returns coefficients c_i such that expr = Σ c_i * basis_vectors[i]
    
    Uses sympy's linear system solving for robust coefficient extraction.
    """
    if expr == S.Zero:
        return [S.Zero] * len(basis_vectors)
    
    # Expand everything
    expr = expand(expr)
    basis_expanded = [expand(basis_vec) for basis_vec in basis_vectors]
    
    # Collect all unique ket terms from expression and basis vectors
    all_kets = set()
    
    def collect_kets(expression):
        if isinstance(expression, Add):
            for term in expression.args:
                if isinstance(term, Ket):
                    all_kets.add(term)
                elif isinstance(term, Mul):
                    for factor in term.args:
                        if isinstance(factor, Ket):
                            all_kets.add(factor)
        elif isinstance(expression, Ket):
            all_kets.add(expression)
        elif isinstance(expression, Mul):
            for factor in expression.args:
                if isinstance(factor, Ket):
                    all_kets.add(factor)
    
    collect_kets(expr)
    for basis_vec in basis_expanded:
        collect_kets(basis_vec)
    
    all_kets = sorted(list(all_kets), key=str)
    
    if not all_kets:
        return [S.Zero] * len(basis_vectors)
    
    # Create coefficient extraction function
    def get_coeff_of_ket(expression, target_ket):
        """Get coefficient of target_ket in expression"""
        if expression == S.Zero:
            return S.Zero
        
        expanded = expand(expression)
        if isinstance(expanded, Add):
            coeff = S.Zero
            for term in expanded.args:
                if isinstance(term, Ket) and term == target_ket:
                    coeff += S.One
                elif isinstance(term, Mul):
                    ket_part = None
                    coeff_part = S.One
                    for factor in term.args:
                        if isinstance(factor, Ket):
                            ket_part = factor
                        else:
                            coeff_part *= factor
                    if ket_part == target_ket:
                        coeff += coeff_part
            return coeff
        elif isinstance(expanded, Ket) and expanded == target_ket:
            return S.One
        elif isinstance(expanded, Mul):
            ket_part = None
            coeff_part = S.One
            for factor in expanded.args:
                if isinstance(factor, Ket):
                    ket_part = factor
                else:
                    coeff_part *= factor
            if ket_part == target_ket:
                return coeff_part
        
        return S.Zero
    
    # Set up linear system: expr = Σ c_i * basis_i
    # For each ket k: coeff_of_k(expr) = Σ c_i * coeff_of_k(basis_i)
    
    n_basis = len(basis_vectors)
    equations = []
    coeffs = symbols([f'c_{i}' for i in range(n_basis)])
    
    for ket in all_kets:
        expr_coeff = get_coeff_of_ket(expr, ket)
        basis_coeffs = [get_coeff_of_ket(basis_expanded[i], ket) for i in range(n_basis)]
        
        # Skip if this ket doesn't appear anywhere
        if expr_coeff == S.Zero and all(bc == S.Zero for bc in basis_coeffs):
            continue
        
        # Build equation: expr_coeff = Σ c_i * basis_coeffs[i]
        lhs = sum(coeffs[i] * basis_coeffs[i] for i in range(n_basis))
        equations.append(Eq(lhs, expr_coeff))
    
    if not equations:
        return [S.Zero] * len(basis_vectors)
    
    # Solve the linear system
    try:
        solution = solve(equations, coeffs)
        if isinstance(solution, dict):
            result = [simplify(solution.get(coeffs[i], S.Zero)) for i in range(n_basis)]
        else:
            # Handle case where solution is a list (parametric solutions)
            result = [S.Zero] * n_basis
            print(f"Warning: Could not find unique solution. Equations: {equations}")
    except Exception as e:
        print(f"Error solving linear system: {e}")
        print(f"Equations: {equations}")
        result = [S.Zero] * n_basis
    
    return result