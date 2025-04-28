from sympy import Add, Mul, S, Matrix
from sympy.physics.quantum import Ket, TensorProduct

import itertools

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


def pretty_ket(expr):
    """
    Transforms an expression made of Ket, TensorProduct, Add or Mul in a single Ket
    e.g. 2*|+> ⊗ |-> ⊗ 3*|+> ⊗ |-> = 6*|+-+->
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
            terms.append(pretty_ket(term))
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
        new_ket = pretty_ket(ket_expr)
        if coeff == 1:
            return new_ket
        else:
            return coeff * new_ket

    else:
        raise TypeError("Input must be a Ket, TensorProduct, Add or Mul of them.")

