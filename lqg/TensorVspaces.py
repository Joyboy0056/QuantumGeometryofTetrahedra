from sympy import pprint, Add, Mul, S, symbols, Matrix, sympify
from sympy.physics.quantum import Ket, TensorProduct

import itertools
import sympy as sp

class Vspace:
    """
    Inizializza uno spazio vettoriale con la sua base.
    :param basis: lista di vettori SymPy che compongono la base.
    """
    def __init__(self, basis):
        self.basis = basis
        self.dimension = len(basis)

        # Genera i vettori non simbolici per la base canonica
        self.vec_basis = [Matrix.eye(self.dimension)[:, j] for j in range(self.dimension)]

        # Registra la relazione in un dizionario
        self.basis_dict = {self.basis[i]: self.vec_basis[i] for i in range(self.dimension)}

        # Dizionario per memorizzare la rappresentazione degli elementi generati
        self.elements_dict = {}

    def get_element(self, *coeffs):
        """Restituisce la combinazione lineare dei vettori di base con i coefficienti dati
            e registra la sua rappresentazione simbolica e vettoriale.
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
    Classe per il prodotto tensoriale di k spazi vettoriali.
    """

    def __init__(self, *Vspaces):
        """
            Ogni Vspace deve essere un oggetto della classe Vspace
        """
        self.Vspaces = Vspaces
        self.k = len(Vspaces)  # Numero di fattori
        self.dimensions = [V.dimension for V in Vspaces]  # lista con le dimensioni dei fattori
        self.total_dimension = 1
        for dim in self.dimensions:
            self.total_dimension *= dim  # dimensione del prodotto

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



    def get_tensor_basis(self):
        res = []
        for tuple in list(itertools.product(*self.bases)):
            res.append(TensorProduct(*tuple))
        return res


def pretty_ket(expr):
    """
    Trasforma un'espressione composta da Ket, TensorProduct, o Add in Ket con etichetta '++-+'
    """
    if isinstance(expr, Ket):
        aux = ''
        for single_ket in expr.args:
            aux += f'{single_ket.label[0]}'
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
            coeff = S.One
            ket_expr = term

            if isinstance(term, Mul):
                coeffs = [f for f in term.args if not isinstance(f, (Ket, TensorProduct))]
                kets = [f for f in term.args if isinstance(f, (Ket, TensorProduct))]
                if coeffs:
                    coeff = Mul(*coeffs)
                if kets:
                    ket_expr = kets[0]

            new_ket = pretty_ket(ket_expr)
            if coeff == 1:
                terms.append(new_ket)
            else:
                terms.append(coeff * new_ket)
        return Add(*terms)

    else:
        raise TypeError("Input must be a Ket, TensorProduct, or Add of them.")
