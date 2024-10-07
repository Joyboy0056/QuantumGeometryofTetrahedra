import numpy as np
import sympy as sp

from ket import Ket2D, Ket3D
ket_2d = Ket2D()
ket_3d = Ket3D()

import scipy.linalg as la

import math

class GeometricOperators:
    def __init__(self):
        # Definiamo le matrici di Pauli
        self.sigma_1 = np.array([[0, 1], [1, 0]], dtype=complex)
        self.sigma_2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
        self.sigma_3 = np.array([[1, 0], [0, -1]], dtype=complex)

        # Definiamo le matrici tau_a
        self.tau_1 = -1j/2 * self.sigma_1
        self.tau_2 = -1j/2 * self.sigma_2
        self.tau_3 = -1j/2 * self.sigma_3

        self.tau = [self.tau_1, self.tau_2, self.tau_3]

    def lie_operator(self, i, a):
        identity = np.eye(2, dtype=complex)
        matrices = [identity, identity, identity, identity]
        matrices[i - 1] = self.tau[a - 1]  # Inserisci la matrice tau corrispondente
        result = matrices[0]
        for mat in matrices[1:]:
            result = np.kron(result, mat)  # Prodotto di Kronecker per costruire la matrice finale
        return result

    def diedral_operator(self, i, j):
        result = np.zeros((16, 16), dtype=complex)
        for a in range(1, 4):
            L_i_a = self.lie_operator(i, a)
            L_j_a = self.lie_operator(j, a)
            result += np.dot(L_i_a, L_j_a)  # Somma dei prodotti L_i_a * L_j_a
        return result

    def apply_operator_to_vector(self, operator_matrix, vector):
        return np.dot(operator_matrix, vector)  # Applica l'operatore al vettore 16x1


    #Qui mettiamo le funzioni per individuare i coefficienti delle matrici diedrali
    def parse_ket_expression(self, ket_str):
        """Parsa una stringa di kets del tipo 'a * |ABCD> + b * |EFGH>'."""
        ket_terms = ket_str.split(' + ')
        parsed_kets = {}

        for term in ket_terms:
            if ' * ' in term:  # Controlla se c'è un coefficiente esplicito
                coeff_str, ket = term.split(' * ')
                coeff = sp.sympify(coeff_str)
            else:
                # Caso del vettore nullo o ket senza coefficiente
                if term.strip() == "0" or term.strip() == "0.+0.j":
                    continue  # Ignora i termini nulli
                else:
                    coeff = sp.sympify(1)  # Se non c'è coefficiente, è implicito 1
                    ket = term

            parsed_kets[ket] = coeff

        if not parsed_kets:  # Se il dizionario è vuoto, ritorna 'null vector'
            return 'null vector'

        return parsed_kets

    def find_combination_coefficients(self, x_str, v1_str, v2_str):
        """Trova i coefficienti a e b tali che x = a*v1 + b*v2."""
        x = self.parse_ket_expression(x_str)

        if x == 'null vector':  # Se x è un vettore nullo, restituisci una stringa appropriata
            return 'null vector'

        v1 = self.parse_ket_expression(v1_str)
        v2 = self.parse_ket_expression(v2_str)

        # Equazioni per risolvere a*v1 + b*v2 = x
        equations = []
        a, b = sp.symbols('a b')

        for ket in x:
            # Ottieni i coefficienti di v1 e v2 per ogni ket
            v1_coeff = v1.get(ket, 0)
            v2_coeff = v2.get(ket, 0)
            x_coeff = x[ket]

            # Imposta l'equazione: a * v1_coeff + b * v2_coeff = x_coeff
            equation = sp.Eq(a * v1_coeff + b * v2_coeff, x_coeff)
            equations.append(equation)

        # Risolvi il sistema di equazioni
        solution = sp.solve(equations, (a, b))

        return solution

    def compose_matrix(self, dict1, dict2):
        """Crea una matrice NumPy che contiene i coefficienti di dict1 e dict2 come colonne."""
        # Chiedi all'utente se desidera contrarre

        # Ottieni i valori da dict1, oppure usa [0] * lunghezza se dict1 è vuoto
        col1 = list(dict1.values()) if dict1 else [0]

        # Ottieni i valori da dict2, oppure usa [0] * lunghezza se dict2 è vuoto
        col2 = list(dict2.values()) if dict2 else [0]

        # Trova la lunghezza massima tra le due colonne
        max_length = max(len(col1), len(col2))

        # Completa le colonne con zeri fino alla lunghezza massima
        col1 += [0] * (max_length - len(col1))
        col2 += [0] * (max_length - len(col2))

        # Componi la matrice come un array NumPy
        matrix = np.array([[col1[i], col2[i]] for i in range(max_length)])

        return matrix

    # Function for nice-displaying the results
    def print_matrix(self, matrix):
        for row in matrix:
            print("  ".join("{:7.3f}".format(float(x)) for x in row))  # Usa float(x) per gestire numeri reali
        print("\n")


    def diedral_matrix(self, indices, basis):
        """Funzione che calcola la matrice di D_ij rispetto alla base data
        :param indices: tupla di indici (i, j)
        :param basis: lista di vettori NON ket nella forma [v1, v2]"""
        i, j = indices[0], indices[1]
        v1, v2 = basis[0], basis[1]
        diedral = self.diedral_operator(i, j)
        D_ijv1, D_ijv2 = self.apply_operator_to_vector(diedral,
                                                                      v1), self.apply_operator_to_vector(
            diedral, v2)
        V1, V2 = ket_2d.vector_to_ket(v1), ket_2d.vector_to_ket(v2)
        D_ijV1, D_ijV2 = ket_2d.vector_to_ket(D_ijv1), ket_2d.vector_to_ket(D_ijv2)

        col1 = self.find_combination_coefficients(D_ijV1, V1, V2)
        col2 = self.find_combination_coefficients(D_ijV2, V1, V2)

        D_ij = self.compose_matrix(col1, col2)

        return D_ij

    def diedral_commutator(self, indices1, indices2, basis):
        """Funzione che calcola il commutatore tra D_ij e D_kl rispetto alla base data
            :param indices1: tupla di indici (i, j)
            :param indices2: tupla di indici (k, l)
            :param basis: lista di vettori NON ket nella forma [v1, v2]"""
        i, j, k, l = indices1[0], indices1[1], indices2[0], indices2[1]
        D_ij = self.diedral_matrix((i, j), basis)
        D_kl = self.diedral_matrix((k, l), basis)

        return np.dot(D_ij, D_kl) - np.dot(D_kl, D_ij)

    def commutatore(self, basis):
        """Funzione che calcola l'operatore volume^2 rispetto ad una data base a meno di una costante 2/9 e un segno
            :param basis: lista di vettori NON ket nella forma [v1, v2]"""

        D_13 = self.diedral_matrix((1, 3), basis)
        D_12 = self.diedral_matrix((1, 2), basis)
        commut = np.dot(D_13, D_12) - np.dot(D_12, D_13)
        #vol = 2/9*np.abs(commut)
        return commut

    def orthogonalize_basis(self, basis):
        """Ortogonalizza una base di vettori utilizzando il processo di Gram-Schmidt.

        :param basis: Lista di vettori (array numpy o SymPy Matrix)
        :return: Base ortogonalizzata come lista di vettori
        """
        orthogonal_basis = []

        for vector in basis:
            # Rimuovi la componente di ciascun vettore già ortogonalizzato
            for ortho_vector in orthogonal_basis:
                projection = np.dot(np.conjugate(ortho_vector), vector) / np.dot(np.conjugate(ortho_vector),
                                                                                 ortho_vector)
                vector -= projection * ortho_vector

            # Aggiungi il nuovo vettore ortogonalizzato alla base, se non nullo
            if np.linalg.norm(vector) > 1e-10:  # Controllo per evitare vettori nulli
                orthogonal_basis.append(vector / np.linalg.norm(vector))  # Normalizza il vettore

        return orthogonal_basis

    def change_basis(self, matrix, indices_list):
        """
        Trasforma una matrice dalla base (v1, v2, ...) alla nuova base (w1, w2, ...)
        dove ogni w_i è una combinazione lineare di v1, v2, ..., vn.

        Qui n = 2

        Args:
        matrix (np.ndarray): La matrice da trasformare, espressa in una certa base.
        indices_list (list): è una lista di due tuple nella forma (a1, b1), (a2, b2).

        Returns:
        np.ndarray: La matrice trasformata nella nuova base (w1, w2).
        """
        v1 = np.array([1, 0])
        v2 = np.array([0, 1])

        #pesco i coefficienti della comb lin per w1 e w2
        a1 = indices_list[0][0]
        b1 = indices_list[0][1]
        a2 = indices_list[1][0]
        b2 = indices_list[1][1]

        w1 = a1 * v1 + b1 * v2
        w2 = a2 * v1 + b2 * v2

        # Costruisce la matrice di cambiamento di base
        change_of_basis_matrix = np.column_stack((w1, w2))

        # Inverte la matrice di cambiamento di base
        change_of_basis_matrix_inv = np.linalg.inv(change_of_basis_matrix)

        # Trasforma la matrice nella nuova base
        #transformed_matrix = change_of_basis_matrix_inv @ matrix @ change_of_basis_matrix
        # @ per il prodotto matriciale non funziona bene in jupyter quindi facciamo così
        transformed_matrix = np.dot(np.dot(change_of_basis_matrix_inv, matrix), change_of_basis_matrix)

        return transformed_matrix



class GeometricOperators1:
    def __init__(self):

        self.tau1 = np.array([[0, -1j * math.sqrt(2) / 2, 0], [-1j * math.sqrt(2) / 2, 0, -1j * math.sqrt(2) / 2],
                         [0, -1j * math.sqrt(2) / 2, 0]], dtype=complex)
        self.tau2 = np.array([[0, -math.sqrt(2) / 2, 0], [math.sqrt(2) / 2, 0, -math.sqrt(2) / 2], [0, math.sqrt(2) / 2, 0]],
                        dtype=complex)
        self.tau3 = np.array([[-1j, 0, 0], [0, 0, 0], [0, 0, 1j]], dtype=complex)

        self.tau = [self.tau1, self.tau2, self.tau3]


    def lie_operator(self, i, a):
        identity = np.eye(3, dtype=complex)
        matrices = [identity, identity, identity, identity]
        matrices[i - 1] = self.tau[a - 1]  # Inserisci la matrice tau corrispondente
        result = matrices[0]
        for mat in matrices[1:]:
            result = np.kron(result, mat)  # Prodotto di Kronecker per costruire la matrice finale
        return result

    def diedral_operator(self, i, j):
        result = np.zeros((81, 81), dtype=complex)
        for a in range(1, 4):
            L_i_a = self.lie_operator(i, a)
            L_j_a = self.lie_operator(j, a)
            result += np.dot(L_i_a, L_j_a)  # Somma dei prodotti L_i_a * L_j_a
        return result

    def apply_operator_to_vector(self, operator_matrix, vector):
        return np.dot(operator_matrix, vector)  # Applica l'operatore al vettore 81x1

    # Qui mettiamo le funzioni per individuare i coefficienti delle matrici diedrali
    def parse_ket_expression(self, ket_str):
        """Parsa una stringa di kets del tipo 'a * |ABCD> + b * |EFGH>'."""
        ket_terms = ket_str.split(' + ')
        parsed_kets = {}

        for term in ket_terms:
            if ' * ' in term:  # Controlla se c'è un coefficiente esplicito
                coeff_str, ket = term.split(' * ')
                coeff = sp.sympify(coeff_str)
            else:
                # Caso del vettore nullo o ket senza coefficiente
                if term.strip() == "0" or term.strip() == "0.+0.j":
                    continue  # Ignora i termini nulli
                else:
                    coeff = sp.sympify(1)  # Se non c'è coefficiente, è implicito 1
                    ket = term

            parsed_kets[ket] = coeff

        if not parsed_kets:  # Se il dizionario è vuoto, ritorna 'null vector'
            return 'null vector'

        return parsed_kets


    def find_combination_coefficients(self, x_str, v1_str, v2_str, v3_str):
        """Trova i coefficienti a, b e c tali che x = a*v1 + b*v2 + c*v3."""
        x = self.parse_ket_expression(x_str)
        #è un dict adesso
        if x == 'null vector':  # Se x è un vettore nullo, restituisci una stringa appropriata
            return 'null vector'

        v1 = self.parse_ket_expression(v1_str)
        v2 = self.parse_ket_expression(v2_str)
        v3 = self.parse_ket_expression(v3_str)

        # Equazioni per risolvere a*v1 + b*v2 = x
        equations = []
        a, b, c = sp.symbols('a b c')

        for ket in x:
            # Ottieni i coefficienti di v1 e v2 per ogni ket
            v1_coeff = v1.get(ket, 0)
            v2_coeff = v2.get(ket, 0)
            v3_coeff = v3.get(ket, 0)
            x_coeff = x[ket]

            # Imposta l'equazione: a * v1_coeff + b * v2_coeff = x_coeff
            equation = sp.Eq(a * v1_coeff + b * v2_coeff + c * v3_coeff, x_coeff)
            equations.append(equation)

        # Risolvi il sistema di equazioni
        solution = sp.solve(equations, (a, b, c))

        return solution

    def compose_matrix(self, dict1, dict2, dict3):
        """Crea una matrice NumPy che contiene i coefficienti di dict1, dict2 e dict3 come colonne."""

        # Ottieni i valori da dict1, oppure usa [0] * lunghezza se dict1 è vuoto
        col1 = list(dict1.values()) if dict1 else [0]

        # Ottieni i valori da dict2, oppure usa [0] * lunghezza se dict2 è vuoto
        col2 = list(dict2.values()) if dict2 else [0]

        col3 = list(dict3.values()) if dict2 else [0]

        # Trova la lunghezza massima tra le due colonne
        max_length = max(len(col1), len(col2), len(col3))

        # Completa le colonne con zeri fino alla lunghezza massima
        col1 += [0] * (max_length - len(col1))
        col2 += [0] * (max_length - len(col2))
        col3 += [0] * (max_length - len(col3))

        # Componi la matrice come un array NumPy
        matrix = np.array([[col1[i], col2[i], col3[i]] for i in range(max_length)])

        return matrix

    # Function for nice-displaying the results
    def print_matrix(self, matrix):
        for row in matrix:
            print("  ".join("{:7.3f}".format(float(x)) for x in row))  # Usa float(x) per gestire numeri reali
        print("\n")

    def diedral_matrix(self, indices, basis):
        """Funzione che calcola la matrice di D_ij rispetto alla base data
        :param indices: tupla di indici (i, j)
        :param basis: lista di vettori NON ket nella forma [v1, v2]"""
        i, j = indices[0], indices[1]
        v1, v2, v3 = basis[0], basis[1], basis[2]
        diedral = self.diedral_operator(i, j)
        D_ijv1, D_ijv2, D_ijv3 = self.apply_operator_to_vector(diedral,
                                                       v1), self.apply_operator_to_vector(
            diedral, v2), self.apply_operator_to_vector(diedral, v3)
        V1, V2, V3 = ket_3d.vector_to_ket(v1), ket_3d.vector_to_ket(v2), ket_3d.vector_to_ket(v3)
        D_ijV1, D_ijV2, D_ijV3 = ket_3d.vector_to_ket(D_ijv1), ket_3d.vector_to_ket(D_ijv2), ket_3d.vector_to_ket(D_ijv3)

        col1 = self.find_combination_coefficients(D_ijV1, V1, V2, V3)
        col2 = self.find_combination_coefficients(D_ijV2, V1, V2, V3)
        col3 = self.find_combination_coefficients(D_ijV3, V1, V2, V3)

        D_ij = self.compose_matrix(col1, col2, col3)

        return D_ij

