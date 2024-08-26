# In this file we try to handle the case of the quantum tetrahedra with 4 spin=1
# which corresponds to the spin network described by the representation {rho^1}^{\otimes4}

import numpy as np

from LieOperators import lie_operator, squared_lie_operator, kronecker_product_4, print_matrix

from DiedralOperators import diedral_operator

# By noticing that the support space of rho^1 is C^3, we shall start by defining

e1 = np.array([1, 0, 0])
e2 = np.array([0, 1, 0])
e3 = np.array([0, 0, 1])

# now compute the 3^4=81 vectors spanning {C^3}^{\otimes4}

basis_vectors = [(e1, e1, e1, e1), (e1, e1, e1, e2), (e1, e1, e1, e3),
                 (e1, e1, e2, e1), (e1, e1, e2, e2), (e1, e1, e2, e3),
                 (e1, e1, e3, e1), (e1, e1, e3, e2), (e1, e1, e3, e3),
                 (e1, e2, e1, e1), (e1, e2, e1, e2), (e1, e2, e1, e3),
                 (e1, e2, e2, e1), (e1, e2, e2, e2), (e1, e2, e2, e3),
                 (e1, e2, e3, e1), (e1, e2, e3, e2), (e1, e2, e3, e3),
                 (e1, e3, e1, e1), (e1, e3, e1, e2), (e1, e3, e1, e3),
                 (e1, e3, e2, e1), (e1, e3, e2, e2), (e1, e3, e2, e3),
                 (e1, e3, e3, e1), (e1, e3, e3, e2), (e1, e3, e3, e3),
                 (e2, e1, e1, e1), (e2, e1, e1, e2), (e2, e1, e1, e3),
                 (e2, e1, e2, e1), (e2, e1, e2, e2), (e2, e1, e2, e3),
                 (e2, e1, e3, e1), (e2, e1, e3, e2), (e2, e1, e3, e3),
                 (e2, e2, e1, e1), (e2, e2, e1, e2), (e2, e2, e1, e3),
                 (e2, e2, e2, e1), (e2, e2, e2, e2), (e2, e2, e2, e3),
                 (e2, e2, e3, e1), (e2, e2, e3, e2), (e2, e2, e3, e3),
                 (e2, e3, e1, e1), (e2, e3, e1, e2), (e2, e3, e1, e3),
                 (e2, e3, e2, e1), (e2, e3, e2, e2), (e2, e3, e2, e3),
                 (e2, e3, e3, e1), (e2, e3, e3, e2), (e2, e3, e3, e3),
                 (e3, e1, e1, e1), (e3, e1, e1, e2), (e3, e1, e1, e3),
                 (e3, e1, e2, e1), (e3, e1, e2, e2), (e3, e1, e2, e3),
                 (e3, e1, e3, e1), (e3, e1, e3, e2), (e3, e1, e3, e3),
                 (e3, e2, e1, e1), (e3, e2, e1, e2), (e3, e2, e1, e3),
                 (e3, e2, e2, e1), (e3, e2, e2, e2), (e3, e2, e2, e3),
                 (e3, e2, e3, e1), (e3, e2, e3, e2), (e3, e2, e3, e3),
                 (e3, e3, e1, e1), (e3, e3, e1, e2), (e3, e3, e1, e3),
                 (e3, e3, e2, e1), (e3, e3, e2, e2), (e3, e3, e2, e3),
                 (e3, e3, e3, e1), (e3, e3, e3, e2), (e3, e3, e3, e3)]

# Calculate the 81 canonical basis vectors of C^3\otimes4
vectors = [kronecker_product_4(*vectors) for vectors in basis_vectors]

# Display the vectors
for i, vector in enumerate(vectors):
    print("Vector {}:\n{}\n".format(i+1, vector))

#ALTERNATIVELY, we can use a generalized method based on the following function
from itertools import product

def generate_combinations(elements, n):
    """
    Generate all possible combinations of length n from the given list of elements.

    :param elements: List of elements to combine
    :param n: Length of each combination
    :return: List of tuples representing all possible combinations of length n
    """
    return list(product(elements, repeat=n))

#through

#e1 = np.array([1, 0, 0])
#e2 = np.array([0, 1, 0])
#e3 = np.array([0, 0, 1])

elements = [e1, e2, e3]
basis_vectors = generate_combinations(elements, 4)

vectors = [kronecker_product_4(*vectors) for vectors in basis_vectors]

print()
# Display the vectors
for i, vector in enumerate(vectors):
    print("Vector {}:\n{}\n".format(i+1, vector))

import math

#Let us construct the framework for Lie (then, diedral) operators in the spin 1 representation
tau1 = np.array([[0,-1j*math.sqrt(2)/2,0],[-1j*math.sqrt(2)/2,0,-1j*math.sqrt(2)/2],[0,-1j*math.sqrt(2)/2,0]], dtype=complex)
tau2 = np.array([[0,-math.sqrt(2)/2,0],[math.sqrt(2)/2,0,-math.sqrt(2)/2],[0,math.sqrt(2)/2,0]], dtype=complex)
tau3 = np.array([[-1j,0,0],[0,0,0],[0,0,1j]], dtype=complex)

tau = [tau1, tau2, tau3]

#Write down the 3-dimensional isotropic subspace of the spin 1 tetrahedron

def find_index(t, basis_vectors):
    """
    Trova l'indice di una tupla t all'interno della lista basis_vectors.

    Args:
    t: tupla da cercare.
    basis_vectors: lista di tuple.

    Returns:
    L'indice della tupla t all'interno di basis_vectors se presente, altrimenti -1.
    """
    try:
        return basis_vectors.index(t)
    except ValueError:
        return -1
#l'idea e' buona, ma possiamo fare a mano facendo direttamente i prodotti

#V1 = (e1,e3,e1,e3) + (e3,e1,e3,e1) - (e1,e3,e2,e2) - (e2,e2,e1,e3) - (e2,e2,e3,e1) - (e3,e1,e2,e2) + (e1,e3,e3,e1) + (e2,e2,e2,e2) + (e3,e1,e1,e3)
#V2 = (e1,e1,e3,e3) + (e3,e3,e1,e1) - (e1,e2,e2,e3) - (e2,e1,e3,e2) - (e2,e3,e1,e2) - (e3,e2,e2,e1) + (e1,e3,e1,e3) + (e3,e1,e3,e1) + (e2,e2,e2,e2)
#V3 = (e1,e2,e2,e3) + (e3,e2,e2,e1) - (e2,e2,e2,e2) - (e2,e1,e2,e3) - (e1,e3,e1,e3) - (e1,e2,e3,e2) - (e2,e3,e2,e1) - (e3,e1,e3,e1) - (e3,e2,e1,e2) + (e3,e1,e2,e2) + (e2,e3,e1,e2) + (e2,e2,e3,e1) + (e2,e2,e1,e3) + (e2,e1,e3,e2) + (e1,e3,e2,e2)

V1 = (kronecker_product_4(e1, e3, e1, e3) + kronecker_product_4(e3, e1, e3, e1) -
      kronecker_product_4(e1, e3, e2, e2) - kronecker_product_4(e2, e2, e1, e3) -
      kronecker_product_4(e2, e2, e3, e1) - kronecker_product_4(e3, e1, e2, e2) +
      kronecker_product_4(e1, e3, e3, e1) + kronecker_product_4(e2, e2, e2, e2) +
      kronecker_product_4(e3, e1, e1, e3))

V2 = (kronecker_product_4(e1, e1, e3, e3) + kronecker_product_4(e3, e3, e1, e1) -
      kronecker_product_4(e1, e2, e2, e3) - kronecker_product_4(e2, e1, e3, e2) -
      kronecker_product_4(e2, e3, e1, e2) - kronecker_product_4(e3, e2, e2, e1) +
      kronecker_product_4(e1, e3, e1, e3) + kronecker_product_4(e3, e1, e3, e1) +
      kronecker_product_4(e2, e2, e2, e2))

V3 = (kronecker_product_4(e1, e2, e2, e3) + kronecker_product_4(e3, e2, e2, e1) -
      kronecker_product_4(e2, e2, e2, e2) - kronecker_product_4(e2, e1, e2, e3) -
      kronecker_product_4(e1, e3, e1, e3) - kronecker_product_4(e1, e2, e3, e2) -
      kronecker_product_4(e2, e3, e2, e1) - kronecker_product_4(e3, e1, e3, e1) -
      kronecker_product_4(e3, e2, e1, e2) + kronecker_product_4(e3, e1, e2, e2) +
      kronecker_product_4(e2, e3, e1, e2) + kronecker_product_4(e2, e2, e3, e1) +
      kronecker_product_4(e2, e2, e1, e3) + kronecker_product_4(e2, e1, e3, e2) +
      kronecker_product_4(e1, e3, e2, e2))

#Now we gotta define Lie and diedral operators adapted to the support space of rho^1\otimes4

def lie_operator1(i, a):
    identity = np.eye(3, dtype=complex)
    matrices = [identity, identity, identity, identity]
    matrices[i - 1] = tau[a - 1]
    result = matrices[0]
    for mat in matrices[1:]:
        result = np.kron(result, mat)
    return result

def diedral_operator1(i, j):
    result = np.zeros((81, 81), dtype=complex)
    for a in range(1, 4):
        L_i_a = lie_operator1(i, a)
        L_j_a = lie_operator1(j, a)
        result += np.dot(L_i_a, L_j_a)
    return result


V = np.array([V1, V2, V3])
print_matrix(V)
print(V1)
print()
print(V2)
print()
print(V3)

print(np.linalg.matrix_rank(V))

def normalize_vector(v):
    """Normalizza un vettore"""
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm

def gram_schmidt(vectors):
    """Ortonormalizza un insieme di vettori usando il processo di Gram-Schmidt"""
    orthonormal_vectors = []
    for v in vectors:
        w = v - sum(np.dot(v, b) * b for b in orthonormal_vectors)
        w_normalized = normalize_vector(w)
        if np.linalg.norm(w_normalized) > 1e-10:  # Aggiunge solo vettori non nulli
            orthonormal_vectors.append(w_normalized)
    return orthonormal_vectors

def gram_schmidt_qr(vectors):
    """
    Ortonormalizza una lista di vettori utilizzando la decomposizione QR di NumPy.

    Parametri:
    vectors - una lista di vettori (ogni vettore e' un array NumPy)

    Restituisce:
    Una matrice NumPy di vettori ortogonali
    """
    # Converte la lista di vettori in una matrice NumPy
    matrix = np.vstack(vectors)

    # Applica la decomposizione QR
    Q, _ = np.linalg.qr(matrix.T)

    # Ritorna la matrice dei vettori ortogonali
    return Q.T


def restrizione_matrice_sottospazio(A, vettori_base):
    """
    Restringe una matrice simmetrica A ad un sottospazio invariante definito dai vettori base.
    Per ottenere un risultato coerente i vettori devono essere ortonogonali, poiche' si sta usando
    che l'inversa di P coincide con la trasposta, ma questo e' vero se e solo se P e' ortogonale

    :param A: np.ndarray, matrice simmetrica n x n
    :param vettori_base: list of np.ndarray, lista di vettori che definiscono il sottospazio invariante
    :return: np.ndarray, matrice ristretta al sottospazio invariante
    """
    # Costruisci la matrice di cambiamento di base P
    P = np.column_stack(vettori_base)

    #Verifica che P abbia rango pieno
    if np.linalg.matrix_rank(P) != P.shape[1]:
        raise ValueError("I vettori base non formano una base indipendente.")

    # Estrai la sottomatrice corrispondente al sottospazio invariante
    k = len(vettori_base)
    A_restr = np.dot(np.dot(P.T, A), P)
    A_restr = A_restr[:k, :k]

    return A_restr


orth_basis = gram_schmidt(V)
print(np.linalg.matrix_rank(orth_basis))

D_12_restr = restrizione_matrice_sottospazio(diedral_operator1(1,2), orth_basis)
print("Matrice D_12 ristretta alla proiezione ortogonale del sottospazio invariante:")
print_matrix(D_12_restr)


D_13_restr = restrizione_matrice_sottospazio(diedral_operator1(1,3), orth_basis)
print("Matrice D_13 ristretta alla proiezione ortogonale del sottospazio invariante:")
print_matrix(D_13_restr)

#orth_basis2 = gram_schmidt(V)

#A_restr = restrizione_matrice_sottospazio(diedral_operator1(1,2), orth_basis2)
#print "Matrice D_12 ristretta alla proiezione ortogonale del sottospazio invariante, altro metodo:"
#print_matrix(A_restr)
#print()
#A_restr = restrizione_matrice_sottospazio(diedral_operator1(1,3), orth_basis2)
#print "Matrice D_13 ristretta alla proiezione ortogonale del sottospazio invariante, altro metodo:"
#print_matrix(A_restr)

#vengono uguali
print('Volume quadrato:')
Volume2 = np.abs(np.dot(D_13_restr,D_12_restr) - np.dot(D_12_restr,D_13_restr))
print_matrix(Volume2)


def print_multiples_of_sqrt3_over_2k_and_sqrt2_over_2k(k_values):
    """
    Stampa i multipli di sqrt(3)/2k e sqrt(2)/2k per i valori di k dati.

    Parametri:
    k_values - una lista di valori k per cui calcolare i multipli
    """
    for k in k_values:
        # Calcola i multipli di sqrt(3)/2k
        sqrt3_over_2k = (k * math.sqrt(3)) / 2

        # Calcola i multipli di sqrt(2)/2k
        sqrt2_over_2k = (k * math.sqrt(2)) / 2

        # Stampa i risultati con approssimazione decimale
        print("Per k =", k)
        print("Multipli di k sqrt(3)/2:", sqrt3_over_2k)
        print("Multipli di k sqrt(2)/2:", sqrt2_over_2k)
        print()


# Valori di k
k_values = [0.16, 0.2, 0.25, 0.33, 0.5, 1, 2, 3, 4, 5, 6, 8]

# Stampa i multipli
print_multiples_of_sqrt3_over_2k_and_sqrt2_over_2k(k_values)


print(float(27*18 + 20*6 + 28*9 + 31*18 + 29*6 + 27*6 + 26*12 + 30*6 + 21*6)/(6*7 + 9*5))
print(float(27.2*110/30))

Volume2_Wolfram = np.abs(np.array([[0, 3/4-1/(4*math.sqrt(2))-math.sqrt(2), 1/4*(3*math.sqrt(3/2))-math.sqrt(3)],[-3/4+1/(4*math.sqrt(2))+math.sqrt(2),0,1/8+1/(4*math.sqrt(2))-math.sqrt(3)/4],[-1/4*(3*math.sqrt(3/2))+math.sqrt(3), -1/8-1/(4*math.sqrt(2))+math.sqrt(3)/4, 0]],dtype=complex))

print_matrix(Volume2_Wolfram)

#print approximation of V^2

D_12_approx = np.array([[2,0,0],[0,-1/4,math.sqrt(3)/4],[0,math.sqrt(3)/4,math.sqrt(2)/4]],dtype=complex)

D_13_approx = np.array([[0,math.sqrt(2)/2,math.sqrt(3)/2],[math.sqrt(2)/2,1/2,1/2],[math.sqrt(3)/2,1/2,1]],dtype=complex)

print_matrix(np.abs(np.dot(D_13_restr,D_12_approx)-np.dot(D_12_approx,D_13_approx)))
