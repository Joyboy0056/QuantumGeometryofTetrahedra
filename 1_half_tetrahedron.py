import numpy as np

from LieOperators import lie_operator, squared_lie_operator, kronecker_product_4, print_matrix

from DiedralOperators import diedral_operator, D_12, D_13

# Define the basis vectors of C^2
e1 = np.array([1, 0])
e2 = np.array([0, 1])


# Function to calculate the Kronecker product of four vectors
def kronecker_product_4(v1, v2, v3, v4):
    return np.kron(np.kron(np.kron(v1, v2), v3), v4)


# All combinations of (v1, v2) repeated 4 times
# it spans a 16 dimensional vector space
basis_vectors = [(e1, e1, e1, e1), (e1, e1, e1, e2), (e1, e1, e2, e1), (e1, e1, e2, e2),
                 (e1, e2, e1, e1), (e1, e2, e1, e2), (e1, e2, e2, e1), (e1, e2, e2, e2),
                 (e2, e1, e1, e1), (e2, e1, e1, e2), (e2, e1, e2, e1), (e2, e1, e2, e2),
                 (e2, e2, e1, e1), (e2, e2, e1, e2), (e2, e2, e2, e1), (e2, e2, e2, e2)]

# Calculate the 16 canonical basis vectors of C^2\otimes4
vectors = [kronecker_product_4(*vectors) for vectors in basis_vectors]

# Display the vectors
# for i, vector in enumerate(vectors):
#    print("Vector {}:\n{}\n".format(i+1, vector))



# Calola la matrice di D_ij rispetto a una base isotropa con

def matrix_coeffs(mat,basis):
    n = len(basis)
    res = np.zeros((n, n),dtype=complex)
    for i in range(0, n):
        for j in range(0, n):
            res[i, j] = np.dot(basis[i], np.dot(mat, basis[j]))
    return res

#ma prima devi cambiare base: Cioe' devi mettere D_ij nella base v1,v2, perche la funzione diedralOps li calcola nella base canonica

def cambio_base(operatore, base_vec):
    # Calcola l'inversa della matrice di cambio di base
    cambio_base_inv = np.linalg.inv(base_vec)

    # Calcola la matrice nell'altro sistema di base
    operatore_nuova_base = np.dot(np.dot(cambio_base_inv, operatore), base_vec)

    return operatore_nuova_base

# Calcola D_{12} ristretto al sottospazio invariante Inv(rho^1/2 \otimes4), che e' bidimensionale generato da

v1 = vectors[5] - vectors[6] - vectors[9] + vectors[10]
v2 = vectors[3] - vectors[5] - vectors[10] + vectors[12]


#-------------------------------------------------------------------------------------------

D_12_v1 = np.dot(D_12, v1)

print(D_12_v1)
print()

D_12_v2 = np.dot(D_12, v2)

print(D_12_v2)
print()
print(D_12_v1 == v1 * 0.75)
print(D_12_v2 == (v1 * -0.5 + v2 * -0.25))
print()

print_matrix(matrix_coeffs(D_12,[v1,v2]))

D_13_v1 = np.dot(D_13, v1)

print(D_13_v1)
print()

D_13_v2 = np.dot(D_13, v2)

print(D_13_v2)
print()
print(D_13_v1 == v1 * 0.25 + v2 * 0.5)
print(D_13_v2 == (v1 * 0.5 + v2 * 0.25))
print()

print_matrix(matrix_coeffs(D_13,[v1,v2]))


# Now: once the matrix of D_12 on (v1,v2) is given, we diagonalize it and we get
e_1 = v1
e_2 = v1 + v2 * 2
# as eigenvectors respectively of 3/4 and -1/4. Hence, the matrix D_12 on (e1,e2) should be diag(3/4,-1/4)

D_12_e1 = np.dot(D_12, e_1)
D_12_e2 = np.dot(D_12, e_2)
print()
print(D_12_e1 == e_1 * 0.75)
print(D_12_e2 == e_2 * -0.25)
print('VERIFICA')

print_matrix(matrix_coeffs(D_12,[e_1,e_2]))

# on the other hand, the matrix of D_13 on the eigenbasis of D_12 should not be diagonal

D_13_e1 = np.dot(D_13, e_1)
D_13_e2 = np.dot(D_13, e_2)
print()

print(D_13_e1 == e_2 * 0.25)
print(D_13_e2 == e_1 * 0.75 + e_2 * 0.5)
print()
print(D_13_e1)
print(D_13_e2)
print()
print(e_1)
print(e_2)

print()

# We can see that the basis v1, v2 is not orthonormal. It can be made orthonormal through

import math

V2 = vectors[5] - vectors[6] - vectors[9] + vectors[10]
V1 = vectors[3] - vectors[5] - vectors[10] + vectors[12]

print(V2)
print(V1)

E1 = V1*0.5
E2 = 1/math.sqrt(3)*V2-1/(2*math.sqrt(3))*V1

print()
print(E2==np.array([0,0,0,1/math.sqrt(3),0,-0.5*1/math.sqrt(3),0,0,-0.5*1/math.sqrt(3),-0.5*1/math.sqrt(3),0,1/math.sqrt(3),0,0,0]))
print()

print(np.dot(diedral_operator(1,2),E1)==(-math.sqrt(3)/4)*E2)
print(np.dot(diedral_operator(1,2),E2)==E2*0.5+E1*-(math.sqrt(3)/4)) #corretto

print()


print(np.dot(diedral_operator(1,3),E1)==-math.sqrt(3)/4*E2)
print(np.dot(diedral_operator(1,3),E2)==E2*0.5+E1*(math.sqrt(3)/4)) #corretto

print()

D_13_E1=np.dot(D_13,E1)
D_13_E2=np.dot(D_13,E2)

D_12_E1=np.dot(D_12,E1)
D_12_E2=np.dot(D_12,E2)


#NOTA! Qui v1 e v2 sono invertiti in V2 e V1

print()
print('Da qui in poi mettiamo a display i risultati facendo i conti direttamente con le matrici')
print()

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

orth_basis = gram_schmidt_qr([v1,v2])

A_restr = restrizione_matrice_sottospazio(diedral_operator(1,2), orth_basis)
print("Matrice D_12 ristretta alla proiezione ortogonale del sottospazio invariante:")
print_matrix(A_restr)
print()
A_restr = restrizione_matrice_sottospazio(diedral_operator(1,3), orth_basis)
print("Matrice D_13 ristretta alla proiezione ortogonale del sottospazio invariante:")
print_matrix(A_restr)


print()
print('Perfezioniamo il discorso:')

v1 = (kronecker_product_4(e1, e2, e1, e2) - kronecker_product_4(e1, e2, e2, e1) - kronecker_product_4(e2, e1, e1, e2) + kronecker_product_4(e2, e1, e2, e1))
v2 = (kronecker_product_4(e1, e1, e2, e2) - kronecker_product_4(e1, e2, e1, e2) - kronecker_product_4(e2, e1, e2, e1) + kronecker_product_4(e2, e2, e1, e1))

orth_basis = gram_schmidt([v1,v2])

D_12_restr = restrizione_matrice_sottospazio(diedral_operator(1,2), orth_basis)
D_13_restr = restrizione_matrice_sottospazio(diedral_operator(1,3),orth_basis)

print("Matrice D_12 ristretta alla proiezione ortogonale del sottospazio invariante:")
print_matrix(D_12_restr)

print("Matrice D_13 ristretta alla proiezione ortogonale del sottospazio invariante:")
print_matrix(D_13_restr)

print("Operatore volume quadrato ristretto ala proiezione ortogonale del sottospazio invariante:")
V_squared = np.dot(D_13_restr, D_12_restr) - np.dot(D_12_restr, D_13_restr)
print_matrix(V_squared)


eigenvectors = [orth_basis[1] - orth_basis[0], orth_basis[0] + orth_basis[1]]

print("Diagonal form of the volume squared in its eigenvectors' basis, up to a constant:")
print_matrix(restrizione_matrice_sottospazio(np.abs(np.dot(D_13,D_12)-np.dot(D_12,D_13)),orth_basis))

print("Non-diagonal form of the non-compatible dihedral angle:")
print_matrix(restrizione_matrice_sottospazio(D_12,eigenvectors))
