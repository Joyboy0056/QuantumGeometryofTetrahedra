import numpy as np

# Let us construct the multilinear framework
# for the ground state tetrahedron of spin (1/2,1/2,1/2,1/2)

# -----------------------------------------------------------

# Define the basis vectors of C^2
e1 = np.array([1, 0])
e2 = np.array([0, 1])

# Function to calculate the Kronecker product of four vectors
def kronecker_product_4(v_1, v_2, v_3, v_4):
    return np.kron(np.kron(np.kron(v_1, v_2), v_3), v_4)


# All combinations of (e1, e2) repeated 4 times
# it spans a 16 dimensional vector space
basis_vectors = [(e1, e1, e1, e1), (e1, e1, e1, e2), (e1, e1, e2, e1), (e1, e1, e2, e2),
                 (e1, e2, e1, e1), (e1, e2, e1, e2), (e1, e2, e2, e1), (e1, e2, e2, e2),
                 (e2, e1, e1, e1), (e2, e1, e1, e2), (e2, e1, e2, e1), (e2, e1, e2, e2),
                 (e2, e2, e1, e1), (e2, e2, e1, e2), (e2, e2, e2, e1), (e2, e2, e2, e2)]

# Compute the 16 canonical basis vectors of C^2\otimes4
vectors = [kronecker_product_4(*vectors) for vectors in basis_vectors]

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

# e.g. vectors = generate_combinations([e1, e2],4)-------------------------------------


# Let us eventually compute Inv(rho^1/2 \otimes4), being 2-dimensional and spanned by

v1 = vectors[5] - vectors[6] - vectors[9] + vectors[10]
v2 = vectors[3] - vectors[5] - vectors[10] + vectors[12]

#-------------------------------------------------------------------------------------
# Once the framework is complete, we can start discussing the case of the maximal
# commuting sub-algebra of the angle D_12 and the four areas from
# the geometric algebra (V^2, D_00, D_11, D_22, D_33, D_12) for this spin tetrahedron
#-------------------------------------------------------------------------------------


# Let us define our fundamental auxiliar functions

# Function computing the Lie operator (L_i)_a, for any i=0,1,2,3 and a=1,2,3
def lie_operator(i, a):
    # Define Pauli matrices
    sigma_1 = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma_3 = np.array([[1, 0], [0, -1]], dtype=complex)

    # Define a different basis of su(2): the tau matrices
    tau_1 = -1j / 2 * sigma_1
    tau_2 = -1j / 2 * sigma_2
    tau_3 = -1j / 2 * sigma_3

    tau = [tau_1, tau_2, tau_3]

    # Lie operator construction
    identity = np.eye(2, dtype=complex)
    matrices = [identity, identity, identity, identity]
    matrices[i - 1] = tau[a - 1]
    result = matrices[0]
    for mat in matrices[1:]:
        result = np.kron(result, mat)
    return result

# Function for nice-displaying the results
def print_matrix(matrix):
    for row in matrix:
        print("  ".join("{:7.3f}".format(x.real) + "{:+7.3f}j".format(x.imag) for x in row))
    print("\n")

# Function computing the dihedral operator D_jj, for any j=0,1,2,3
def diedral_operator(i, j):
    result = np.zeros((16, 16), dtype=complex)
    for a in range(1, 4):
        L_i_a = lie_operator(i, a)
        L_j_a = lie_operator(j, a)
        result += np.dot(L_i_a, L_j_a)
    return result

# Functions for orthonormalize a set of vectors


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

#-------------------------------------------------------------------------------------------

# First, being our operators symmetric 2x2 matrices, we orthogonalize the basis (v1, v2)

orth_basis = gram_schmidt([v1, v2]) # The output is a list containing E1, E2 - from the theory.


# Up to now, dihedral ops (as well as Lie ones) are defined in the whole product space,
# support of rho, hence they result as 16 dimensional square matrices. For that, we must
# select the correct restriction of them to the 2-dimensional isotropic subspace spanned by E1, E2.

def matrix_to_subspace(A, basis_vectors):
    """
    It restricts a symmetric matrix A to an invariant subspace spanned by basis_vectors.
    In order to get a coherent result, such vectors must be orthogonal in R^3, the property of
    the inverse of P coincide with the transpose being used.

    :param A: np.ndarray, symmetric matrix n x n
    :param basis_vectors: list of np.ndarray, vectors defining the invariant subspace stored in a list
    :return: np.ndarray, matrix restricted to the invariant subspace
    """
    # Construct the matrix of change of basis
    P = np.column_stack(basis_vectors)

    # Check that P has max rank
    if np.linalg.matrix_rank(P) != P.shape[1]:
        raise ValueError("Basis vectors do not provide an independent basis.")

    # Extract the sub-matrix in correspondence of the invariant subspace
    k = len(basis_vectors)
    res = np.dot(np.dot(P.T, A), P)
    res = res[:k, :k]

    return res


# LET'S START DO OUR COMPUTATIONS
# We compute our operators D_12, D_00,..., D_33 and we show that they simultaneously diagonalize in "orth_basis"
# while the volume [D_13, D_12] (up some constant) does not.

print("\nWe are showing that dihedral operators (D_12, D_00,..., D_33) form a maximal commuting sub-algebra of "
      "(D_12, D_00, D_11, D_22, D_33, V^2),\nin the sense that they simultaneously diagonalize (meaning, with respect to"
      " the same basis), leaving the volume non diagonal.\n")

D12 = matrix_to_subspace(diedral_operator(1,2),orth_basis)

print("Angle operator D_12 restricted to the basis (E1, E2):\n")
print_matrix(D12)

D00 = matrix_to_subspace(diedral_operator(0,0),orth_basis)

print("Area operator D_00 restricted to the basis (E1, E2):\n")
print_matrix(D00)

D11 = matrix_to_subspace(diedral_operator(1,1),orth_basis)

print("Area operator D_11 restricted to the basis (E1, E2):\n")
print_matrix(D11)

D22 = matrix_to_subspace(diedral_operator(2,2),orth_basis)

print("Area operator D_22 restricted to the basis (E1, E2):\n")
print_matrix(D22)

D33 = matrix_to_subspace(diedral_operator(3,3),orth_basis)

print("Area operator D_33 restricted to the basis (E1, E2):\n")
print_matrix(D33)

#aux
D13 = matrix_to_subspace(diedral_operator(1,3), orth_basis)
Volume2 = np.abs(np.dot(D13,D12) - np.dot(D12,D13))


print("Volume operator V^2 restricted to the basis (E1, E2):\n")
print_matrix(Volume2)

print("CVD\n")