import numpy as np

from LieOperators import lie_operator, squared_lie_operator, print_matrix

# Definisci le matrici di Pauli
sigma_1 = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_3 = np.array([[1, 0], [0, -1]], dtype=complex)

# Definisci le matrici tau_a
tau_1 = -1j/2 * sigma_1
tau_2 = -1j/2 * sigma_2
tau_3 = -1j/2 * sigma_3

# Lista delle matrici tau
tau = [tau_1, tau_2, tau_3]

# Funzione per calcolare l'operatore di Lie L_{i,a}

#def lie_operator(i, a):
#    identity = np.eye(2, dtype=complex)
 #   matrices = [identity, identity, identity, identity]
  #  matrices[i - 1] = tau[a - 1]
   # result = matrices[0]
    #for mat in matrices[1:]:
     #   result = np.kron(result, mat)
    #return result


#diedral operators

def diedral_operator(i, j):
    result = np.zeros((16, 16), dtype=complex)
    for a in range(1, 4):
        L_i_a = lie_operator(i, a)
        L_j_a = lie_operator(j, a)
        result += np.dot(L_i_a, L_j_a)
    return result

# Esempio di uso: calcola D_{1, 2}
D_12 = diedral_operator(1, 2)
print("Dimensione di D_{1, 2} =", D_12.shape)
print("D_{12} =")
print_matrix(D_12)

D_13=diedral_operator(1,3)


print()


#Notice that, the squared area operators D_00, D_11, D_22, D_33 are computed as squared Lie operators

print(squared_lie_operator(4)==diedral_operator(4,4))

#hence, actually, the funcion "squared_lie" is useless

print_matrix(diedral_operator(1,2))

print_matrix(diedral_operator(3,3))

# Try through orthogonalisation


