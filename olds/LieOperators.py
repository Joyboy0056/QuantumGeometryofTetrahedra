import numpy as np

#if __name__=='__main__':

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
def lie_operator(i, a):
    identity = np.eye(2, dtype=complex)
    matrices = [identity, identity, identity, identity]
    matrices[i - 1] = tau[a - 1]
    result = matrices[0]
    for mat in matrices[1:]:
        result = np.kron(result, mat)
    return result

# Funzione per stampare le matrici in modo leggibile
def print_matrix(matrix):
    for row in matrix:
        print("  ".join("{:7.3f}".format(x.real) + "{:+7.3f}j".format(x.imag) for x in row))
    print("\n")

## Calcola e stampa le matrici L_{i,a}
for i in range(1, 4):
    for a in range(1, 4):
        L_i_a = lie_operator(i, a)
        print("L_{{{},{}}} =".format(i, a))
        print_matrix(lie_operator(i,a))

#print(np.kron(np.kron(np.eye(2),np.eye(2)),np.kron(np.eye(2),np.eye(2))))


# UN'APPLICAZIONE

# Define the basis vectors of C^2
e1 = np.array([1, 0])
e2 = np.array([0, 1])

# Function to calculate the Kronecker product of four vectors
def kronecker_product_4(v1, v2, v3, v4):
    return np.kron(np.kron(np.kron(v1, v2), v3), v4)

# All combinations of (v1, v2) repeated 4 times
basis_vectors = [(e1, e1, e1, e1), (e1, e1, e1, e2), (e1, e1, e2, e1), (e1, e1, e2, e2),
                 (e1, e2, e1, e1), (e1, e2, e1, e2), (e1, e2, e2, e1), (e1, e2, e2, e2),
                 (e2, e1, e1, e1), (e2, e1, e1, e2), (e2, e1, e2, e1), (e2, e1, e2, e2),
                 (e2, e2, e1, e1), (e2, e2, e1, e2), (e2, e2, e2, e1), (e2, e2, e2, e2)]

# Calculate the 16 canonical basis vectors of C^2\otimes4
vectors = [kronecker_product_4(*vectors) for vectors in basis_vectors]

# Display the vectors
#for i, vector in enumerate(vectors):
#    print("Vector {}:\n{}\n".format(i+1, vector))

# Calcola D_{12} ristretto al sottospazio invariante Inv(rho)
v1 = vectors[5]-vectors[6]-vectors[9]+vectors[10]
v2 = vectors[3]-vectors[5]-vectors[10]+vectors[12]

L_2_1_v1 = np.dot(lie_operator(2,1),v1)
L_2_2_v1 = np.dot(lie_operator(2,2),v1)
L_2_3_v1 = np.dot(lie_operator(2,3),v1)

D_12_v1 = np.dot(lie_operator(1,1),L_2_1_v1)+np.dot(lie_operator(1,2),L_2_2_v1)+np.dot(lie_operator(1,3),L_2_3_v1)

L_2_1_v2 = np.dot(lie_operator(2,1),v2)
L_2_2_v2 = np.dot(lie_operator(2,2),v2)
L_2_3_v2 = np.dot(lie_operator(2,3),v2)

D_12_v2 = np.dot(lie_operator(1,1),L_2_1_v2)+np.dot(lie_operator(1,2),L_2_2_v2)+np.dot(lie_operator(1,3),L_2_3_v2)

L_3_1_v1 = np.dot(lie_operator(3,1),v1)
L_3_2_v1 = np.dot(lie_operator(3,2),v1)
L_3_3_v1 = np.dot(lie_operator(3,3),v1)

D_13_v1=np.dot(lie_operator(1,1),L_3_1_v1)+np.dot(lie_operator(1,2),L_3_2_v1)+np.dot(lie_operator(1,3),L_3_3_v1)

L_3_1_v2 = np.dot(lie_operator(3,1),v2)
L_3_2_v2 = np.dot(lie_operator(3,2),v2)
L_3_3_v2 = np.dot(lie_operator(3,3),v2)

D_13_v2=np.dot(lie_operator(1,1),L_3_1_v2)+np.dot(lie_operator(1,2),L_3_2_v2)+np.dot(lie_operator(1,3),L_3_3_v2)

print(v1)
print(v2)
# print()
# print(D_12_v1)
# print()
# print(D_12_v2)

# print(v1*0.75)
print(D_12_v1==v1*0.75)
print(D_12_v2==(v1*-0.5+v2*-0.25))
print()
print(D_13_v1==v1*0.25+v2*0.5)
print(D_13_v2==v1*0.5+v2*0.25)


#Compute (L_i)^2=L_i,a\delta^ab L_i,b

def squared_lie_operator(i):
    result = np.zeros((16, 16), dtype=complex)
    for a in range(1, 4):
        L_i_a = lie_operator(i, a)
        result += np.dot(L_i_a, L_i_a)
    return result

print()
print("The Lie squared operator L_i^2 is -3/4 the identity, for any i=1,2,3:")
for i in range(1,4):
    L_i_squared = squared_lie_operator(i)
    print("(L_{})^2 =".format(i))
    print_matrix(L_i_squared)


#L_0,a=-(L_1,a+L_2,a+L_3,a)



