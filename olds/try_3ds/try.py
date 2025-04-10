import sympy as sp
import itertools

class Ket2D:
    def __init__(self):
        # Genera tutte le combinazioni di |ABCD> con A, B, C, D in base alla dimensione 2 e la epsilon
        self.kets = [sp.symbols(f"|{''.join(combination)}>") for combination in itertools.product('+-', repeat=4)]
        self.C2 = ['+', '-']
        self.epsilon = {
            ('+', '+'): 0,
            ('-', '-'): 0,
            ('+', '-'): 1,
            ('-', '+'): -1
        }

        self.ket1 = sum(
            self.epsilon[(A, B)]*self.epsilon[(C, D)]*sp.symbols(f'|{A}{B}{C}{D}>')
            for A in self.C2 for B in self.C2 for C in self.C2 for D in self.C2
        )
        self.ket2 = sum(
            self.epsilon[(A, D)] * self.epsilon[(B, C)] * sp.symbols(f'|{A}{B}{C}{D}>')
            for A in self.C2 for B in self.C2 for C in self.C2 for D in self.C2
        )
        self.ket3 = sum(
            self.epsilon[(A, C)] * self.epsilon[(B, D)] * sp.symbols(f'|{A}{B}{C}{D}>')
            for A in self.C2 for B in self.C2 for C in self.C2 for D in self.C2
        ) # ket3 = ket1 + ket2



class Ket3D:
    def __init__(self):
        self.C2 = ['+', '-']
        self.kets = [sp.symbols(f"|{A1}{A2}|{B1}{B2}|{C1}{C2}|{D1}{D2}>")
                     for A1 in self.C2 for A2 in self.C2 if A1 <= A2
                     for B1 in self.C2 for B2 in self.C2 if B1 <= B2
                     for C1 in self.C2 for C2 in self.C2 if C1 <= C2
                     for D1 in self.C2 for D2 in self.C2 if D1 <= D2
                     ]
        self.epsilon = {
            ('+', '+'): 0,
            ('-', '-'): 0,
            ('+', '-'): 1,
            ('-', '+'): -1
        }
        self.isomorphism = {
            ('+', '+'): (1, 1),  # |1> con coefficiente 1
            ('+', '-'): (0, sp.sqrt(2)/2),  # |0> con coefficiente 1/√2
            ('-', '+'): (0, sp.sqrt(2)/2),  # |0> con coefficiente 1/√2
            ('-', '-'): (-1, 1)  # |-1> con coefficiente 1
        }
    def generate_symmetric_ket1(self):
        """Costruisce il ket imponendo un vincolo di simmetria sugli indici."""
        ket_symmetric = 0

        for A1 in self.C2:
            for A2 in self.C2:
                A_pair = tuple(sorted([A1, A2]))  # Rende la coppia simmetrica
                for B1 in self.C2:
                    for B2 in self.C2:
                        B_pair = tuple(sorted([B1, B2]))
                        for C1 in self.C2:
                            for C2 in self.C2:
                                C_pair = tuple(sorted([C1, C2]))
                                for D1 in self.C2:
                                    for D2 in self.C2:
                                        D_pair = tuple(sorted([D1, D2]))

                                        coeff = (self.epsilon[(A1, B1)] * self.epsilon[(A2, B2)] *
                                                 self.epsilon[(C1, D1)] * self.epsilon[(C2, D2)])

                                        ket = sp.Symbol(
                                            f"|{A_pair[0]}{A_pair[1]}|{B_pair[0]}{B_pair[1]}|{C_pair[0]}{C_pair[1]}|{D_pair[0]}{D_pair[1]}>")
                                        ket_symmetric += coeff * ket

        return ket_symmetric


    def generate_symmetric_ket2(self):
        """Costruisce il ket imponendo un vincolo di simmetria sugli indici."""
        ket_symmetric = 0

        for A1 in self.C2:
            for A2 in self.C2:
                A_pair = tuple(sorted([A1, A2]))  # Rende la coppia simmetrica
                for B1 in self.C2:
                    for B2 in self.C2:
                        B_pair = tuple(sorted([B1, B2]))
                        for C1 in self.C2:
                            for C2 in self.C2:
                                C_pair = tuple(sorted([C1, C2]))
                                for D1 in self.C2:
                                    for D2 in self.C2:
                                        D_pair = tuple(sorted([D1, D2]))

                                        coeff = (self.epsilon[(A1, D1)] * self.epsilon[(A2, D2)] *
                                                 self.epsilon[(B1, C1)] * self.epsilon[(B2, D2)])

                                        ket = sp.Symbol(
                                            f"|{A_pair[0]}{A_pair[1]}|{B_pair[0]}{B_pair[1]}|{C_pair[0]}{C_pair[1]}|{D_pair[0]}{D_pair[1]}>")
                                        ket_symmetric += coeff * ket

        return ket_symmetric


    def generate_symmetric_ket3(self):
        """Costruisce il ket imponendo un vincolo di simmetria sugli indici."""
        ket_symmetric = 0

        for A1 in self.C2:
            for A2 in self.C2:
                A_pair = tuple(sorted([A1, A2]))  # Rende la coppia simmetrica
                for B1 in self.C2:
                    for B2 in self.C2:
                        B_pair = tuple(sorted([B1, B2]))
                        for C1 in self.C2:
                            for C2 in self.C2:
                                C_pair = tuple(sorted([C1, C2]))
                                for D1 in self.C2:
                                    for D2 in self.C2:
                                        D_pair = tuple(sorted([D1, D2]))

                                        coeff = (self.epsilon[(A1, D1)] * self.epsilon[(A2, B2)] *
                                                 self.epsilon[(B1, C1)] * self.epsilon[(C2, D2)])

                                        ket = sp.Symbol(
                                            f"|{A_pair[0]}{A_pair[1]}|{B_pair[0]}{B_pair[1]}|{C_pair[0]}{C_pair[1]}|{D_pair[0]}{D_pair[1]}>")
                                        ket_symmetric += coeff * ket

        return sp.simplify(ket_symmetric)

    def get_ket2_C3(self):
        """Calcola la somma applicando l'isomorfismo"""
        ket_transformed = 0

        for A1 in self.C2:
            for A2 in self.C2:
                for B1 in self.C2:
                    for B2 in self.C2:
                        for C1 in self.C2:
                            for C2 in self.C2:
                                for D1 in self.C2:
                                    for D2 in self.C2:
                                        # Calcola il coefficiente epsilon
                                        coeff = (self.epsilon[(A1, D1)] * self.epsilon[(A2, D2)] *
                                                 self.epsilon[(C1, B1)] * self.epsilon[(C2, B2)])

                                        if coeff != 0:  # Evita termini nulli
                                            # Applica l'isomorfismo a ogni coppia
                                            (state_A, coeff_A) = self.isomorphism[(A1, A2)]
                                            (state_B, coeff_B) = self.isomorphism[(B1, B2)]
                                            (state_C, coeff_C) = self.isomorphism[(C1, C2)]
                                            (state_D, coeff_D) = self.isomorphism[(D1, D2)]

                                            # Calcola il coefficiente totale
                                            total_coeff = coeff * coeff_A * coeff_B * coeff_C * coeff_D

                                            # Crea il nuovo ket con stati trasformati
                                            ket = sp.Symbol(f"|{state_A}|{state_B}|{state_C}|{state_D}>")

                                            # Somma al risultato finale
                                            ket_transformed += total_coeff * ket

        return ket_transformed


    def get_ket3_C3(self):
        """Calcola la somma applicando l'isomorfismo"""
        ket_transformed = 0

        for A1 in self.C2:
            for A2 in self.C2:
                for B1 in self.C2:
                    for B2 in self.C2:
                        for C1 in self.C2:
                            for C2 in self.C2:
                                for D1 in self.C2:
                                    for D2 in self.C2:
                                        # Calcola il coefficiente epsilon
                                        coeff = (self.epsilon[(A1, D1)] * self.epsilon[(A2, B2)] *
                                                 self.epsilon[(C1, B1)] * self.epsilon[(C2, D2)])

                                        if coeff != 0:  # Evita termini nulli
                                            # Applica l'isomorfismo a ogni coppia
                                            (state_A, coeff_A) = self.isomorphism[(A1, A2)]
                                            (state_B, coeff_B) = self.isomorphism[(B1, B2)]
                                            (state_C, coeff_C) = self.isomorphism[(C1, C2)]
                                            (state_D, coeff_D) = self.isomorphism[(D1, D2)]

                                            # Calcola il coefficiente totale
                                            total_coeff = coeff * coeff_A * coeff_B * coeff_C * coeff_D

                                            # Crea il nuovo ket con stati trasformati
                                            ket = sp.Symbol(f"|{state_A}|{state_B}|{state_C}|{state_D}>")

                                            # Somma al risultato finale
                                            ket_transformed += total_coeff * ket

        return ket_transformed


    def get_ket1_C3(self):
        """Calcola la somma applicando l'isomorfismo"""
        ket_transformed = 0

        for A1 in self.C2:
            for A2 in self.C2:
                for B1 in self.C2:
                    for B2 in self.C2:
                        for C1 in self.C2:
                            for C2 in self.C2:
                                for D1 in self.C2:
                                    for D2 in self.C2:
                                        # Calcola il coefficiente epsilon
                                        coeff = (self.epsilon[(A1, B1)] * self.epsilon[(A2, B2)] *
                                                 self.epsilon[(C1, D1)] * self.epsilon[(C2, D2)])

                                        if coeff != 0:  # Evita termini nulli
                                            # Applica l'isomorfismo a ogni coppia
                                            (state_A, coeff_A) = self.isomorphism[(A1, A2)]
                                            (state_B, coeff_B) = self.isomorphism[(B1, B2)]
                                            (state_C, coeff_C) = self.isomorphism[(C1, C2)]
                                            (state_D, coeff_D) = self.isomorphism[(D1, D2)]

                                            # Calcola il coefficiente totale
                                            total_coeff = coeff * coeff_A * coeff_B * coeff_C * coeff_D

                                            # Crea il nuovo ket con stati trasformati
                                            ket = sp.Symbol(f"|{state_A}|{state_B}|{state_C}|{state_D}>")

                                            # Somma al risultato finale
                                            ket_transformed += total_coeff * ket

        return ket_transformed


print(Ket3D().generate_symmetric_ket1())
print(Ket3D().generate_symmetric_ket2())
print(Ket3D().generate_symmetric_ket3())

print(Ket2D().kets)
print(len(Ket2D().kets))
print(Ket3D().kets)
print(len(Ket3D().kets))

print('ket: ', f"|{Ket3D().isomorphism[('+','-')][0]}>")
print('coeff: ', Ket3D().isomorphism[('+','-')][1])


print('\nv1:')
sp.pprint(Ket3D().get_ket1_C3())
print('\nv2:')
sp.pprint(Ket3D().get_ket2_C3())
print('\nv3:')
sp.pprint(Ket3D().get_ket3_C3())

import sympy as sp

import sympy as sp

# Definiamo le matrici di Pauli in SymPy
sigma_1 = sp.Matrix([[0, 1], [1, 0]])
sigma_2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
sigma_3 = sp.Matrix([[1, 0], [0, -1]])

# Definiamo le matrici tau_a
tau_1 = -sp.I/2 * sigma_1
tau_2 = -sp.I/2 * sigma_2
tau_3 = -sp.I/2 * sigma_3

tau = [tau_1, tau_2, tau_3]

def lie_operator(i, a):
    identity = sp.eye(2)  # Matrice identità in SymPy
    matrices = [identity, identity, identity, identity]
    matrices[i - 1] = tau[a - 1]  # Inserisci la matrice tau corrispondente
    result = matrices[0]
    for mat in matrices[1:]:
        result = sp.tensorproduct(result, mat)  # Prodotto di Kronecker per costruire la matrice finale
    return sp.Matrix(result.reshape(2**4,2**4))


# Funzione per l'operatore diedrale
def diedral_operator(i, j):
    result = sp.MutableDenseNDimArray(sp.zeros(16, 16))  # Matrice iniziale

    for a in range(1, 4):
        L_i_a = lie_operator(i, a)
        L_j_a = lie_operator(j, a)

        tensor_product = sp.tensorproduct(L_i_a, L_j_a).as_mutable().reshape(16, 16)

        result += tensor_product  # Sommiamo il tensore alla matrice

    return result



sp.pprint(lie_operator(1,2))
#sp.pprint(diedral_operator(1, 2))



# Definiamo i vettori base di C^3
e1 = sp.Matrix([1, 0])
e2 = sp.Matrix([0, 1])
e3 = sp.Matrix([0, 0, 1])

# Lista dei vettori base
basis_vectors = [e1, e2]

# Generiamo tutte le combinazioni dei 4 indici
n = len(basis_vectors)
base_states = []
for indices in itertools.product(range(n), repeat=4):
    v1, v2, v3, v4 = [basis_vectors[i] for i in indices]
    tensor_product = sp.tensorproduct(sp.tensorproduct(sp.tensorproduct(v1, v2), v3), v4)
    base_states.append(tensor_product)

for i in range(len(base_states)):
    base_states[i] = sp.Matrix(base_states[i].reshape(n**4,1))

# qui base_states è una lista di 81 vettori di C^81

support = ['+', '-']
kets = [sp.symbols(f"|{A}{B}{C}{D}>")
                    for A in support
                    for B in support
                    for C in support
                    for D in support
                     ]

# qui kets sono i 16 ket di C^2\otimes4


# Stampa un esempio di vettore base
#sp.pprint(base_states[0])  # Primo vettore base
print(f"Numero totale di vettori base: {len(base_states)}")  # Dovrebbe essere 3^4 = 81

#sp.pprint(sp.tensorproduct(lie_operator(1,2), base_states[0]))

print(lie_operator(1,2).shape)
print(base_states[3].shape)

Lie_12 = sp.Matrix(lie_operator(1,2).reshape(n**4,n**4))
L_12_v1 = Lie_12*base_states[0]

sp.pprint(L_12_v1)

ket_to_vector = {kets[i]: base_states[i] for i in range(len(kets))}

coefs = [sp.symbols(f'a{j+1}') for j in range(len(base_states))]
sp.pprint(coefs)

sp.pprint(coefs[0]*base_states[0] + coefs[1]*base_states[1])

#sp.pprint(sum(coefs[j]*base_states[j] for j in range(len(base_states))))


#--------------------------------------------------------

#definiamo l'operatore di Lie solo come azione su un singolo |A>

class SymbolicKet:
    def __init__(self, elements):
        """
        Inizializza il ket con una lista di elementi.
        Esempio: elements = ['+', '-', '+', '-']
        """
        self.elements = list(elements)  # Memorizza i simboli come lista immutabile

    def __repr__(self):
        """Rappresentazione testuale del ket"""
        return f"|{''.join(map(str, self.elements))}>"

    def __getitem__(self, index):
        """Permette di accedere a un elemento specifico del ket"""
        return self.elements[index]

    def __setitem__(self, index, value):
        """Permette di modificare un elemento del ket"""
        self.elements[index] = value




def is_multiple(v1, v2):
    """
    Controlla se v1 è un multiplo di v2 e restituisce il coefficiente se vero.

    Args:
        v1 (sp.Matrix): Primo vettore.
        v2 (sp.Matrix): Secondo vettore.

    Returns:
        tuple (bool, coeff): (True, coeff) se v1 è un multiplo di v2, altrimenti (False, None).
    """
    if v1.shape != v2.shape:
        return False, None  # Devono avere la stessa dimensione

    # Trova gli indici delle componenti non nulle di v2
    indices = [i for i in range(v2.shape[0]) if v2[i] != 0]

    if not indices:
        return False, None  # Evita divisioni per zero

    # Calcola il rapporto tra le componenti corrispondenti di v1 e v2
    coeffs = [v1[i] / v2[i] for i in indices]

    # Se tutti i coefficienti sono uguali, restituisci il primo come fattore comune
    if all(sp.simplify(c - coeffs[0]) == 0 for c in coeffs):
        return sp.simplify(coeffs[0])
    return None


def lie(i, a, ket):
    # Definiamo le matrici di Pauli in SymPy
    sigma_1 = sp.Matrix([[0, 1], [1, 0]])
    sigma_2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sigma_3 = sp.Matrix([[1, 0], [0, -1]])

    ket_dict = {
        '+': sp.Matrix([1, 0]),
        '-': sp.Matrix([0, 1])
    }

    # Definiamo le matrici tau_a
    tau = [-sp.I / 2 * sigma_1, -sp.I / 2 * sigma_2, -sp.I / 2 * sigma_3]

    new_ket = ket
    aux = tau[a] * ket_dict[ket[i]] #così è un vettore però
    for key, value in ket_dict.items():
        mul = is_multiple(aux, value)
        if mul != 0:
            new_ket[i] = key
            new_ket = mul * sp.symbols(f'{ket}')

    return new_ket


# Esempio di utilizzo
kets = [
        SymbolicKet(f"{A}{B}{C}{D}")
        for A in ['+', '-'] for B in ['+', '-']
        for C in ['+', '-'] for D in ['+', '-']
    ]


print(kets)
ket1 = kets[0]

for j in range(len(kets)):
    sp.pprint(lie(1,2, kets[j]))


def LieOperator(i, a, ket):
    """Estende l'operatre di Lie a una somma di ket"""

