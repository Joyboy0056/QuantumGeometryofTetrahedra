import sympy as sp
from itertools import product

class C:
    def __init__(self, dimension):
        self.dim = dimension
        self.vec_basis = [sp.eye(self.dim)[:, j] for j in range(self.dim)]

        if self.dim == 2:
            self.ket_basis = [sp.symbols(f'|{A}>') for A in ['+', '-']]

            self.sigma = [
                sp.Matrix([[0, 1], [1, 0]]),
                sp.Matrix([[0, -sp.I], [sp.I, 0]]),
                sp.Matrix([[1, 0], [0, -1]])
            ]
            self.tau = [-sp.I / 2 * self.sigma[0],
                        -sp.I / 2 * self.sigma[1],
                        -sp.I / 2 * self.sigma[2]]

            self.symbol_to_index = {'+': 0, '-': 1}
            self.index_to_symbol = {0: '+', 1: '-'}

    def apply_tau(self, a, ket_symbol):
        if self.dim != 2:
            raise ValueError("tau è definito solo per C^2")

        ket_str = str(ket_symbol)
        ket_char = ket_str[1]
        ket_index = self.symbol_to_index[ket_char]

        transformed_ket = self.tau[a] * self.vec_basis[ket_index]

        result = 0
        for i in range(self.dim):
            if transformed_ket[i] != 0:
                new_ket_symbol = sp.symbols(f"|{self.index_to_symbol[i]}>")
                result += transformed_ket[i] * new_ket_symbol

        return result

    def __repr__(self):
        return f"C^{self.dim}"


class ProductRepr:
    def __init__(self, *spaces):
        self.spaces = spaces
        self.total_dim = [space.dim for space in spaces]
        index_sets = [space.ket_basis for space in spaces]
        self.ket_basis = {idx: sp.symbols(f"|{''.join(str(k)[1] for k in idx)}>")
                          for idx in product(*index_sets)}

        if self.spaces[0].dim == 2:
            self.epsilon = {
                ('+', '+'): 0,
                ('-', '-'): 0,
                ('+', '-'): 1,
                ('-', '+'): -1
            }

            self.ket1 = sum(
                self.epsilon[(A, B)] * self.epsilon[(C, D)] * sp.symbols(f'|{A}{B}{C}{D}>')
                for A in ['+', '-'] for B in ['+', '-'] for C in ['+', '-'] for D in ['+', '-']
            )
            self.ket2 = sum(
                self.epsilon[(A, D)] * self.epsilon[(B, C)] * sp.symbols(f'|{A}{B}{C}{D}>')
                for A in ['+', '-'] for B in ['+', '-'] for C in ['+', '-'] for D in ['+', '-']
            )
            self.ket3 = sum(
                self.epsilon[(A, C)] * self.epsilon[(B, D)] * sp.symbols(f'|{A}{B}{C}{D}>')
                for A in ['+', '-'] for B in ['+', '-'] for C in ['+', '-'] for D in ['+', '-']
            )

    def __repr__(self):
        dims = " ⊗ ".join([f"C^{s.dim}" for s in self.spaces])
        return f"TensorC({dims})"

    def lie_operator(self, a, i, ket_expr):
        """
        Applica tau_a su una combinazione lineare di kets.

        Args:
            a (int): Indice dell'operatore tau.
            i (int): is the space_idx Indice dello spazio su cui agisce tau.
            ket_expr (sympy.Expr): Espressione simbolica di kets.

        Returns:
            Espressione simbolica della combinazione lineare risultante.
        """
        if not isinstance(ket_expr, sp.Expr):
            ket_expr = sp.sympify(ket_expr)

        result = 0
        for term in ket_expr.as_ordered_terms():
            coeff, ket_symbol = term.as_coeff_Mul()

            ket_str = str(ket_symbol)[1:-1]
            ket_chars = list(ket_str)

            if self.spaces[i].dim != 2:
                raise ValueError("tau è definito solo per C^2")

            transformed_ket = self.spaces[i].apply_tau(a, sp.symbols(f"|{ket_chars[i]}>"))

            for subterm in transformed_ket.as_ordered_terms():
                sub_coeff, new_ket = subterm.as_coeff_Mul()
                new_ket_str = ket_chars[:]
                new_ket_str[i] = str(new_ket)[1]
                new_ket_symbol = sp.symbols(f"|{''.join(new_ket_str)}>")

                result += coeff * sub_coeff * new_ket_symbol

        return result


# Test delle classi
C2 = C(2)
repr1half = ProductRepr(C2, C2, C2, C2)

# Test con combinazione lineare di stati
for a in range(3):
    print(f"τ{a+1}[ {repr1half.ket1} ] =", repr1half.lie_operator(a, 3, repr1half.ket1))

import sympy as sp
import sympy.physics.quantum as spq


# Funzione Pauli generalizzata
def Pauli(ind, ket, n):
    res = spq.Ket('')  # inizializzo un ket vuoto

    if ind == 1:  # Pauli-X
        if ket == spq.Ket('+'):
            res = spq.Ket('-')
        elif ket == spq.Ket('-'):
            res = spq.Ket('+')
        else:
            res = apply_pauli_x(ket, n)

    elif ind == 2:  # Pauli-Y
        if ket == spq.Ket('+'):
            res = sp.I * spq.Ket('-')
        elif ket == spq.Ket('-'):
            res = -sp.I * spq.Ket('+')
        else:
            res = apply_pauli_y(ket, n)

    elif ind == 3:  # Pauli-Z
        if ket == spq.Ket('+'):
            res = spq.Ket('+')
        elif ket == spq.Ket('-'):
            res = -spq.Ket('-')
        else:
            res = apply_pauli_z(ket, n)

    else:
        raise ValueError('No Pauli matrix detected for this index.')

    return res


# Funzione per l'azione di Pauli-X su ket in uno spazio di dimensione n
def apply_pauli_x(ket, n):
    # Esegui azioni generali su base |0>, |1>, ..., |n-1> per Pauli-X
    print(f"Pauli-X action on {ket}")
    return ket  # Placeholder, implementa la logica per n>2


# Funzione per l'azione di Pauli-Y su ket in uno spazio di dimensione n
def apply_pauli_y(ket, n):
    print(f"Pauli-Y action on {ket}")
    return ket  # Placeholder


# Funzione per l'azione di Pauli-Z su ket in uno spazio di dimensione n
def apply_pauli_z(ket, n):
    print(f"Pauli-Z action on {ket}")
    return ket  # Placeholder


# Test con un ket di base
ket_plus = spq.Ket('+')
ket_minus = spq.Ket('-')

# Test con Pauli-X per j = 1/2 (dimensione 2)
print("Testing Pauli-X on Ket(+):")
result = Pauli(1, ket_plus, 2)
print(result)

print("Testing Pauli-X on Ket(-):")
result = Pauli(1, ket_minus, 2)
print(result)

# Test con Pauli-Y per j = 1/2 (dimensione 2)
print("Testing Pauli-Y on Ket(+):")
result = Pauli(2, ket_plus, 2)
print(result)

print("Testing Pauli-Y on Ket(-):")
result = Pauli(2, ket_minus, 2)
print(result)

# Test con Pauli-Z per j = 1/2 (dimensione 2)
print("Testing Pauli-Z on Ket(+):")
result = Pauli(3, ket_plus, 2)
print(result)

print("Testing Pauli-Z on Ket(-):")
result = Pauli(3, ket_minus, 2)
print(result)

# Test con Pauli-X su uno spazio di dimensione 3 (ad esempio, j = 3/2)
print("Testing Pauli-X on Ket(0) in 3D space:")
ket_3d = spq.Ket('0')  # o un altro stato di base
result = Pauli(1, ket_3d, 3)
print(result)


def su2_generators(j):
    """
    Restituisce i generatori dell'algebra su(2) per spin j,
    come matrici di dimensione (2j+1)x(2j+1).
    """
    dim = int(2 * j + 1)  # Dimensione dello spazio

    # Gestiamo il caso j semintero
    if j == int(j) + 1/2:

        m_values = [m/2 for m in range(-int(2*j), int(2*j)+1, 1)]
        Jz = sp.diag(*m_values)

    # Matrice Jz (diagonale con autovalori m) per j intero
    Jz = sp.diag(*[m for m in range(-int(j), int(j) + 1)])

    # Operatori di innalzamento e abbassamento J+ e J-
    Jp = sp.zeros(dim, dim)
    Jm = sp.zeros(dim, dim)

    for m in range(dim - 1):
        m_val = -j + m  # Valore m dell'autovalore attuale
        coeff = sp.sqrt((j - m_val) * (j + m_val + 1))
        Jp[m, m + 1] = coeff
        Jm[m + 1, m] = coeff  # J_- è la trasposta di J_+

    # Calcolo di Jx e Jy
    Jx = (Jp + Jm)
    Jy = (Jp - Jm) / sp.I

    return [sp.simplify(Jx), sp.simplify(Jy), sp.simplify(Jz)]


# Esempio per j = 1/2 e j = 1
J_half = su2_generators(1/2)
J_one = su2_generators(1)
J_3half = su2_generators(3/2)

print("Generatori per j=1/2 (Matrici di Pauli):")
for J in J_half:
    sp.pprint(J)
    print()

print("Generatori per j=1:")
for J in J_one:
    sp.pprint(J)
    print()

print("Generatori per j=3/2:")
for J in J_3half:
    sp.pprint(J)
    print()

j = 1/2
for m in range(-int(2*j), int(2*j)+1, 1):
    print(sp.simplify(m/2))
