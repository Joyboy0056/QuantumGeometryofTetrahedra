import sympy as sp
import sympy.physics.quantum as spq

ket1 = spq.Ket('+')
ket2 = spq.Ket('-')


sp.pprint(ket1+2*ket1+ket2)

from sympy.physics.paulialgebra import Pauli, evaluate_pauli_product
from sympy import I
sp.pprint(evaluate_pauli_product(I*Pauli(1)*Pauli(2)))

sp.pprint(Pauli(1))
sp.pprint(Pauli(2))
sp.pprint(evaluate_pauli_product(Pauli(3)*ket1))

def Pauli(ind, ket):
    """Funzione che modella l'azione delle matrici di pauli σ = [σ₁, σ₂, σ₃] su C2
        :param ind: è un int 1, 2 o 3
        :param ket: è un ket di sympy.physics.quantum.Ket
    """
    res = spq.Ket('')  # inizializzo un ket vuoto
    if ind == 1:
        if ket == spq.Ket('+'):
            res = spq.Ket('-')
        elif ket == spq.Ket('-'):
            res = spq.Ket('+')

    elif ind == 2:
        if ket == spq.Ket('+'):
            res = sp.I * spq.Ket('-')
        elif ket == spq.Ket('-'):
            res = -sp.I * spq.Ket('+')

    elif ind == 3:
        if ket == spq.Ket('+'):
            res = spq.Ket('+')
        elif ket == spq.Ket('-'):
            res = -spq.Ket('-')

    else:
        raise('No Pauli matrix detected for this index.')

    return res




sp.pprint(Pauli(1, ket1))
sp.pprint(Pauli(1, ket2))
print()
sp.pprint(Pauli(2, ket1))
sp.pprint(Pauli(2, ket2))
print()
sp.pprint(Pauli(3, ket1))
sp.pprint(Pauli(3, ket2))


epsilon = {
                ('+', '+'): 0,
                ('-', '-'): 0,
                ('+', '-'): 1,
                ('-', '+'): -1
            }

v1 = sum(
    epsilon[(A, B)] * epsilon[(C, D)] * spq.Ket(f'{A}{B}{C}{D}')
    for A in ['+', '-'] for B in ['+', '-'] for C in ['+', '-'] for D in ['+', '-']
)
v2 = sum(
    epsilon[(A, D)] * epsilon[(B, C)] * spq.Ket(f'{A}{B}{C}{D}')
    for A in ['+', '-'] for B in ['+', '-'] for C in ['+', '-'] for D in ['+', '-']
)

sp.pprint(v1)
print()
sp.pprint(v2)


def parse_ket(term):
    """Parsa un termine della forma coeff * |ABCD> e restituisce [coeff, (|A>, |B>, |C>, |D>)]"""
    if isinstance(term, sp.Mul):  # Se è una moltiplicazione (coefficiente * ket)
        coefficient = term.args[0]
        ket_term = term.args[1]
    else:  # Se non ha coefficiente numerico esplicito
        coefficient = 1
        ket_term = term
    # Se il ket è un prodotto tensoriale, otteniamo le componenti
    components = tuple(ket_term.args) if isinstance(ket_term, spq.TensorProduct) else (ket_term,)

    return [coefficient, components]



from sympy.physics.quantum import TensorProduct, Ket


def lie(i, a, ket_tp):
    """Applica l'operatore di Lie (i, a) ad un ket o una somma di kets."""

    if isinstance(ket_tp, sp.Add):  # Se è una somma di termini
        return sp.Add(*(lie(i, a, term) for term in ket_tp.args))

    elif isinstance(ket_tp, sp.Mul):  # Se è un prodotto (coefficiente * ket)
        coeffs, kets = sp.Mul.make_args(ket_tp), []
        for term in coeffs:
            if isinstance(term, spq.TensorProduct):  # Trova il primo TensorProduct
                kets.append(term)
        coeff = sp.Mul(*(set(coeffs) - set(kets)))  # Estrai il coefficiente

        if not kets:  # Se non ci sono kets, restituisci solo il coefficiente
            return coeff

        new_ket = lie(i, a, kets[0])  # Applica l'operatore al ket
        return coeff * new_ket if new_ket != 0 else 0  # Evita prodotti con 0

    elif isinstance(ket_tp, spq.TensorProduct):  # Se è un prodotto tensoriale di kets
        components = list(ket_tp.args)
        components[i] = (-sp.I / 2) * Pauli(a, components[i])  # Applica tau_a
        return spq.TensorProduct(*components)

    else:
        return ket_tp  # Mantieni eventuali scalari separati



# Esempio di utilizzo

Inv = [
    sum(
    epsilon[(A, B)] * epsilon[(C, D)] * spq.TensorProduct(spq.Ket(f'{A}'), spq.Ket(f'{B}'), spq.Ket(f'{C}'), spq.Ket(f'{D}'))
    for A in ['+', '-'] for B in ['+', '-'] for C in ['+', '-'] for D in ['+', '-']
        ),
    sum(
    epsilon[(A, D)] * epsilon[(B, C)] * spq.TensorProduct(spq.Ket(f'{A}'), spq.Ket(f'{B}'), spq.Ket(f'{C}'), spq.Ket(f'{D}'))
    for A in ['+', '-'] for B in ['+', '-'] for C in ['+', '-'] for D in ['+', '-']
    )
]

sp.pprint(Inv[0])
sp.pprint(lie(0, 2, Inv[1]))
print()


def diedral(alpha, beta, ket):
    res = spq.Ket('')
    for a in range(1, 4):
        aux = lie(beta, a, ket)
        res += lie(alpha, a, aux)
    return sp.simplify(res)

sp.pprint(diedral(1,3, Inv[0]))