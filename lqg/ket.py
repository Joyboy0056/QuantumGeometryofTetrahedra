import sympy as sp
import itertools
from collections import defaultdict

class Ket2D:
    def __init__(self):
        # Genera tutte le combinazioni di |ABCD> con A, B, C, D in base alla dimensione 2
        self.kets = [f"|{''.join(combination)}>" for combination in itertools.product('+-', repeat=4)]

    def __str__(self):
        return ', '.join(self.kets)

    def contract_with_epsilon(self, index_pair_A_B, index_pair_C_D):
        # Dizionario epsilon per la dimensione 2
        epsilon = {
            ('+', '+'): 0,
            ('-', '-'): 0,
            ('+', '-'): 1,
            ('-', '+'): -1
        }

        contraction_terms = []

        # Itera su tutte le combinazioni di kets
        for ket in self.kets:
            labels = ket[1:5]  # Ottiene i caratteri dalla stringa ket

            # Calcola il prodotto delle epsilon
            product_eps = epsilon[(labels[index_pair_A_B[0]], labels[index_pair_A_B[1]])] * \
                          epsilon[(labels[index_pair_C_D[0]], labels[index_pair_C_D[1]])]

            if product_eps != 0:
                contraction_terms.append(f"{product_eps} * {ket}")

        return " + ".join(contraction_terms)



class Ket3D:
    def __init__(self):
        # Genera tutte le combinazioni di |A1 A2|B1 B2|C1 C2|D1 D2>
        self.kets = [f"|{''.join(combination)}>" for combination in itertools.product('+-', repeat=8)]

    def __str__(self):
        return ', '.join(self.kets)

    def contract_with_epsilon(self, index_pairs):
        # Dizionario epsilon per la dimensione 2
        epsilon = {
            ('+', '+'): 0,
            ('-', '-'): 0,
            ('+', '-'): 1,
            ('-', '+'): -1
        }

        contraction_terms = defaultdict(int)

        # Itera su tutte le combinazioni di kets
        for ket in self.kets:
            labels = ket[1:9]  # Ottiene i caratteri dalla stringa ket

            product_eps = 1
            for pair in index_pairs:
                product_eps *= epsilon[(labels[pair[0]], labels[pair[1]])]

            if product_eps != 0:
                contraction_terms[ket] += product_eps  # Somma i coefficienti per gli stessi ket

        return contraction_terms

    def map_to_C3(self, contraction_terms):
        # Funzione per mappare i termini contratti nello spazio C^3
        mapped_terms = []

        for ket, coefficient in contraction_terms.items():
            # Estrai i simboli A1, A2, B1, B2, C1, C2, D1, D2 dalla rappresentazione |A1 A2|B1 B2|C1 C2|D1 D2>
            ket_part = ket  # Ottiene la parte del ket |A1 A2|B1 B2|C1 C2|D1 D2>
            labels = ket_part[1:9]  # Prende i caratteri effettivi

            # Raggruppiamo le coppie (A1, A2), (B1, B2), (C1, C2), (D1, D2)
            pairs = [(labels[i], labels[i + 1]) for i in range(0, 8, 2)]

            # Applichiamo l'isomorfismo alle coppie
            mapped_ket = []
            for pair in pairs:
                if pair == ('+', '+'):
                    mapped_ket.append('|1>')
                elif pair in {('+', '-'), ('-', '+')}:
                    mapped_ket.append('|0>')
                elif pair == ('-', '-'):
                    mapped_ket.append('|-1>')

            # Mappiamo i coefficienti e costruiamo il termine finale
            mapped_terms.append(f"{coefficient}/4 * {' '.join(mapped_ket)}")

        # Raggruppa i termini simili sommando i coefficienti
        grouped_terms = defaultdict(int)
        for term in mapped_terms:
            coefficient, ket = term.split(' * ')
            coefficient_value = sp.sympify(coefficient)  # Usa sympy per convertire la frazione
            grouped_terms[ket] += float(coefficient_value)

        # Costruisce la stringa finale
        final_terms = []
        for ket, coeff in grouped_terms.items():
            if coeff != 0:  # Ignora i termini con coefficiente 0
                final_terms.append(f"{coeff} * {ket}")

        return " + ".join(final_terms)


