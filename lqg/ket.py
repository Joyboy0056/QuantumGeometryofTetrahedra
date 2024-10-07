from itertools import repeat

import sympy as sp
import itertools
from collections import defaultdict
import numpy as np
import math



class Ket2D:
    def __init__(self):
        # Genera tutte le combinazioni di |ABCD> con A, B, C, D in base alla dimensione 2 e la epsilon
        self.kets = [f"|{''.join(combination)}>" for combination in itertools.product('+-', repeat=4)]
        self.epsilon = {
            ('+', '+'): 0,
            ('-', '-'): 0,
            ('+', '-'): 1,
            ('-', '+'): -1
        }
    def __str__(self):
        return ', '.join(self.kets)

    def compute_ket(self, index_pair_A_B, index_pair_C_D):
        """Funzione che calcola un vettore nello spazio C^2 ⊗ C^2 ⊗ C^2 ⊗ C^2 in base a contrazioni con epsilon,
        nella base dei ket.

        :param index_pair_A_B: tupla che rappresenta la prima coppia di indici per la contrazione.
        :param index_pair_C_D: tupla che rappresenta la seconda coppia di indici per la contrazione.
        :return: Un stringa che rappresenta il vettore ket. """

        contraction_terms = []

        # Itera su tutte le combinazioni di kets
        for ket in self.kets:
            labels = ket[1:5]  # Ottiene i caratteri dalla stringa ket

            # Calcola il prodotto delle epsilon
            product_eps = self.epsilon[(labels[index_pair_A_B[0]], labels[index_pair_A_B[1]])] * \
                          self.epsilon[(labels[index_pair_C_D[0]], labels[index_pair_C_D[1]])]

            if product_eps != 0:
                contraction_terms.append(f"{product_eps} * {ket}")

        return " + ".join(contraction_terms)

    def compute_vector(self, index_pairs):
        """
        Calcola un vettore nello spazio C^2 ⊗ C^2 ⊗ C^2 ⊗ C^2 in base a contrazioni con epsilon.

        :param index_pairs: Lista di tuple che rappresentano le coppie di indici per le contrazioni.
        :return: Un vettore come array numpy.
        """
        vector = np.zeros(len(self.kets), dtype=complex)  # Vettore di dimensione 16 come array numpy

        for ket in self.kets:
            labels = ket[1:5]  # Ottieni i caratteri dalla stringa ket
            product_eps = 1

            # Calcola il prodotto delle epsilon per ogni coppia di indici
            for (i, j) in index_pairs:
                product_eps *= self.epsilon[(labels[i], labels[j])]

            # Se il prodotto è diverso da zero, aggiungi il contributo al vettore
            if product_eps != 0:
                idx = self.kets.index(ket)
                vector[idx] += product_eps  # Aggiungi al vettore

        return vector


    #la prossima funzione servità per gesire le ortogonalizzazioni
    def map_to_sqrt(self, coeff):
        """Mappa un coefficiente ai multipli di sqrt(2) o sqrt(3), se possibile."""
        sqrt_multiples = {
            math.sqrt(3): "sqrt(3)",
            math.sqrt(3) / 2: "sqrt(3)/2",
            math.sqrt(3) / 3: "sqrt(3)/3",
            math.sqrt(3) / 4: "sqrt(3)/4",
            math.sqrt(3) / 5: "sqrt(3)/5",
            math.sqrt(3) / 6: "sqrt(3)/6",
            math.sqrt(2): "sqrt(2)",
            math.sqrt(2) / 2: "sqrt(2)/2",
            math.sqrt(2) / 3: "sqrt(2)/3",
            math.sqrt(2) / 4: "sqrt(2)/4",
            math.sqrt(2) / 5: "sqrt(2)/5"
        }

        # Tolleranza per il confronto tra float
        tolerance = 1e-6

        # Verifica se il coefficiente è un multiplo di sqrt(2) o sqrt(3)
        for value, symbol in sqrt_multiples.items():
            if abs(coeff - value) < tolerance:
                return symbol
            elif abs(coeff + value) < tolerance:
                return f"-{symbol}"

        # Se non corrisponde a nessun multiplo noto, restituisce il coeff originale
        return coeff

    def vector_to_ket(self, vector):
        """Mappa un vector nella sua rappresentazione nella base dei kets"""

        # Controlla se tutti i coefficienti sono nulli
        if all(coeff == 0 for coeff in vector):
            return "null vector"
        ket_representation = []
        for i, coeff in enumerate(vector):
            if coeff != 0:
                possibly_sqrt_mult_coeff = self.map_to_sqrt(coeff.real)
                ket_representation.append(f"{possibly_sqrt_mult_coeff} * {self.kets[i]}")
        return " + ".join(ket_representation)


   #Towards ket_to_vector: first some aux functions
    def kronecker_product_4(self, v1, v2, v3, v4):
        """Function to calculate the Kronecker product of four vectors"""
        return np.kron(np.kron(np.kron(v1, v2), v3), v4)

    def create_ket_map(self):
        """Questa funzione implementa in un dict l'isomorfismo tra le rappr vector e ket in C^2\otimes4"""
        # Define the basis vectors of C^2
        e1 = np.array([1, 0]) # |+>
        e2 = np.array([0, 1]) # |->
        # Create all combinations of (e1, e2) repeated 4 times
        # it spans a 16 dimensional vector space
        #basis_vectors = [
         #   (e1, e1, e1, e1), (e1, e1, e1, e2), (e1, e1, e2, e1), (e1, e1, e2, e2),
          #  (e1, e2, e1, e1), (e1, e2, e1, e2), (e1, e2, e2, e1), (e1, e2, e2, e2),
           # (e2, e1, e1, e1), (e2, e1, e1, e2), (e2, e1, e2, e1), (e2, e1, e2, e2),
            #(e2, e2, e1, e1), (e2, e2, e1, e2), (e2, e2, e2, e1), (e2, e2, e2, e2)]
        #meglio farlo con itertools
        basis_vectors = list(itertools.product((e1, e2), repeat=4))
        # Calculate the 16 canonical basis vectors of C^2\otimes4
        vectors = [self.kronecker_product_4(*vec) for vec in basis_vectors]

        # Create the ket_map
        ket_map = {}
        ket_labels = [
            '|++++>', '|+++->', '|++-+>', '|++-->',
            '|+-++>', '|+-+->', '|+--+>', '|+--->',
            '|-+++>', '|-++->', '|-+-+>', '|-+-->',
            '|--++>', '|--+->', '|---+>', '|---->'
        ]

        # Map each ket label to its corresponding vector
        for label, vector in zip(ket_labels, vectors):
            ket_map[label] = vector

        return ket_map

    def ket_to_vector(self, ket_str):
        """
        Converte un ket espresso come stringa nella sua rappresentazione vettoriale.

        Args:
        ket_str (str): La stringa del ket da convertire.
        ket_map (dict): Un dizionario che mappa i kets a vettori.

        Returns:
        np.ndarray: Il vettore corrispondente al ket.
        """
        ket_map = self.create_ket_map()

        # Inizializza il vettore
        vector = np.zeros(len(ket_map), dtype=complex)

        # Dividi il ket in termini
        terms = ket_str.split(' + ')

        for term in terms:
            # Controlla se c'è un coefficiente
            if ' * ' in term:
                coeff_str, ket = term.split(' * ')
                coeff = np.complex128(eval(coeff_str))  # Usa eval per convertire il coefficiente in complesso
            else:
                coeff = 1  # Coefficiente implicito di 1
                ket = term.strip()

            # Aggiungi il coefficiente al vettore corrispondente al ket
            if ket in ket_map:
                vector += coeff * ket_map[ket]
            else:
                raise ValueError(f"Ket '{ket}' non riconosciuto nella ket_map.")

        return vector

    # Highlighting the vector structure of C^2\otimes4 also for the ket repr
    def ket_sum_vec_repr(self, vector1, vector2):
        """
        Somma due vettori ket.

        :param vector1: Il primo vettore ket come array numpy.
        :param vector2: Il secondo vettore ket come array numpy.
        :return: Un nuovo vettore ket che rappresenta la somma.
        """
        if vector1.shape != vector2.shape:
            raise ValueError("Vector must be of the same dimension to be summed.")
        return self.vector_to_ket(vector1 + vector2)

    def ket_sum(self, ket_expr1, ket_expr2):
        """
        Somma due espressioni di kets date come stringhe e restituisce la nuova espressione.

        Args:
        ket_expr1 (str): La prima espressione di kets.
        ket_expr2 (str): La seconda espressione di kets.

        Returns:
        str: La somma delle due espressioni di kets.
        """

        # Suddivide le due espressioni in termini di ket
        ket_terms1 = ket_expr1.split(' + ')
        ket_terms2 = ket_expr2.split(' + ')

        # Combina i termini in una lista
        combined_terms = ket_terms1 + ket_terms2

        # Crea un dizionario per gestire i coefficienti
        ket_dict = {}

        for term in combined_terms:
            if ' * ' in term:  # Se c'è un coefficiente esplicito
                coeff_str, ket = term.split(' * ')
                coeff = sp.sympify(coeff_str)  # Converte il coefficiente in un oggetto sympy
            else:
                # Se non c'è coefficiente esplicito, è implicito 1
                coeff = sp.sympify(1)
                ket = term

            # Somma i coefficienti per il ket
            if ket in ket_dict:
                ket_dict[ket] += coeff
            else:
                ket_dict[ket] = coeff

        # Costruisce la nuova espressione a partire dal dizionario
        new_terms = []
        for ket, coeff in ket_dict.items():
            if coeff == 1:
                new_terms.append(f"{ket}")
            elif coeff == -1:
                new_terms.append(f"-{ket}")
            else:
                new_terms.append(f"{coeff} * {ket}")

        # Combina i nuovi termini in una stringa
        return ' + '.join(new_terms)


    def ket_scalar_mult(self, ket_expr, scalar_str):
        """
        Moltiplica un'espressione di ket per uno scalare dato come stringa e restituisce la nuova espressione.

        Args:
        ket_expr (str): Un'espressione di kets nella forma 'a * |ket1> + b * |ket2> + ...'
        scalar_str (str): Lo scalare dato come stringa, es. '0.5', 'sqrt(2)'.

        Returns:
        str: La nuova espressione dei kets moltiplicata per lo scalare.
        """

        scalar = sp.sympify(scalar_str)  # Converte lo scalare in un oggetto sympy
        ket_terms = ket_expr.split(' + ')  # Suddivide l'espressione in termini di ket
        new_terms = []

        for term in ket_terms:
            if ' * ' in term:  # Se c'è un coefficiente esplicito
                coeff_str, ket = term.split(' * ')
                coeff = sp.sympify(coeff_str)  # Converte il coefficiente in un oggetto sympy
            else:
                # Se non c'è coefficiente esplicito, è implicito 1
                coeff = sp.sympify(1)
                ket = term

            # Moltiplica il coefficiente per lo scalare
            new_coeff = scalar * coeff

            # Se il coefficiente è 1, ometti il coefficiente
            if new_coeff == 1:
                new_terms.append(f"{ket}")
            # Se il coefficiente è -1, ometti il coefficiente ma mantieni il segno
            elif new_coeff == -1:
                new_terms.append(f"-{ket}")
            else:
                # Aggiungi il nuovo coefficiente e il ket
                new_terms.append(f"{new_coeff} * {ket}")

        # Combina i nuovi termini in una stringa
        return ' + '.join(new_terms)


class Ket3D:
    def __init__(self):
        # Genera tutte le combinazioni di |A1 A2|B1 B2|C1 C2|D1 D2>
        self.kets = [f"|{''.join(combination)}>" for combination in itertools.product('+-', repeat=8)]

        self.epsilon = {
            ('+', '+'): 0,
            ('-', '-'): 0,
            ('+', '-'): 1,
            ('-', '+'): -1
        }
        # Definisco l'isomorfismo tra coppie di bit e il mapping in C^3
        self.isomorphism = {
            '++': (1, 1),  # |1> con coefficiente 1
            '+-': (0, np.sqrt(2)),  # |0> con coefficiente √2
            '-+': (0, np.sqrt(2)),  # |0> con coefficiente √2
            '--': (-1, 1)  # |-1> con coefficiente 1
        }

    def __str__(self):
        return ', '.join(self.kets)

    def eps_contraction(self, index_pairs):
        contraction_terms = defaultdict(int)

        # Itera su tutte le combinazioni di kets
        for ket in self.kets:
            labels = ket[1:9]  # Ottiene i caratteri dalla stringa ket

            product_eps = 1
            for pair in index_pairs:
                product_eps *= self.epsilon[(labels[pair[0]], labels[pair[1]])]

            if product_eps != 0:
                contraction_terms[ket] += product_eps
        return contraction_terms

    def C3_mapping(self, contraction_terms):
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
                    mapped_ket.append('1')
                elif pair in {('+', '-'), ('-', '+')}:
                    mapped_ket.append('0')
                elif pair == ('-', '-'):
                    mapped_ket.append('-1')

            # Mappiamo i coefficienti e costruiamo il termine finale
            mapped_terms.append(f"{coefficient}/4 * |{'|'.join(mapped_ket)}>")

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


    def compute_ket(self, index_pairs):
        """:param: lista di 4 tuple
        :return: stringa con vettore ket contratto"""

        contraction_terms = defaultdict(int)

        # Itera su tutte le combinazioni di kets
        for ket in self.kets:
            labels = ket[1:9]  # Ottiene i caratteri dalla stringa ket

            product_eps = 1
            for pair in index_pairs:
                product_eps *= self.epsilon[(labels[pair[0]], labels[pair[1]])]

            if product_eps != 0:
                contraction_terms[ket] += product_eps  # Somma i coefficienti per gli stessi ket

        #a questo punto i termini contratti sono stati archiviati nel dict, partiamo con il mapping
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
                    mapped_ket.append('1')
                elif pair in {('+', '-'), ('-', '+')}:
                    mapped_ket.append('0')
                elif pair == ('-', '-'):
                    mapped_ket.append('-1')

            # Mappiamo i coefficienti e costruiamo il termine finale
            mapped_terms.append(f"{coefficient}/4 * |{'|'.join(mapped_ket)}>")

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

    def compute_vector(self, index_pairs):
        """
        Calcola un vettore nello spazio C^3 ⊗ C^3 ⊗ C^3 ⊗ C^3 in base a contrazioni con epsilon.

        :param index_pairs: Lista di tuple che rappresentano le coppie di indici per le contrazioni.
        :return: Un vettore come array numpy.
        """
        # Vettore di dimensione 81 come array numpy
        vector = np.zeros((3, 3, 3, 3), dtype=complex)  # Dimensione corretta

        for ket in self.kets:
            labels = ket[1:9]  # Ottieni i caratteri dalla stringa ket |1 2|3 4|5 6|7 8>
            product_eps = 1

            # Calcola il prodotto delle epsilon per ogni coppia di indici
            for (i, j) in index_pairs:
                product_eps *= self.epsilon[(labels[i], labels[j])]

            # Se il prodotto è diverso da zero, aggiungi il contributo al vettore
            if product_eps != 0:
                idx = self.kets.index(ket)

                # Calcola la posizione nel vettore C^3
                # Assicurati che l'indice non superi i limiti
                c3_idx = (idx // 27, (idx // 9) % 3, (idx // 3) % 3, idx % 3)

                # Controllo per gli indici
                if all(0 <= idx < 3 for idx in c3_idx):
                    vector[c3_idx] += product_eps  # Aggiungi al vettore
                else:
                    pass
                    #print(f"Index out of bounds for ket {ket}: c3_idx = {c3_idx}")
        # Restituisci il vettore come array 81x1
        return vector.flatten()  # Concatenare in un vettore 81x1

    # la prossima funzione servità per gesire le ortogonalizzazioni
    def map_to_sqrt1(self, coeff):
        """Mappa un coefficiente ai multipli di sqrt(2) o sqrt(3), se possibile."""
        sqrt_multiples = {
            math.sqrt(3): "sqrt(3)",
            math.sqrt(3) / 2: "sqrt(3)/2",
            math.sqrt(3) / 3: "sqrt(3)/3",
            math.sqrt(3) / 4: "sqrt(3)/4",
            math.sqrt(3) / 5: "sqrt(3)/5",
            math.sqrt(3) / 6: "sqrt(3)/6",
            math.sqrt(2): "sqrt(2)",
            math.sqrt(2) / 2: "sqrt(2)/2",
            math.sqrt(2) / 3: "sqrt(2)/3",
            math.sqrt(2) / 4: "sqrt(2)/4",
            math.sqrt(2) / 5: "sqrt(2)/5",
            2 * math.sqrt(2) : "2 * sqrt(2)",
            2 : "2"
        }

        # Tolleranza per il confronto tra float
        tolerance = 1e-6

        # Verifica se il coefficiente è un multiplo di sqrt(2) o sqrt(3)
        for value, symbol in sqrt_multiples.items():
            if abs(coeff - value) < tolerance:
                return symbol
            elif abs(coeff + value) < tolerance:#questo gestisce i segni
                return f"-{symbol}"

        # Se non corrisponde a nessun multiplo noto, restituisce il coeff originale
        return coeff

    def map_to_sqrt2(self, coeff):
        sqrt2 = sp.sqrt(2)
        sqrt3 = sp.sqrt(3)

        # Creiamo espressioni con multipli di sqrt(2) e sqrt(3)
        res = sp.nsimplify(coeff, [sqrt2, sqrt3], tolerance=1e-10)
        return res

    def map_ket_to_c3(self, ket):
        """
        Mappa un ket del tipo |++|-+|-+|++> nel corrispondente |A|B|C|D> con A,B,C,D ∈ {1, 0, -1}
        secondo l'isomorfismo tra C^2 ⊗ C^2 e C^3, includendo i coefficienti corretti.
        """
        # Dividi il ket in coppie di simboli
        pairs = [ket[i:i + 2] for i in range(1, 9, 2)]  # Ottiene solo i caratteri interni come coppie

        # Mappa le coppie in C^3 usando l'isomorfismo e raccoglie i coefficienti
        c3_representation = []
        total_coeff = 1  # Inizializza il coefficiente totale a 1
        for pair in pairs:
            state, coeff = self.isomorphism[pair]
            c3_representation.append(str(state))  # Aggiungi lo stato |A|B|C|D>
            total_coeff *= coeff  # Moltiplica i coefficienti di normalizzazione

        # Ritorna il vettore ket nella forma |A|B|C|D> insieme al coefficiente complessivo
        return f"|{'|'.join(c3_representation)}>", total_coeff

    def vector_to_ket(self, vector):
        """Mappa un vector nella sua rappresentazione nella base dei kets con il mapping C^2 ⊗ C^2 -> C^3."""
        # Controlla se tutti i coefficienti sono nulli
        if all(coeff == 0 for coeff in vector):
            return "null vector"

        ket_representation = []
        for i, coeff in enumerate(vector):
            if coeff != 0:
                c3_ket, normalization_coeff = self.map_ket_to_c3(self.kets[i])  # Mappa il ket nello spazio C^3
                final_coeff = coeff.real * normalization_coeff  # Applica il coefficiente di normalizzazione
                final_coeff = self.map_to_sqrt2(final_coeff)
                ket_representation.append(f"{final_coeff} * {c3_ket}") #è una lista di stringhe

        return " + ".join(ket_representation)




    def ket_sum(self, vector1, vector2):
        """
            Somma due vettori ket.

            :param vector1: Il primo vettore ket come array numpy.
            :param vector2: Il secondo vettore ket come array numpy.
            :return: Un nuovo vettore ket che rappresenta la somma.
            """
        if vector1.shape != vector2.shape:
            raise ValueError("Vector must be of the same dimension to be summed.")
        return self.vector_to_ket(vector1 + vector2)