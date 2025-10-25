# 🐍 Specchietto Metodi Speciali Python (“Dunder Methods”)

> Guida rapida ai principali metodi speciali Python con esempi di utilizzo pratici.
> (*dunder* = “double underscore”, es. `__init__`, `__len__`...)

---

### 🧱 Costruzione e ciclo di vita

| Metodo                         | Quando viene usato                           | Esempio               |
| ------------------------------ | -------------------------------------------- | --------------------- |
| `__new__(cls, *a, **kw)`       | Creazione dell’istanza (prima di `__init__`) | `obj = MyClass()`     |
| `__init__(self, *a, **kw)`     | Inizializzazione dell’istanza                | `obj = MyClass(x, y)` |
| `__del__(self)`                | Distruttore (garbage collection)             | `del obj`             |
| `__init_subclass__(cls)`       | Quando nasce una sottoclasse                 | —                     |
| `__class_getitem__(cls, item)` | Supporto per `MyClass[T]`                    | —                     |

---

### 🪞 Rappresentazione e conversione

| Metodo                   | Usato in...                        | Esempio        |
| ------------------------ | ---------------------------------- | -------------- |
| `__repr__(self)`         | Rappresentazione ufficiale (debug) | `repr(obj)`    |
| `__str__(self)`          | Rappresentazione leggibile         | `print(obj)`   |
| `__format__(self, spec)` | Formattazione avanzata             | `f"{obj:.2f}"` |
| `__bytes__(self)`        | Conversione in bytes               | `bytes(obj)`   |

---

### 🔍 Attributi e accesso dinamico

| Metodo                         | Quando                                    | Esempio            |
| ------------------------------ | ----------------------------------------- | ------------------ |
| `__getattr__(self, name)`      | Attributo inesistente                     | `obj.missing_attr` |
| `__getattribute__(self, name)` | Ogni accesso (⚠️ attenzione a ricorsione) | —                  |
| `__setattr__(self, name, val)` | Assegnazione attributi                    | `obj.x = 1`        |
| `__delattr__(self, name)`      | Cancellazione attributo                   | `del obj.x`        |
| `__dir__(self)`                | Lista attributi                           | `dir(obj)`         |

---

### ⚙️ Protocollo Descriptor

| Metodo                            | Usato in...                   | Esempio        |
| --------------------------------- | ----------------------------- | -------------- |
| `__get__(self, inst, owner)`      | Lettura                       | `obj.attr`     |
| `__set__(self, inst, val)`        | Scrittura                     | `obj.attr = 1` |
| `__delete__(self, inst)`          | Cancellazione                 | `del obj.attr` |
| `__set_name__(self, owner, name)` | Alla definizione della classe | —              |

---

### 📦 Container e sequenze

| Metodo                        | Usato in...                   | Esempio          |
| ----------------------------- | ----------------------------- | ---------------- |
| `__len__(self)`               | Lunghezza                     | `len(obj)`       |
| `__getitem__(self, key)`      | Accesso                       | `obj[key]`       |
| `__setitem__(self, key, val)` | Assegnazione                  | `obj[key] = val` |
| `__delitem__(self, key)`      | Cancellazione                 | `del obj[key]`   |
| `__iter__(self)`              | Iterazione                    | `for x in obj:`  |
| `__reversed__(self)`          | Iterazione inversa            | `reversed(obj)`  |
| `__contains__(self, x)`       | Appartenenza                  | `x in obj`       |
| `__missing__(self, key)`      | Chiave mancante (dict custom) | —                |

---

### 🔁 Iteratori e generatori

| Metodo           | Usato in...                                                  | Esempio         |
| ---------------- | ------------------------------------------------------------ | --------------- |
| `__iter__(self)` | Inizia l’iterazione                                          | `for x in obj:` |
| `__next__(self)` | Restituisce elemento successivo                              | `next(obj)`     |
| *(Generatori)*   | Usano `yield` e supportano `.send()`, `.throw()`, `.close()` | —               |

---

### ☎️ Callability & Context Manager

| Metodo                           | Usato in...      | Esempio     |
| -------------------------------- | ---------------- | ----------- |
| `__call__(self, *a, **kw)`       | Oggetto callable | `obj()`     |
| `__enter__(self)`                | Inizio contesto  | `with obj:` |
| `__exit__(self, exc_t, exc, tb)` | Uscita contesto  | —           |

---

### ⚖️ Confronti e hash

| Metodo     | Operatore | Esempio     |
| ---------- | --------- | ----------- |
| `__lt__`   | `<`       | `a < b`     |
| `__le__`   | `<=`      | `a <= b`    |
| `__eq__`   | `==`      | `a == b`    |
| `__ne__`   | `!=`      | `a != b`    |
| `__gt__`   | `>`       | `a > b`     |
| `__ge__`   | `>=`      | `a >= b`    |
| `__hash__` | Hashing   | `hash(obj)` |

---

### ➕ Operatori aritmetici e bitwise

| Categoria | Metodi                                                                                             | Esempi                                                 |
| --------- | -------------------------------------------------------------------------------------------------- | ------------------------------------------------------ |
| Unari     | `__neg__`, `__pos__`, `__abs__`, `__invert__`                                                      | `-obj`, `+obj`, `abs(obj)`, `~obj`                     |
| Binari    | `__add__`, `__sub__`, `__mul__`, `__matmul__`, `__truediv__`, `__floordiv__`, `__mod__`, `__pow__` | `a + b`, `a @ b`, `a ** b`                             |
| Riflessi  | `__radd__`, `__rsub__`, ecc.                                                                       | Usati se l’operando sinistro non supporta l’operazione |
| In-place  | `__iadd__`, `__imul__`, ecc.                                                                       | `a += b`, `a *= b`                                     |
| Bitwise   | `__and__`, `__or__`, `__xor__`, `__lshift__`, `__rshift__` (+ riflessi/in-place)                   | `a & b`, `a << 2`                                      |

---

### 💾 Pickling e serializzazione

| Metodo                                | Scopo                                |
| ------------------------------------- | ------------------------------------ |
| `__getstate__`, `__setstate__`        | Stato durante pickling               |
| `__getnewargs__`, `__getnewargs_ex__` | Argomenti per `__new__` nel pickling |
| `__reduce__`, `__reduce_ex__`         | Controllo totale sul pickling        |

---

### 🧬 Classi e metaclassi

| Metodo                                   | Descrizione                                |
| ---------------------------------------- | ------------------------------------------ |
| `__prepare__`                            | Prepara il namespace per la classe         |
| `__mro_entries__`                        | Supporto PEP 560 (typing)                  |
| `__instancecheck__`, `__subclasscheck__` | Personalizzano `isinstance` / `issubclass` |
| `__call__` *(su metaclass)*              | Controlla l’istanziazione della classe     |

---

### ⚡ Async / Await

| Metodo            | Usato in...              | Esempio               |
| ----------------- | ------------------------ | --------------------- |
| `__await__(self)` | Supporto per `await obj` | `await obj`           |
| `__aiter__(self)` | Async iterator           | `async for x in obj:` |
| `__anext__(self)` | Prossimo elemento async  | `await anext(obj)`    |

---

### 🧩 Altri hook utili

| Metodo                     | Usato in...               | Esempio              |
| -------------------------- | ------------------------- | -------------------- |
| `__bool__(self)`           | Verità logica             | `if obj:`            |
| `__sizeof__(self)`         | Dimensione in byte        | `sys.getsizeof(obj)` |
| `__copy__`, `__deepcopy__` | Supporto al modulo `copy` | —                    |

---
