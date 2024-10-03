# Quantum Geometry of Tetrahedra

You can find here a little pdf file where I tried to describe the scenario of quantum geometry of space, offered by the Loop Quantisation of Gravity (LQG), well discussed in  my master's thesis, as user-friendly as possible.

Thereafter, the quantum geometry of space is modeled by the theory in a discrete lattice made of spin tetrahedra and one is able to describe the quantum geometric properties of space as observables of the theory: they reduce to some self-adjoint operators among finite dimensional vector spaces with real spectra containing the possible outcomes of the measurements, in a perfect quantum flavor.

Everything lies within the framework of multilinear algebra and here you can find Python code computing such matricial geometric observables of the gravitational field.



### 1. **Theoretical Context**
   - **Quantum Theory**: Gravitational field with gauge symmetries $\text{SU}(2)$ and diffeomorphisms.
   - **Dynamic Variable**: Discretized connection on a lattice $\Gamma$ with  N nodes and L links, mapped to $\text{SU}(2)^L$.
   - **Quantum States**: Functionals of discrete connections, expressed in the Peter-Weyl product basis in terms of the fundamental representation $\rho^{1/2}$ of $\text{SU}(2)$.

### 2. **Definition of Representations**
   - Product representation $\rho = \bigotimes_{i=1}^{k} (\rho^{j_i})$ identifying the lattice.
   - Use of the theorem $\rho^{j} = \text{Sym}^{2j}$, being $\text{Sym}^k:=\bigodot_{i=1}^k\rho^{1/2}:\text{SU}(2)\to\text{GL}(\mathbb{C}^{k+1})$.

### 3. **Lie Operators**
   - **Construction**: Utilizing the Pauli matrices $\sigma_a$ and defining $\tau_a = -\frac{i}{2} \sigma_a$.
   - **Form**: Lie operators ${L_{(node, link)}}_a$ as tensor products of identities and $\tau_a$.

### 4. **Spin Networks**
   - An element of a Hilbert space, $\text{Inv}(\rho) = \{ v \in V \|\ \rho(U)v = v,\ \forall U \in \text{SU}(2)\}$.
   - Dimension of $\text{Inv}(\rho)$ for the ground state quantum tetrahedron: generated by vectors $v_1$ and $v_2$.

### 5. **Gauge Invariant Operators**
   - Defined as $D_{\alpha\beta} = {L_{\alpha}}_a \delta^{ab} {L_{\beta}}_b$.
   - Restriction to the subspace $\text{span}(v_1, v_2)$.

### 6. **Goal of Generalization**
   - **Generalize** to tetrahedra with an arbitrary number of legs.
   - **Develop** a program to calculate gauge invariant operators and Lie operators in a multilinear algebra context.

### 7. **Steps for Program Development**
   - **Define Data Structures**: Classes for the lattice, representations, and operators.
   - **Construct Operators**: Implement logic to compute Lie operators and gauge invariant operators.
   - **Calculate States**: Generate and manipulate quantum states in the Peter-Weyl basis.
   - **Multilinear Algebra**: Use libraries like NumPy or SciPy for matrix and tensor operations.
   - **Testing and Validation**: Write unit tests and utilize visualization methods to verify results.
