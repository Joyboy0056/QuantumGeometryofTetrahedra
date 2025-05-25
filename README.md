# Quantum Geometry of Tetrahedra

You can find here a [little pdf](https://github.com/Joyboy0056/QuantumGeometryofTetrahedra/docs/Loop_Quantisation_of_Space.pdf) file where I tried to describe the scenario of quantum geometry of space, offered by the Loop Quantisation of Gravity (LQG), well discussed in  my master's thesis, as user-friendly as possible.

Thereafter, the quantum geometry of space is modeled by the theory in a discrete lattice made of spin tetrahedra and one is able to describe the quantum geometric properties of space as observables of the theory: they reduce to some self-adjoint operators among finite dimensional vector spaces with real spectra containing the possible outcomes of the measurements, in a perfect quantum flavor.

Everything lies within the framework of multilinear algebra and here you can find Python code computing such matricial geometric observables of the gravitational field.



### 1. **Theoretical Context**
   - **Quantum Theory**: Gravitational field with gauge symmetries $\text{SU}(2)$ and diffeomorphisms.
   - **Dynamic Variable**: Discretized connection on a lattice $\Gamma$ with  $N$ nodes and $L$ links, i.e. a map in $\text{SU}(2)^L$.
   - **Quantum States**: Functionals of discrete connections, expressed in the Peter-Weyl product basis in terms of the fundamental representation $\rho^{1/2}$ of $\text{SU}(2)$.

### 2. **Definition of Representations**
   - Product representation $\rho = \bigotimes_{i=1}^{k} (\rho^{j_i})$ identifying the lattice $\Gamma\leftrightarrow(j_1,...,j_k)$ with $1$ node and $k$ links.
   - Use of the theorem $\rho^{j} = \text{Sym}^{2j}$, being $\text{Sym}^k:=\bigodot_{i=1}^k\rho^{1/2}: \text{SU}(2)\to\text{GL}(\mathbb{C}^{k+1})$.

### 3. **Lie Operators**
   - **Construction**: Utilizing the Pauli matrices $\sigma_a$ and defining $T\rho(\tau_a) = -\frac{i}{2} T\rho(\sigma_a)$, $T$ being the tangent map.
   - **Form**: Lie operators ${L_{(node, link)}}_a$ as tensor products of identities and $\tau_a$ in the $T\rho$ representation.

### 4. **Spin Networks**
   - Spinnets are elements of a Hilbert space $H_{(j_1,...,j_k)}$ corresponding to an isotropic subspace $\text{Inv}(\rho) = {v \in V \|\ \rho(U)v = v,\ \forall U \in \text{SU}(2)}$, being $V\cong\bigotimes_{i=1}^k\mathbb{C^{2j_i+1}}$ the support vector space of $\rho$.
   - Dimension of $\text{Inv}\left({\rho^{1/2}}^{\otimes4}\right)$ for the ground state quantum tetrahedron is 2.
   - Dimension of $\text{Inv}\left({\rho^1}^{\otimes4}\right)$ is instead 3.
   - In general, dimension of $\text{Inv}\left({\rho^j}^{\otimes4}\right)$ turns out to be $2j+1$.

### 5. **Gauge Invariant Operators**
   - Also called geometric operators, are defined as $D_{\alpha\beta} = (L_\alpha)^a\delta_{ab}(L_\beta)^b$.
   - Restriction to the subspace $\text{Inv}(\rho)$.

### 6. **Goal of Generalization**
   - **Generalize** to tetrahedra with an arbitrary number of legs.
   - **Develop** a program to calculate gauge invariant operators and Lie operators in a multilinear algebra context.

## Fruibility
You can mainly refer to the folder [lqg](https://github.com/Joyboy0056/QuantumGeometryofTetrahedra/tree/main/lqg), that's where you can find a comprehensible workflow of the theory of quantum spin $(j,j,j,j)$ tetrahedra. In particular, there you will spot my libraries on Kets and Quantum Operators for both the cases $j=\frac{1}{2}$ and $j=1$.
