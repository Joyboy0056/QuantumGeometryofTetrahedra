from typing import List, Optional, Tuple
from sympy import Symbol, Matrix

from functools import reduce
from operator import mul
import re

def _format_dimension_superscript(dim: int) -> str:
    """Convert dimension to superscript notation"""
    superscripts = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    return str(dim).translate(superscripts)

class TensorProduct(Symbol):
    """Custom symbol for tensor products that displays with ⊗"""
    
    def __new__(cls, *args):
        # Check if all factors are ket states to use compact notation
        if all(cls._is_ket(arg) for arg in args):
            # Compact notation: |+> ⊗ |-> ⊗ |+> becomes |+−+>
            ket_content = ''.join(cls._extract_ket_content(arg) for arg in args)
            compact_name = f"|{ket_content}>"
            obj = Symbol.__new__(cls, compact_name)
        else:
            # Standard notation with ⊗
            obj = Symbol.__new__(cls, '⊗'.join(str(arg) for arg in args))
        
        obj.factors = args
        obj.is_compact_ket = all(cls._is_ket(arg) for arg in args)
        return obj
    
    @staticmethod
    def _is_ket(expr) -> bool:
        """Check if expression is a ket state like |+> or |0>"""
        s = str(expr)
        # Match patterns like |+>, |−>, |0>, |1>, etc.
        return bool(re.match(r'^\|.+\>$', s))
    
    @staticmethod
    def _extract_ket_content(expr) -> str:
        """Extract content from ket: |+> -> +"""
        s = str(expr)
        match = re.match(r'^\|(.+)\>$', s)
        if match:
            return match.group(1)
        return str(expr)
    
    def __repr__(self):
        if self.is_compact_ket:
            return str(self)
        return ' ⊗ '.join(str(f) for f in self.factors)
    
    def __str__(self):
        if self.is_compact_ket:
            ket_content = ''.join(self._extract_ket_content(f) for f in self.factors)
            return f"|{ket_content}>"
        return ' ⊗ '.join(str(f) for f in self.factors)
    
    def _latex(self, printer=None):
        if self.is_compact_ket:
            ket_content = ''.join(self._extract_ket_content(f) for f in self.factors)
            return f"|{ket_content}\\rangle"
        return r' \otimes '.join(str(f) for f in self.factors)

class VectorSpace:
    """A class for a vector space"""
    def __init__(
            self,
            basis: List[Symbol], 
            name: Optional[str] = None, 
            validate: bool = True,
            show_dimension: bool = True
    ):
        """
            Initialize a vector space with its basis.
            :param basis: list of SymPy vectors for the basis
            :param name: optional name for the vector space
            :param validate: whether to validate the basis
        """
        self.basis = basis
        self.dimension = len(basis)
        self.base_name = name or "V"
        self.show_dimension = show_dimension

        # Create full name with superscript dimension
        if show_dimension and self.dimension > 0:
            self.name = f"{self.base_name}{_format_dimension_superscript(self.dimension)}"
        else:
            self.name = self.base_name

        if validate and self.dimension > 0:
                if len(set(basis)) != self.dimension:
                    raise ValueError("Basis vectors must be distinct")

        self.vector_basis = [
            Matrix.eye(self.dimension)[:,j] for j in range(self.dimension)
        ]
        self.basis_dict = {
            self.basis[j]: self.vector_basis[j] for j in range(self.dimension)
        }
        self.elements_dict = {}

    # Standard methods
    def _compute_symbolic(self, coeffs: Tuple) -> Symbol:
        """Compute symbolic representation"""
        return sum(coeff * svec for coeff, svec in zip(coeffs, self.basis))
    
    def _compute_vectorial(self, coeffs: Tuple) -> Matrix:
        """Compute vectorial representation"""
        return Matrix([coeffs]).T
    
    def get_vector_repr(self, element):
        """Get the vector representation of a symbolic element"""
        return self.elements_dict.get(element)

    def get_element(self, *coeffs) -> Symbol:
        """
            Returns the linear combination of the basis vectors with the inpunt coefficients
            and records both its symbolic and vectorial representation in the "elements_dict"
        """
        if len(coeffs) != self.dimension:
            raise ValueError(f"Expected {self.dimension} coefficients, got {len(coeffs)}")
        
        element = self._compute_symbolic(tuple(coeffs))
        
        if element not in self.elements_dict:
            vector_element = self._compute_vectorial(tuple(coeffs))
            self.elements_dict[element] = vector_element

        return element
    

    # ------- Dunder methods -------
    def __repr__(self):
        return f"VectorSpace(name='{self.name}')"
    
    def __str__(self):
        basis_str = ', '.join(str(b) for b in self.basis)
        return f"{self.name}: basis=[{basis_str}]"
    
    def __len__(self):
        return self.dimension
    
    def __getitem__(self, index):
        """Access basis vectors by index"""
        return self.basis[index]
    
    def __contains__(self, element):
        """Check if element is in the basis"""
        return element in self.basis

    def __eq__(self, other):
        if not isinstance(other, VectorSpace):
            return False
        return self.dimension == other.dimension and set(self.basis) == set(other.basis)
    
    def __iter__(self):
        """Iterate over basis vectors"""
        return iter(self.basis)
    
    def __matmul__(self, other: 'VectorSpace') -> 'TensorProductSpace':
        """
        Tensor product using @ operator
        V1 @ V2 creates the tensor product space
        """
        return TensorProductSpace(self, other)
    

class TensorProductSpace(VectorSpace):
    """Tensor product of vector spaces"""
    
    def __init__(self, *spaces: VectorSpace):
        """
        Create tensor product of vector spaces.
        :param spaces: VectorSpace instances to tensor together
        """
        if len(spaces) < 2:
            raise ValueError("Need at least 2 spaces for tensor product")
        
        # Flatten nested tensor products
        self.factor_spaces: List[VectorSpace] = self._flatten_spaces(spaces)
        self.num_factors = len(self.factor_spaces)
        
        # Generate tensor product basis using SymPy's TensorProduct
        tensor_basis = self._generate_tensor_basis()
        
        # Compute dimension (product of individual dimensions)
        dimension = reduce(mul, (space.dimension for space in self.factor_spaces))
        
        # Create name: C² ⊗ C³
        base_name = ' ⊗ '.join(space.name for space in self.factor_spaces)
        
        # Set attributes manually
        self.basis = tensor_basis
        self.dimension = dimension
        self.base_name = base_name
        self.name = base_name
        self.show_dimension = False
        
        # Create correct vector_basis for tensor product space
        self.vector_basis = [
            Matrix.eye(self.dimension)[:,j] for j in range(self.dimension)
        ]
        
        # Create basis_dict and populate elements_dict with basis vectors
        self.basis_dict = {
            self.basis[j]: self.vector_basis[j] for j in range(self.dimension)
        }
        self.elements_dict = {
            self.basis[j]: self.vector_basis[j] for j in range(self.dimension)
        }
        
        # Store factor dimensions for later use
        self.factor_dimensions = [space.dimension for space in self.factor_spaces]
    
    def _flatten_spaces(self, spaces):
        """Flatten nested TensorProductSpace into list of basic VectorSpaces"""
        flattened = []
        for space in spaces:
            if isinstance(space, TensorProductSpace):
                flattened.extend(space.factor_spaces)
            else:
                flattened.append(space)
        return flattened
    
    def _generate_tensor_basis(self) -> List[TensorProduct]:
        """Generate basis for tensor product space using SymPy's TensorProduct"""
        # Start with first space basis
        current_basis = [[b] for b in self.factor_spaces[0].basis]
        
        # Iteratively tensor with each subsequent space
        for space in self.factor_spaces[1:]:
            new_basis = []
            for existing in current_basis:
                for new_b in space.basis:
                    new_basis.append(existing + [new_b])
            current_basis = new_basis
        
        # Create SymPy TensorProduct objects
        tensor_basis = []
        for basis_combo in current_basis:
            tensor_symbol = TensorProduct(*basis_combo)
            tensor_basis.append(tensor_symbol)
        
        return tensor_basis
    
    def get_tensor_element(self, *elements):
        """
        Create tensor product of elements from factor spaces.
        :param elements: one element from each factor space
        """
        if len(elements) != self.num_factors:
            raise ValueError(f"Expected {self.num_factors} elements, got {len(elements)}")
        
        # Get vector representations from factor spaces
        vectors = []
        for elem, space in zip(elements, self.factor_spaces):
            vec = space.get_vector_repr(elem)
            if vec is None:
                raise ValueError(f"Element {elem} not found in space {space.name}. "
                               f"Did you create it with space.get_element(...)?")
            vectors.append(vec)
        
        # Compute Kronecker product of vectors
        tensor_vector = self._kronecker_product(vectors)
        
        # Create symbolic representation using SymPy's TensorProduct
        tensor_symbol = TensorProduct(*elements)
        
        # Store in elements dict
        self.elements_dict[tensor_symbol] = tensor_vector
        
        return tensor_symbol
    
    def _kronecker_product(self, vectors: List[Matrix]) -> Matrix:
        """Compute Kronecker product of matrices/vectors"""
        result = vectors[0]
        for vec in vectors[1:]:
            # Kronecker product
            new_result = Matrix.zeros(result.rows * vec.rows, 1)
            for i in range(result.rows):
                for j in range(vec.rows):
                    new_result[i * vec.rows + j] = result[i] * vec[j]
            result = new_result
        return result
    
    def __repr__(self):
        return f"TensorProductSpace({self.name})"
    
    def __str__(self):
        basis_str = ', '.join(str(b) for b in self.basis)
        
        # Show total dimension
        total_dim_str = _format_dimension_superscript(self.dimension)
        return f"{self.name} (≅ {self.factor_spaces[0].base_name}{total_dim_str}): basis=[{basis_str}]"
