from sympy import pprint
from sympy.physics.quantum import Ket
from basis_spaces.vector_spaces import VectorSpace

# Example usage
if __name__ == "__main__":
    # Create quantum states using SymPy's Ket
    plus, minus = Ket("+"), Ket("-")
    
    # Create 2D Hilbert space
    C2 = VectorSpace([plus, minus], name="C")
    
    print("Single space:", C2)
    print()
    
    # Create tensor product of 4 spaces
    H4 = C2 @ C2 @ C2 @ C2
    
    print("Tensor product:", H4)
    print()
    
    # Create elements
    state1 = C2.get_element(1, 2)  # |+> + 2|−>
    state2 = C2.get_element(2, 1)  # 2|+> + |−>
    state3 = C2.get_element(1, 0)  # |+>
    state4 = C2.get_element(0, 1)  # |−>
    
    for j, state in enumerate([state1, state2, state3, state4]):
        print(f"State {j}:")
        pprint(state)
    print()
    
    # Tensor product of states
    tensor_state = H4.get_tensor_element(state1, state2, state3, state4)
    print("Tensor product of states:")
    pprint(tensor_state)
    print()
    
    # Show some basis states - now works with pprint!
    print("Basis states examples:")
    for i in [0, 1, 5, 10, 15]:
        print(f"Basis[{i}]:")
        print(H4.basis[i])
        print(f"Vector: {H4.elements_dict[H4.basis[i]].T}")
        print()


#  Mostra il branch Git corrente
# function parse_git_branch() {
#   git branch 2>/dev/null | grep '^*'$ß
# }

# function venv_name() {
#   if [ -n "$VIRTUAL_ENV" ]; then
#     if [ -f ".venv-name" ]; then
#       cat .venv-name
#     else
#       echo "${VIRTUAL_ENV:t}"
#     fi
#   fi
# }

# # Prompt colorato
# autoload -Uz colors && colors
# setopt PROMPT_SUBST
# export PS1="$(venv_name) %F{green}%n$