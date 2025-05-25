
from src.QuantumGeometry import SpinTetrahedron
from src.utilia import bcolors
from sympy import pprint, sqrt, simplify

s0 = SpinTetrahedron(1/2, 1/2, 1/2, 1/2)
# ! class authomatically prints Inv(ρ)-spinnet basis

print()
print(f"{bcolors.BLACK}={bcolors.ENDC}"*60)
print(f"{bcolors.BLACK}VERIFYING SECTION 5 `Quantum Geometry of Space` from `docs`{bcolors.ENDC}")
print(f"{bcolors.BLACK}={bcolors.ENDC}"*60)


print(f"\n{bcolors.BOLD}{bcolors.BLACK}Section 5.1.{bcolors.ENDC}")
print(f"{bcolors.BLACK}-{bcolors.ENDC}"*12)

# By the theory E1=1/2 v1 and E2=sqrt(3)/6 v1 + sqrt(3)/3 v2
E1, E2 = 1/2*s0.Inv[0], simplify(sqrt(3)/6 * s0.Inv[0] + sqrt(3)/3 * s0.Inv[1])


print(f'\n{bcolors.ITALIC}Orthonormal spinnet basis:{bcolors.ENDC}\n\nE1 =')
pprint(E1)
print('\nE2 =')
pprint(E2)

print(f'\nDihedral angle {bcolors.BOLD}{bcolors.ORANGE}D_{{1,2}}{bcolors.ENDC} in (E1, E2)\n')
D_12_orth = s0.dihedral_matrix(0,1, basis=[E1, E2])
pprint(D_12_orth)

for alpha in range(4):
    print(f'\nSquared area {bcolors.BOLD}{bcolors.RED}D_{{{alpha},{alpha}}}{bcolors.ENDC} in (E1, E2)\n')
    pprint(s0.dihedral_matrix(alpha, alpha, basis=[E1, E2]))

print(f'\nSquared volume {bcolors.BOLD}{bcolors.GREEN}i|[D_{{1,3}},D_{{1,2}}]|{bcolors.ENDC} in (E1, E2):\n')
pprint(s0.volume_squared(coeff=1, basis=[E1, E2]))
print(f'\n{bcolors.ITALIC}{bcolors.BLACK}fuzzy{bcolors.ENDC} ✅')

#--

print(f"\n\n{bcolors.BOLD}{bcolors.BLACK}Section 5.2. Canonical LQG{bcolors.ENDC}")
print(f"{bcolors.BLACK}-{bcolors.ENDC}"*26)

e0, e1 = simplify(E2 - E1), simplify(E1 + E2)

print(f'\n{bcolors.ITALIC}Eigenbasis for the squared volume:{bcolors.ENDC}\n\ne0 = E2 - E1 =')
pprint(e0)
print('\ne1 = E1 + E2 =\n')
pprint(e1)

# Prima calcolo vol sulla base (E1, E2)
vol = s0.volume_squared(coeff=1, basis=[E1,E2])

# Poi cambio base vs la sua eigenbasis
vol_diag = s0.change_basis_matrix(vol, [E1, E2], [e0,e1])

print(f'\nSquared volume {bcolors.BOLD}{bcolors.GREEN}i|[D_{{1,3}},D_{{1,2}}]|{bcolors.ENDC} in [e0, e1]\n')
pprint(vol_diag)

for alpha in range(4):
    print(f'\nSquared area {bcolors.BOLD}{bcolors.RED}D_{{{alpha},{alpha}}}{bcolors.ENDC} in [e0, e1]\n')
    pprint(s0.dihedral_matrix(alpha, alpha, basis=[e0,e1]))

print(f'\nDihedral angle {bcolors.BOLD}{bcolors.ORANGE}D_{{1,2}}{bcolors.ENDC} in [e0, e1]:\n')
pprint(s0.dihedral_matrix(0,1, basis=[e0, e1]))
print(f'\n{bcolors.ITALIC}{bcolors.BLACK}fuzzy{bcolors.ENDC} ✅')


print('\n'+f"{bcolors.BLACK}={bcolors.ENDC}"*30)
print(f"{bcolors.BLACK}TEST COMPLETED SUCCESSFULLY!{bcolors.ENDC}")
print(f"{bcolors.BLACK}={bcolors.ENDC}"*30)