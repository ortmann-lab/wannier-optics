#!/usr/bin/env python
"""
Analyses the first peaks of the spectrum
"""

import numpy as np
import matplotlib.pyplot as plt
# from scipy.signal import savgol_filter


matrix = np.zeros((4,4), dtype=np.complex128)

# density-density interaction
matrix[0,0] = 1.0
matrix[1,1] = -3.0

# excitonic TIs
matrix[0,2] += 2.5 + 1.j
matrix[2,0] += 2.5 - 1.j

matrix[1,2] += -0.5 - 2.1j
matrix[2,1] += -0.5 + 2.1j

matrix[0,2] += 2.5 + 1.j
matrix[2,0] += 2.5 - 1.j

matrix[0,2] += .5
matrix[2,0] += .5

matrix[3,2] += -2.1j
matrix[2,3] += 2.1j

matrix[3,0] += 5 - 3.2j
matrix[0,3] += 5 + 3.2j


matrix = -matrix   # Coulomb matrix has a minus sign in the hamiltonian

for i in range(4):
    matrix[i,i] += 2.0  # onsite energies

print(matrix)

eval, evec = np.linalg.eigh(matrix)

print(eval)


# different S,S'
matrix2 = np.zeros((2,2), dtype=np.complex128)

matrix2[0,1] = -6.6 +1.j
matrix2[1,0] = -6.6 -1.j

matrix2 = -matrix2   # Coulomb matrix has a minus sign in the hamiltonian

for i in range(2):
    matrix2[i,i] += 2.0  # onsite energies

print(matrix2)
eval2, evec = np.linalg.eigh(matrix2)
print(eval2)


# different S=S'
matrix3 = np.zeros((4,4), dtype=np.complex128)

matrix3[0,2] += 2.0 -3.2j
matrix3[2,0] += 2.0 +3.2j

matrix3[0,2] += 1.0
matrix3[2,0] += 1.0

matrix3[1,0] += 2.0 -3.1j
matrix3[0,1] += 2.0 +3.1j

matrix3 = -matrix3   # Coulomb matrix has a minus sign in the hamiltonian

for i in range(4):
    matrix3[i,i] += 2.0  # onsite energies

print(matrix3)
eval3, evec = np.linalg.eigh(matrix3)
print(eval3)





dos = np.loadtxt("output/T020Kdensity.dat")
gamma = 1.0

# rescale energy
#dos[:,0] = dos[:,0]*gamma



# plot
plt.semilogy(dos[:,0], dos[:,1])
plt.plot([2,2], [1e-4,1e2], color='C1')  # onsite energies
for e in eval:
    plt.plot([e,e], [1e-4,1e2], color='k')

for e in eval2:
    plt.plot([e,e], [1e-4,1e2], color='C2')

for e in eval3:
    plt.plot([e,e], [1e-4,1e2], color='C3')

plt.xlabel("Energy (eV)")
plt.ylabel("Exciton DOS (arb.)")
plt.tight_layout()
plt.show()

