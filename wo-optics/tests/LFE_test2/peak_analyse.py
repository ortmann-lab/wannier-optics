#!/usr/bin/env python
"""
Analyses the first peaks of the spectrum
"""

import numpy as np
import matplotlib.pyplot as plt
# from scipy.signal import savgol_filter

# the hamiltonian has block diagonal form
# 1st block
matrix = np.zeros((4,4), dtype=np.complex128)
matrix[0,1] = 1.0
matrix[1,0] = 1.0
matrix[0,2] = 0.10+1.j*1.0
matrix[2,0] = 0.10-1.j*1.0
matrix[0,3] = -0.10+1.j*2.0
matrix[3,0] = -0.10-1.j*2.0
matrix[1,2] = 2.0 - 1.j*0.1
matrix[2,1] = 2.0 + 1.j*0.1
matrix[2,3] = 2.5 - 1.j*0.5
matrix[3,2] = 2.5 + 1.j*0.5

matrix[0,0] = 3.1

matrix *=2.  # LFE have prefactor 2 in the Hamiltonian

# onsite energies
for i in range(4):
    matrix[i,i] += 2.0

# transfer integrals
matrix[0,2] +=  0.5
matrix[2,0] +=  0.5
matrix[1,3] +=  0.5
matrix[3,1] +=  0.5

matrix[0,1] +=  0.3
matrix[1,0] +=  0.3
matrix[2,3] +=  0.3
matrix[3,2] +=  0.3

print(matrix)
eval, evec = np.linalg.eigh(matrix)
print(eval)

# 2nd block
matrix2 = np.zeros((4,4), dtype=np.complex128)
matrix2[0,2] +=  0.5  # single particle TIs
matrix2[2,0] +=  0.5
matrix2[1,3] +=  0.5
matrix2[3,1] +=  0.5

matrix2[0,1] +=  0.3
matrix2[1,0] +=  0.3
matrix2[2,3] +=  0.3
matrix2[3,2] +=  0.3

matrix2[1,1] +=  -1.7 *2 # LFE
matrix2[1,3] +=  (1.1 + 1.j*2.0) *2
matrix2[3,1] +=  (1.1 - 1.j*2.0) *2

# onsite energies
for i in range(4):
    matrix2[i,i] += 2.0

print(matrix2)
eval2, evec = np.linalg.eigh(matrix2)
print(eval2)


# other blocks
matrix3 = np.zeros((4,4), dtype=np.complex128)
matrix3[0,2] +=  0.5
matrix3[2,0] +=  0.5
matrix3[1,3] +=  0.5
matrix3[3,1] +=  0.5

matrix3[0,1] +=  0.3
matrix3[1,0] +=  0.3
matrix3[2,3] +=  0.3
matrix3[3,2] +=  0.3

# onsite energies
for i in range(4):
    matrix3[i,i] += 2.0

print(matrix3)
eval3, evec = np.linalg.eigh(matrix3)
print(eval3)


dos = np.loadtxt("output/T020Kdensity.dat")
gamma = 0.5

# rescale energy
#dos[:,0] = dos[:,0]*gamma

# plot
plt.semilogy(dos[:,0], dos[:,1])

# plt.plot([2,2], [1e-4,1e2], color='C1')  # onsite energy
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

