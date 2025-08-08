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

matrix[1,1] += -1.5
matrix[2,2] += 1.36

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

dos = np.loadtxt("output/T020Kdensity.dat")
gamma = 1.0

# rescale energy
#dos[:,0] = dos[:,0]*gamma



# plot
plt.semilogy(dos[:,0], dos[:,1])

for e in eval:
    plt.plot([e,e], [1e-4,1e2], color='k')


plt.xlabel("Energy (eV)")
plt.ylabel("Exciton DOS (arb.)")
plt.tight_layout()
plt.show()

