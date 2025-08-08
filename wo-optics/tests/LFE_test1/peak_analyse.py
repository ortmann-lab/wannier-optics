#!/usr/bin/env python
"""
Analyses the first peaks of the spectrum
"""

import numpy as np
import matplotlib.pyplot as plt
# from scipy.signal import savgol_filter


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

matrix[1,1] = 3.1

# matrix[0,1] = 1.0
# matrix[1,0] = 1.0
# matrix[0,2] = 1.0
# matrix[2,0] = 1.0
# matrix[0,3] = 1.0
# matrix[3,0] = 1.0
# matrix[1,2] = 1.0
# matrix[2,1] = 1.0
# matrix[2,3] = 1.0
# matrix[3,2] = 1.0

matrix *=2.  # LFE have prefactor 2 in the Hamiltonian

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

plt.plot([2,2], [1e-4,1e2], color='C1')  # onsite energy
plt.plot([3.1415*2+2.,3.1415*2+2.], [1e-4,1e2], color='C1')  # onsite energy (via LFE)
for e in eval:
    plt.plot([e,e], [1e-4,1e2], color='k')


plt.xlabel("Energy (eV)")
plt.ylabel("Exciton DOS (arb.)")
plt.tight_layout()
plt.show()

