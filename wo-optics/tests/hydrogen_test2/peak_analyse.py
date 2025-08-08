#!/usr/bin/env python
"""
Analyses the first peaks of the spectrum
"""

import numpy as np
import matplotlib.pyplot as plt
# from scipy.signal import savgol_filter

ARB_TO_EV = 14.39964535

dos = np.loadtxt("output/T000Kdensity.dat")
t = 1.0
a = 15.0
WW = 1.2
epsilon = 1./WW
e0 = 13.0

gamma = t
bandedge = -12*t+e0


rydberg = ARB_TO_EV**2 / (8.*t* a**2 * epsilon**2) #* 2
bohr_radius = 4*t*a**2*epsilon / ARB_TO_EV
corrections = -rydberg**2/(24.* t) * 5

print("Exciton-Rydberg-Energy = ", rydberg)
print("Exciton Bohr radius = ", bohr_radius, bohr_radius/a)
print("k^4 corrections (n=1, l=0) = ", corrections)

# rescale energy
#dos[:,0] = dos[:,0]*gamma

# get first peak
dos_prim = (dos[1:,1] - dos[:-1,1]) / dos[:-1,1]
dos_ = dos[:-1,:]
first_peak = dos_[dos_prim<0,0][0]

# first_peak = first_peak*4

plt.semilogy(dos[:,0], dos[:,1])

# for n in np.arange(1,7):
#     peak = first_peak/(n**2)
#     print(n, peak)
#     plt.plot([peak, peak], [0,1.])

# rydberg = -first_peak + bandedge #+ (first_peak**2/(12.* t) * (4*1/(0.5) - 3))
# for n in np.arange(1,5):
#     for l in np.arange(0,n):
#         correction = rydberg**2/(n**4)/(24.* t) * (4*n/(l+0.5) - 3)
#         peak = bandedge - (rydberg/(n**2) + correction)
#         print(n,l, peak)
#         plt.plot([peak, peak], [0,1.])


for n in np.arange(1,7):
    for l in np.arange(0,1): # 0 to n !!!!
        correction = 0.0 #rydberg**2/(n**4)/(24.* t) * (4*n/(l+0.5) - 3)
        peak = bandedge - (rydberg/(n**2) - correction)
        print(n,l, peak)
        plt.plot([peak, peak], [0,1.], color='k', ls='--')


# for n in np.arange(1,3):
#     peak = bandedge - rydberg/(n**2)
#     print(n, peak)
#     plt.plot([peak, peak], [0,1.])


# plt.plot([bandedge, bandedge], [0,1.], color='k')
# plt.plot([-bandedge, -bandedge], [0,1.], color='k')

plt.ylim([1e-7,1e-2])
plt.xlim([bandedge-rydberg-0.01, bandedge+0.1])
plt.xlabel("Energy (eV)")
plt.ylabel("Exciton DOS (arb.)")
plt.tight_layout()
plt.savefig("Exciton_DOS.png", dpi=300)
plt.show()

