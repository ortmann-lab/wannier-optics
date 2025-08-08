#!/usr/bin/env python
"""
Plots the exciton density of states (DOS) from the output of wo-optics.x and compares
it with the analytical result.
"""
import os.path as osp
import numpy as np
import matplotlib.pyplot as plt
ARB_TO_EV = 14.39964535 # e^2/(4pi epsilon_0) in eV*Angström


# model parameters
t = 1.0             # hopping parameter (transfer integral) of the tight-binding models in eV
e0 = 13.0           # onsite-energy of tight-binding model for conduction bands in eV
                    # (determines bandgap)
a = 15.0            # lattice constant in Angström
epsilon = 1.        # relative dielectric constant

file_DOS = "output/T000Kdensity.dat"  # output file from wo-optics.x


# read exciton DOS
if not osp.exists(file_DOS):
    print(f"No exciton DOS found. Please run wo-optics.x and make sure that {file_DOS} exists.")
    exit(1)
dos = np.loadtxt(file_DOS)


# calculate exciton energies analytically
bandedge = -12*t+e0
rydberg = ARB_TO_EV**2 / (8.*t* a**2 * epsilon**2)
bohr_radius = 4*t*a**2*epsilon / ARB_TO_EV
corrections = -rydberg**2/(24.* t) * 5

# print summary
print("Model parameters:")
print(f"  Hopping parameter = {t} eV")
print(f"  Onsite energy     = {e0} eV")
print(f"  Lattice constant  = {a} Angström")
print(f"  Rel. dielectric constant = {epsilon}")
print("\nResults:")
print(f"  Bandgap                = {bandedge:.3f} eV")
print(f"  Exciton-Rydberg-Energy = {rydberg:.3f} eV")
print(f"  Exciton Bohr radius    = {bohr_radius:.3f} Angström")
print(f"                         = {bohr_radius/a:.3f} Lattice constants")
print(f"  k^4 corrections (n=1, l=0)  = {corrections:.3f} eV")


# plot exciton DOS and compare with analytical result (Rydberg series)
plt.semilogy(dos[:,0], dos[:,1])
print("\nExcition energies:")
for n in np.arange(1,5):  # quantum number n
    for l in np.arange(0,n): # angular momentum quantum number from 0 to n
        correction = rydberg**2/(n**4)/(24.* t) * (4*n/(l+0.5) - 3)
        peak = bandedge - (rydberg/(n**2) - correction)
        print(f"  n={n}, l={l}, energy={peak:.5f} eV")
        plt.plot([peak, peak], [0,1.], color='k', ls='--')

plt.ylim([1e-7,1e-2])
plt.xlim([bandedge-rydberg-0.01, bandedge+0.1])
plt.xlabel("Energy (eV)")
plt.ylabel("Exciton DOS (arb.)")
plt.tight_layout()
plt.show()
