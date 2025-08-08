#!/usr/bin/env python
"""
Compares the spectrum to the experiment
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams["savefig.directory"] = os.path.dirname(os.getcwd())
matplotlib.rcParams.update({'font.size': 16})


print("Reading sprectrum.dat and Si_exp.csv")
data = np.loadtxt("spectrum.dat")
data_exp = np.loadtxt("Si_exp.csv", usecols=[0,1], delimiter=',')

qp_shift = 0.9
print(f"Apply quasi-particle shift of {qp_shift} eV to the spectrum.")

# create figure
fig = plt.figure(figsize=(4.8,6))
ax1 = fig.add_subplot(111)
ax1.plot(data_exp[:,0], data_exp[:,1], '--',lw=2, color="C1", label="Experiment")
ax1.plot(data[:,0]+qp_shift, data[:,1], lw=2, color="k", label="MLWF (exciton)")
ax1.set_xlabel("Energy (eV)")
ax1.set_ylabel("Im $ϵ(\omega)$")
ax1.set_xlim([2.5,5.5])
ax1.set_ylim([0,75])
plt.legend(fontsize=13)
plt.tight_layout()
plt.savefig("comparison_Si.png")
plt.show()