#!/usr/bin/env python

import numpy as np

ARB_TO_EV = 14.39964535
alpha = 1.0
N = 20
L = 11.

prefactor = ARB_TO_EV / (2*np.pi**2)
print(prefactor)

value = 0.0
for x in np.arange(-N,N+1):
    for y in np.arange(-N,N+1):
        for z in np.arange(-N,N+1):

            if x==0 and y==0 and z==0:
                continue

            # print(x,y,z)

            qx = 2*np.pi / L * x
            qy = 2*np.pi / L * y
            qz = 2*np.pi / L * z

            # print(qx,qy,qz, ARB_TO_EV / (qx**2+qy**2+qz**2))
            # print(qx*L/(2*np.pi),qy*L/(2*np.pi),qz*L/(2*np.pi))

    
            value += np.exp(- (qx**2+qy**2+qz**2)/(2*alpha)) / (qx**2+qy**2+qz**2)

value = value * prefactor

print(value)