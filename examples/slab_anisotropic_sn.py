# Analytical Benchmark Test Set for Criticality Code and Verification
# Problem No.:32 [PUa-1-1-SL]
# Exact keff = 1.

from scarabee import *
import numpy as np
import matplotlib.pyplot as plt

Et = np.array([1.])
Ea = np.array([1. - 0.733333])
Es = np.array([[[0.733333]],
              [[0.333333]]])
Ef = np.array([0.266667])
vEf = np.array([2.5 * 0.266667])
chi = np.array([1.])

dtr = np.array([0.333333])

fuel = CrossSection(Et, dtr, Ea, Es, Ef, vEf, chi)

# Define Cells

thickness =  0.79606
Nx = 5000
# slab
dx = thickness / Nx

sn = ReflectorSN(Nx*[fuel], Nx*[dx], 128, True)
print(sn.anisotropic, sn.max_legendre_order)
sn.solve()

x = np.arange(start=0., stop=thickness, step=dx)
flux = np.zeros(Nx)
for i in range(Nx):
  flux[i] = sn.flux(i, 0)

plt.plot(x, flux)
plt.show()
