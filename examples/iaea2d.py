from scarabee import *
import numpy as np
import matplotlib.pyplot as plt

Et = np.array([0.222222, 0.833333])
D = 1. / (3. * Et)
Ea = np.array([0.010120, 0.080032])
vEf = np.array([0., 0.135])
Ef = np.array([0., 0.135])
chi = np.array([1., 0.])
Es = np.array([[0.00, 0.02],
               [0.00, 0.00]])
a1 = DiffusionCrossSection(D, Ea, Es, Ef, vEf, chi, "A1")

Et = np.array([0.222222, 0.833333])
D = 1. / (3. * Et)
Ea = np.array([0.010120, 0.085032])
vEf = np.array([0., 0.135])
Ef = np.array([0., 0.135])
chi = np.array([1., 0.])
Es = np.array([[0.00, 0.02],
               [0.00, 0.00]])
a2 = DiffusionCrossSection(D, Ea, Es, Ef, vEf, chi, "A2")

Et = np.array([0.222222, 0.833333])
D = 1. / (3. * Et)
Ea = np.array([0.010120, 0.130032])
vEf = np.array([0., 0.135])
Ef = np.array([0., 0.135])
chi = np.array([1., 0.])
Es = np.array([[0.00, 0.02],
               [0.00, 0.00]])
a3 = DiffusionCrossSection(D, Ea, Es, Ef, vEf, chi, "A3")

Et = np.array([0.166667, 1.111111])
D = 1. / (3. * Et)
Ea = np.array([0.000160, 0.010024])
Es = np.array([[0.00, 0.04],
               [0.00, 0.00]])
a4 = DiffusionCrossSection(D, Ea, Es, "A4")

#        Half Assemblies down this column
#        |
#        V
tiles = [a3, a2, a2, a2, a3, a2, a2, a1, a4, # <- Half Assemblies in along this row
         a2, a2, a2, a2, a2, a2, a2, a1, a4,
         a2, a2, a2, a2, a2, a2, a1, a1, a4,
         a2, a2, a2, a2, a2, a2, a1, a4, a4,
         a3, a2, a2, a2, a3, a1, a1, a4, 0.,
         a2, a2, a2, a2, a1, a1, a4, a4, 0.,
         a2, a2, a1, a1, a1, a4, a4, 0., 0.,
         a1, a1, a1, a4, a4, a4, 0., 0., 0.,
         a4, a4, a4, a4, 0., 0., 0., 0., 0.]

dx = np.array([10., 20., 20., 20., 20., 20., 20., 20., 20.])
nx = np.array([9,   17,  17,  17,  17,  17,  17,  17,  17])

dy = np.array([20., 20., 20., 20., 20., 20., 20., 20., 10.])
ny = np.array([17,  17,  17,  17,  17,  17,  17,  17,  9])

geom = DiffusionGeometry(tiles, dx, nx, dy, ny, 1., 0., 0., 1.)
solver = FDDiffusionDriver(geom)
solver.solve()
flux, x, y, _ = solver.flux()

plt.pcolormesh(y, x, flux[0,:,:], cmap='jet')
plt.title("Group 1")
plt.show()

plt.pcolormesh(y, x, flux[1,:,:], cmap='jet')
plt.title("Group 2")
plt.show()

power, x, y, _ = solver.power()
plt.pcolormesh(y, x, power, cmap='jet')
plt.title("Power")
plt.show()
