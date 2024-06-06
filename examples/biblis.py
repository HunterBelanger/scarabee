from scarabee import *
import numpy as np
import matplotlib.pyplot as plt

Et = np.array([0.2315941, 0.9057971])
D = 1. / (3.*Et)
Ea = np.array([0.0102940, 0.0905100])
vEf = np.array([0.0064285, 0.1091100])
Ef = np.array([0.0064285, 0.1091100])
chi = np.array([1., 0.])
Es = np.array([[0., 0.0170270],
               [0., 0.]])
a8 = DiffusionCrossSection(D, Ea, Es, Ef, vEf, chi, "A8")

Et = np.array([0.2316584, 0.9060433])
D = 1. / (3.*Et)
Ea = np.array([0.0101650, 0.0880240])
vEf = np.array([0.0061908, 0.1035800])
Ef = np.array([0.0061908, 0.1035800])
chi = np.array([1., 0.])
Es = np.array([[0., 0.0171250],
               [0., 0.]])
a7 = DiffusionCrossSection(D, Ea, Es, Ef, vEf, chi, "A7")

Et = np.array([0.2317229, 0.9095043])
D = 1. / (3.*Et)
Ea = np.array([0.0101320, 0.0873140])
vEf = np.array([0.0064285, 0.1091100])
Ef = np.array([0.0064285, 0.1091100])
chi = np.array([1., 0.])
Es = np.array([[0., 0.0171920],
               [0., 0.]])
a6 = DiffusionCrossSection(D, Ea, Es, Ef, vEf, chi, "A6")

Et = np.array([0.2317873, 0.9095043])
D = 1. / (3.*Et)
Ea = np.array([0.0100030, 0.0848280])
vEf = np.array([0.0061908, 0.1035800])
Ef = np.array([0.0061908, 0.1035800])
chi = np.array([1., 0.])
Es = np.array([[0., 0.0172900],
               [0., 0.]])
a5 = DiffusionCrossSection(D, Ea, Es, Ef, vEf, chi, "A5")

Et = np.array([0.2316584, 0.9162544])
D = 1. / (3.*Et)
Ea = np.array([0.0103630, 0.0914080])
vEf = np.array([0.0074527, 0.1323600])
Ef = np.array([0.0074527, 0.1323600])
chi = np.array([1., 0.])
Es = np.array([[0., 0.0171010],
               [0., 0.]])
a4 = DiffusionCrossSection(D, Ea, Es, Ef, vEf, chi, "A4")

Et = np.array([0.2525253, 1.2025012])
D = 1. / (3.*Et)
Ea = np.array([0.0026562, 0.0715960])
Es = np.array([[0., 0.0231060],
               [0., 0.]])
a3 = DiffusionCrossSection(D, Ea, Es, "A3")

Et = np.array([0.2320293, 0.9167583])
D = 1. / (3.*Et)
Ea = np.array([0.0096785, 0.0784360])
vEf = np.array([0.0061908, 0.1035800])
Ef = np.array([0.0061908, 0.1035800])
chi = np.array([1., 0.])
Es = np.array([[0., 0.0176210],
               [0., 0.]])
a2 = DiffusionCrossSection(D, Ea, Es, Ef, vEf, chi, "A2")

Et = np.array([0.2321263, 0.9170105])
D = 1. / (3.*Et)
Ea = np.array([0.0095042, 0.0750580])
vEf = np.array([0.0058708, 0.0960670])
Ef = np.array([0.0058708, 0.0960670])
chi = np.array([1., 0.])
Es = np.array([[0., 0.0177540],
               [0., 0.]])
a1 = DiffusionCrossSection(D, Ea, Es, Ef, vEf, chi, "A1")

#        Half Assemblies down this column
#        |
#        V
tiles = [a1,  a8,  a2,  a6,  a1,  a7,  a1,  a4,  a3, # <- Half Assemblies in along this row
         a8,  a1,  a8,  a2,  a8,  a1,  a1,  a4,  a3,
         a2,  a8,  a1,  a8,  a2,  a7,  a1,  a4,  a3,
         a6,  a2,  a8,  a2,  a8,  a1,  a8,  a4,  a3,
         a1,  a8,  a2,  a8,  a2,  a5,  a4,  a3,  a3,
         a7,  a1,  a7,  a1,  a5,  a4,  a4,  a3,  0.,
         a1,  a1,  a1,  a8,  a4,  a4,  a3,  a3,  0.,
         a4,  a4,  a4,  a4,  a3,  a3,  a3,  0.,  0.,
         a3,  a3,  a3,  a3,  a3,  0.,  0.,  0.,  0.]

dx = np.array([11.5613, 23.1226, 23.1226, 23.1226, 23.1226, 23.1226, 23.1226, 23.1226, 23.1226])
nx = np.array([5,       10,      10,      10,      10,      10,      10,      10,      10])

dy = np.array([23.1226, 23.1226, 23.1226, 23.1226, 23.1226, 23.1226, 23.1226, 23.1226, 11.5613])
ny = np.array([10,      10,      10,      10,      10,      10,      10,      10,      5])

geom = DiffusionGeometry(tiles, dx, nx, dy, ny, 1., 0., 0., 1.)
solver = FDDiffusionDriver(geom)
solver.solve()
flux, x, y, _ = solver.flux()

plt.pcolormesh(y, x, flux[0,:,:])
plt.title("Group 1")
plt.show()

plt.pcolormesh(y, x, flux[1,:,:])
plt.title("Group 2")
plt.show()
