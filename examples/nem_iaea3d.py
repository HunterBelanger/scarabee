from scarabee import *
import numpy as np
import matplotlib.pyplot as plt

# Outer Fuel
Et = np.array([0.222222, 0.833333])
D = 1. / (3. * Et)
Ea = np.array([0.010, 0.080])
vEf = np.array([0.0, 0.135])
Ef = vEf
chi = np.array([1., 0.])
Es = np.array([[0.1922, 0.020],
               [0.0,    0.7533]])
a1 = DiffusionCrossSection(D, Ea, Es, Ef, vEf, chi, "A1")

# Inner Fuel
Et = np.array([0.222222, 0.833333])
D = 1. / (3. * Et)
Ea = np.array([0.010, 0.085])
vEf = np.array([0.0, 0.135])
Ef = vEf
Es = np.array([[0.1922, 0.020],
               [0.0,    0.7483]])
a2 = DiffusionCrossSection(D, Ea, Es, Ef, vEf, chi, "A2")

# Inner Fuel + Control Rods
Et = np.array([0.222222, 0.833333])
D = 1. / (3. * Et)
Ea = np.array([0.0100, 0.1300])
vEf = np.array([0.0, 0.135])
Ef = vEf
Es = np.array([[0.1922, 0.020],
               [0.0,    0.7033]])
a3 = DiffusionCrossSection(D, Ea, Es, Ef, vEf, chi, "A3")

# Reflector
Et = np.array([0.166667, 1.111111])
D = 1. / (3. * Et)
Ea = np.array([0.0, 0.01])
Es = np.array([[0.1267, 0.04],
               [0.0,    1.1011]])
a4 = DiffusionCrossSection(D, Ea, Es, "A4")

# Reflector + Control Rods
Et = np.array([0.166667, 1.111111])
D = 1. / (3. * Et)
Ea = np.array([0.0, 0.055])
Es = np.array([[0.0, 0.040],
               [0.0, 0.000]])
a5 = DiffusionCrossSection(D, Ea, Es, "A5")


#        Half Assemblies down this column
#        |
#        V
tiles = [# Top Reflector
         a5, a4, a4, a4, a5, a4, a4, a4, a4, # <- Half Assemblies
         a4, a4, a4, a4, a4, a4, a4, a4, a4,
         a4, a4, a5, a4, a4, a4, a4, a4, a4,
         a4, a4, a4, a4, a4, a4, a4, a4, a4,
         a5, a4, a4, a4, a5, a4, a4, a4, 0.,
         a4, a4, a4, a4, a4, a4, a4, a4, 0.,
         a4, a4, a4, a4, a4, a4, a4, 0., 0.,
         a4, a4, a4, a4, a4, a4, 0., 0., 0.,
         a4, a4, a4, a4, 0., 0., 0., 0., 0.,

         # Fuel + Control Rods
         a3, a2, a2, a2, a3, a2, a2, a1, a4, # <- Half Assemblies
         a2, a2, a2, a2, a2, a2, a2, a1, a4,
         a2, a2, a3, a2, a2, a2, a1, a1, a4,
         a2, a2, a2, a2, a2, a2, a1, a4, a4,
         a3, a2, a2, a2, a3, a1, a1, a4, 0.,
         a2, a2, a2, a2, a1, a1, a4, a4, 0.,
         a2, a2, a1, a1, a1, a4, a4, 0., 0.,
         a1, a1, a1, a4, a4, a4, 0., 0., 0.,
         a4, a4, a4, a4, 0., 0., 0., 0., 0.,

         # Fuel
         a3, a2, a2, a2, a3, a2, a2, a1, a4, # <- Half Assemblies
         a2, a2, a2, a2, a2, a2, a2, a1, a4,
         a2, a2, a2, a2, a2, a2, a1, a1, a4,
         a2, a2, a2, a2, a2, a2, a1, a4, a4,
         a3, a2, a2, a2, a3, a1, a1, a4, 0.,
         a2, a2, a2, a2, a1, a1, a4, a4, 0.,
         a2, a2, a1, a1, a1, a4, a4, 0., 0.,
         a1, a1, a1, a4, a4, a4, 0., 0., 0.,
         a4, a4, a4, a4, 0., 0., 0., 0., 0.,
        
         # Bottom Reflector
         a4, a4, a4, a4, a4, a4, a4, a4, a4, # <- Half Assemblies
         a4, a4, a4, a4, a4, a4, a4, a4, a4,
         a4, a4, a4, a4, a4, a4, a4, a4, a4,
         a4, a4, a4, a4, a4, a4, a4, a4, a4,
         a4, a4, a4, a4, a4, a4, a4, a4, 0.,
         a4, a4, a4, a4, a4, a4, a4, a4, 0.,
         a4, a4, a4, a4, a4, a4, a4, 0., 0.,
         a4, a4, a4, a4, a4, a4, 0., 0., 0.,
         a4, a4, a4, a4, 0., 0., 0., 0., 0.
        ]

dx = np.array([10., 20., 20., 20., 20., 20., 20., 20., 20.])
nx = np.array([1,   2,   2,   2,   2,   2,   2,   2,   2])

dy = np.array([20., 20., 20., 20., 20., 20., 20., 20., 10.])
ny = np.array([2,   2,   2,   2,   2,   2,   2,   2,   1])

dz = np.array([20., 13*20., 4*20., 20.])
nz = np.array([1,   13,  4,   1])

geom = DiffusionGeometry(tiles, dx, nx, dy, ny, dz, nz, 1., 0., 0., 1., 0., 0.)
solver = NEMDiffusionDriver(geom)
solver.solve()

x_max = np.sum(dx)
dx = x_max / 250.
x = np.arange(start=0., stop=x_max+dx, step=dx)

y_max = np.sum(dy)
dy = y_max / 250.
y = np.arange(start=0., stop=y_max+dy, step=dy)

z_max = np.sum(dz)
dz = z_max / 100.
z = np.arange(start=0., stop=z_max+dz, step=dz)

flux = solver.flux(x, y, z)

power = solver.power(x, y, z)

m = int(flux.shape[-1]/2)

plt.pcolormesh(y, x, flux[0,:,:,m], cmap='jet')
plt.title("Group 1")
plt.show()

plt.pcolormesh(y, x, flux[1,:,:,m], cmap='jet')
plt.title("Group 2")
plt.show()

plt.pcolormesh(y, x, power[:,:,m], cmap='jet')
plt.title("Power")
plt.show()

m = int(flux.shape[-2]/2)
plt.pcolormesh(z, x, flux[0,:,m,:], cmap='jet')
plt.title("Group 1")
plt.show()

plt.pcolormesh(z, x, flux[1,:,m,:], cmap='jet')
plt.title("Group 2")
plt.show()

plt.pcolormesh(z, x, power[:,m,:], cmap='jet')
plt.title("Power")
plt.show()
