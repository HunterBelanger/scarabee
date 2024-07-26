from scarabee import *
import matplotlib.pyplot as plt

# Get diffusion cross sections
a1 = load_diffusion_xs("F16_0.npy")
a2 = load_diffusion_xs("F24_0.npy")
a3 = load_diffusion_xs("F31_0.npy")
rf = load_diffusion_xs("reflector.npy")

# Define nodal geometry
#            R   P   N   M   L   K   J   H   G   F   E   D   C   B   A
tiles = [0., 0., 0., 0., rf, rf, rf, rf, rf, rf, rf, rf, rf, 0., 0., 0., 0.,
         0., 0., rf, rf, rf, a3, a3, a3, a3, a3, a3, a3, rf, rf, rf, 0., 0., #  1 
         0., rf, rf, a3, a3, a3, a1, a3, a1, a3, a1, a3, a3, a3, rf, rf, 0., #  2 
         0., rf, a3, a3, a2, a1, a2, a1, a2, a1, a2, a1, a2, a3, a3, rf, 0., #  3
         rf, rf, a3, a2, a2, a2, a1, a2, a1, a2, a1, a2, a2, a2, a3, rf, rf, #  4
         rf, a3, a3, a1, a2, a1, a2, a1, a2, a1, a2, a1, a2, a1, a3, a3, rf, #  5
         rf, a3, a1, a2, a1, a2, a1, a2, a1, a2, a1, a2, a1, a2, a1, a3, rf, #  6
         rf, a3, a3, a1, a2, a1, a2, a1, a2, a1, a2, a1, a2, a1, a3, a3, rf, #  7
         rf, a3, a1, a2, a1, a2, a1, a2, a1, a2, a1, a2, a1, a2, a1, a3, rf, #  8
         rf, a3, a3, a1, a2, a1, a2, a1, a2, a1, a2, a1, a2, a1, a3, a3, rf, #  9
         rf, a3, a1, a2, a1, a2, a1, a2, a1, a2, a1, a2, a1, a2, a1, a3, rf, # 10
         rf, a3, a3, a1, a2, a1, a2, a1, a2, a1, a2, a1, a2, a1, a3, a3, rf, # 11
         rf, rf, a3, a2, a2, a2, a1, a2, a1, a2, a1, a2, a2, a2, a3, rf, rf, # 12
         0., rf, a3, a3, a2, a1, a2, a1, a2, a1, a2, a1, a2, a3, a3, rf, 0., # 13
         0., rf, rf, a3, a3, a3, a1, a3, a1, a3, a1, a3, a3, a3, rf, rf, 0., # 14
         0., 0., rf, rf, rf, a3, a3, a3, a3, a3, a3, a3, rf, rf, rf, 0., 0., # 15
         0., 0., 0., 0., rf, rf, rf, rf, rf, rf, rf, rf, rf, 0., 0., 0., 0.,
        ]

dx = np.array(17*[21.50364])
nx = 2*np.array(17*[1])

dy = np.array(17*[21.50364])
ny = 2*np.array(17*[1])

dz = np.array([20.])
nz = np.array([1])

geom = DiffusionGeometry(tiles, dx, nx, dy, ny, dz, nz, 0., 0., 0., 0., 1., 1.)
solver = NEMDiffusionDriver(geom)
solver.solve()

x_max = np.sum(dx)
dx = x_max / 500.
x = np.arange(start=0., stop=x_max+dx, step=dx)

y_max = np.sum(dy)
dy = y_max / 500.
y = np.arange(start=0., stop=y_max+dy, step=dy)

z = np.array([0.])

flux = solver.flux(x, y, z)[:,:,:,0]

plt.pcolormesh(y, x, flux[0,:,:], cmap='jet')
plt.title("Group 1")
plt.show()

plt.pcolormesh(y, x, flux[1,:,:], cmap='jet')
plt.title("Group 2")
plt.show()

power = solver.power(x, y, z)[:,:,0]
plt.pcolormesh(y, x, power, cmap='jet')
plt.title("Power")
plt.show()
