from scarabee import *
import matplotlib.pyplot as plt

# Get diffusion cross sections
a1_00 = load_diffusion_xs("F16_0.npy")
a2_00 = load_diffusion_xs("F24_0.npy")
a3_00 = load_diffusion_xs("F31_0.npy")
a2_12 = load_diffusion_xs("F24_12.npy")
a2_16 = load_diffusion_xs("F24_16.npy")
a3_06 = load_diffusion_xs("F31_6.npy")
a3_15 = load_diffusion_xs("F31_15.npy")
a3_16 = load_diffusion_xs("F31_16.npy")
a3_20 = load_diffusion_xs("F31_20.npy")
rf___  = load_diffusion_xs("reflector.npy")

# Define nodal geometry
#                 R      P      N      M      L      K      J      H      G      F      E      D      C      B      A
tiles = [0.   , 0.   , 0.   , 0.   , rf___, rf___, rf___, rf___, rf___, rf___, rf___, rf___, rf___, 0.   , 0.   , 0.   , 0.,

         0.   , 0.   , rf___, rf___, rf___, a3_00, a3_06, a3_00, a3_06, a3_00, a3_06, a3_00, rf___, rf___, rf___, 0.   , 0.,    #  1 

         0.   , rf___, rf___, a3_00, a3_00, a3_16, a1_00, a3_20, a1_00, a3_20, a1_00, a3_16, a3_00, a3_00, rf___, rf___, 0.,    #  2 

         0.   , rf___, a3_00, a3_15, a2_16, a1_00, a2_16, a1_00, a2_16, a1_00, a2_16, a1_00, a2_16, a3_15, a3_00, rf___, 0.,    #  3

         rf___, rf___, a3_00, a2_16, a2_00, a2_16, a1_00, a2_12, a1_00, a2_12, a1_00, a2_16, a2_00, a2_16, a3_00, rf___, rf___, #  4

         rf___, a3_00, a3_16, a1_00, a2_16, a1_00, a2_12, a1_00, a2_12, a1_00, a2_12, a1_00, a2_16, a1_00, a3_16, a3_00, rf___, #  5

         rf___, a3_06, a1_00, a2_16, a1_00, a2_12, a1_00, a2_12, a1_00, a2_12, a1_00, a2_12, a1_00, a2_16, a1_00, a3_06, rf___, #  6

         rf___, a3_00, a3_20, a1_00, a2_12, a1_00, a2_12, a1_00, a2_16, a1_00, a2_12, a1_00, a2_12, a1_00, a3_20, a3_00, rf___, #  7

         rf___, a3_06, a1_00, a2_16, a1_00, a2_12, a1_00, a2_16, a1_00, a2_16, a1_00, a2_12, a1_00, a2_16, a1_00, a3_06, rf___, #  8

         rf___, a3_00, a3_20, a1_00, a2_12, a1_00, a2_12, a1_00, a2_16, a1_00, a2_12, a1_00, a2_12, a1_00, a3_20, a3_00, rf___, #  9

         rf___, a3_06, a1_00, a2_16, a1_00, a2_12, a1_00, a2_12, a1_00, a2_12, a1_00, a2_12, a1_00, a2_16, a1_00, a3_06, rf___, # 10

         rf___, a3_00, a3_16, a1_00, a2_16, a1_00, a2_12, a1_00, a2_12, a1_00, a2_12, a1_00, a2_16, a1_00, a3_16, a3_00, rf___, # 11

         rf___, rf___, a3_00, a2_16, a2_00, a2_16, a1_00, a2_12, a1_00, a2_12, a1_00, a2_16, a2_00, a2_16, a3_00, rf___, rf___, # 12

         0.   , rf___, a3_00, a3_15, a2_16, a1_00, a2_16, a1_00, a2_16, a1_00, a2_16, a1_00, a2_16, a3_15, a3_00, rf___, 0.,    # 13

         0.   , rf___, rf___, a3_00, a3_00, a3_16, a1_00, a3_20, a1_00, a3_20, a1_00, a3_16, a3_00, a3_00, rf___, rf___, 0.,    # 14

         0.   , 0.   , rf___, rf___, rf___, a3_00, a3_06, a3_00, a3_06, a3_00, a3_06, a3_00, rf___, rf___, rf___, 0.   , 0.,    # 15

         0.   , 0.   , 0.   , 0.   , rf___, rf___, rf___, rf___, rf___, rf___, rf___, rf___, rf___, 0.   , 0.   , 0.   , 0.,
        ]

# Define assembly dimension and the number of nodes per assembly
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

