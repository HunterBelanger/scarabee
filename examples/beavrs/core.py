from scarabee import *
import matplotlib.pyplot as plt
import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(linewidth=sys.maxsize)

# Get diffusion cross sections
a1_00_  = DiffusionData.load("F16_0.npz")
a2_00_  = DiffusionData.load("F24_0.npz")
a3_00_  = DiffusionData.load("F31_0.npz")
a2_12_  = DiffusionData.load("F24_12.npz")
a2_16_  = DiffusionData.load("F24_16.npz")

a3_06U  = DiffusionData.load("F31_6U.npz")

a3_06R  = DiffusionData.load("F31_6U.npz")
a3_06R.rotate_clockwise()

a3_06D  = DiffusionData.load("F31_6U.npz")
a3_06D.rotate_clockwise()
a3_06D.rotate_clockwise()

a3_06L  = DiffusionData.load("F31_6U.npz")
a3_06L.rotate_clockwise()
a3_06L.rotate_clockwise()
a3_06L.rotate_clockwise()

a3_152  = DiffusionData.load("F31_15II.npz")

a3_151  = DiffusionData.load("F31_15II.npz")
a3_151.rotate_clockwise()

a3_154  = DiffusionData.load("F31_15II.npz")
a3_154.rotate_clockwise()
a3_154.rotate_clockwise()

a3_153  = DiffusionData.load("F31_15II.npz")
a3_153.rotate_clockwise()
a3_153.rotate_clockwise()
a3_153.rotate_clockwise()

a3_16_  = DiffusionData.load("F31_16.npz")
a3_20_  = DiffusionData.load("F31_20.npz")
rf____  = DiffusionData.load("reflector.npz")
rf____.adf = np.array([[1.25, 1.25, 1.25, 1.25],
                       [0.3,  0.3,  0.3,  0.3]])


# Define nodal geometry
#                R       P       N       M       L       K       J       H       G       F       E       D       C       B       A
tiles = [0.    , 0.    , 0.    , 0.    , rf____, rf____, rf____, rf____, rf____, rf____, rf____, rf____, rf____, 0.    , 0.    , 0.    , 0.,

         0.    , 0.    , rf____, rf____, rf____, a3_00_, a3_06D, a3_00_, a3_06D, a3_00_, a3_06D, a3_00_, rf____, rf____, rf____, 0.    , 0.,     #  1 

         0.    , rf____, rf____, a3_00_, a3_00_, a3_16_, a1_00_, a3_20_, a1_00_, a3_20_, a1_00_, a3_16_, a3_00_, a3_00_, rf____, rf____, 0.,     #  2 

         0.    , rf____, a3_00_, a3_154, a2_16_, a1_00_, a2_16_, a1_00_, a2_16_, a1_00_, a2_16_, a1_00_, a2_16_, a3_153, a3_00_, rf____, 0.,     #  3

         rf____, rf____, a3_00_, a2_16_, a2_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_16_, a2_00_, a2_16_, a3_00_, rf____, rf____, #  4

         rf____, a3_00_, a3_16_, a1_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_16_, a1_00_, a3_16_, a3_00_, rf____, #  5

         rf____, a3_06R, a1_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_16_, a1_00_, a3_06L, rf____, #  6

         rf____, a3_00_, a3_20_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a3_20_, a3_00_, rf____, #  7

         rf____, a3_06R, a1_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_16_, a1_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_16_, a1_00_, a3_06L, rf____, #  8

         rf____, a3_00_, a3_20_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a3_20_, a3_00_, rf____, #  9

         rf____, a3_06R, a1_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_16_, a1_00_, a3_06L, rf____, # 10

         rf____, a3_00_, a3_16_, a1_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_16_, a1_00_, a3_16_, a3_00_, rf____, # 11

         rf____, rf____, a3_00_, a2_16_, a2_00_, a2_16_, a1_00_, a2_12_, a1_00_, a2_12_, a1_00_, a2_16_, a2_00_, a2_16_, a3_00_, rf____, rf____, # 12

         0.    , rf____, a3_00_, a3_151, a2_16_, a1_00_, a2_16_, a1_00_, a2_16_, a1_00_, a2_16_, a1_00_, a2_16_, a3_152, a3_00_, rf____, 0.,     # 13

         0.    , rf____, rf____, a3_00_, a3_00_, a3_16_, a1_00_, a3_20_, a1_00_, a3_20_, a1_00_, a3_16_, a3_00_, a3_00_, rf____, rf____, 0.,     # 14

         0.    , 0.    , rf____, rf____, rf____, a3_00_, a3_06U, a3_00_, a3_06U, a3_00_, a3_06U, a3_00_, rf____, rf____, rf____, 0.    , 0.,     # 15

         0.    , 0.    , 0.    , 0.    , rf____, rf____, rf____, rf____, rf____, rf____, rf____, rf____, rf____, 0.    , 0.    , 0.    , 0.,
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

pin_power, x, y = solver.pin_power(z)
msk = np.where(pin_power == 0.)
nmsk = np.where(pin_power != 0.)
pin_power[msk] = np.nan

avg = np.mean(pin_power[nmsk])
pin_power /= avg
print(np.max(pin_power[nmsk]))
print(np.min(pin_power[nmsk]))

plt.pcolormesh(y, x, pin_power[:,:,0], cmap='jet')
plt.title("Pin Power")
plt.show()

