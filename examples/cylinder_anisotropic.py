# Analytical Benchmark Test Set for Criticality Code and Verification
# Problem No.: 36 [Ua-1-1-CY]
# Exact keff = 1.

from scarabee import *
import numpy as np
import matplotlib.pyplot as plt


Et = np.array([0.32640])
Ea = np.array([0.32640-0.248064])
Es = np.array([[[0.248064]],
               [[0.042432]]])
Ef = np.array([0.065280])
vEf = np.array([2.70 * 0.065280])
chi = np.array([1.])
dtr = np.array([0.])
fuel = CrossSection(Et, dtr, Ea, Es, Ef, vEf, chi)

# fictitious material
Ea = np.array([1E-4])
Et = Ea
Es = np.array([[[0.]],
               [[0.]]])
Ef = np.array([0.])
vEf = np.array([0.])
chi = np.array([0.])
dtr = np.array([0.])
fictitious_mat = CrossSection(Et, dtr, Ea, Es, Ef, vEf, chi)

# Creating Cell

# max length in x and y
radius = 5.514296811
length_x = 2 * radius
length_y = 2 * radius
# number of elements in x and y
Nx = 100
Ny = 100
# delta x and y for element
dx = length_x / Nx
dy = length_y / Ny
# low points
x0 = -length_x * 0.5
y0 = -length_y * 0.5

# fuel-cell
F = EmptyCell(fuel, dx, dy)
# ficititious material cell
N = EmptyCell(fictitious_mat, dx, dy)

# creating the cylinder with small cubes
c2d_dx = [dx] * Nx
c2d_dy = [dy] * Ny

root_univ = Cartesian2D(c2d_dx, c2d_dy)

lattice = []
for j in range(0, Ny):
    # mid-y in the element
    y = y0 + (j + 0.5) * dy

    for i in range(0, Nx):
        # mid-x in the element
        x = x0 + (i + 0.5) * dx
        # check, if the mid point is inside the cylinder
        if (np.sqrt(x*x + y*y) < radius):
            lattice.append(F)
        else:
            lattice.append(N)
            
root_univ.set_tiles(lattice)

moc = MOCDriver(root_univ, anisotropic = True)
moc.x_min_bc = BoundaryCondition.Vacuum
moc.x_max_bc = BoundaryCondition.Vacuum
moc.y_min_bc = BoundaryCondition.Vacuum
moc.y_max_bc = BoundaryCondition.Vacuum
moc.generate_tracks(128, 0.01, YamamotoTabuchi6())
moc.flux_tolerance = 1.E-5

# moc.plot()
moc.solve()

#flux, x, y = moc.rasterize_flux(100, 100)


