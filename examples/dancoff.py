from pyScarabee import *
import numpy as np
import matplotlib.pyplot as plt

#Et = np.array([0.4336])
#Ea = np.array([0.4336])
Et = np.array([1.E5])
Ea = np.array([1.E5])
Es = np.array([[0.]])
UO2 = TransportXS(Et, Ea, Es)

Et = np.array([0.2943])
Ea = np.array([0.2943])
Es = np.array([[0.]])
CLD = TransportXS(Et, Ea, Es)

Et = np.array([1.035])
Ea = np.array([1.035])
Es = np.array([[0.]])
H2O = TransportXS(Et, Ea, Es)

# Define Cells
pitch = 1.26
fuel_rad = 0.4025
clad_rad = 0.475

radii = [fuel_rad, clad_rad]
mats =  [UO2, CLD, H2O]
P = SimplePinCell(radii, mats, pitch, pitch)

W = EmptyCell(H2O, pitch, pitch)

dx = [pitch]*4
c2d = Cartesian2D(dx, dx)
c2d.set_tiles([P, P, P, P,
               P, P, P, P,
               W, P, P, P,
               W, W, P, P])

moc = MOCDriver(c2d)
moc.generate_tracks(64, 0.05, YamamotoTabuchi6())

# Set the source
u = Direction(1.,0.)
for i in range(4):
  for j in range(4):
    r = Vector((i+0.5)*pitch + moc.x_min, (j+0.5)*pitch+ moc.y_min)
    if i == 0 and j in [0, 1]:
      moc.set_extern_src(r, u, 0, H2O.Et(0))
    elif i == 1 and j == 1:
      moc.set_extern_src(r, u, 0, H2O.Et(0))
    else:
      moc.set_extern_src(r, u, 0, UO2.Et(0))
      moc.set_extern_src(r+Vector(fuel_rad,0.), u, 0, CLD.Et(0))
      moc.set_extern_src(r+Vector(clad_rad,0.), u, 0, H2O.Et(0))

moc.flux_tolerance = 1.E-5
moc.sim_mode = SimulationMode.FixedSource
moc.solve()

flux, x, y = moc.rasterize_flux(1000, 1000)

for g in range(moc.ngroups):
  plt.title("Flux in group {}".format(g+1))
  plt.pcolormesh(x, y, flux[g,:,:], cmap='jet')
  plt.show()
