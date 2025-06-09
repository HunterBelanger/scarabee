from scarabee import *
import numpy as np
import matplotlib.pyplot as plt

Et = np.array([1.59206E-01, 4.12970E-01, 5.90310E-01, 5.84350E-01, 7.18000E-01, 1.25445E+00, 2.65038E+00])
Ea = np.array([6.01050E-04, 1.57930E-05, 3.37160E-04, 1.94060E-03, 5.74160E-03, 1.50010E-02, 3.72390E-02])
Es = np.array([[4.44777E-02, 1.13400E-01, 7.23470E-04, 3.74990E-06, 5.31840E-08, 0.00000E+00, 0.00000E+00],
               [0.00000E+00, 2.82334E-01, 1.29940E-01, 6.23400E-04, 4.80020E-05, 7.44860E-06, 1.04550E-06],
               [0.00000E+00, 0.00000E+00, 3.45256E-01, 2.24570E-01, 1.69990E-02, 2.64430E-03, 5.03440E-04],
               [0.00000E+00, 0.00000E+00, 0.00000E+00, 9.10284E-02, 4.15510E-01, 6.37320E-02, 1.21390E-02],
               [0.00000E+00, 0.00000E+00, 0.00000E+00, 7.14370E-05, 1.39138E-01, 5.11820E-01, 6.12290E-02],
               [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.21570E-03, 6.99913E-01, 5.37320E-01],
               [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.32440E-01, 2.48070E+00]])
H2Oxs = CrossSection(Et, Ea, Es)

# Define Cells
pitch = 1.26

NW = 6 # Number of empty water cells in a pin cell
dx = [pitch / NW] * NW
WC = EmptyCell(H2Oxs, pitch / NW, pitch / NW)
WT = Cartesian2D(dx, dx)
WT.set_tiles([WC]*NW*NW)

# Water assembly
dx = [pitch]*17
WAS = Cartesian2D(dx, dx)
WAS.set_tiles([WT]*17*17)

# YOU MUST CREATE A COARSE MESH THAT FITS PROPERLY WITH FSR BOUNDARIES
# OTHERWISE BAD THINGS HAPPEN
dxc = [pitch/2]*17*2

moc = MOCDriver(WAS)
moc.x_min_bc = BoundaryCondition.Vacuum
moc.x_max_bc = BoundaryCondition.Vacuum
moc.y_min_bc = BoundaryCondition.Vacuum
moc.y_max_bc = BoundaryCondition.Vacuum
moc.cmfd = CMFD(dxc,dxc,[[0,0], [1,1], [2,2], [3,3], [4,4], [5,5], [6,6]])
# Works with several damping factors, but may become unstable
# with larger values
moc.cmfd.set_damping(0.2)
moc.cmfd.flux_limiting = True
moc.generate_tracks(32, 0.05, YamamotoTabuchi6())
moc.set_extern_src(Vector(0.,0.), Direction(0.,1.), 0, 1.)
moc.flux_tolerance = 1.E-5
moc.sim_mode = SimulationMode.FixedSource
moc.solve()

flux, x, y = moc.rasterize_flux(1000, 1000)

for g in range(moc.ngroups):
  plt.title("Flux in group {}".format(g+1))
  plt.pcolormesh(x, y, flux[g,:,:], cmap='jet')
  plt.show()

