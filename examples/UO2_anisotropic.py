from scarabee import *
import numpy as np
import matplotlib.pyplot as plt


UO2xs = CrossSection.load("UO2_L3_CAS7")


H2Oxs = CrossSection.load("H2O_L3_CAS7")


# Define Cells
pitch = 1.26

radii = [0.4,   0.54,  0.63]
mats =  [UO2xs, UO2xs, H2Oxs, H2Oxs]
U2 = PinCell(radii, mats, pitch, pitch)

dx = [pitch]*17
UO2 = Cartesian2D(dx, dx)
UO2.set_tiles([U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,
               U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,
               U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,
               U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,
               U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,
               U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,
               U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,
               U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,
               U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,
               U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,
               U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,
               U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,
               U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,
               U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,
               U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,
               U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,
               U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2,U2])

NW = 4 # Number of empty water cells in a pin cell
dx = [pitch / NW] * NW
WC = EmptyCell(H2Oxs, pitch / NW, pitch / NW)
WT = Cartesian2D(dx, dx)
WT.set_tiles([WC]*NW*NW)

# Water assembly
dx = [pitch]*17
WAS = Cartesian2D(dx, dx)
WAS.set_tiles([WT]*17*17)

# Core assembly
dx = [pitch*17, pitch*17 , pitch*17, pitch*17]
core = Cartesian2D(dx, dx)
core.set_tiles([UO2, WAS, WAS, WAS,
                WAS, WAS, WAS, WAS,
                WAS, WAS, WAS, WAS,
                WAS, WAS, WAS, WAS])

#If you don't use pin-level cells CMFD doesn't converge
dx_cmfd = [pitch]*17*4

moc = MOCDriver(core, anisotropic = True)
moc.x_max_bc = BoundaryCondition.Vacuum
moc.y_min_bc = BoundaryCondition.Vacuum

moc.cmfd = CMFD(dx_cmfd, dx_cmfd, [[0,1],[2,4],[5,6]])
moc.cmfd.damping = 0.6
moc.cmfd.flux_limiting = False
moc.cmfd.larsen_correction = True
moc.generate_tracks(64, 0.05, YamamotoTabuchi6())

moc.keff_tolerance = 1.E-5
moc.flux_tolerance = 1.E-5
moc.solve()

'''
flux, x, y = moc.rasterize_flux(1000, 1000)

for g in range(moc.ngroups):
  plt.title("Flux in group {}".format(g+1))
  plt.pcolormesh(x, y, flux[g,:,:], cmap='jet')
  plt.savefig(f"C5G7_plot_g{g}.png")
'''