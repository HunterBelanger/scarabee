# This is problem 1 from the following paper:
# N. SUGIMURA and A. YAMAMOTO, “Evaluation of Dancoff Factors in Complicated
# Geometry using the Method of Characteristics,” J. Nucl. Sci. Technol.,
# vol. 43, no. 10, pp. 1182–1187, 2006, doi: 10.1080/18811248.2006.9711210.

from pyScarabee import *
import numpy as np

Et = np.array([1.E5])
Ea = np.array([1.E5])
Es = np.array([[0.]])
Fuel = TransportXS(Et, Ea, Es, "Fuel")

Mod_pot_xs = 0.746157
Et = np.array([Mod_pot_xs])
Ea = np.array([Mod_pot_xs])
Mod = TransportXS(Et, Ea, Es, "Mod")

# Define Cells
pitch = 1.27
iso_pitch = 100.*pitch
fuel_rad = 0.635

# Simulate isolated geometry
radii = [fuel_rad]
mats =  [Fuel, Mod]
isolated_pin = SimplePinCell(radii, mats, iso_pitch, iso_pitch)
iso_c2d = Cartesian2D([iso_pitch], [iso_pitch])
iso_c2d.set_tiles([isolated_pin])
iso_moc = MOCDriver(iso_c2d)
for i in range(iso_moc.nfsr):
  i_xs = iso_moc.xs(i)
  if i_xs.name == "Fuel":
    iso_moc.set_extern_src(i, 0, 0.)
  else:
    iso_moc.set_extern_src(i, 0, Mod_pot_xs)
iso_moc.sim_mode = SimulationMode.FixedSource
#iso_moc.x_min_bc = BoundaryCondition.Vacuum
#iso_moc.x_max_bc = BoundaryCondition.Vacuum
#iso_moc.y_min_bc = BoundaryCondition.Vacuum
#iso_moc.y_max_bc = BoundaryCondition.Vacuum
iso_moc.generate_tracks(64, 0.005, Legendre4())
iso_moc.flux_tolerance = 1.E-6
iso_moc.solve()

# Simulate real geomemtry
radii = [fuel_rad]
mats =  [Fuel, Mod]
P = SimplePinCell(radii, mats, pitch, pitch)
c2d = Cartesian2D([pitch], [pitch])
c2d.set_tiles([P])
moc = MOCDriver(c2d)
# Set the source
for i in range(moc.nfsr):
  i_xs = moc.xs(i)
  if i_xs.name == "Fuel":
    moc.set_extern_src(i, 0, 0.)
  else:
    moc.set_extern_src(i, 0, Mod_pot_xs)
moc.generate_tracks(64, 0.005, Legendre4())
moc.flux_tolerance = 1.E-6
moc.sim_mode = SimulationMode.FixedSource
moc.solve()

# Calculate Dancoff correction factor
iso_flux = iso_moc.flux(Vector(0.,0.), Direction(1.,0.), 0)
flux = moc.flux(Vector(0.,0.), Direction(1.,0.), 0)
C = (iso_flux - flux) / iso_flux
D = 1. - C

print()
print("Dancoff Correction: {}".format(C))
print("Dancoff Factor:     {}".format(D))