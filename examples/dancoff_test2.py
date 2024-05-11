# This is problem 2 from the following paper:
# N. SUGIMURA and A. YAMAMOTO, “Evaluation of Dancoff Factors in Complicated
# Geometry using the Method of Characteristics,” J. Nucl. Sci. Technol.,
# vol. 43, no. 10, pp. 1182–1187, 2006, doi: 10.1080/18811248.2006.9711210.

from pyScarabee import *
import numpy as np
np.set_printoptions(linewidth=200)

Et = np.array([1.E5])
Ea = np.array([1.E5])
Es = np.array([[0.]])
Fuel = TransportXS(Et, Ea, Es, "Fuel")

Mod_pot_xs = 1.05
Et = np.array([Mod_pot_xs])
Ea = np.array([Mod_pot_xs])
Mod = TransportXS(Et, Ea, Es, "Mod")

Clad_pot_xs = 0.30
Et = np.array([Clad_pot_xs])
Ea = np.array([Clad_pot_xs])
Clad = TransportXS(Et, Ea, Es, "Clad")

Thmb_pot_xs = 0.30
Et = np.array([Thmb_pot_xs])
Ea = np.array([Thmb_pot_xs])
Thmb = TransportXS(Et, Ea, Es, "Thmb")

# Define Cells
pitch = 1.26
iso_pitch = 7.*pitch
fuel_rad = 0.41
clad_rad = 0.48
mod_rad = 0.57
thmb_rad = 0.61

# Simulate isolated geometry
radii = [fuel_rad, clad_rad]
mats =  [Fuel, Clad, Mod]
isolated_pin = SimplePinCell(radii, mats, iso_pitch, iso_pitch)
iso_c2d = Cartesian2D([iso_pitch], [iso_pitch])
iso_c2d.set_tiles([isolated_pin])
iso_moc = MOCDriver(iso_c2d)
for i in range(iso_moc.nfsr):
  i_xs = iso_moc.xs(i)
  if i_xs.name == "Fuel":
    iso_moc.set_extern_src(i, 0, 0.)
  elif i_xs.name == "Clad":
    iso_moc.set_extern_src(i, 0, Clad_pot_xs)
  else:
    iso_moc.set_extern_src(i, 0, Mod_pot_xs)
iso_moc.sim_mode = SimulationMode.FixedSource
iso_moc.generate_tracks(64, 0.01, Legendre12())
iso_moc.flux_tolerance = 1.E-6
iso_moc.solve()

# Simulate real geomemtry
radii = [fuel_rad, clad_rad]
mats =  [Fuel, Clad, Mod]
F = SimplePinCell(radii, mats, pitch, pitch)
radii = [mod_rad, thmb_rad]
mats = [Mod, Thmb, Mod]
G = SimplePinCell(radii, mats, pitch, pitch)
c2d = Cartesian2D(17*[pitch], 17*[pitch])
c2d.set_tiles([F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,
               F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,
               F,F,F,F,F,G,F,F,G,F,F,G,F,F,F,F,F,
               F,F,F,G,F,F,F,F,F,F,F,F,F,G,F,F,F,
               F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,
               F,F,G,F,F,G,F,F,G,F,F,G,F,F,G,F,F,
               F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,
               F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,
               F,F,G,F,F,G,F,F,G,F,F,G,F,F,G,F,F,
               F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,
               F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,
               F,F,G,F,F,G,F,F,G,F,F,G,F,F,G,F,F,
               F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,
               F,F,F,G,F,F,F,F,F,F,F,F,F,G,F,F,F,
               F,F,F,F,F,G,F,F,G,F,F,G,F,F,F,F,F,
               F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,
               F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F])
moc = MOCDriver(c2d)
# Set the source
for i in range(moc.nfsr):
  i_xs = moc.xs(i)
  if i_xs.name == "Fuel":
    moc.set_extern_src(i, 0, 0.)
  elif i_xs.name == "Clad":
    moc.set_extern_src(i, 0, Clad_pot_xs)
  elif i_xs.name == "Thmb":
    moc.set_extern_src(i, 0, Thmb_pot_xs)
  else:
    moc.set_extern_src(i, 0, Mod_pot_xs)
moc.generate_tracks(64, 0.01, Legendre12())
moc.flux_tolerance = 1.E-6
moc.sim_mode = SimulationMode.FixedSource
moc.solve()

# Calculate Dancoff correction factor
iso_flux = iso_moc.flux(Vector(0.,0.), Direction(1.,0.), 0)

D = np.zeros((17,17))
for j in range(17):
  y = moc.y_max - (j+0.5)*pitch
  for i in range(17):
    x = moc.x_min + (i+0.5)*pitch
    xs = moc.xs(Vector(x,y), Direction(1.,0.))
    if xs.name == "Fuel":
      flux = moc.flux(Vector(x,y), Direction(1.,0.), 0)
      C = (iso_flux - flux) / iso_flux
      D[j,i] = 1. - C

print()
for i in range(9):
  print(D[i+8, 8:i+9])