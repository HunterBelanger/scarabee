import openmc
import numpy as np

# Keff = 1.32669 +/- 0.00030

# create a model to tie together geometry, materials, settings, and tallies
model = openmc.Model()

groups = openmc.mgxs.EnergyGroups(np.logspace(-5, 7, 8))

uo2_data = openmc.XSdata('uo2', groups)
uo2_data.order = 0
uo2_data.set_total([1.77949E-01, 3.29805E-01, 4.80388E-01, 5.54367E-01, 3.11801E-01, 3.95168E-01, 5.64406E-01], temperature=294.)
uo2_data.set_absorption([8.02480E-03, 3.71740E-03, 2.67690E-02, 9.62360E-02, 3.00200E-02, 1.11260E-01, 2.82780E-01], temperature=294.)
uo2_data.set_fission([7.21206E-03, 8.19301E-04, 6.45320E-03, 1.85648E-02, 1.78084E-02, 8.30348E-02, 2.16004E-01], temperature=294.)
uo2_data.set_nu_fission([2.78145*7.21206E-03, 2.47443*8.19301E-04, 2.43383*6.45320E-03, 2.43380*1.85648E-02, 2.43380*1.78084E-02, 2.43380*8.30348E-02, 2.43380*2.16004E-01], temperature=294.)
uo2_data.set_chi([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0., 0., 0.], temperature=294.)
scat_mat = [[[1.27537E-01, 4.23780E-02, 9.43740E-06, 5.51630E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00],
             [0.00000E+00, 3.24456E-01, 1.63140E-03, 3.14270E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00],
             [0.00000E+00, 0.00000E+00, 4.50940E-01, 2.67920E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00],
             [0.00000E+00, 0.00000E+00, 0.00000E+00, 4.52565E-01, 5.56640E-03, 0.00000E+00, 0.00000E+00],
             [0.00000E+00, 0.00000E+00, 0.00000E+00, 1.25250E-04, 2.71401E-01, 1.02550E-02, 1.00210E-08],
             [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.29680E-03, 2.65802E-01, 1.68090E-02],
             [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 8.54580E-03, 2.73080E-01]]]
scat_mat = np.array(scat_mat)
scat_mat = np.rollaxis(scat_mat, 0, 3)
uo2_data.set_scatter_matrix(scat_mat, temperature=294.)

h2o_data = openmc.XSdata('h2o', groups)
h2o_data.order = 0
h2o_data.set_total([1.59206E-01, 4.12970E-01, 5.90310E-01, 5.84350E-01, 7.18000E-01, 1.25445E+00, 2.65038E+00], temperature=294.)
h2o_data.set_absorption([6.01050E-04, 1.57930E-05, 3.37160E-04, 1.94060E-03, 5.74160E-03, 1.50010E-02, 3.72390E-02], temperature=294.)
scat_mat = [[[4.44777E-02, 1.13400E-01, 7.23470E-04, 3.74990E-06, 5.31840E-08, 0.00000E+00, 0.00000E+00],
             [0.00000E+00, 2.82334E-01, 1.29940E-01, 6.23400E-04, 4.80020E-05, 7.44860E-06, 1.04550E-06],
             [0.00000E+00, 0.00000E+00, 3.45256E-01, 2.24570E-01, 1.69990E-02, 2.64430E-03, 5.03440E-04],
             [0.00000E+00, 0.00000E+00, 0.00000E+00, 9.10284E-02, 4.15510E-01, 6.37320E-02, 1.21390E-02],
             [0.00000E+00, 0.00000E+00, 0.00000E+00, 7.14370E-05, 1.39138E-01, 5.11820E-01, 6.12290E-02],
             [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.21570E-03, 6.99913E-01, 5.37320E-01],
             [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.32440E-01, 2.48070E+00]]]
scat_mat = np.array(scat_mat)
scat_mat = np.rollaxis(scat_mat, 0, 3)
h2o_data.set_scatter_matrix(scat_mat, temperature=294.)

# Initialize the library
mg_cross_sections_file = openmc.MGXSLibrary(groups)

# Add the UO2 data to it
mg_cross_sections_file.add_xsdata(uo2_data)
mg_cross_sections_file.add_xsdata(h2o_data)

# And write to disk
mg_cross_sections_file.export_to_hdf5('mgxs.h5')

# For every cross section data set in the library, assign an openmc.Macroscopic object to a material
material_dict = {}
for xs in ['uo2', 'h2o']:
  material_dict[xs] = openmc.Material(name=xs)
  material_dict[xs].set_density('macro', 1.)
  material_dict[xs].add_macroscopic(xs)

# Instantiate a Materials collection, register all Materials, and export to XML
materials = openmc.Materials(material_dict.values())

# Set the location of the cross sections file to our pre-written set
materials.cross_sections = 'mgxs.h5'

model.materials = materials

fuel_rad = openmc.ZCylinder(r=0.54)
UO2 = openmc.Cell(region=-fuel_rad, fill=material_dict['uo2'])

cell_rad = openmc.ZCylinder(r=1.25/np.sqrt(np.pi), boundary_type='white')
WTR = openmc.Cell(region=+fuel_rad & -cell_rad, fill=material_dict['h2o'])

cells = [UO2, WTR]

root_uni = openmc.Universe(cells=cells)
# Create Geometry and set root Universe
geometry = openmc.Geometry(root_uni)
model.geometry = geometry

# OpenMC simulation parameters
batches = 6000
inactive = 100
particles = 100000

# Instantiate a Settings object
settings = openmc.Settings()
settings.batches = batches
settings.inactive = inactive
settings.particles = particles

# Tell OpenMC this is a multi-group problem
settings.energy_mode = 'multi-group'

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-0.54, -1., -1., 0.54, 1., 1.]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings.source = openmc.Source(space=uniform_dist)

# Tell OpenMC we want to run in eigenvalue mode
settings.run_mode = 'eigenvalue'

model.settings = settings

model.run()
    