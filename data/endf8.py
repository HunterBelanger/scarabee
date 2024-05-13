import frendy as fdy
import h5py

#base = "/mnt/c/Users/hunte/Documents/nuclear_data/ENDF-B-VIII.0/neutrons/"
#tslbase = "/mnt/c/Users/hunte/Documents/nuclear_data/ENDF-B-VIII.0/thermal_scatt/"
base = "/mnt/c/Users/BELANH2/Documents/nuclear_data/ENDF-VIII.0/endf/neutrons/"
tslbase = "/mnt/c/Users/BELANH2/Documents/nuclear_data/ENDF-VIII.0/endf/thermal_scatt/"
lib_name = "ENDF/B-VIII.0"
#temps = [293., 500., 600., 800., 1000., 1500., 2000.]
temps = [293.6]
dil_u238 = [1.E1, 2.E1, 5.E1, 1.E2, 3.E2, 1.E3, 1.E4, 1.E5, 1.E6, 1.E8]
dil_hvy = [1.E1, 1.E2, 3.E2, 1.E3, 3.E3, 1.E4, 3.E4, 1.E5, 1.E6, 1.E8]
dil_oth = [1.E0, 1.E1, 1.E2, 1.E3, 1.E4, 1.E5, 1.E6, 1.E8]

# Set the default group strucutre
fdy.set_default_group_structure("SHEM-281")

print(fdy.get_default_group_structure())

h5 = h5py.File('endf8_shem281.h5', 'w')
h5.attrs['group-structure'] = fdy.get_default_group_structure()
h5.attrs['group-bounds'] = fdy.get_default_group_bounds()
h5.attrs['ngroups'] = len(fdy.get_default_group_bounds()) - 1
h5.attrs['library'] = lib_name

# Process TSL based evaluations
N = fdy.FrendyMG()
N.name = "H1_H2O"
N.endf_file = base + "n-001_H_001.endf"
N.tsl_file = tslbase + "tsl-HinH2O.endf"
N.tsl_type = "hh2o"
N.label = "H1 in H2O from ENDF/B-8.0"
#N.temps = [283.6, 293.6, 300., 323.6, 350., 373.6, 400., 423.6, 450., 473.6, 500., 523.6, 550., 573.6, 600., 623.6, 650., 800.]
N.temps = [293.6]
N.process(h5)

# Free Gas Nuclides
N = fdy.FrendyMG()
N.name = "H1"
N.endf_file = base + "n-001_H_001.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "O16"
N.endf_file = base + "n-008_O_016.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Zr90"
N.endf_file = base + "n-040_Zr_090.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Zr91"
N.endf_file = base + "n-040_Zr_091.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Zr92"
N.endf_file = base + "n-040_Zr_092.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Zr94"
N.endf_file = base + "n-040_Zr_094.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Zr96"
N.endf_file = base + "n-040_Zr_096.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "U235"
N.endf_file = base + "n-092_U_235.endf"
N.label = "U235 from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "U238"
N.endf_file = base + "n-092_U_238.endf"
N.label = "U238 from ENDF/B-8.0"
N.temps = temps
N.process(h5)

h5.close()