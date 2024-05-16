import frendy as fdy
import h5py

base = "/mnt/c/Users/hunte/Documents/nuclear_data/ENDF-B-VIII.0/neutrons/"
tslbase = "/mnt/c/Users/hunte/Documents/nuclear_data/ENDF-B-VIII.0/thermal_scatt/"
#base = "/mnt/c/Users/BELANH2/Documents/nuclear_data/ENDF-VIII.0/endf/neutrons/"
#tslbase = "/mnt/c/Users/BELANH2/Documents/nuclear_data/ENDF-VIII.0/endf/thermal_scatt/"
lib_name = "ENDF/B-VIII.0"
temps = [293., 500., 600., 800., 1000., 1500., 2000.]
#temps = [293.6]

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
N.temps = [283.6, 293.6, 300., 323.6, 350., 373.6, 400., 423.6, 450., 473.6, 500., 523.6, 550., 573.6, 600., 623.6, 650., 800.]
N.process(h5)

# Free Gas Nuclides
N = fdy.FrendyMG()
N.name = "H1"
N.endf_file = base + "n-001_H_001.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "He3"
N.endf_file = base + "n-002_He_003.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "He4"
N.endf_file = base + "n-002_He_004.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "B10"
N.endf_file = base + "n-005_B_010.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "B11"
N.endf_file = base + "n-005_B_011.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "C12"
N.endf_file = base + "n-006_C_012.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "C13"
N.endf_file = base + "n-006_C_013.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "N14"
N.endf_file = base + "n-007_N_014.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "N15"
N.endf_file = base + "n-007_N_015.endf"
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
N.name = "O17"
N.endf_file = base + "n-008_O_017.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "O18"
N.endf_file = base + "n-008_O_018.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Al27"
N.endf_file = base + "n-013_Al_027.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Si28"
N.endf_file = base + "n-014_Si_028.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Si29"
N.endf_file = base + "n-014_Si_029.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Si30"
N.endf_file = base + "n-014_Si_030.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ar36"
N.endf_file = base + "n-018_Ar_036.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ar38"
N.endf_file = base + "n-018_Ar_038.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ar40"
N.endf_file = base + "n-018_Ar_040.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Cr50"
N.endf_file = base + "n-024_Cr_050.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Cr52"
N.endf_file = base + "n-024_Cr_052.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Cr53"
N.endf_file = base + "n-024_Cr_053.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Cr54"
N.endf_file = base + "n-024_Cr_054.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Mn55"
N.endf_file = base + "n-025_Mn_055.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Fe54"
N.endf_file = base + "n-026_Fe_054.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Fe55"
N.endf_file = base + "n-026_Fe_055.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Fe56"
N.endf_file = base + "n-026_Fe_056.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Fe57"
N.endf_file = base + "n-026_Fe_057.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Fe58"
N.endf_file = base + "n-026_Fe_058.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Co59"
N.endf_file = base + "n-027_Co_059.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ni58"
N.endf_file = base + "n-028_Ni_058.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ni60"
N.endf_file = base + "n-028_Ni_060.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ni61"
N.endf_file = base + "n-028_Ni_061.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ni62"
N.endf_file = base + "n-028_Ni_062.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ni64"
N.endf_file = base + "n-028_Ni_064.endf"
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
N.name = "Sn112"
N.endf_file = base + "n-050_Sn_112.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Sn114"
N.endf_file = base + "n-050_Sn_114.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Sn115"
N.endf_file = base + "n-050_Sn_115.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Sn116"
N.endf_file = base + "n-050_Sn_116.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Sn117"
N.endf_file = base + "n-050_Sn_117.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Sn118"
N.endf_file = base + "n-050_Sn_118.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Sn119"
N.endf_file = base + "n-050_Sn_119.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Sn120"
N.endf_file = base + "n-050_Sn_120.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Sn122"
N.endf_file = base + "n-050_Sn_122.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Sn124"
N.endf_file = base + "n-050_Sn_124.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "I135"
N.endf_file = base + "n-053_I_135.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Xe134"
N.endf_file = base + "n-054_Xe_134.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Xe135"
N.endf_file = base + "n-054_Xe_135.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Sm149"
N.endf_file = base + "n-062_Sm_149.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "U234"
N.endf_file = base + "n-092_U_234.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "U235"
N.endf_file = base + "n-092_U_235.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "U238"
N.endf_file = base + "n-092_U_238.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

h5.close()