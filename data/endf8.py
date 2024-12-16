import frendy as fdy
import h5py
import numpy as np

#base = "/mnt/c/Users/hunte/Documents/nuclear_data/ENDF-B-VIII.0/neutrons/"
#tslbase = "/mnt/c/Users/hunte/Documents/nuclear_data/ENDF-B-VIII.0/thermal_scatt/"
base = "/mnt/c/Users/BELANH2/Documents/nuclear_data/ENDF-VIII.0/endf/neutrons/"
tslbase = "/mnt/c/Users/BELANH2/Documents/nuclear_data/ENDF-VIII.0/endf/thermal_scatt/"
lib_name = "ENDF/B-VIII.0"
temps = [293., 500., 600., 800., 1000., 1500., 2000.]
#temps = [293.6]

# Condensation Schemes for SHEM-281 group structure
# 25 groups (similar to CASMO-25 group structure)
cond_spec = np.array([[0, 8], [9, 11], [12, 13], [14, 17], [18, 22], [23, 25],
             [26, 36], [37, 55], [56, 58], [59, 79], [80, 120],
             [121, 142], [143, 195], [196, 216], [217, 229], [230, 233],
             [234, 236], [237, 240], [241, 246], [247, 253], [254, 256],
             [257, 262], [263, 268], [269, 272], [273, 280]])

few_grp_cond_spec = np.array([[0, 18], [19, 24]])

ref_few_grp_cond_spec = np.array([[0, 246], [247, 280]])

# Set the default group strucutre
fdy.set_default_group_structure("SHEM-281")

fdy.set_default_max_legendre_moments(3)

print(fdy.get_default_group_structure())

h5 = h5py.File('endf8_shem281.h5', 'w')
#h5 = h5py.File('endf8_xmas172.h5', 'w')
h5.attrs['group-structure'] = fdy.get_default_group_structure()
h5.attrs['group-bounds'] = fdy.get_default_group_bounds()
h5.attrs['ngroups'] = len(fdy.get_default_group_bounds()) - 1
h5.attrs['library'] = lib_name

if fdy.get_default_group_structure() == "SHEM-281":
    h5.attrs['macro-group-condensation-scheme'] = cond_spec
    h5.attrs['few-group-condensation-scheme'] = few_grp_cond_spec
    h5.attrs['reflector-few-group-condensation-scheme'] = ref_few_grp_cond_spec

# We process U235 first so that we can steal its fission spectrum for
# performing transport correciton calculations
N = fdy.FrendyMG()
N.name = "U235"
N.endf_file = base + "n-092_U_235.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

chi = N.chi[:]

# Process TSL based evaluations
N = fdy.FrendyMG()
N.name = "H1_H2O"
N.endf_file = base + "n-001_H_001.endf"
N.tsl_file = tslbase + "tsl-HinH2O.endf"
N.tsl_type = "hh2o"
N.label = "H1 in H2O from ENDF/B-8.0"
N.temps = [283.6, 293.6, 300., 323.6, 350., 373.6, 400., 423.6, 450., 473.6, 500., 523.6, 550., 573.6, 600., 623.6, 650., 800.]
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "H1_H2O_TC"
N.endf_file = base + "n-001_H_001.endf"
N.tsl_file = tslbase + "tsl-HinH2O.endf"
N.tsl_type = "hh2o"
N.label = "H1 in H2O from ENDF/B-8.0"
N.temps = [283.6, 293.6, 300., 323.6, 350., 373.6, 400., 423.6, 450., 473.6, 500., 523.6, 550., 573.6, 600., 623.6, 650., 800.]
N.dilutions = [1.E10]
N.process(h5, chi)

N = fdy.FrendyMG()
N.name = "H2_D2O"
N.endf_file = base + "n-001_H_002.endf"
N.tsl_file = tslbase + "tsl-DinD2O.endf"
N.tsl_type = "dd2o"
N.label = "H2 in D2O from ENDF/B-8.0"
N.temps = [283.6, 293.6, 300., 323.6, 350., 373.6, 400., 423.6, 450., 473.6, 500., 523.6, 550., 573.6, 600., 623.6, 650.]
N.dilutions = [1.E10]
N.process(h5)

# Free Gas Nuclides
N = fdy.FrendyMG()
N.name = "H1"
N.endf_file = base + "n-001_H_001.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "H2"
N.endf_file = base + "n-001_H_002.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "He3"
N.endf_file = base + "n-002_He_003.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "He4"
N.endf_file = base + "n-002_He_004.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Li6"
N.endf_file = base + "n-003_Li_006.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Li7"
N.endf_file = base + "n-003_Li_007.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Be9"
N.endf_file = base + "n-003_Be_009.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "B10"
N.endf_file = base + "n-005_B_010.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "B11"
N.endf_file = base + "n-005_B_011.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "C12"
N.endf_file = base + "n-006_C_012.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "C13"
N.endf_file = base + "n-006_C_013.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "N14"
N.endf_file = base + "n-007_N_014.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "N15"
N.endf_file = base + "n-007_N_015.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "O16"
N.endf_file = base + "n-008_O_016.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "O17"
N.endf_file = base + "n-008_O_017.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "O18"
N.endf_file = base + "n-008_O_018.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Na23"
N.endf_file = base + "n-011_Na_023.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Mg24"
N.endf_file = base + "n-012_Mg_024.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Mg25"
N.endf_file = base + "n-012_Mg_025.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Mg26"
N.endf_file = base + "n-012_Mg_026.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Al27"
N.endf_file = base + "n-013_Al_027.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Si28"
N.endf_file = base + "n-014_Si_028.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Si29"
N.endf_file = base + "n-014_Si_029.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Si30"
N.endf_file = base + "n-014_Si_030.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ar36"
N.endf_file = base + "n-018_Ar_036.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ar38"
N.endf_file = base + "n-018_Ar_038.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ar40"
N.endf_file = base + "n-018_Ar_040.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Cr50"
N.endf_file = base + "n-024_Cr_050.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
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
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Mn55"
N.endf_file = base + "n-025_Mn_055.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
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
N.dilutions = [1.E10]
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
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Fe58"
N.endf_file = base + "n-026_Fe_058.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Co59"
N.endf_file = base + "n-027_Co_059.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ni58"
N.endf_file = base + "n-028_Ni_058.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ni60"
N.endf_file = base + "n-028_Ni_060.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ni61"
N.endf_file = base + "n-028_Ni_061.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ni62"
N.endf_file = base + "n-028_Ni_062.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ni64"
N.endf_file = base + "n-028_Ni_064.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Cu63"
N.endf_file = base + "n-029_Cu_063.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
N.process(h5)

N = fdy.FrendyMG()
N.name = "Cu65"
N.endf_file = base + "n-029_Cu_065.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
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
N.name = "Mo92"
N.endf_file = base + "n-042_Mo_092.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Mo94"
N.endf_file = base + "n-042_Mo_094.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Mo95"
N.endf_file = base + "n-042_Mo_095.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Mo96"
N.endf_file = base + "n-042_Mo_096.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Mo97"
N.endf_file = base + "n-042_Mo_097.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Mo98"
N.endf_file = base + "n-042_Mo_098.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Mo100"
N.endf_file = base + "n-042_Mo_100.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ag107"
N.endf_file = base + "n-047_Ag_107.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ag109"
N.endf_file = base + "n-047_Ag_109.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Ag111"
N.endf_file = base + "n-047_Ag_111.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Cd106"
N.endf_file = base + "n-048_Cd_106.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Cd108"
N.endf_file = base + "n-048_Cd_108.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Cd110"
N.endf_file = base + "n-048_Cd_110.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Cd111"
N.endf_file = base + "n-048_Cd_111.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Cd112"
N.endf_file = base + "n-048_Cd_112.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Cd113"
N.endf_file = base + "n-048_Cd_113.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Cd114"
N.endf_file = base + "n-048_Cd_114.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Cd116"
N.endf_file = base + "n-048_Cd_116.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "In113"
N.endf_file = base + "n-049_In_113.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "In115"
N.endf_file = base + "n-049_In_115.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Sn112"
N.endf_file = base + "n-050_Sn_112.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = [1.E10]
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
N.name = "Gd152"
N.endf_file = base + "n-064_Gd_152.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Gd154"
N.endf_file = base + "n-064_Gd_154.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Gd155"
N.endf_file = base + "n-064_Gd_155.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Gd156"
N.endf_file = base + "n-064_Gd_156.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Gd157"
N.endf_file = base + "n-064_Gd_157.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Gd158"
N.endf_file = base + "n-064_Gd_158.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Gd160"
N.endf_file = base + "n-064_Gd_160.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Th232"
N.endf_file = base + "n-090_Th_232.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "U233"
N.endf_file = base + "n-092_U_233.endf"
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
N.name = "U236"
N.endf_file = base + "n-092_U_236.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "U238"
N.endf_file = base + "n-092_U_238.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Pu236"
N.endf_file = base + "n-094_Pu_236.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Pu237"
N.endf_file = base + "n-094_Pu_237.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Pu238"
N.endf_file = base + "n-094_Pu_238.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Pu239"
N.endf_file = base + "n-094_Pu_239.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Pu240"
N.endf_file = base + "n-094_Pu_240.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Pu241"
N.endf_file = base + "n-094_Pu_241.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Pu242"
N.endf_file = base + "n-094_Pu_242.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Pu243"
N.endf_file = base + "n-094_Pu_243.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = fdy.FrendyMG()
N.name = "Pu244"
N.endf_file = base + "n-094_Pu_244.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)


h5.close()
