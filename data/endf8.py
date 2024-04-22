from frendy import FrendyMG
import h5py

base = "/mnt/c/Users/Hunter/Documents/nuclear_data/ENDF-B-VIII.0_neutrons/"
tslbase = "/mnt/c/Users/Hunter/Documents/nuclear_data/ENDF-B-VIII.0_thermal_scatt/"
temps = [293., 500., 600., 800., 1000., 1500., 2000.]
dil_u238 = [1.E1, 2.E1, 5.E1, 1.E2, 3.E2, 1.E3, 1.E4, 1.E5, 1.E6, 1.E8]
dil_hvy = [1.E1, 1.E2, 3.E2, 1.E3, 3.E3, 1.E4, 3.E4, 1.E5, 1.E6, 1.E8]
dil_oth = [1.E0, 1.E1, 1.E2, 1.E3, 1.E4, 1.E5, 1.E6, 1.E8]

h5 = h5py.File('endf8.h5', 'w')
h5.attrs['ngroups'] = 281
h5.attrs['group-structure'] = 'SHEM-281'
h5.attrs['library'] = 'ENDF/B-VIII.0'

# Process TSL based evaluations
N = FrendyMG()
N.name = "H1_H2O"
N.endf_file = base + "n-001_H_001.endf"
N.tsl_file = tslbase + "tsl-HinH2O.endf"
N.label = "H1 in H2O from ENDF/B-8.0"
N.temps = [283.6, 293.6, 300., 323.6, 350., 373.6, 400., 423.6, 450., 473.6, 500., 523.6, 550., 573.6, 600., 623.6, 650., 800.]
N.process(h5)

N = FrendyMG()
N.name = "H1_CH2"
N.endf_file = base + "n-001_H_001.endf"
N.tsl_file = tslbase + "tsl-HinCH2.endf"
N.label = "H1 in Poly (CH2) from ENDF/B-8.0"
N.temps = [77., 196., 233., 293.6, 300., 303., 313., 323., 333., 343., 350.]
N.process(h5)

N = FrendyMG()
N.name = "H2_D2O"
N.endf_file = base + "n-001_H_002.endf"
N.tsl_file = tslbase + "tsl-DinD2O.endf"
N.label = "H2 in D2O from ENDF/B-8.0"
N.temps = [283.6, 293.6, 300., 323.6, 350., 373.6, 400., 423.6, 450., 473.6, 500., 523.6, 550., 573.6, 600., 623.6, 650.]
N.process(h5)

N = FrendyMG()
N.name = "O16_D2O"
N.endf_file = base + "n-008_O_016.endf"
N.tsl_file = tslbase + "tsl-OinD2O.endf"
N.label = "O16 in D2O from ENDF/B-8.0"
N.temps = [283.6, 293.6, 300., 323.6, 350., 373.6, 400., 423.6, 450., 473.6, 500., 523.6, 550., 573.6, 600., 623.6, 650.]
N.process(h5)

N = FrendyMG()
N.name = "O17_D2O"
N.endf_file = base + "n-008_O_017.endf"
N.tsl_file = tslbase + "tsl-OinD2O.endf"
N.label = "O17 in D2O from ENDF/B-8.0"
N.temps = [283.6, 293.6, 300., 323.6, 350., 373.6, 400., 423.6, 450., 473.6, 500., 523.6, 550., 573.6, 600., 623.6, 650.]
N.process(h5)

N = FrendyMG()
N.name = "O18_D2O"
N.endf_file = base + "n-008_O_018.endf"
N.tsl_file = tslbase + "tsl-OinD2O.endf"
N.label = "O18 in D2O from ENDF/B-8.0"
N.temps = [283.6, 293.6, 300., 323.6, 350., 373.6, 400., 423.6, 450., 473.6, 500., 523.6, 550., 573.6, 600., 623.6, 650.]
N.process(h5)

N = FrendyMG()
N.name = "C12_Graphite"
N.endf_file = base + "n-006_C_012.endf"
N.tsl_file = tslbase + "tsl-crystalline-graphite.endf"
N.label = "C12 in Graphite from ENDF/B-8.0"
N.temps = [296., 400., 500., 600., 700., 800., 1000., 1200., 1600., 2000.]
N.process(h5)

N = FrendyMG()
N.name = "C13_Graphite"
N.endf_file = base + "n-006_C_013.endf"
N.tsl_file = tslbase + "tsl-crystalline-graphite.endf"
N.label = "C13 in Graphite from ENDF/B-8.0"
N.temps = [296., 400., 500., 600., 700., 800., 1000., 1200., 1600., 2000.]
N.process(h5)

N = FrendyMG()
N.name = "O16_UO2"
N.endf_file = base + "n-008_O_016.endf"
N.tsl_file = tslbase + "tsl-OinUO2.endf"
N.label = "O16 in UO2 from ENDF/B-8.0"
N.temps = [296., 400., 500., 600., 700., 800., 1000., 1200.]
N.process(h5)

N = FrendyMG()
N.name = "O17_UO2"
N.endf_file = base + "n-008_O_017.endf"
N.tsl_file = tslbase + "tsl-OinUO2.endf"
N.label = "O17 in UO2 from ENDF/B-8.0"
N.temps = [296., 400., 500., 600., 700., 800., 1000., 1200.]
N.process(h5)

N = FrendyMG()
N.name = "O18_UO2"
N.endf_file = base + "n-008_O_018.endf"
N.tsl_file = tslbase + "tsl-OinUO2.endf"
N.label = "O18 in UO2 from ENDF/B-8.0"
N.temps = [296., 400., 500., 600., 700., 800., 1000., 1200.]
N.process(h5)

N = FrendyMG
N.name = "U234_UO2"
N.endf_file = base + "n-092_U_234.endf"
N.tsl_file = tslbase + "tsl-UinUO2.endf"
N.label = "U234 in UO2 from ENDF/B-8.0"
N.temps = [296., 400., 500., 600., 700., 800., 1000., 1200.]
N.process(h5)

N = FrendyMG
N.name = "U235_UO2"
N.endf_file = base + "n-092_U_235.endf"
N.tsl_file = tslbase + "tsl-UinUO2.endf"
N.label = "U235 in UO2 from ENDF/B-8.0"
N.temps = [296., 400., 500., 600., 700., 800., 1000., 1200.]
N.process(h5)

N = FrendyMG
N.name = "U238_UO2"
N.endf_file = base + "n-092_U_238.endf"
N.tsl_file = tslbase + "tsl-UinUO2.endf"
N.label = "U238 in UO2 from ENDF/B-8.0"
N.temps = [296., 400., 500., 600., 700., 800., 1000., 1200.]
N.process(h5)

# Free Gas Nuclides
N = FrendyMG()
N.name = "H1"
N.endf_file = base + "n-001_H_001.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)
N.add_to_hdf5(h5)

N = FrendyMG()
N.name = "H2"
N.endf_file = base + "n-001_H_002.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "H3"
N.endf_file = base + "n-001_H_003.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "He3"
N.endf_file = base + "n-002_He_003.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "He4"
N.endf_file = base + "n-002_He_004.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Li6"
N.endf_file = base + "n-003_Li_006.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Li7"
N.endf_file = base + "n-003_Li_007.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Be7"
N.endf_file = base + "n-004_Be_007.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "B10"
N.endf_file = base + "n-005_B_010.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "B11"
N.endf_file = base + "n-005_B_011.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "C12"
N.endf_file = base + "n-006_C_012.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "C13"
N.endf_file = base + "n-006_C_013.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "N14"
N.endf_file = base + "n-007_N_014.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "N15"
N.endf_file = base + "n-007_N_015.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "O16"
N.endf_file = base + "n-008_O_016.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "O17"
N.endf_file = base + "n-008_O_017.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "O18"
N.endf_file = base + "n-008_O_018.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "F19"
N.endf_file = base + "n-009_F_019.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Na23"
N.endf_file = base + "n-011_Na_023.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Mg24"
N.endf_file = base + "n-012_Mg_024.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Mg25"
N.endf_file = base + "n-012_Mg_025.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Mg26"
N.endf_file = base + "n-012_Mg_026.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Al27"
N.endf_file = base + "n-013_Al_027.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Si28"
N.endf_file = base + "n-014_Si_028.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Si29"
N.endf_file = base + "n-014_Si_029.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Si30"
N.endf_file = base + "n-014_Si_030.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "P31"
N.endf_file = base + "n-015_P_031.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "S32"
N.endf_file = base + "n-016_S_032.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "S33"
N.endf_file = base + "n-016_S_033.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "S34"
N.endf_file = base + "n-016_S_034.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "S36"
N.endf_file = base + "n-016_S_036.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Cl35"
N.endf_file = base + "n-017_Cl_035.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Cl37"
N.endf_file = base + "n-017_Cl_037.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Ar36"
N.endf_file = base + "n-018_Ar_036.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Ar38"
N.endf_file = base + "n-018_Ar_038.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Ar40"
N.endf_file = base + "n-018_Ar_040.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "K39"
N.endf_file = base + "n-019_K_039.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "K40"
N.endf_file = base + "n-019_K_040.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "K41"
N.endf_file = base + "n-019_K_041.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Cr50"
N.endf_file = base + "n-024_Cr_050.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Cr52"
N.endf_file = base + "n-024_Cr_052.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Cr53"
N.endf_file = base + "n-024_Cr_053.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Cr54"
N.endf_file = base + "n-024_Cr_054.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Mn55"
N.endf_file = base + "n-025_Mn_055.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Fe54"
N.endf_file = base + "n-026_Fe_054.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Fe56"
N.endf_file = base + "n-026_Fe_056.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = dil_oth
N.process(h5)

N = FrendyMG()
N.name = "Fe57"
N.endf_file = base + "n-026_Fe_057.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Fe58"
N.endf_file = base + "n-026_Fe_058.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Ni58"
N.endf_file = base + "n-028_Ni_058.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Ni60"
N.endf_file = base + "n-028_Ni_060.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Ni61"
N.endf_file = base + "n-028_Ni_061.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Ni62"
N.endf_file = base + "n-028_Ni_062.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Ni64"
N.endf_file = base + "n-028_Ni_064.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Zr90"
N.endf_file = base + "n-040_Zr_090.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = dil_oth
N.process(h5)

N = FrendyMG()
N.name = "Zr91"
N.endf_file = base + "n-040_Zr_091.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = dil_oth
N.process(h5)

N = FrendyMG()
N.name = "Zr92"
N.endf_file = base + "n-040_Zr_092.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = dil_oth
N.process(h5)

N = FrendyMG()
N.name = "Zr94"
N.endf_file = base + "n-040_Zr_094.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Zr96"
N.endf_file = base + "n-040_Zr_096.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Ag107"
N.endf_file = base + "n-047_Ag_107.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = dil_oth
N.process(h5)

N = FrendyMG()
N.name = "Ag109"
N.endf_file = base + "n-047_Ag_109.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = dil_oth
N.process(h5)

N = FrendyMG()
N.name = "Cd106"
N.endf_file = base + "n-048_Cd_106.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Cd108"
N.endf_file = base + "n-048_Cd_108.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Cd110"
N.endf_file = base + "n-048_Cd_110.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = dil_oth
N.process(h5)

N = FrendyMG()
N.name = "Cd111"
N.endf_file = base + "n-048_Cd_111.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Cd112"
N.endf_file = base + "n-048_Cd_112.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Cd113"
N.endf_file = base + "n-048_Cd_113.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Cd114"
N.endf_file = base + "n-048_Cd_114.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Cd116"
N.endf_file = base + "n-048_Cd_116.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "In113"
N.endf_file = base + "n-049_In_113.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "In115"
N.endf_file = base + "n-049_In_115.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.dilutions = dil_oth
N.process(h5)

N = FrendyMG()
N.name = "Sn112"
N.endf_file = base + "n-050_Sn_112.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Sn114"
N.endf_file = base + "n-050_Sn_114.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Sn115"
N.endf_file = base + "n-050_Sn_115.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Sn116"
N.endf_file = base + "n-050_Sn_116.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Sn117"
N.endf_file = base + "n-050_Sn_117.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Sn118"
N.endf_file = base + "n-050_Sn_118.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Sn119"
N.endf_file = base + "n-050_Sn_119.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Sn120"
N.endf_file = base + "n-050_Sn_120.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Sn122"
N.endf_file = base + "n-050_Sn_122.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG()
N.name = "Sn124"
N.endf_file = base + "n-050_Sn_124.endf"
N.label = N.name + " from ENDF/B-8.0"
N.temps = temps
N.process(h5)

N = FrendyMG
N.name = "U234"
N.endf_file = base + "n-092_U_234.endf"
N.label = "U234 from ENDF/B-8.0"
N.temps = temps
N.dilutions = dil_hvy
N.process(h5)

N = FrendyMG
N.name = "U235"
N.endf_file = base + "n-092_U_235.endf"
N.label = "U235 from ENDF/B-8.0"
N.temps = temps
N.dilutions = dil_hvy
N.process(h5)

N = FrendyMG
N.name = "U238"
N.endf_file = base + "n-092_U_238.endf"
N.label = "U238 from ENDF/B-8.0"
N.temps = temps
N.dilutions = dil_u238
N.process(h5)

N = FrendyMG
N.name = "Pu238"
N.endf_file = base + "n-094_Pu_238.endf"
N.label = "Pu238 from ENDF/B-8.0"
N.temps = temps
N.dilutions = dil_hvy
N.process(h5)

N = FrendyMG
N.name = "Pu239"
N.endf_file = base + "n-094_Pu_239.endf"
N.label = "Pu239 from ENDF/B-8.0"
N.temps = temps
N.dilutions = dil_hvy
N.process(h5)

N = FrendyMG
N.name = "Pu240"
N.endf_file = base + "n-094_Pu_240.endf"
N.label = "Pu240 from ENDF/B-8.0"
N.temps = temps
N.dilutions = dil_hvy
N.process(h5)

N = FrendyMG
N.name = "Pu241"
N.endf_file = base + "n-094_Pu_241.endf"
N.label = "Pu241 from ENDF/B-8.0"
N.temps = temps
N.dilutions = dil_hvy
N.process(h5)

N = FrendyMG
N.name = "Pu242"
N.endf_file = base + "n-094_Pu_242.endf"
N.label = "Pu242 from ENDF/B-8.0"
N.temps = temps
N.dilutions = dil_hvy
N.process(h5)

h5.close()