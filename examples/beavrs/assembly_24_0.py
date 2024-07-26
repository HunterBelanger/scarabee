from scarabee import *

name = "F24_0"

set_output_file(name+"_out.txt")

ndl = NDLibrary('/mnt/c/Users/BELANH2/Documents/nuclear_data/endf8_shem281.h5')

cond_spec = [[0, 3], [4, 8], [9, 11], [12, 13], [14, 17], [18, 22], [23, 25], [26, 29],
             [30, 32], [33, 36], [37, 39], [40, 42], [43, 48], [49, 52], [53, 55],
             [56, 58], [59, 60], [61, 64], [65, 67], [68, 70], [71, 73], [74, 76],
             [77, 80], [81, 82], [83, 84], [85, 88], [89, 92], [93, 128], [129, 148],
             [149, 195], [196, 249], [250, 264], [265, 280]]

few_grp_cond_spec = [[0, 30], [31, 32]]

# Define all Materials
Fuel24Comp = MaterialComposition()
Fuel24Comp.fractions = Fraction.Atoms
Fuel24Comp.add_nuclide("O16",  4.5830e-02)
Fuel24Comp.add_nuclide("O17",  1.7411e-05)
Fuel24Comp.add_nuclide("O18",  9.1898e-05)
Fuel24Comp.add_nuclide("U234", 4.4842e-06)
Fuel24Comp.add_nuclide("U235", 5.5814e-04)
Fuel24Comp.add_nuclide("U238", 2.2407e-02)
Fuel24 = Material(Fuel24Comp, 800., ndl)

CladComp = MaterialComposition()
CladComp.fractions = Fraction.Atoms
CladComp.add_nuclide("Cr50",  3.2962e-06)
CladComp.add_nuclide("Cr52",  6.3564e-05)
CladComp.add_nuclide("Cr53",  7.2076e-06)
CladComp.add_nuclide("Cr54",  1.7941e-06)
CladComp.add_nuclide("Fe54",  8.6698e-06)
CladComp.add_nuclide("Fe56",  1.3610e-04)
CladComp.add_nuclide("Fe57",  3.1431e-06)
CladComp.add_nuclide("Fe58",  4.1829e-07)
CladComp.add_nuclide("O16",   3.0744e-04)
CladComp.add_nuclide("O17",   1.1680e-07)
CladComp.add_nuclide("O18",   6.1648e-07)
CladComp.add_nuclide("Sn112", 4.6735e-06)
CladComp.add_nuclide("Sn114", 3.1799e-06)
CladComp.add_nuclide("Sn115", 1.6381e-06)
CladComp.add_nuclide("Sn116", 7.0055e-05)
CladComp.add_nuclide("Sn117", 3.7003e-05)
CladComp.add_nuclide("Sn118", 1.1669e-04)
CladComp.add_nuclide("Sn119", 4.1387e-05)
CladComp.add_nuclide("Sn120", 1.5697e-04)
CladComp.add_nuclide("Sn122", 2.2308e-05)
CladComp.add_nuclide("Sn124", 2.7897e-05)
CladComp.add_nuclide("Zr90",  2.1828e-02)
CladComp.add_nuclide("Zr91",  4.7601e-03)
CladComp.add_nuclide("Zr92",  7.2759e-03)
CladComp.add_nuclide("Zr94",  7.3734e-03)
CladComp.add_nuclide("Zr96",  1.1879e-03)
Clad = Material(CladComp, 575., ndl)

HeComp = MaterialComposition()
HeComp.add_nuclide("He3", 4.8089e-10)
HeComp.add_nuclide("He4", 2.4044e-04)
He = Material(HeComp, 575., ndl) 

BH2O = 14
RATIO_LW = (4.9456e-02)/((4.9456e-02)+(7.7035e-06))
RATIO_HW = 1. - RATIO_LW
WaterComp = MaterialComposition()
WaterComp.add_nuclide("B10",     7.9714e-06)
WaterComp.add_nuclide("B11",     3.2247e-05)
WaterComp.add_nuclide("H1_H2O",  4.9456e-02)
WaterComp.add_nuclide("O16",     RATIO_LW*2.4673e-02)
WaterComp.add_nuclide("O17",     RATIO_LW*9.3734e-06)
WaterComp.add_nuclide("O18",     RATIO_LW*4.9474e-05)
#WaterComp.add_nuclide("H2_D2O",  7.7035e-06)
#WaterComp.add_nuclide("O16_D2O", RATIO_HW*2.4673e-02)
#WaterComp.add_nuclide("O17_D2O", RATIO_HW*9.3734e-06)
#WaterComp.add_nuclide("O18_D2O", RATIO_HW*4.9474e-05)
WaterComp.add_nuclide("H1",  7.7035e-06)
WaterComp.add_nuclide("O16", RATIO_HW*2.4673e-02)
WaterComp.add_nuclide("O17", RATIO_HW*9.3734e-06)
WaterComp.add_nuclide("O18", RATIO_HW*4.9474e-05)
Water = Material(WaterComp, 575., ndl)

# Define a guide tube
gt = GuideTube(inner_radius=0.50419, outer_radius=0.54610, clad=Clad)

# Define fuel pin
fp = FuelPin(fuel=Fuel24, fuel_radius=0.39218, gap=He, gap_width=0.40005-0.39218, clad=Clad, clad_width=0.45720-0.40005)

# Define assembly
asmbly = PWRAssembly(pitch=1.25984, moderator=Water, shape=(17, 17), ndl=ndl)
asmbly.condensation_scheme = cond_spec
asmbly.few_group_condensation_scheme = few_grp_cond_spec
asmbly.num_azimuthal_angles = 64
asmbly.track_spacing = 0.01
asmbly.pins = [fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
               fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
               fp, fp, fp, fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp, fp, fp, fp,
               fp, fp, fp, gt, fp, fp, fp, fp, fp, fp, fp, fp, fp, gt, fp, fp, fp,
               fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
               fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp,
               fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
               fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
               fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp,
               fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
               fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
               fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp,
               fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
               fp, fp, fp, gt, fp, fp, fp, fp, fp, fp, fp, fp, fp, gt, fp, fp, fp,
               fp, fp, fp, fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp, fp, fp, fp,
               fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
               fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp]
asmbly.solve()
save_diffusion_xs(name+".npy", asmbly.diffusion_xs)
