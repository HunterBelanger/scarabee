from scarabee import *
import numpy as np
import matplotlib.pyplot as plt

name = "F16_0"

set_output_file(name+"_out.txt")

ndl = NDLibrary()

# 25 groups (similar to CASMO-25 group structure)
cond_spec = [[0, 8], [9, 11], [12, 13], [14, 17], [18, 22], [23, 25],
             [26, 36], [37, 55], [56, 58], [59, 79], [80, 120],
             [121, 142], [143, 195], [196, 216], [217, 229], [230, 233],
             [234, 236], [237, 240], [241, 246], [247, 253], [254, 256],
             [257, 262], [263, 268], [269, 272], [273, 280]]
few_grp_cond_spec = [[0, 18], [19, 24]]

# Define all Materials
Fuel16Comp = MaterialComposition()
Fuel16Comp.fractions = Fraction.Atoms
Fuel16Comp.add_nuclide("O16",  4.5897e-02)
Fuel16Comp.add_nuclide("O17",  1.7436e-05)
Fuel16Comp.add_nuclide("O18",  9.2032e-05)
Fuel16Comp.add_nuclide("U234", 3.0131e-06)
Fuel16Comp.add_nuclide("U235", 3.7503e-04)
Fuel16Comp.add_nuclide("U238", 2.2625e-02)
Fuel16 = Material(Fuel16Comp, 575., ndl)

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

WaterComp = MaterialComposition()
WaterComp.add_nuclide("B10",     7.9714e-06)
WaterComp.add_nuclide("B11",     3.2247e-05)
WaterComp.add_nuclide("H1_H2O_TC",  4.9456e-02 + 7.7035e-06)
WaterComp.add_nuclide("O16",     2.4673e-02)
WaterComp.add_nuclide("O17",     9.3734e-06)
WaterComp.add_nuclide("O18",     4.9474e-05)
Water = Material(WaterComp, 575., ndl)

# Define a guide tube
gt = GuideTube(inner_radius=0.56134, outer_radius=0.60198, clad=Clad)

# Define fuel pin
fp = FuelPin(fuel=Fuel16, fuel_radius=0.39218, gap=He, gap_radius=0.40005,
             clad=Clad, clad_radius=0.45720)

# Define assembly
asmbly = PWRAssembly(pitch=1.25984, moderator=Water, shape=(17, 17), ndl=ndl)
asmbly.condensation_scheme = cond_spec
asmbly.few_group_condensation_scheme = few_grp_cond_spec
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
asmbly.save_diffusion_data(name+".npz")
