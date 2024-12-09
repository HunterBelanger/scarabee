from scarabee import *

name = "F24_16"

set_output_file(name+"_out.txt")

ndl = NDLibrary()

# Define all Materials
Fuel24Comp = MaterialComposition()
Fuel24Comp.fractions = Fraction.Atoms
Fuel24Comp.add_nuclide("O16",  4.5830e-02)
Fuel24Comp.add_nuclide("O17",  1.7411e-05)
Fuel24Comp.add_nuclide("O18",  9.1898e-05)
Fuel24Comp.add_nuclide("U234", 4.4842e-06)
Fuel24Comp.add_nuclide("U235", 5.5814e-04)
Fuel24Comp.add_nuclide("U238", 2.2407e-02)
Fuel24 = Material(Fuel24Comp, 575., ndl)

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

# SS304
SS304Comp = MaterialComposition()
SS304Comp.add_nuclide("Cr50", 7.6778e-04)
SS304Comp.add_nuclide("Cr52", 1.4806e-02)
SS304Comp.add_nuclide("Cr53", 1.6789e-03)
SS304Comp.add_nuclide("Cr54", 4.1791e-04)
SS304Comp.add_nuclide("Fe54", 3.4620e-03)
SS304Comp.add_nuclide("Fe56", 5.4345e-02)
SS304Comp.add_nuclide("Fe57", 1.2551e-03)
SS304Comp.add_nuclide("Fe58", 1.6703e-04)
SS304Comp.add_nuclide("Mn55", 1.7604e-03)
SS304Comp.add_nuclide("Ni58", 5.6089e-03)
SS304Comp.add_nuclide("Ni60", 2.1605e-03)
SS304Comp.add_nuclide("Ni61", 9.3917e-05)
SS304Comp.add_nuclide("Ni62", 2.9945e-04)
SS304Comp.add_nuclide("Ni64", 7.6261e-05)
SS304Comp.add_nuclide("Si28", 9.5281e-04)
SS304Comp.add_nuclide("Si29", 4.8381e-05)
SS304Comp.add_nuclide("Si30", 3.1893e-05)
SS304 = Material(SS304Comp, 575., ndl)

BSiGlassComp = MaterialComposition()
BSiGlassComp.add_nuclide("Al27", 1.7352e-03)
BSiGlassComp.add_nuclide("B10",  9.6506e-04)
BSiGlassComp.add_nuclide("B11",  3.9189e-03)
BSiGlassComp.add_nuclide("O16",  4.6514e-02)
BSiGlassComp.add_nuclide("O17",  1.7671e-05)
BSiGlassComp.add_nuclide("O18",  9.3268e-05)
BSiGlassComp.add_nuclide("Si28", 1.6926e-02)
BSiGlassComp.add_nuclide("Si29", 8.5944e-04)
BSiGlassComp.add_nuclide("Si30", 5.6654e-04)
BSiGlass = Material(BSiGlassComp, 575., ndl)

AirComp = MaterialComposition()
AirComp.add_nuclide("Ar36", 7.8730e-09)
AirComp.add_nuclide("Ar38", 1.4844e-09)
AirComp.add_nuclide("Ar40", 2.3506e-06)
AirComp.add_nuclide("C12",  6.7539e-08)
AirComp.add_nuclide("C13",  7.5658e-10)
AirComp.add_nuclide("N14",  1.9680e-04)
AirComp.add_nuclide("N15",  7.2354e-07)
AirComp.add_nuclide("O16",  5.2866e-05)
AirComp.add_nuclide("O17",  2.0084e-08)
AirComp.add_nuclide("O18",  1.0601e-07)
Air = Material(AirComp, 575., ndl)

# Define a guide tube
gt = GuideTube(inner_radius=0.56134, outer_radius=0.60198, clad=Clad)

# Define fuel pin
fp = FuelPin(fuel=Fuel24, fuel_radius=0.39218, gap=He, gap_radius=0.40005,
             clad=Clad, clad_radius=0.45720)

# Define a poison pin
bp = BurnablePoisonPin(center=Air, center_radius=0.214, poison_clad=SS304,
                       inner_poison_clad_radius=0.23051, gap=He,
                       inner_gap_radius=0.24130, outer_gap_radius=0.43688,
                       poison=BSiGlass, poison_radius=0.42672,
                       outer_poison_clad_radius=0.48387, guide_tube_clad=Clad,
                       inner_moderator_radius=0.56134,
                       guide_tube_radius=0.60198)

# Define assembly
asmbly = PWRAssembly(pitch=1.25984, moderator=Water, shape=(17, 17), ndl=ndl)
asmbly.pins = [fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
               fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
               fp, fp, fp, fp, fp, bp, fp, fp, bp, fp, fp, bp, fp, fp, fp, fp, fp,
               fp, fp, fp, bp, fp, fp, fp, fp, fp, fp, fp, fp, fp, bp, fp, fp, fp,
               fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
               fp, fp, bp, fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp, bp, fp, fp,
               fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
               fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
               fp, fp, bp, fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp, bp, fp, fp,
               fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
               fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
               fp, fp, bp, fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp, bp, fp, fp,
               fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
               fp, fp, fp, bp, fp, fp, fp, fp, fp, fp, fp, fp, fp, bp, fp, fp, fp,
               fp, fp, fp, fp, fp, bp, fp, fp, bp, fp, fp, bp, fp, fp, fp, fp, fp,
               fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
               fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp]
asmbly.solve()
asmbly.save_diffusion_data(name+".npz")
