from scarabee import *

name = "reflector"

set_output_file(name+"_out.txt")

ndl = NDLibrary('/home/hunter/projects/scarabee/data/endf8_shem281.h5')

cond_spec = [[0, 3], [4, 8], [9, 11], [12, 13], [14, 17], [18, 22], [23, 25], [26, 29],
             [30, 32], [33, 36], [37, 39], [40, 42], [43, 48], [49, 52], [53, 55],
             [56, 58], [59, 60], [61, 64], [65, 67], [68, 70], [71, 73], [74, 76],
             [77, 80], [81, 82], [83, 84], [85, 88], [89, 92], [93, 128], [129, 148],
             [149, 195], [196, 249], [250, 264], [265, 280]]

few_grp_cond_spec = [[0, 30], [31, 32]]

# Define all Materials
Fuel31Comp = MaterialComposition()
Fuel31Comp.fractions = Fraction.Atoms
Fuel31Comp.add_nuclide("O16",  4.5853e-02)
Fuel31Comp.add_nuclide("O17",  1.7420e-05)
Fuel31Comp.add_nuclide("O18",  9.1942e-05)
Fuel31Comp.add_nuclide("U234", 5.7987e-06)
Fuel31Comp.add_nuclide("U235", 7.2175e-04)
Fuel31Comp.add_nuclide("U238", 2.2253e-02)
Fuel31 = Material(Fuel31Comp, 575., ndl)

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

Bafflecomp = MaterialComposition()
Bafflecomp.fractions = Fraction.Atoms
Bafflecomp.add_nuclide("Cr50", 7.6778E-4)
Bafflecomp.add_nuclide("Cr52", 1.4806E-2)
Bafflecomp.add_nuclide("Cr53", 1.6789E-3)
Bafflecomp.add_nuclide("Cr54", 4.1791E-4)
Bafflecomp.add_nuclide("Fe54", 3.4620E-3)
Bafflecomp.add_nuclide("Fe56", 5.4345E-2)
Bafflecomp.add_nuclide("Fe57", 1.2551E-3)
Bafflecomp.add_nuclide("Fe58", 1.6703E-4)
Bafflecomp.add_nuclide("Mn55", 1.7604E-3)
Bafflecomp.add_nuclide("Ni58", 5.6089E-3)
Bafflecomp.add_nuclide("Ni60", 2.1605E-3)
Bafflecomp.add_nuclide("Ni61", 9.3917E-5)
Bafflecomp.add_nuclide("Ni62", 2.9945E-4)
Bafflecomp.add_nuclide("Ni64", 7.6261E-5)
Bafflecomp.add_nuclide("Si28", 9.5281E-4)
Bafflecomp.add_nuclide("Si29", 4.8381E-5)
Bafflecomp.add_nuclide("Si30", 3.1893E-5)
Baffle = Material(Bafflecomp, 575., ndl)

# Define a guide tube
gt = GuideTube(inner_radius=0.56134, outer_radius=0.60198, clad=Clad)

# Define fuel pin
fp = FuelPin(fuel=Fuel31, fuel_radius=0.39218, gap=He, gap_radius=0.40005, clad=Clad, clad_radius=0.45720)

# Define assembly
refl = PWRReflector(pitch=1.25984, moderator=Water, shape=(17, 17), gap_width= 0.1627, baffle_width=2.2225, baffle=Baffle, ndl=ndl)
refl.condensation_scheme = cond_spec
refl.few_group_condensation_scheme = few_grp_cond_spec
refl.num_azimuthal_angles = 64
refl.track_spacing = 0.01
refl.plot_assembly = True
refl.pins = [fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp,
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
refl.solve()
refl.save_diffusion_data(name+".npz")