from scarabee import *

name = "F24_16"

set_output_file(name+"_out.txt")

ndl = NDLibrary()

# Define all Materials
Fuel24Comp = MaterialComposition(name="Fuel 2.4%")
Fuel24Comp.add_nuclide("O16",  4.5830e-02)
Fuel24Comp.add_nuclide("O17",  1.7411e-05)
Fuel24Comp.add_nuclide("O18",  9.1898e-05)
Fuel24Comp.add_nuclide("U234", 4.4842e-06)
Fuel24Comp.add_nuclide("U235", 5.5814e-04)
Fuel24Comp.add_nuclide("U238", 2.2407e-02)
Fuel24 = Material(Fuel24Comp, 575., ndl)

CladComp = MaterialComposition(Fraction.Weight, name="Zircaloy 4")
CladComp.add_element('O', 0.00125)
CladComp.add_element('Cr', 0.0010)
CladComp.add_element('Fe', 0.0021)
CladComp.add_element('Zr', 0.98115)
CladComp.add_element('Sn', 0.0145)
Clad = Material(CladComp, 575., 6.55, DensityUnits.g_cm3 , ndl)

HeComp = MaterialComposition(Fraction.Atoms, name="He Gas")
HeComp.add_element("He", 1.)
He = Material(HeComp, 575., 0.0015981, DensityUnits.g_cm3, ndl) 

WaterComp = MaterialComposition(Fraction.Atoms, name="Water")
WaterComp.add_nuclide("H1_H2O_TC",  4.9456e-02 + 7.7035e-06)
WaterComp.add_element("B", 7.9714e-06 + 3.2247e-05)
WaterComp.add_element("O", 2.4673e-02 + 9.3734e-06 + 4.9474e-05)
Water = Material(WaterComp, 575., ndl)

SS304Comp = MaterialComposition(Fraction.Weight, name="SS304")
SS304Comp.add_element('Si', 0.0060)
SS304Comp.add_element('Cr', 0.1900)
SS304Comp.add_element('Mn', 0.0200)
SS304Comp.add_element('Fe', 0.6840)
SS304Comp.add_element('Ni', 0.1000)
SS304 = Material(SS304Comp, 575., 8.03, DensityUnits.g_cm3, ndl)

BSiGlassComp = MaterialComposition(Fraction.Weight, name="Borosilicate Glass")
BSiGlassComp.add_element('O',   0.5481)
BSiGlassComp.add_element('Si',  0.3787)
BSiGlassComp.add_element('Al',  0.0344)
BSiGlassComp.add_nuclide('B10', 0.0071)
BSiGlassComp.add_nuclide('B11', 0.0317)
BSiGlass = Material(BSiGlassComp, 575., 2.26, DensityUnits.g_cm3, ndl)

AirComp = MaterialComposition(Fraction.Atoms, name="Air")
AirComp.add_element('O',  0.2095)
AirComp.add_element('N',  0.7809)
AirComp.add_element('Ar', 0.00933)
AirComp.add_element('C',  0.00027)
Air = Material(AirComp, 575., 0.000616, DensityUnits.g_cm3, ndl)

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
asmbly.save_diffusion_data(name+".bin")
