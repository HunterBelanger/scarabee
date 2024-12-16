from scarabee import *

name = "F31_0"

set_output_file(name+"_reflector_out.txt")

ndl = NDLibrary()

# Define all Materials
Fuel31Comp = MaterialComposition(name="Fuel 3.1%")
Fuel31Comp.add_nuclide("O16",  4.5853e-02)
Fuel31Comp.add_nuclide("O17",  1.7420e-05)
Fuel31Comp.add_nuclide("O18",  9.1942e-05)
Fuel31Comp.add_nuclide("U234", 5.7987e-06)
Fuel31Comp.add_nuclide("U235", 7.2175e-04)
Fuel31Comp.add_nuclide("U238", 2.2253e-02)
Fuel31 = Material(Fuel31Comp, 575., ndl)
Fuel31.max_legendre_order = 3

CladComp = MaterialComposition(Fraction.Weight, name="Zircaloy 4")
CladComp.add_element('O', 0.00125)
CladComp.add_element('Cr', 0.0010)
CladComp.add_element('Fe', 0.0021)
CladComp.add_element('Zr', 0.98115)
CladComp.add_element('Sn', 0.0145)
Clad = Material(CladComp, 575., 6.55, DensityUnits.g_cm3 , ndl)
Clad.max_legendre_order = 3

HeComp = MaterialComposition(Fraction.Atoms, name="He Gas")
HeComp.add_element("He", 1.)
He = Material(HeComp, 575., 0.0015981, DensityUnits.g_cm3, ndl) 
He.max_legendre_order = 3

WaterComp = MaterialComposition(Fraction.Atoms, name="Water")
WaterComp.add_nuclide("H1_H2O_TC",  4.9456e-02 + 7.7035e-06)
WaterComp.add_element("B", 7.9714e-06 + 3.2247e-05)
WaterComp.add_element("O", 2.4673e-02 + 9.3734e-06 + 4.9474e-05)
Water = Material(WaterComp, 575., ndl)
Water.max_legendre_order = 3

SS304Comp = MaterialComposition(Fraction.Weight, name="SS304")
SS304Comp.add_element('Si', 0.0060)
SS304Comp.add_element('Cr', 0.1900)
SS304Comp.add_element('Mn', 0.0200)
SS304Comp.add_element('Fe', 0.6840)
SS304Comp.add_element('Ni', 0.1000)
SS304 = Material(SS304Comp, 575., 8.03, DensityUnits.g_cm3, ndl)
SS304.max_legendre_order = 3

# Define a guide tube
gt = GuideTube(inner_radius=0.56134, outer_radius=0.60198, clad=Clad)

# Define fuel pin
fp = FuelPin(fuel=Fuel31, fuel_radius=0.39218, gap=He, gap_radius=0.40005, clad=Clad, clad_radius=0.45720)

# Define assembly
asmbly = PWRAssembly(pitch=1.25984, moderator=Water, shape=(17, 17), ndl=ndl)
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
asmbly.save_diffusion_data(name+".bin")

# Solve reflector
print()
refl = Reflector(asmbly.average_fuel_pin, moderator=asmbly.moderator_xs,
                 assembly_width=21.50364, gap_width=0.1627,
                 baffle_width=2.2225, baffle=SS304, ndl=ndl)
refl.solve()
refl.save_diffusion_data("reflector.bin")
