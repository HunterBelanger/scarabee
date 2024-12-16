from scarabee import *
import numpy as np
import matplotlib.pyplot as plt

name = "F16_0"

set_output_file(name+"_out.txt")

ndl = NDLibrary()

# Define all Materials
Fuel16Comp = MaterialComposition(Fraction.Atoms, name="Fuel 1.6%")
Fuel16Comp.add_nuclide("O16",  4.5897e-02)
Fuel16Comp.add_nuclide("O17",  1.7436e-05)
Fuel16Comp.add_nuclide("O18",  9.2032e-05)
Fuel16Comp.add_nuclide("U234", 3.0131e-06)
Fuel16Comp.add_nuclide("U235", 3.7503e-04)
Fuel16Comp.add_nuclide("U238", 2.2625e-02)
Fuel16 = Material(Fuel16Comp, 575., ndl)

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

# Define a guide tube
gt = GuideTube(inner_radius=0.56134, outer_radius=0.60198, clad=Clad)

# Define fuel pin
fp = FuelPin(fuel=Fuel16, fuel_radius=0.39218, gap=He, gap_radius=0.40005,
             clad=Clad, clad_radius=0.45720)

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
