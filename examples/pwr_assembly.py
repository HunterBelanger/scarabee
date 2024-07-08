from scarabee import *
import numpy as np
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=np.inf)
np.set_printoptions(4)

ndl = NDLibrary('/mnt/c/Users/hunte/Documents/nuclear_data/scarabee/endf8_shem281.h5')
#ndl = NDLibrary('/mnt/c/Users/hunte/Documents/nuclear_data/scarabee/endf8_xmas172.h5')

UO2comp = MaterialComposition()
UO2comp.fractions = Fraction.Atoms
UO2comp.add_nuclide("U235", 7.0803E-4)
UO2comp.add_nuclide("U238", 2.2604E-2)
UO2comp.add_nuclide("O16",  4.6624E-2)
UO2 = Material(UO2comp, 293.6, ndl)

Zirccomp = MaterialComposition()
Zirccomp.fractions = Fraction.Atoms
Zirccomp.add_nuclide("Zr90", 2.2200E-2)
Zirccomp.add_nuclide("Zr91", 4.8280E-3)
Zirccomp.add_nuclide("Zr92", 7.3713E-3)
Zirccomp.add_nuclide("Zr94", 7.5006E-3)
Zirccomp.add_nuclide("Zr96", 1.2070E-3)
Zirc = Material(Zirccomp, 293.6, ndl)

Watercomp = MaterialComposition()
Watercomp.fractions = Fraction.Atoms
Watercomp.add_nuclide("H1_H2O", 6.6630E-2)
Watercomp.add_nuclide("O16",    3.3315E-2)
Water = Material(Watercomp, 293.6, ndl)

fp = FuelPin(fuel=UO2, fuel_radius=0.4095, clad=Zirc, clad_width=(0.475-0.3095))
gt = GuideTube(inner_radius=0.57, outer_radius=0.61, clad=Zirc)

asmbly = PWRAssembly(pitch=1.26, gap=0., moderator=Water, shape=(17, 17))
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

asmbly.get_fuel_dancoff_corrections()
asmbly.get_clad_dancoff_corrections()

tmp = np.array(asmbly.fuel_dancoff_corrections)
tmp = tmp.reshape((17, 17))
print(tmp)
print()

#tmp = np.array(asmbly.clad_dancoff_corrections)
#tmp = tmp.reshape((17, 17))
#print(tmp)
#print()

