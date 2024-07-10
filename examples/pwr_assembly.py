from scarabee import *
import numpy as np
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=np.inf)
np.set_printoptions(4)

ndl = NDLibrary('/mnt/c/Users/hunte/Documents/nuclear_data/endf8_shem281.h5')
#ndl = NDLibrary('C:\\Users\\hunte\\Documents\\nuclear_data\\endf8_shem281.h5')

cond_spec = [[0, 3], [4, 8], [9, 11], [12, 13], [14, 17], [18, 22], [23, 25], [26, 29],
             [30, 32], [33, 36], [37, 39], [40, 42], [43, 48], [49, 52], [53, 55],
             [56, 58], [59, 60], [61, 64], [65, 67], [68, 70], [71, 73], [74, 76],
             [77, 80], [81, 82], [83, 84], [85, 88], [89, 92], [93, 128], [129, 148],
             [149, 195], [196, 249], [250, 264], [265, 280]]

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

fp = FuelPin(fuel=UO2, fuel_radius=0.4095, clad=Zirc, clad_width=(0.475-0.4095))
gt = GuideTube(inner_radius=0.57, outer_radius=0.61, clad=Zirc)

asmbly = PWRAssembly(pitch=1.26, gap=0., moderator=Water, shape=(17, 17), ndl=ndl)
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

asmbly.condensation_scheme = cond_spec

asmbly.solve()

flux, x, y = asmbly.moc.rasterize_flux(500, 500)
for g in range(asmbly.moc.ngroups):
  plt.title("Flux in group {}".format(g+1))
  plt.pcolormesh(x, y, flux[g,:,:], cmap='jet')
  plt.show()

