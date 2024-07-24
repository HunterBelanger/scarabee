from scarabee import *
import matplotlib.pyplot as plt
import numpy as np
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=np.inf)
np.set_printoptions(4)

ndl = NDLibrary('/mnt/c/Users/BELANH2/Documents/nuclear_data/endf8_shem281.h5')
#ndl = NDLibrary('C:\\Users\\BELANH2\\Documents\\nuclear_data\\endf8_shem281.h5')

cond_spec = [[0, 3], [4, 8], [9, 11], [12, 13], [14, 17], [18, 22], [23, 25], [26, 29],
             [30, 32], [33, 36], [37, 39], [40, 42], [43, 48], [49, 52], [53, 55],
             [56, 58], [59, 60], [61, 64], [65, 67], [68, 70], [71, 73], [74, 76],
             [77, 80], [81, 82], [83, 84], [85, 88], [89, 92], [93, 128], [129, 148],
             [149, 195], [196, 249], [250, 264], [265, 280]]

few_grp_cond_spec = [[0, 30], [31, 32]]

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

fp = FuelPin(fuel=UO2, fuel_radius=0.4095, clad=Zirc, clad_width=(0.475-0.4095), fuel_rings=8)
gt = GuideTube(inner_radius=0.57, outer_radius=0.61, clad=Zirc)

asmbly = PWRAssembly(pitch=1.26, moderator=Water, shape=(17, 17), ndl=ndl)
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
asmbly.few_group_condensation_scheme = few_grp_cond_spec
asmbly.num_azimuthal_angles = 128
asmbly.track_spacing = 0.01
asmbly.plot_assembly = True

asmbly.solve()

flux, x, y = asmbly.moc.rasterize_flux(500, 500)
for g in range(asmbly.moc.ngroups):
  plt.title("Flux in group {}".format(g+1))
  plt.pcolormesh(x, y, flux[g,:,:], cmap='jet')
  plt.show()

print()
print("Performing Reflector Calculation")
print()

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
Baffle = Material(Bafflecomp, 293.6, ndl)

refl = Reflector(asmbly.average_fuel_pin, moderator=asmbly.moderator_xs, assembly_width=21.50364, gap_width=0.1627, baffle_width=2.2225, baffle=Baffle, ndl=ndl)
refl.condensation_scheme = [[0, 249], [250, 280]]
refl.solve()
