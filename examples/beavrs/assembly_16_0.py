from scarabee import (
    NDLibrary,
    MaterialComposition,
    Material,
    Fraction,
    DensityUnits,
    set_output_file,
)
from scarabee.reseau import FuelPin, GuideTube, PWRAssembly, Symmetry

name = "F16_0"

set_output_file(name + "_out.txt")

ndl = NDLibrary()

# Define all Materials
Fuel16Comp = MaterialComposition(Fraction.Atoms, name="Fuel 1.6%")
Fuel16Comp.add_leu(1.6, 1.0)
Fuel16Comp.add_element("O", 2.0)
Fuel16 = Material(Fuel16Comp, 575.0, 10.31341, DensityUnits.g_cm3, ndl)

CladComp = MaterialComposition(Fraction.Weight, name="Zircaloy 4")
CladComp.add_element("O", 0.00125)
CladComp.add_element("Cr", 0.0010)
CladComp.add_element("Fe", 0.0021)
CladComp.add_element("Zr", 0.98115)
CladComp.add_element("Sn", 0.0145)
Clad = Material(CladComp, 575.0, 6.55, DensityUnits.g_cm3, ndl)

HeComp = MaterialComposition(Fraction.Atoms, name="He Gas")
HeComp.add_element("He", 1.0)
He = Material(HeComp, 575.0, 0.0015981, DensityUnits.g_cm3, ndl)

# Define a guide tube
gt = GuideTube(inner_radius=0.56134, outer_radius=0.60198, clad=Clad)

# Define fuel pin
fp = FuelPin(
    fuel=Fuel16,
    fuel_radius=0.39218,
    gap=He,
    gap_radius=0.40005,
    clad=Clad,
    clad_radius=0.45720,
)

cells = [
    [fp, fp, fp, fp, fp, fp, fp, fp, fp],
    [fp, fp, fp, fp, fp, fp, fp, fp, fp],
    [gt, fp, fp, gt, fp, fp, fp, fp, fp],
    [fp, fp, fp, fp, fp, gt, fp, fp, fp],
    [fp, fp, fp, fp, fp, fp, fp, fp, fp],
    [gt, fp, fp, gt, fp, fp, gt, fp, fp],
    [fp, fp, fp, fp, fp, fp, fp, fp, fp],
    [fp, fp, fp, fp, fp, fp, fp, fp, fp],
    [gt, fp, fp, gt, fp, fp, gt, fp, fp],
]

# Define assembly
asmbly = PWRAssembly(
    pitch=1.25984,
    assembly_pitch=21.50364,
    shape=(17, 17),
    symmetry=Symmetry.Quarter,
    moderator_pressure=15.5132,
    moderator_temp=575.0,
    boron_ppm=975.0,
    cells=cells,
    ndl=ndl,
)

asmbly.solve()

diffusion_data = asmbly._compute_diffusion_data()
diffusion_data.save(name + ".bin")
