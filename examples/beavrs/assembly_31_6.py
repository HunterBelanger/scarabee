from scarabee import (
    NDLibrary,
    MaterialComposition,
    Material,
    Fraction,
    DensityUnits,
    set_output_file,
)
from scarabee.reseau import FuelPin, GuideTube, BurnablePoisonRod, PWRAssembly, Symmetry

name = "F31_6R"

set_output_file(name + "_out.txt")

ndl = NDLibrary()

# Define all Materials
Fuel31Comp = MaterialComposition(Fraction.Atoms, name="Fuel 3.1%")
Fuel31Comp.add_leu(3.1, 1.0)
Fuel31Comp.add_element("O", 2.0)
Fuel31 = Material(Fuel31Comp, 575.0, 10.30166, DensityUnits.g_cm3, ndl)

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

SS304Comp = MaterialComposition(Fraction.Weight, name="SS304")
SS304Comp.add_element("Si", 0.0060)
SS304Comp.add_element("Cr", 0.1900)
SS304Comp.add_element("Mn", 0.0200)
SS304Comp.add_element("Fe", 0.6840)
SS304Comp.add_element("Ni", 0.1000)
SS304 = Material(SS304Comp, 575.0, 8.03, DensityUnits.g_cm3, ndl)

BSiGlassComp = MaterialComposition(Fraction.Weight, name="Borosilicate Glass")
BSiGlassComp.add_element("O", 0.5481)
BSiGlassComp.add_element("Si", 0.3787)
BSiGlassComp.add_element("Al", 0.0344)
BSiGlassComp.add_nuclide("B10", 0.0071)
BSiGlassComp.add_nuclide("B11", 0.0317)
BSiGlass = Material(BSiGlassComp, 575.0, 2.26, DensityUnits.g_cm3, ndl)

AirComp = MaterialComposition(Fraction.Atoms, name="Air")
AirComp.add_element("O", 0.2095)
AirComp.add_element("N", 0.7809)
AirComp.add_element("Ar", 0.00933)
AirComp.add_element("C", 0.00027)
Air = Material(AirComp, 575.0, 0.000616, DensityUnits.g_cm3, ndl)

# Define a guide tube
gt = GuideTube(inner_radius=0.56134, outer_radius=0.60198, clad=Clad)

# Define guide tube with a burnable poison rod
bpr = BurnablePoisonRod(
    center=Air,
    clad=SS304,
    gap=He,
    poison=BSiGlass,
    center_radius=0.214,
    inner_clad_radius=0.23051,
    inner_gap_radius=0.24130,
    poison_radius=0.42672,
    outer_gap_radius=0.43688,
    outer_clad_radius=0.48387,
)

bp = GuideTube(inner_radius=0.56134, outer_radius=0.60198, clad=Clad, fill=bpr)

# Define fuel pin
fp = FuelPin(
    fuel=Fuel31,
    fuel_radius=0.39218,
    gap=He,
    gap_radius=0.40005,
    clad=Clad,
    clad_radius=0.45720,
)

cells = [
    [fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp],
    [fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp],
    [fp, fp, fp, fp, fp, gt, fp, fp, gt, fp, fp, bp, fp, fp, fp, fp, fp],
    [fp, fp, fp, gt, fp, fp, fp, fp, fp, fp, fp, fp, fp, bp, fp, fp, fp],
    [fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp],
    [fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp, bp, fp, fp],
    [fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp],
    [fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp, fp],
    [fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp, gt, fp, fp],
]

# Define assembly
asmbly = PWRAssembly(
    pitch=1.25984,
    assembly_pitch=21.50364,
    shape=(17, 17),
    symmetry=Symmetry.Half,
    moderator_pressure=15.5132,
    moderator_temp=575.0,
    boron_ppm=975.0,
    cells=cells,
    ndl=ndl,
)

asmbly.solve()
asmbly.diffusion_data.save(name + ".bin")
