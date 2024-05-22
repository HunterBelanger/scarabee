from scarabee import *
import numpy as np
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=np.inf)

ndl = NDLibrary("endf8_shem281.h5")

UO2_16_c = MaterialComposition()
UO2_16_c.fractions = Fraction.Atoms
UO2_16_c.add_nuclide("O16",  4.5897e-02)
UO2_16_c.add_nuclide("O17",  1.7436e-05)
UO2_16_c.add_nuclide("O18",  9.2032e-05)
UO2_16_c.add_nuclide("U234", 3.0131e-06)
UO2_16_c.add_nuclide("U235", 3.7503e-04)
UO2_16_c.add_nuclide("U238", 2.2625e-02)
UO2_16 = Material(UO2_16_c, 800., ndl)

He_c = MaterialComposition()
He_c.fractions = Fraction.Atoms
He_c.add_nuclide("He3", 4.8089e-10)
He_c.add_nuclide("He4", 2.4044e-04)
He = Material(He_c, 575., ndl)

Z4_c = MaterialComposition()
Z4_c.fractions = Fraction.Atoms
Z4_c.add_nuclide("Cr50",  3.2962e-06)
Z4_c.add_nuclide("Cr52",  6.3564e-05)
Z4_c.add_nuclide("Cr53",  7.2076e-06)
Z4_c.add_nuclide("Cr54",  1.7941e-06)
Z4_c.add_nuclide("Fe54",  8.6698e-06)
Z4_c.add_nuclide("Fe56",  1.3610e-04)
Z4_c.add_nuclide("Fe57",  3.1431e-06)
Z4_c.add_nuclide("Fe58",  4.1829e-07)
Z4_c.add_nuclide("O16",   3.0744e-04)
Z4_c.add_nuclide("O17",   1.1680e-07)
Z4_c.add_nuclide("O18",   6.1648e-07)
Z4_c.add_nuclide("Sn112", 4.6735e-06)
Z4_c.add_nuclide("Sn114", 3.1799e-06)
Z4_c.add_nuclide("Sn115", 1.6381e-06)
Z4_c.add_nuclide("Sn116", 7.0055e-05)
Z4_c.add_nuclide("Sn117", 3.7003e-05)
Z4_c.add_nuclide("Sn118", 1.1669e-04)
Z4_c.add_nuclide("Sn119", 4.1387e-05)
Z4_c.add_nuclide("Sn120", 1.5697e-04)
Z4_c.add_nuclide("Sn122", 2.2308e-05)
Z4_c.add_nuclide("Sn124", 2.7897e-05)
Z4_c.add_nuclide("Zr90",  2.1828e-02)
Z4_c.add_nuclide("Zr91",  4.7601e-03)
Z4_c.add_nuclide("Zr92",  7.2759e-03)
Z4_c.add_nuclide("Zr94",  7.3734e-03)
Z4_c.add_nuclide("Zr96",  1.1879e-03)
Z4 = Material(Z4_c, 575., ndl)

RATIO_LW = (4.9456e-02)/((4.9456e-02)+(7.7035e-06))
RATIO_HW = 1. - RATIO_LW

H2O_c = MaterialComposition()
H2O_c.fractions= Fraction.Atoms
H2O_c.add_nuclide("B10",     7.9714e-06)
H2O_c.add_nuclide("B11",     3.2247e-05)
H2O_c.add_nuclide("H1_H2O",  4.9456e-02 + 7.7035e-06)
H2O_c.add_nuclide("O16",     2.4673e-02)
H2O_c.add_nuclide("O17",     9.3734e-06)
H2O_c.add_nuclide("O18",     4.9474e-05)
#H2O_c.add_nuclide("O16",     RATIO_LW*2.4673e-02)
#H2O_c.add_nuclide("O17",     RATIO_LW*9.3734e-06)
#H2O_c.add_nuclide("O18",     RATIO_LW*4.9474e-05)
#H2O_c.add_nuclide("H2_D2O",  7.7035e-06)
#H2O_c.add_nuclide("O16_D2O", RATIO_HW*2.4673e-02)
#H2O_c.add_nuclide("O17_D2O", RATIO_HW*9.3734e-06)
#H2O_c.add_nuclide("O18_D2O", RATIO_HW*4.9474e-05)
H2O = Material(H2O_c, 575., ndl)

fuel_rad  = 0.39218
he_rad    = 0.40005
clad_rad  = 0.45720

gt_wtr_rad = 0.50419
gt_cld_rad = 0.54610

pin_pitch = 1.25984
asm_pitch = 21.50364
gap_width = 0.5 * (asm_pitch - 17.*pin_pitch)

def calculate_fuel_dancoff_correction():
    gap_width = 0.5 * (asm_pitch - 17.*pin_pitch)

    # Make material cross sections from potential cross sections 
    Et = np.array([1.E5])
    Ea = np.array([1.E5])
    Es = np.array([[0.]])
    Fuel = CrossSection(Et, Ea, Es, "Fuel")

    Et = np.array([He.potential_xs])
    Ea = np.array([He.potential_xs])
    Gas = CrossSection(Et, Ea, Es, "Gas")

    Et = np.array([Z4.potential_xs])
    Ea = np.array([Z4.potential_xs])
    Clad = CrossSection(Et, Ea, Es, "Clad")

    Et = np.array([H2O.potential_xs])
    Ea = np.array([H2O.potential_xs])
    Mod = CrossSection(Et, Ea, Es, "Mod")

    # First we do transport on an isolated fuel pin
    print("Solving Isolate Fuel Pin")
    IFP = SimplePinCell([fuel_rad, he_rad, clad_rad], [Fuel, Gas, Clad, Mod], 20.*pin_pitch, 20.*pin_pitch)
    igeom = Cartesian2D([20.*pin_pitch], [20.*pin_pitch])
    igeom.set_tiles([IFP])
    imoc = MOCDriver(igeom)
    for i in range(imoc.nregions):
        xs = imoc.xs(i)
        if xs.name == "Fuel":
            imoc.set_extern_src(i, 0, 0.)
        elif xs.name == "Gap":
            imoc.set_extern_src(i, 0, He.potential_xs)
        elif xs.name == "Clad":
            imoc.set_extern_src(i, 0, Z4.potential_xs)
        else:
            imoc.set_extern_src(i, 0, H2O.potential_xs)
    imoc.sim_mode = SimulationMode.FixedSource
    imoc.generate_tracks(128, 0.02, YamamotoTabuchi6())
    imoc.flux_tolerance = 1.E-5
    imoc.solve()
    iso_flux = imoc.flux(0, 0)
    print()

    # First, we need to make the 4 corner pieces of the gap
    CG = EmptyCell(Mod, gap_width, gap_width) # Corner Gap

    ng = 17*5
    gap_cell_size = (asm_pitch-2.*gap_width) / ng
    
    # Vertical gap tiles
    VG_cell = EmptyCell(Mod, gap_width, gap_cell_size)
    VG = Cartesian2D([gap_width], ng*[gap_cell_size]) # Vertical Gap
    VG.set_tiles(ng*[VG_cell])

    # Horizontal gap tiles
    HG_cell = EmptyCell(Mod, gap_cell_size, gap_width)
    HG = Cartesian2D(ng*[gap_cell_size], [gap_width]) # Horizontal Gap
    HG.set_tiles(ng*[HG_cell])

    # Make the basic pin cells
    FP = SimplePinCell([fuel_rad, he_rad, clad_rad], [Fuel, Gas, Clad, Mod], pin_pitch, pin_pitch)
    GT = SimplePinCell([gt_wtr_rad, gt_cld_rad], [Mod, Clad, Mod], pin_pitch, pin_pitch)
    Pins = Cartesian2D(17*[pin_pitch], 17*[pin_pitch])
    Pins.set_tiles([FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP,
                    FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP,
                    FP, FP, FP, FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP, FP, FP, FP,
                    FP, FP, FP, GT, FP, FP, FP, FP, FP, FP, FP, FP, FP, GT, FP, FP, FP,
                    FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP,
                    FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP,
                    FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP,
                    FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP,
                    FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP,
                    FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP,
                    FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP,
                    FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP,
                    FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP,
                    FP, FP, FP, GT, FP, FP, FP, FP, FP, FP, FP, FP, FP, GT, FP, FP, FP,
                    FP, FP, FP, FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP, FP, FP, FP,
                    FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP,
                    FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP])
    
    geom = Cartesian2D([gap_width, 17.*pin_pitch, gap_width], [gap_width, 17.*pin_pitch, gap_width])
    geom.set_tiles([CG,  HG,   CG,
                    VG,  Pins, VG,
                    CG,  HG,   CG])

    print("Solving Fuel Assembly")
    moc = MOCDriver(geom)
    for i in range(moc.nregions):
        xs = moc.xs(i)
        if xs.name == "Fuel":
            moc.set_extern_src(i, 0, 0.)
        elif xs.name == "Gap":
            moc.set_extern_src(i, 0, He.potential_xs)
        elif xs.name == "Clad":
            moc.set_extern_src(i, 0, Z4.potential_xs)
        else:
            moc.set_extern_src(i, 0, H2O.potential_xs)
    moc.sim_mode = SimulationMode.FixedSource
    moc.generate_tracks(128, 0.02, YamamotoTabuchi6())
    moc.flux_tolerance = 1.E-5
    moc.solve()
    print()

    C = np.zeros((17, 17))
    x = moc.x_min + gap_width + 0.5*pin_pitch
    for i in range(17):
        y = moc.y_min + gap_width + 0.5*pin_pitch
        for j in range(17):
            xs = moc.xs(Vector(x, y), Direction(1., 0.))
            if xs.name == "Fuel":
                flux = moc.flux(Vector(x, y), Direction(1., 0.), 0)
                C[i,j] = (iso_flux - flux) / iso_flux
            y += pin_pitch
        x += pin_pitch
    
    return C

def calculate_clad_dancoff_correction():
    # Make material cross sections from potential cross sections 
    Et = np.array([UO2_16.potential_xs])
    Ea = np.array([UO2_16.potential_xs])
    Es = np.array([[0.]])
    Fuel = CrossSection(Et, Ea, Es, "Fuel")

    Et = np.array([He.potential_xs])
    Ea = np.array([He.potential_xs])
    Gas = CrossSection(Et, Ea, Es, "Gas")

    Et = np.array([1.E5])
    Ea = np.array([1.E5])
    Clad = CrossSection(Et, Ea, Es, "Clad")

    Et = np.array([H2O.potential_xs])
    Ea = np.array([H2O.potential_xs])
    Mod = CrossSection(Et, Ea, Es, "Mod")

    # First we do transport on an isolated fuel pin
    print("Solving Isolate Fuel Pin")
    IFP = SimplePinCell([fuel_rad, he_rad, clad_rad], [Fuel, Gas, Clad, Mod], 20.*pin_pitch, 20.*pin_pitch)
    igeom = Cartesian2D([20.*pin_pitch], [20.*pin_pitch])
    igeom.set_tiles([IFP])
    imoc = MOCDriver(igeom)
    for i in range(imoc.nregions):
        xs = imoc.xs(i)
        if xs.name == "Fuel":
            imoc.set_extern_src(i, 0, UO2_16.potential_xs)
        elif xs.name == "Gap":
            imoc.set_extern_src(i, 0, He.potential_xs)
        elif xs.name == "Clad":
            imoc.set_extern_src(i, 0, 0.)
        else:
            imoc.set_extern_src(i, 0, H2O.potential_xs)
    imoc.sim_mode = SimulationMode.FixedSource
    imoc.generate_tracks(128, 0.02, YamamotoTabuchi6())
    imoc.flux_tolerance = 1.E-5
    imoc.solve()
    iso_flux = imoc.flux(2, 0)
    print()

    # First, we need to make the 4 corner pieces of the gap
    CG = EmptyCell(Mod, gap_width, gap_width) # Corner Gap

    ng = 17*5
    gap_cell_size = (asm_pitch-2.*gap_width) / ng
    
    # Vertical gap tiles
    VG_cell = EmptyCell(Mod, gap_width, gap_cell_size)
    VG = Cartesian2D([gap_width], ng*[gap_cell_size]) # Vertical Gap
    VG.set_tiles(ng*[VG_cell])

    # Horizontal gap tiles
    HG_cell = EmptyCell(Mod, gap_cell_size, gap_width)
    HG = Cartesian2D(ng*[gap_cell_size], [gap_width]) # Horizontal Gap
    HG.set_tiles(ng*[HG_cell])

    # Make the basic pin cells
    FP = SimplePinCell([fuel_rad, he_rad, clad_rad], [Fuel, Gas, Clad, Mod], pin_pitch, pin_pitch)
    GT = SimplePinCell([gt_wtr_rad, gt_cld_rad], [Mod, Clad, Mod], pin_pitch, pin_pitch)
    Pins = Cartesian2D(17*[pin_pitch], 17*[pin_pitch])
    Pins.set_tiles([FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP,
                    FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP,
                    FP, FP, FP, FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP, FP, FP, FP,
                    FP, FP, FP, GT, FP, FP, FP, FP, FP, FP, FP, FP, FP, GT, FP, FP, FP,
                    FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP,
                    FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP,
                    FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP,
                    FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP,
                    FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP,
                    FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP,
                    FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP,
                    FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP,
                    FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP,
                    FP, FP, FP, GT, FP, FP, FP, FP, FP, FP, FP, FP, FP, GT, FP, FP, FP,
                    FP, FP, FP, FP, FP, GT, FP, FP, GT, FP, FP, GT, FP, FP, FP, FP, FP,
                    FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP,
                    FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP, FP])
    
    geom = Cartesian2D([gap_width, 17.*pin_pitch, gap_width], [gap_width, 17.*pin_pitch, gap_width])
    geom.set_tiles([CG,  HG,   CG,
                    VG,  Pins, VG,
                    CG,  HG,   CG])

    print("Solving Fuel Assembly")
    moc = MOCDriver(geom)
    for i in range(moc.nregions):
        xs = moc.xs(i)
        if xs.name == "Fuel":
            moc.set_extern_src(i, 0, UO2_16.potential_xs)
        elif xs.name == "Gap":
            moc.set_extern_src(i, 0, He.potential_xs)
        elif xs.name == "Clad":
            moc.set_extern_src(i, 0, 0.)
        else:
            moc.set_extern_src(i, 0, H2O.potential_xs)
    moc.sim_mode = SimulationMode.FixedSource
    moc.generate_tracks(128, 0.02, YamamotoTabuchi6())
    moc.flux_tolerance = 1.E-5
    moc.solve()
    print()

    C = np.zeros((17, 17))
    x = moc.x_min + gap_width + 0.5*pin_pitch
    for i in range(17):
        y = moc.y_min + gap_width + 0.5*pin_pitch
        for j in range(17):
            r = Vector(x+0.5*(he_rad + clad_rad), y)
            xs = moc.xs(r, Direction(1., 0.))
            if xs.name == "Clad":
                flux = moc.flux(r, Direction(1., 0.), 0)
                C[i,j] = (iso_flux - flux) / flux
            y += pin_pitch
        x += pin_pitch

    # Use local averaging for the dancoff correction in the guide tubes 
    for i in range(17):
        for j in range(17):
            if C[i,j] == 0.:
                C[i,j] = 0.25 * (C[i-1,j] + C[i+1,j] + C[i,j-1] + C[i,j+1])
    
    return C

def run_assembly(C_fuel, C_clad):
    Mod = H2O.dilution_xs(H2O.size*[1.E10], ndl)
    Gas = He.dilution_xs(He.size*[1.E10], ndl)

    # First, we need to make the 4 corner pieces of the gap
    CG = EmptyCell(Mod, gap_width, gap_width) # Corner Gap

    ng = 17*5
    gap_cell_size = (asm_pitch-2.*gap_width) / ng
    
    # Vertical gap tiles
    VG_cell = EmptyCell(Mod, gap_width, gap_cell_size)
    VG = Cartesian2D([gap_width], ng*[gap_cell_size]) # Vertical Gap
    VG.set_tiles(ng*[VG_cell])

    # Horizontal gap tiles
    HG_cell = EmptyCell(Mod, gap_cell_size, gap_width)
    HG = Cartesian2D(ng*[gap_cell_size], [gap_width]) # Horizontal Gap
    HG.set_tiles(ng*[HG_cell])

    # Make all the pins
    Pins = Cartesian2D(17*[pin_pitch], 17*[pin_pitch])

    gt_pins = [       (5,2), (8,2), (11,2),
                  (3,3),                 (13,3),
               (2,5), (5,5), (8,5), (11,5), (14,5),
               (2,8), (5,8), (8,8), (11,8), (14,8),
               (2,11), (5,11), (8,11), (11,11), (14,11),
                  (3,13),                 (13,13),
                        (5,14), (8,14), (11,14)]

    pins_list = np.array(17*[17*[None]])
    for j in range(9):
        for i in range(j, 9):
            print("Making pin ({},{})".format(i, j))
            if (i, j) not in gt_pins:
                Ee = 1. / (2. * fuel_rad)
                Fuel = UO2_16.carlvik_xs(C_fuel[i,j], Ee, ndl)

                Ee = 1. / (2. * (clad_rad-he_rad))
                Clad = Z4.roman_xs(C_clad[i,j], Ee, ndl)

                FP = PinCell([0.8*fuel_rad, fuel_rad, he_rad, clad_rad, 0.5*pin_pitch], [Fuel, Fuel, Gas, Clad, Mod, Mod], pin_pitch, pin_pitch)
                pins_list[j,i] = FP
            else:
                Ee = 1. / (2. * (gt_cld_rad - gt_wtr_rad))
                Clad = Z4.roman_xs(C_clad[i,j], Ee, ndl)

                GT = PinCell([0.5*gt_wtr_rad, gt_wtr_rad, gt_cld_rad, 0.5*pin_pitch], [Mod, Mod, Clad, Mod, Mod], pin_pitch, pin_pitch)
                pins_list[j,i] = GT

            # Make symetric copies
            pins_list[i,j]    = pins_list[j,i]
            pins_list[16-i,j] = pins_list[j,i]
            pins_list[i,16-j] = pins_list[j,i]
            pins_list[16-j,i] = pins_list[j,i]
            pins_list[j,16-i] = pins_list[j,i]
            pins_list[16-i,16-j] = pins_list[j,i]
            pins_list[16-j,16-i] = pins_list[j,i]

    Pins.set_tiles(list(pins_list.flatten()))

    geom = Cartesian2D([gap_width, 17.*pin_pitch, gap_width], [gap_width, 17.*pin_pitch, gap_width])
    geom.set_tiles([CG,  HG,   CG,
                    VG,  Pins, VG,
                    CG,  HG,   CG])

    moc = MOCDriver(geom)
    moc.generate_tracks(128, 0.01, YamamotoTabuchi6())
    moc.flux_tolerance = 1.E-5
    moc.solve()

def main():
    C_fuel = calculate_fuel_dancoff_correction()
    C_clad = calculate_clad_dancoff_correction()
    
    print("Fuel Dancoff Correction Factors:")
    print(C_fuel)
    print()

    print("Clad Dancoff Correction Factors:")
    print(C_clad)
    print()

    run_assembly(C_fuel, C_clad)


if __name__ == "__main__":
    main()