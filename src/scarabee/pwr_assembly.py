from _scarabee import *
import numpy as np
from typing import Tuple, List, Optional
from enum import Enum

class FuelPin:
    def __init__(self, fuel: Material, fuel_radius: float, clad: Material, clad_width: float, gap: Optional[Material] = None, gap_width: Optional[float] = None):
        self.fuel = fuel
        self.fuel_radius = fuel_radius
        self.clad = clad
        self.clad_width = clad_width
        self.gap = gap
        self.gap_width = gap_width

        if self.gap is not None and self.gap_width is None:
            raise RuntimeError('Fuel gap material is provided, but no gap width.')
        elif self.gap_width is not None and self.gap is None:
            raise RuntimeError('Fuel gap width provided, but no gap material.')

    def make_fuel_dancoff_cell(self, pitch: float, moderator: Material):
        # We first determine all the radii
        radii = []
        radii.append(self.fuel_radius)
        if self.gap is not None:
            radii.append(radii[-1] + self.gap_width)
        radii.append(radii[-1] + self.clad_width)

        mats = []

        # This returns a cell for calculating the fuel pin Dancoff factor.
        # As such, the fuel XS has infinite values.
        Et = np.array([1.E5])
        Ea = np.array([1.E5])
        Es = np.array([[0.]])
        Fuel = CrossSection(Et, Ea, Es, "Fuel")
        mats.append(Fuel)

        if self.gap is not None:
            Et[0] = self.gap.potential_xs
            Ea[0] = self.gap.potential_xs
            gap = CrossSection(Et, Ea, Es, "Gap")
            mats.append(Gap)

        Et[0] = self.clad.potential_xs
        Ea[0] = self.clad.potential_xs
        Clad = CrossSection(Et, Ea, Es, "Clad")
        mats.append(Clad)

        Et[0] = moderator.potential_xs
        Ea[0] = moderator.potential_xs
        Mod = CrossSection(Et, Ea, Es, "Moderator")
        mats.append(Mod)

        return SimplePinCell(radii, mats, pitch, pitch)

    def make_clad_dancoff_cell(self, pitch: float, moderator: Material):
        # We first determine all the radii
        radii = []
        radii.append(self.fuel_radius)
        if self.gap is not None:
            radii.append(radii[-1] + self.gap_width)
        radii.append(radii[-1] + self.clad_width)

        mats = []
 
        Et = np.array([self.fuel.potential_xs])
        Ea = np.array([self.fuel.potential_xs])
        Es = np.array([[0.]])
        Fuel = CrossSection(Et, Ea, Es, "Fuel")
        mats.append(Fuel)

        if self.gap is not None:
            Et[0] = self.gap.potential_xs
            Ea[0] = self.gap.potential_xs
            gap = CrossSection(Et, Ea, Es, "Gap")
            mats.append(Gap)
        
        # This returns a cell for calculating the fuel pin Dancoff factor.
        # As such, the clad XS has infinite values.
        Et[0] = 1.E5
        Ea[0] = 1.E5
        Clad = CrossSection(Et, Ea, Es, "Clad")
        mats.append(Clad)

        Et[0] = moderator.potential_xs
        Ea[0] = moderator.potential_xs
        Mod = CrossSection(Et, Ea, Es, "Moderator")
        mats.append(Mod)

        return SimplePinCell(radii, mats, pitch, pitch)

    def make_cylindrical_cell(self, pitch: float, dancoff_fuel: float, moderator: Material, ndl: NDLibrary, dancoff_clad: Optional[float] = None, clad_dilution = 1.E10):
        # We first determine all the radii
        radii = []
        radii.append(self.fuel_radius)
        if self.gap is not None:
            radii.append(radii[-1] + self.gap_width)
        radii.append(radii[-1] + self.clad_width)
        radii.append(np.sqrt(pitch*pitch / np.pi))
        
        # Next, we determine all the materials.
        # This requires applying self shielding to the fuel and cladding
        mats = []

        # First, treat the fuel
        Ee = 1. / (2. * self.fuel_radius) # Fuel escape xs
        mats.append(self.fuel.carlvik_xs(dancoff_fuel, Ee, ndl))

        # Next, add the gap (if present)
        if self.gap is not None:
            mats.append(self.gap.dilution_xs(self.gap.size*[1.E10], ndl))

        # Add the cladding
        if dancoff_clad is not None:
            Ee = 1. / (2. * self.clad_width)
            mats.append(self.clad.roman_xs(dancoff_clad, Ee, ndl))
        else:
            mats.append(self.clad.dilution_xs(self.clad.size*[clad_dilution], ndl))

        # Finally, add moderator
        mats.append(moderator.dilution_xs(moderator.size*[1.E10], ndl))

        return CylindricalCell(radii, mats)

class GuideTube:
    def __init__(self, inner_radius: float, outer_radius: float, clad: Material):
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.clad = clad

        if self.outer_radius <= self.inner_radius:
            raise RuntimeError("Outer radius must be > inner radius.")

    def make_clad_dancoff_cell(self, pitch: float, moderator: Material):
        # We first determine all the radii
        radii = []
        radii.append(self.inner_radius)
        radii.append(self.outer_radius)

        mats = []
 
        Et = np.array([moderator.potential_xs])
        Ea = np.array([moderator.potential_xs])
        Es = np.array([[0.]])
        InnerMod = CrossSection(Et, Ea, Es, "Inner Moderator")
        mats.append(InnerMod)

        # This returns a cell for calculating the fuel pin Dancoff factor.
        # As such, the clad XS has infinite values.
        Et[0] = 1.E5
        Ea[0] = 1.E5
        Clad = CrossSection(Et, Ea, Es, "Clad")
        mats.append(Clad)

        Et[0] = moderator.potential_xs
        Ea[0] = moderator.potential_xs
        OuterMod = CrossSection(Et, Ea, Es, "Outer Moderator")
        mats.append(OuterMod)

        return SimplePinCell(radii, mats, pitch, pitch)

    def make_fuel_dancoff_cell(self, pitch: float, moderator: Material):
        # We first determine all the radii
        radii = []
        radii.append(self.inner_radius)
        radii.append(self.outer_radius)

        mats = []
 
        Et = np.array([moderator.potential_xs])
        Ea = np.array([moderator.potential_xs])
        Es = np.array([[0.]])
        InnerMod = CrossSection(Et, Ea, Es, "Inner Moderator")
        mats.append(InnerMod)

        Et[0] = self.clad.potential_xs
        Ea[0] = self.clad.potential_xs
        Clad = CrossSection(Et, Ea, Es, "Clad")
        mats.append(Clad)

        Et[0] = moderator.potential_xs
        Ea[0] = moderator.potential_xs
        OuterMod = CrossSection(Et, Ea, Es, "Outer Moderator")
        mats.append(OuterMod)

        return SimplePinCell(radii, mats, pitch, pitch)

    def make_cylindrical_cell(self, pitch: float, moderator: Material, buffer_radius: float, buffer: CrossSection, ndl: NDLibrary, dancoff_clad: Optional[float] = None, clad_dilution: float = 1.E10):
        # We first determine all the radii
        radii = []
        radii.append(self.fuel_radius)
        if self.gap is not None:
            radii.append(radii[-1] + self.gap_width)
        radii.append(radii[-1] + self.clad_width)
        radii.append(np.sqrt(pitch*pitch / np.pi))
        if radii[-1] >= buffer_radius:
            raise RuntimeError("Buffer radius is smaller than the radius of the cell.")
        radii.append(buffer_radius)

        # Next, we determine all the materials.
        # This requires applying self shielding to the fuel and cladding
        mats = []

        mod_xs = moderator.dilution_xs(moderator.size*[1.E10], ndl)
        mats.append(mod_xs)

        # Add the cladding
        if dancoff_clad is not None:
            Ee = 1. / (2. * (self.outer_radius - self.inner_radius))
            mats.append(self.clad.roman_xs(dancoff_clad, Ee, ndl))
        else:
            mats.append(self.clad.dilution_xs(self.clad.size*[clad_dilution], ndl))

        # Add outer moderator
        mats.append(mod_xs)

        # Add the buffer
        mats.append(buffer)

        return CylindricalCell(radii, mats)

class PWRAssembly:
    def __init__(self, pitch: float, gap: float, moderator: Material, shape: Tuple[int, int]):
        self.pitch = pitch
        self.gap = gap
        self.moderator = moderator
        self.shape = shape
        self.pins = []
        
        # MOC parameters for computing dancoff corrections
        self.dt_dancoff = 0.05
        self.na_dancoff = 64
        self.iso_dancoff = 20.0

    def get_fuel_dancoff_corrections(self):
        # We first make the system for an isolated fuel pin.
        # We isolate it by multiplying the pitch by 20.
        isolated_fp = None
        for pin in self.pins:
            if isinstance(pin, FuelPin):
                isolated_fp = pin.make_fuel_dancoff_cell(pitch=self.iso_dancoff*self.pitch, moderator=self.moderator)
        if isolated_fp is None:
            raise RuntimeError("No FuelPin type found in pins.")
        iso_geom = Cartesian2D([self.iso_dancoff*self.pitch], [self.iso_dancoff*self.pitch])  
        iso_geom.set_tiles([isolated_fp])
        iso_moc = MOCDriver(iso_geom)

        # Set the source
        for i in range(iso_moc.nfsr):
            i_xs = iso_moc.xs(i)
            if i_xs.name != "Fuel":
                # If we aren't in the fuel, set the source to be the value of
                # the potential xs (should be Et).
                iso_moc.set_extern_src(i, 0, i_xs.Et(0))

        # Solve the isolated pin problem
        iso_moc.x_min_bc = BoundaryCondition.Vacuum
        iso_moc.x_max_bc = BoundaryCondition.Vacuum
        iso_moc.y_min_bc = BoundaryCondition.Vacuum
        iso_moc.y_max_bc = BoundaryCondition.Vacuum
        iso_moc.generate_tracks(self.na_dancoff, self.dt_dancoff, YamamotoTabuchi6())
        iso_moc.sim_mode = SimulationMode.FixedSource
        iso_moc.flux_tolerance = 1.E-5
        iso_moc.solve()
        iso_flux = iso_moc.flux(0, 0)

        # Now we setup the lattice problem
        fuel_df_pins = []
        for pin in self.pins:
            fuel_df_pins.append(pin.make_fuel_dancoff_cell(pitch=self.pitch, moderator=self.moderator))
        geom = Cartesian2D(self.shape[0]*[self.pitch], self.shape[1]*[self.pitch]) 
        geom.set_tiles(fuel_df_pins)
        moc = MOCDriver(geom)

        # Set the source
        for i in range(moc.nfsr):
            i_xs = moc.xs(i)
            if i_xs.name != "Fuel":
                # If we aren't in the fuel, set the source to be the value of
                # the potential xs (should be Et).
                moc.set_extern_src(i, 0, i_xs.Et(0))

        # Solve the lattice problem
        moc.generate_tracks(self.na_dancoff, self.dt_dancoff, YamamotoTabuchi6())
        moc.sim_mode = SimulationMode.FixedSource
        moc.flux_tolerance = 1.E-5
        moc.solve()

        # Now we need to calculate the dancoff correction for each pin
        self.fuel_dancoff_corrections = []
        u = Direction(1., 0.)
        for j in range(self.shape[1]):
            y = moc.y_max - self.gap - (j+0.5)*self.pitch
            for i in range(self.shape[0]):
                x = moc.x_min + self.gap + (i+0.5)*self.pitch
                r = Vector(x, y)
                xs = moc.xs(r, u)
                if xs.name == "Fuel":
                    flux = moc.flux(r, u, 0)
                    C = (iso_flux - flux) / iso_flux
                    self.fuel_dancoff_corrections.append(C)
                else:
                    self.fuel_dancoff_corrections.append(0.)

    def get_clad_dancoff_corrections(self):
        pass

