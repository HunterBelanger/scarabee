from _scarabee import *
import numpy as np
from typing import Optional


class FuelPin:
    """
    Represents a fuel pin, contianing its geometric and material properties.

    Parameters
    ----------
    fuel : Material
        Material representing the fuel composition, temperature, and density.
    fuel_radius : float
        Radius of the fuel pellet.
    gap : Material, optional
        Material representing the gap composition, temperature, and density.
    gap_radius: float, optional
        Outer radius of the gap between the fuel pellet and the cladding
        (if modeled).
    clad: Material
        Material representing the cladding composition, temperature, and
        density.
    clad_radius : float
        Outer radius of the cladding.
    fuel_rings: int
        The number of rings into which the fuel pellet will be divided.
        Default value is 1.
    """

    def __init__(
        self,
        fuel: Material,
        fuel_radius: float,
        gap: Optional[Material],
        gap_radius: Optional[float],
        clad: Material,
        clad_radius: float,
        fuel_rings: int = 1,
    ):
        self.fuel = fuel
        self.fuel_radius = fuel_radius
        self.clad = clad
        self.clad_radius = clad_radius
        self.gap = gap
        self.gap_radius = gap_radius
        self.fuel_rings = fuel_rings
        self.condensed_xs = []

        if self.gap is not None and self.gap_radius is None:
            raise RuntimeError("Fuel gap material is provided, but no gap radius.")
        elif self.gap_radius is not None and self.gap is None:
            raise RuntimeError("Fuel gap radius provided, but no gap material.")

        if fuel_rings <= 0:
            raise RuntimeError("The number of fuel rings must be >= 1.")

    def clad_offset(self):
        if self.gap_radius is not None:
            return Vector(self.gap_radius + 0.5 * (self.clad_radius - self.gap_radius), 0.0)
        else:
            return Vector(self.fuel_radius + 0.5 * (self.clad_radius - self.fuel_radius), 0.0)

    def make_fuel_dancoff_cell(self, pitch: float, moderator: Material):
        # We first determine all the radii
        radii = []
        radii.append(self.fuel_radius)
        if self.gap is not None:
            radii.append(self.gap_radius)
        radii.append(self.clad_radius)

        mats = []

        # This returns a cell for calculating the fuel pin Dancoff factor.
        # As such, the fuel XS has infinite values.
        Et = np.array([1.0e5])
        Ea = np.array([Et[0]])
        Es = np.array([[0.0]])
        Fuel = CrossSection(Et, Ea, Es, "Fuel")
        mats.append(Fuel)

        if self.gap is not None:
            Et[0] = self.gap.potential_xs
            Ea[0] = Et[0]
            gap = CrossSection(Et, Ea, Es, "Gap")
            mats.append(gap)

        Et[0] = self.clad.potential_xs
        Ea[0] = Et[0]
        Clad = CrossSection(Et, Ea, Es, "Clad")
        mats.append(Clad)

        Et[0] = moderator.potential_xs
        Ea[0] = Et[0]
        Mod = CrossSection(Et, Ea, Es, "Moderator")
        mats.append(Mod)

        return SimplePinCell(radii, mats, pitch, pitch)

    def make_clad_dancoff_cell(self, pitch: float, moderator: Material):
        # We first determine all the radii
        radii = []
        radii.append(self.fuel_radius)
        if self.gap is not None:
            radii.append(self.gap_radius)
        radii.append(self.clad_radius)

        mats = []

        Et = np.array([self.fuel.potential_xs])
        Ea = np.array([Et[0]])
        Es = np.array([[0.0]])
        Fuel = CrossSection(Et, Ea, Es, "Fuel")
        mats.append(Fuel)

        if self.gap is not None:
            Et[0] = self.gap.potential_xs
            Ea[0] = Et[0]
            gap = CrossSection(Et, Ea, Es, "Gap")
            mats.append(gap)

        # This returns a cell for calculating the fuel pin Dancoff factor.
        # As such, the clad XS has infinite values.
        Et[0] = 1.0e5
        Ea[0] = Et[0]
        Clad = CrossSection(Et, Ea, Es, "Clad")
        mats.append(Clad)

        Et[0] = moderator.potential_xs
        Ea[0] = Et[0]
        Mod = CrossSection(Et, Ea, Es, "Moderator")
        mats.append(Mod)

        return SimplePinCell(radii, mats, pitch, pitch)

    def make_cylindrical_cell(
        self,
        pitch: float,
        dancoff_fuel: float,
        moderator: CrossSection,
        ndl: NDLibrary,
        dancoff_clad: Optional[float] = None,
        clad_dilution=300.,
    ):
        # We first determine all the radii
        radii = []

        # Add the fuel radius
        if self.fuel_rings == 1:
            radii.append(self.fuel_radius)
        else:
            # We need to subdivide the pellet into rings. We start be getting
            # the fuel pellet volume
            V = np.pi * self.fuel_radius * self.fuel_radius
            Vring = V / float(self.fuel_rings)
            for r in range(self.fuel_rings):
                Rin = 0.0
                if r > 0:
                    Rin = radii[-1]
                Rout = np.sqrt((Vring + np.pi * Rin * Rin) / np.pi)
                if Rout > self.fuel_radius:
                    Rout = self.fuel_radius
                radii.append(Rout)

        # Add gap radius
        if self.gap is not None:
            radii.append(self.gap_radius)

        # Add clad radius
        radii.append(self.clad_radius)

        # Add water radius
        radii.append(np.sqrt(pitch * pitch / np.pi))

        # Next, we determine all the materials.
        # This requires applying self shielding to the fuel and cladding
        mats = []

        # First, treat the fuel
        if self.fuel_rings == 1:
            # For a single pellet, use the standard Carlvik two-term
            Ee = 1.0 / (2.0 * self.fuel_radius)  # Fuel escape xs
            mats.append(self.fuel.carlvik_xs(dancoff_fuel, Ee, ndl))
            mats[-1].name = "Fuel"
        else:
            # We need to apply spatial self shielding.
            for r in range(self.fuel_rings):
                Rin = 0.0
                if r > 0:
                    Rin = radii[r - 1]
                Rout = radii[r]
                mats.append(
                    self.fuel.ring_carlvik_xs(
                        dancoff_fuel, self.fuel_radius, Rin, Rout, ndl
                    )
                )
                mats[-1].name = "Fuel"

        # Next, add the gap (if present)
        if self.gap is not None:
            mats.append(self.gap.dilution_xs(self.gap.size * [1.0e10], ndl))
            mats[-1].name = "Gap"

        # Add the cladding
        if dancoff_clad is not None:
            if self.gap_radius is not None:
                Ee = 1.0 / (2.0 * (self.clad_radius - self.gap_radius))
            else:
                Ee = 1.0 / (2.0 * (self.clad_radius - self.fuel_radius))
            mats.append(self.clad.roman_xs(dancoff_clad, Ee, ndl))
        else:
            mats.append(self.clad.dilution_xs(self.clad.size * [clad_dilution], ndl))
        mats[-1].name = "Clad"

        # Finally, add moderator
        mats.append(moderator)

        return CylindricalCell(radii, mats)

    def make_moc_cell(self, pitch: float):
        radii = []

        if self.fuel_rings == 1:
            radii.append(self.fuel_radius)
        else:
            # We need to subdivide the pellet into rings. We start be getting
            # the fuel pellet volume
            V = np.pi * self.fuel_radius * self.fuel_radius
            Vring = V / float(self.fuel_rings)
            for r in range(self.fuel_rings):
                Rin = 0.0
                if r > 0:
                    Rin = radii[-1]
                Rout = np.sqrt((Vring + np.pi * Rin * Rin) / np.pi)
                radii.append(Rout)

        if self.gap is not None:
            radii.append(self.gap_radius)

        radii.append(self.clad_radius)

        mod_width = 0.5 * pitch - radii[-1]
        radii.append(radii[-1] + 0.8 * mod_width)

        mats = self.condensed_xs.copy()
        mats.append(self.condensed_xs[-1])

        return PinCell(radii, mats, pitch, pitch)


class GuideTube:
    """
    Represents an empty guide tube, contianing its geometric and material
    properties.

    Parameters
    ----------
    inner_radius : float
        Inside radius of the guide tube.
    outer_radius : float
        Outside radius of the guide tube.
    clad : Material
        Material representing the guide tube composition, temperature, and
        density (assumed to be the same material as the fuel pin cladding).
    """

    def __init__(self, inner_radius: float, outer_radius: float, clad: Material):
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.clad = clad
        self.condensed_xs = []

        if self.outer_radius <= self.inner_radius:
            raise RuntimeError("Outer radius must be > inner radius.")

    def clad_offset(self):
        return Vector(0.5 * (self.inner_radius + self.outer_radius), 0.0)

    def make_clad_dancoff_cell(self, pitch: float, moderator: Material):
        # We first determine all the radii
        radii = []
        radii.append(self.inner_radius)
        radii.append(self.outer_radius)

        mats = []

        Et = np.array([moderator.potential_xs])
        Ea = np.array([Et[0]])
        Es = np.array([[0.0]])
        Mod = CrossSection(Et, Ea, Es, "Moderator")
        mats.append(Mod)

        # This returns a cell for calculating the fuel pin Dancoff factor.
        # As such, the clad XS has infinite values.
        Et[0] = 1.0e5
        Ea[0] = Et[0]
        Clad = CrossSection(Et, Ea, Es, "Clad")
        mats.append(Clad)

        mats.append(Mod)

        return SimplePinCell(radii, mats, pitch, pitch)

    def make_fuel_dancoff_cell(self, pitch: float, moderator: Material):
        # We first determine all the radii
        radii = []
        radii.append(self.inner_radius)
        radii.append(self.outer_radius)

        mats = []

        Et = np.array([moderator.potential_xs])
        Ea = np.array([Et[0]])
        Es = np.array([[0.0]])
        Mod = CrossSection(Et, Ea, Es, "Moderator")
        mats.append(Mod)

        Et[0] = self.clad.potential_xs
        Ea[0] = Et[0]
        Clad = CrossSection(Et, Ea, Es, "Clad")
        mats.append(Clad)

        mats.append(Mod)

        return SimplePinCell(radii, mats, pitch, pitch)

    def make_cylindrical_cell(
        self,
        pitch: float,
        moderator: CrossSection,
        buffer_radius: float,
        buffer: CrossSection,
        ndl: NDLibrary,
        dancoff_clad: Optional[float] = None,
        clad_dilution: float = 300.,
    ):
        # We first determine all the radii
        radii = []
        radii.append(self.inner_radius)
        radii.append(self.outer_radius)
        radii.append(np.sqrt(pitch * pitch / np.pi))
        if radii[-1] >= buffer_radius:
            raise RuntimeError("Buffer radius is smaller than the radius of the cell.")
        radii.append(buffer_radius)

        # Next, we determine all the materials.
        # This requires applying self shielding to the fuel and cladding
        mats = []

        mats.append(moderator)

        # Add the cladding
        if dancoff_clad is not None:
            Ee = 1.0 / (2.0 * (self.outer_radius - self.inner_radius))
            mats.append(self.clad.roman_xs(dancoff_clad, Ee, ndl))
        else:
            mats.append(self.clad.dilution_xs(self.clad.size * [clad_dilution], ndl))
        mats[-1].name = "Clad"

        # Add outer moderator
        mats.append(moderator)

        # Add the buffer
        mats.append(buffer)

        return CylindricalCell(radii, mats)

    def make_moc_cell(self, pitch: float):
        r_inner_inner_mod = np.sqrt(0.5 * self.inner_radius * self.inner_radius)
        radii = [r_inner_inner_mod, self.inner_radius, self.outer_radius]
        mats = [self.condensed_xs[0]] + self.condensed_xs.copy()
        return PinCell(radii, mats, pitch, pitch)


class BurnablePoisonPin:
    """
    Represents a burnable poison bin which could be of borosilicate glass
    (BSG), or a wet annular burnable absorber (WABA).

    Parameters
    ----------
    center: Material
        Material at the center of the burnable poison tube.
    center_radius: float
        Radius of the center material in the burnable poison tube.
    poison_clad: Material
        Cladding material for the burnable poison tube.
    inner_poison_clad_radius: float
        Radius of the inner cladding of the burnable poison tube.
    gap: Material, optional
        Material of the gap between the burnable poison and the cladding
        (if present).
    inner_gap_radius: float, optional
        Radius of the inner gap between poison cladding and the poison.
    poison: Material
        Material of the burnable poison.
    poison_radius: float
        Outer radius of the burnable poison.
    outer_gap_radius: float, optional
        Outer radius of the gap between the burnable poison and cladding.
    outer_poison_clad_radius: float
        Outer radius of the cladding of the burnable poison tube.
    inner_moderator_radius: float
        The outer radius of the moderator between the burnable poison tube
        and the guide tube.
    guide_tube_clad: Material
        Material for the guide tube.
    guide_tube_radius: float
        Outer radius of the guide tube.
    """

    def __init__(
        self,
        center: Material,
        center_radius: float,
        poison_clad: Material,
        inner_poison_clad_radius: float,
        gap: Optional[Material],
        inner_gap_radius: Optional[float],
        poison: Material,
        poison_radius: float,
        outer_gap_radius: Optional[float],
        outer_poison_clad_radius: float,
        inner_moderator_radius: float, 
        guide_tube_clad: Material,
        guide_tube_radius: float
    ):
        self.center = center
        self.center_radius = center_radius

        self.poison_clad = poison_clad
        self.inner_poison_clad_radius = inner_poison_clad_radius
        
        self.gap = gap
        self.inner_gap_radius = inner_gap_radius

        self.poison = poison
        self.poison_radius = poison_radius

        self.outer_gap_radius = outer_gap_radius

        self.outer_poison_clad_radius = outer_poison_clad_radius

        self.inner_moderator_radius = inner_moderator_radius

        self.guide_tube_clad = guide_tube_clad
        self.guide_tube_radius = guide_tube_radius

        self.condensed_xs = []

        gap_nones = 0
        if self.gap is None:
            gap_nones += 1
        if self.inner_gap_radius is None:
            gap_nones += 1
        if self.outer_gap_radius is None:
            gap_nones += 1
        if gap_nones != 0 and gap_nones != 3:
            raise RuntimeError("Must provide all three gap parameters.")

    def clad_offset(self):
        return Vector(self.inner_moderator_radius + 0.5*(self.guide_tube_radius - self.inner_moderator_radius), 0.0)

    def make_fuel_dancoff_cell(self, pitch: float, moderator: Material):
        radii = []
        mats = []

        Et = np.array([self.center.potential_xs])
        Ea = np.array([Et[0]])
        Es = np.array([[0.0]])
        radii.append(self.center_radius)
        mats.append(CrossSection(Et, Ea, Es, "BP Center"))

        Et[0] = self.poison_clad.potential_xs
        Ea[0] = Et[0]
        radii.append(self.inner_poison_clad_radius)
        PC = CrossSection(Et, Ea, Es, "Poison Clad") 
        mats.append(PC)

        if self.gap:
            Et[0] = self.gap.potential_xs
            Ea[0] = Et[0]
            radii.append(self.inner_gap_radius)
            Gap = CrossSection(Et, Ea, Es, "Gap")
            mats.append(Gap)

        Et[0] = self.poison.potential_xs
        Ea[0] = Et[0]
        radii.append(self.poison_radius)
        mats.append(CrossSection(Et, Ea, Es, "Poison"))
        
        if self.gap:
            Et[0] = self.gap.potential_xs
            Ea[0] = Et[0]
            radii.append(self.outer_gap_radius)
            mats.append(Gap)

        Et[0] = self.poison_clad.potential_xs
        Ea[0] = Et[0]
        radii.append(self.outer_poison_clad_radius)
        mats.append(PC)

        Et[0] = moderator.potential_xs
        Ea[0] = Et[0]
        radii.append(self.inner_moderator_radius)
        Mod = CrossSection(Et, Ea, Es, "Moderator")
        mats.append(Mod)

        Et[0] = self.guide_tube_clad.potential_xs
        Ea[0] = Et[0]
        radii.append(self.guide_tube_radius)
        mats.append(CrossSection(Et, Ea, Es, "Clad"))

        mats.append(Mod) 

        return SimplePinCell(radii, mats, pitch, pitch)

    def make_clad_dancoff_cell(self, pitch: float, moderator: Material):
        radii = []
        mats = []

        Et = np.array([self.center.potential_xs])
        Ea = np.array([Et[0]])
        Es = np.array([[0.0]])
        radii.append(self.center_radius)
        mats.append(CrossSection(Et, Ea, Es, "BP Center"))

        Et[0] = self.poison_clad.potential_xs
        Ea[0] = Et[0]
        radii.append(self.inner_poison_clad_radius)
        PC = CrossSection(Et, Ea, Es, "Poison Clad") 
        mats.append(PC)

        if self.gap:
            Et[0] = self.gap.potential_xs
            Ea[0] = Et[0]
            radii.append(self.inner_gap_radius)
            Gap = CrossSection(Et, Ea, Es, "Gap")
            mats.append(Gap)

        Et[0] = self.poison.potential_xs
        Ea[0] = Et[0]
        radii.append(self.poison_radius)
        mats.append(CrossSection(Et, Ea, Es, "Poison"))
        
        if self.gap:
            Et[0] = self.gap.potential_xs
            Ea[0] = Et[0]
            radii.append(self.outer_gap_radius)
            mats.append(Gap)

        Et[0] = self.poison_clad.potential_xs
        Ea[0] = Et[0]
        radii.append(self.outer_poison_clad_radius)
        mats.append(PC)

        Et[0] = moderator.potential_xs
        Ea[0] = Et[0]
        radii.append(self.inner_moderator_radius)
        Mod = CrossSection(Et, Ea, Es, "Moderator")
        mats.append(Mod)

        Et[0] = 1.0E5
        Ea[0] = Et[0]
        radii.append(self.guide_tube_radius)
        mats.append(CrossSection(Et, Ea, Es, "Clad"))

        mats.append(Mod) 

        return SimplePinCell(radii, mats, pitch, pitch)

    def make_cylindrical_cell(
        self,
        pitch: float,
        moderator: CrossSection,
        buffer_radius: float,
        buffer: CrossSection,
        ndl: NDLibrary,
        dancoff_clad: Optional[float] = None,
        clad_dilution: float = 300,
        poison_clad_dilution: float = 300
    ):
        radii = []
        mats = []

        radii.append(self.center_radius)
        mats.append(self.center.dilution_xs(self.center.size * [1.0E10], ndl))
        mats[-1].name = "Center"

        radii.append(self.inner_poison_clad_radius)
        PC = self.poison_clad.dilution_xs(self.poison_clad.size * [poison_clad_dilution], ndl) 
        mats.append(PC)
        mats[-1].name = "Poison Clad"

        if self.gap:
            radii.append(self.inner_gap_radius)
            Gap = self.gap.dilution_xs(self.gap.size * [1.0E10], ndl)
            mats.append(Gap)
            mats[-1].name = "Gap"

        radii.append(self.poison_radius)
        mats.append(self.poison.dilution_xs(self.poison.size * [1.0E10], ndl))
        mats[-1].name = "Poison"
        
        if self.gap:
            radii.append(self.outer_gap_radius)
            mats.append(Gap)

        radii.append(self.outer_poison_clad_radius)
        mats.append(PC)

        radii.append(self.inner_moderator_radius)
        mats.append(moderator)

        radii.append(self.guide_tube_radius)
        if dancoff_clad is not None:
            Ee = 1.0 / (2.0 * (self.guide_tube_radius - self.inner_moderator_radius))
            mats.append(self.guide_tube_clad.roman_xs(dancoff_clad, Ee, ndl))
        else:
            mats.append(self.guide_tube_clad.dilution_xs(self.guide_tube_clad.size * [clad_dilution], ndl))
        mats[-1].name = "Clad"
        
        radii.append(np.sqrt(pitch * pitch / np.pi))
        mats.append(moderator) 

        # Add the buffer
        radii.append(buffer_radius)
        mats.append(buffer)

        return CylindricalCell(radii, mats)

    def make_moc_cell(self, pitch: float):
        radii = []
        radii.append(self.center_radius)
        radii.append(self.inner_poison_clad_radius)
        if self.gap:
            radii.append(self.inner_gap_radius)
        radii.append(self.poison_radius)
        if self.gap:
            radii.append(self.outer_gap_radius)
        radii.append(self.outer_poison_clad_radius)
        radii.append(self.inner_moderator_radius)
        radii.append(self.guide_tube_radius)

        return PinCell(radii, self.condensed_xs, pitch, pitch)
