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
    clad: Material
        Material representing the cladding composition, temperature, and
        density.
    clad_width : float
        Thickness of the cladding.
    gap : Material, optional
        Material representing the gap composition, temperature, and density.
    gap_width: float, optional
        Thickness of the gap between the fuel pellet and the cladding
        (if modeled).
    fuel_rings: int
        The number of rings into which the fuel pellet will be divided.
        Default value is 1.
    """

    def __init__(
        self,
        fuel: Material,
        fuel_radius: float,
        clad: Material,
        clad_width: float,
        gap: Optional[Material] = None,
        gap_width: Optional[float] = None,
        fuel_rings: int = 1
    ):
        self.fuel = fuel
        self.fuel_radius = fuel_radius
        self.clad = clad
        self.clad_width = clad_width
        self.gap = gap
        self.gap_width = gap_width
        self.fuel_rings = fuel_rings
        self.condensed_xs = []

        if self.gap is not None and self.gap_width is None:
            raise RuntimeError("Fuel gap material is provided, but no gap width.")
        elif self.gap_width is not None and self.gap is None:
            raise RuntimeError("Fuel gap width provided, but no gap material.")

        if fuel_rings <= 0:
            raise RuntimeError("The number of fuel rings must be >= 1.")

    def clad_offset(self):
        if self.gap_width is not None:
            return Vector(
                self.fuel_radius + self.gap_width + 0.5 * self.clad_width, 0.0
            )
        else:
            return Vector(self.fuel_radius + 0.5 * self.clad_width, 0.0)

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
        Et = np.array([1.0e5])
        Ea = np.array([1.0e5])
        Es = np.array([[0.0]])
        Fuel = CrossSection(Et, Ea, Es, "Fuel")
        mats.append(Fuel)

        if self.gap is not None:
            Et[0] = self.gap.potential_xs
            Ea[0] = self.gap.potential_xs
            gap = CrossSection(Et, Ea, Es, "Gap")
            mats.append(gap)

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
        Es = np.array([[0.0]])
        Fuel = CrossSection(Et, Ea, Es, "Fuel")
        mats.append(Fuel)

        if self.gap is not None:
            Et[0] = self.gap.potential_xs
            Ea[0] = self.gap.potential_xs
            gap = CrossSection(Et, Ea, Es, "Gap")
            mats.append(gap)

        # This returns a cell for calculating the fuel pin Dancoff factor.
        # As such, the clad XS has infinite values.
        Et[0] = 1.0e5
        Ea[0] = 1.0e5
        Clad = CrossSection(Et, Ea, Es, "Clad")
        mats.append(Clad)

        Et[0] = moderator.potential_xs
        Ea[0] = moderator.potential_xs
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
        clad_dilution=1.0e10,
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
                Rin = 0.
                if r > 0:
                    Rin = radii[-1]
                Rout = np.sqrt((Vring + np.pi*Rin*Rin)/np.pi)
                if Rout > self.fuel_radius:
                  Rout = self.fuel_radius
                radii.append(Rout)

        # Add gap radius
        if self.gap is not None:
            radii.append(radii[-1] + self.gap_width)

        # Add clad radius
        radii.append(radii[-1] + self.clad_width)

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
                Rin = 0.
                if r > 0:
                    Rin = radii[r-1]
                Rout = radii[r]
                mats.append(self.fuel.ring_carlvik_xs(dancoff_fuel, self.fuel_radius, Rin, Rout, ndl))
                mats[-1].name = "Fuel"

        # Next, add the gap (if present)
        if self.gap is not None:
            mats.append(self.gap.dilution_xs(self.gap.size * [1.0e10], ndl))
            mats[-1].name = "Gap"

        # Add the cladding
        if dancoff_clad is not None:
            Ee = 1.0 / (2.0 * self.clad_width)
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
                Rin = 0.
                if r > 0:
                    Rin = radii[-1]
                Rout = np.sqrt((Vring + np.pi*Rin*Rin)/np.pi)
                radii.append(Rout)

        if self.gap is not None:
            radii.append(radii[-1] + self.gap_width)

        radii.append(radii[-1] + self.clad_width)

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
        Ea = np.array([moderator.potential_xs])
        Es = np.array([[0.0]])
        Mod = CrossSection(Et, Ea, Es, "Moderator")
        mats.append(Mod)

        # This returns a cell for calculating the fuel pin Dancoff factor.
        # As such, the clad XS has infinite values.
        Et[0] = 1.0e5
        Ea[0] = 1.0e5
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
        Ea = np.array([moderator.potential_xs])
        Es = np.array([[0.0]])
        Mod = CrossSection(Et, Ea, Es, "Moderator")
        mats.append(Mod)

        Et[0] = self.clad.potential_xs
        Ea[0] = self.clad.potential_xs
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
        clad_dilution: float = 1.0e10,
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

