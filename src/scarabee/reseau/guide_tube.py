from .._scarabee import (
    NDLibrary,
    Material,
    CrossSection,
    SimplePinCell,
    PinCell,
)
import numpy as np
from typing import Optional, List, Tuple
import copy


class GuideTube:
    """
    Represents an empty guide tube for a PWR.

    Parameters
    ----------
    clad : Material
        Material which describes the cladding composition, density, and
        temperature.
    inner_radius : float
        Inner radius of the guide tube.
    outer_radius : float
        Outer radius of the guide tube.

    Attributes
    ----------
    clad : Material
        Material which describes the cladding composition, density, and
        temperature.
    inner_radius : float
        Inner radius of the guide tube.
    outer_radius : float
        Outer radius of the guide tube.
    clad_dancoff_factors : list of float
        Dancoff factors to be used when self-shielding the cladding at each
        depletion time step.
    """
    def __init__(self, clad: Material, inner_radius: float, outer_radius: float):
        if inner_radius >= outer_radius:
            raise ValueError("Inner radius must be > outer radius.")

        self._clad = copy.deepcopy(clad)
        self._inner_radius = inner_radius
        self._outer_radius = outer_radius

        # Initialize empty list of Dancoff factors for the cladding
        self._clad_dancoff_factors: List[float] = []
        
        # Initialize empty variables for Dancoff factor calculations.
        # These are all kept private.
        self._clad_dancoff_xs: CrossSection = CrossSection(
            np.array([self.clad.potential_xs]),
            np.array([self.clad.potential_xs]),
            np.array([[0.0]]),
            "Clad",
        )

        # Holds the CrossSection object used for the real transport calc.
        self._clad_xs: Optional[CrossSection] = None

    @property
    def clad(self):
        return self._clad

    @property
    def inner_radius(self):
        return self._inner_radius

    @property
    def outer_radius(self):
        return self._outer_radius

    @property
    def clad_dancoff_factors(self):
        return self._clad_dancoff_factors

    def load_nuclides(self, ndl: NDLibrary) -> None:
        """
        Loads all the nuclides for all current materials into the data library.

        Parameters
        ----------
        ndl : NDLibrary
            Nuclear data library which should load the nuclides.
        """
        self.clad.load_nuclides(ndl)

    def set_xs_for_fuel_dancoff_calculation(self) -> None:
        """
        Sets the 1-group cross sections to calculate the fuel Dancoff factors.
        """
        self._clad_dancoff_xs.set(
            CrossSection(
                np.array([self.clad.potential_xs]),
                np.array([self.clad.potential_xs]),
                np.array([[0.0]]),
                "Clad",
            )
        )

    def set_xs_for_clad_dancoff_calculation(self) -> None:
        """
        Sets the 1-group cross sections to calculate the clad Dancoff factors.
        """
        self._clad_dancoff_xs.set(
            CrossSection(
                np.array([1.0e5]),
                np.array([1.0e5]),
                np.array([[0.0]]),
                "Clad",
            )
        )

    def make_dancoff_moc_cell(
        self, moderator_xs: CrossSection, pitch: float
    ) -> Tuple[SimplePinCell, None, int]:
        """
        Makes a simplified cell suitable for performing Dancoff factor
        calculations.

        Parameters
        ----------
        moderator_xs : CrossSection
            One group cross sections for the moderator. Total should equal
            absorption (i.e. no scattering) and should be equal to the
            macroscopic potential cross section.
        pitch : float
            Width of the cell.

        Returns
        -------
        SimplifiedPinCell, None, int
            Pin cell object for MOC calculation, None, ID of the clad flat
            source region.
        """
        if pitch < 2.0 * self.clad_radius:
            raise ValueError(
                "The fuel pin pitch must be > the diameter of the cladding."
            )

        # First we create list of radii and materials
        radii = []
        xs = []

        radii.append(self.inner_radius) 
        xs.append(moderator_xs)

        radii.append(self.outer_radius)
        xs.append(self._clad_dancoff_xs)

        xs.append(moderator_xs)

        # Make the simple pin cell.
        cell = SimplePinCell(radii, xs, pitch, pitch)

        # Get the FSR IDs for the regions of interest
        cell_fsr_ids = list(cell.get_all_fsr_ids())
        cell_fsr_ids.sort()

        clad_fsr_id = cell_fsr_ids[1]

        return cell, None, clad_fsr_id

    def append_clad_dancoff_factor(self, D) -> None:
        """
        Saves new Dancoff factor for the cladding that will be used for all
        subsequent cross section updates.

        Parameters
        ----------
        D : float
            New Dancoff factor.
        """
        if D < 0.0 or D > 1.0:
            raise ValueError("Dancoff factor must be in range [0, 1].")
        self._clad_dancoff_factors.append(D)

    def set_clad_xs_for_depletion_step(self, t: int, ndl: NDLibrary) -> None:
        """
        Constructs the CrossSection object for the cladding of the guide tube
        at the specified depletion step. The depletion step only changes the
        Dancoff factor, not the cladding composition.

        Parameters
        ----------
        t : int
            Index for the depletion step.
        ndl : NDLibrary
            Nuclear data library to use for cross sections.
        """
        # Compute escape xs
        Ee = 1. / (2. * (self.outer_radius - self.inner_radius))

        # Get / set the xs
        if self._clad_xs is None:
            self._clad_xs = self.clad.roman_xs(self._clad_dancoff_factors[t], Ee, ndl)
        else:
            self._clad_xs.set(self.clad.roman_xs(self._clad_dancoff_factors[t], Ee, ndl))

        if self._clad_xs.name == "":
            self._clad_xs.set_name("Clad")

    def make_moc_cell(self, moderator_xs: CrossSection, pitch: float) -> PinCell:
        """
        Constructs the pin cell object used in for the global MOC simulation.

        Parameters
        ----------
        moderator_xs : CrossSection
            Cross sections to use for the moderator surrounding the fuel pin.
        pitch : float
            Spacing between fuel pins. Must be larger than the outer diameter
            of the cladding.
        """
        if self._clad_xs is None:
            raise RuntimeError("Clad cross section has not yet been built.")
        if pitch < 2.*self.outer_radius:
            raise RuntimeError("Pitch must be >= twice the outer cladding radius.")

        # Create list of inner moderator radii
        radii = []
        V = np.pi * self.inner_radius * self.inner_radius
        Vr = V / 3.
        for ri in range(3):
            Rin = 0.
            if ri > 0:
                Rin = radii[-1]
            Rout = np.sqrt((Vr + np.pi * Rin * Rin) / np.pi)
            if Rout > self.inner_radius:
                Rout = self.inner_radius
            radii.append(Rout)

        # Initialize the cross section lists with the inner moderator xs
        xss = [moderator_xs] * len(radii)

        # Add cladding
        radii.append(self.outer_radius)
        xss.append(self._clad_xs)

        # Add another ring of moderator if possible
        if pitch > 2.*self.outer_radius:
            radii.append(0.5*pitch)
            xss.append(moderator_xs)

        # Add moderator to the end of materials
        xss.append(moderator_xs)

        return PinCell(radii, xss, pitch, pitch)        