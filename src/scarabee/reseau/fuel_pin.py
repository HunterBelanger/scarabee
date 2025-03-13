from .._scarabee import NDLibrary, Material, CrossSection, SimplePinCell, MOCDriver
import numpy as np
from typing import Optional, List, Tuple
import copy


class FuelPin:
    """
    Represents a generic fuel pin. This could be a fuel pin in a PWR, or a fuel
    pin in a BWR next to the channel box or in the channel corner.

    Parameters
    ----------
    fuel : Material
        Material which describes the fuel composition, density, and temperature.
    fuel_radius : float
        Outer radius of the fuel pellet region.
    gap : Material, optional
        Material which describes the composition, density, and temperature of
        the gap between the fuel pellet and the cladding, if present.
    gap_radius : float, optional
        Outer radius of the gap material, if present.
    clad : Material
        Material which describes the cladding composition, density, and
        temperature.
    clad_radius : float
        Outer radius of the cladding.
    num_fuel_rings : int, default 1
        Number of rings which should be used to discretize the fuel material.
        Each ring will be self-shielded and depleted separately.

    Attributes
    ----------
    fuel_radius : float
        Outer radius of the fuel pellet region.
    fuel_ring_materials : list of list of Material
        Contains the Material in each fuel ring for each depletion time step.
    fuel_ring_flux_spectra : list of list of ndarray
        Contains the average flux spectrum in each fuel ring for each depletion
        time step.
    fuel_dancoff_factors : list of float
        Dancoff factors to be used when self-shielding the fuel at each
        depletion time step.
    gap : Material, optional
        Material which describes the composition, density, and temperature of
        the gap between the fuel pellet and the cladding, if present.
    gap_radius : float, optional
        Outer radius of the gap material, if present.
    clad : Material
        Material which describes the cladding composition, density, and
        temperature.
    clad_radius : float
        Outer radius of the cladding.
    clad_dancoff_factors : list of float
        Dancoff factors to be used when self-shielding the cladding at each
        depletion time step.
    num_fuel_rings : int, default 1
        Number of rings which should be used to discretize the fuel material.
        Each ring will be self-shielded and depleted separately.
    """

    def __init__(
        self,
        fuel: Material,
        fuel_radius: float,
        clad: Material,
        clad_radius,
        gap: Optional[Material],
        gap_radius: Optional[float],
        num_fuel_rings: int = 1,
    ):
        if fuel_radius <= 0.0:
            raise ValueError("Fuel radius must be > 0.")
        self._fuel_radius = fuel_radius

        if num_fuel_rings <= 0:
            raise ValueError("Number of fuel rings must be >= 1.")
        self._num_fuel_rings = num_fuel_rings

        # Initialize array of compositions for the fuel. This holds the
        # composition for each fuel ring and for each depletion step.
        self._fuel_ring_materials: List[List[Material]] = []
        for r in range(self.num_fuel_rings):
            # All rings initially start with the same composition
            self._fuel_ring_materials.append([copy.deepcopy(fuel)])

        # Initialize an array to hold the flux spectrum for each fuel ring and
        # for each depletion step.
        self._fuel_ring_flux_spectra: List[List[np.ndarray]] = []
        for r in range(self.num_fuel_rings):
            # All rings initially start with empty flux spectrum list
            self._fuel_ring_flux_spectra.append([])

        # Initialize empty list of Dancoff factors for the fuel
        self._fuel_dancoff_factors: List[float] = []

        # Get gap related parameters
        if gap is None and gap_radius is not None:
            raise ValueError("Gap material is None but gap radius is defined.")
        elif gap is not None and gap_radius is None:
            raise ValueError("Gap material is defined but gap radius is None.")

        if gap_radius is not None and gap_radius <= fuel_radius:
            raise ValueError("Gap radius must be > fuel radius.")

        self._gap = gap
        self._gap_radius = gap_radius

        # Get cladding related parameters
        self._clad = clad

        if clad_radius <= 0.0 or clad_radius <= fuel_radius:
            raise ValueError("Clad radius must be > fuel radius.")
        elif gap_radius is not None and clad_radius <= gap_radius:
            raise ValueError("Clad radius must be > gap radius.")
        self._clad_radius = clad_radius

        # Initialize empty list of Dancoff factors for the cladding
        self._clad_dancoff_factors: List[float] = []

        # Initialize empty variables for Dancoff factor calculations.
        # These are all kept private.
        self._fuel_dancoff_xs: CrossSection = CrossSection(
            np.array([1.0e5]), np.array([1.0e5]), np.array([[0.0]]), "Fuel"
        )
        self._gap_dancoff_xs: Optional[CrossSection] = None
        if self.gap is not None:
            self._gap_dancoff_xs = CrossSection(
                np.array([self.gap.potential_xs]),
                np.array([self.gap.potential_xs]),
                np.array([[0.0]]),
                "Gap",
            )
        self._clad_dancoff_xs: CrossSection = CrossSection(
            np.array([self.clad.potential_xs]),
            np.array([self.clad.potential_xs]),
            np.array([[0.0]]),
            "Clad",
        )

        # Lists of the FSR IDs for each fuel ring, used to homogenize flux
        # spectra for depletion. These will be filled by make_moc_cell.
        self._fuel_ring_fsr_ids: List[List[int]] = []
        for r in range(self.num_fuel_rings):
            self._fuel_ring_fsr_ids.append([])

    @property
    def fuel_radius(self) -> float:
        return self._fuel_radius

    @property
    def num_fuel_rings(self) -> int:
        return self._num_fuel_rings

    @property
    def fuel_ring_materials(self) -> List[List[Material]]:
        return self._fuel_ring_materials

    @property
    def fuel_ring_flux_spectra(self) -> List[List[Material]]:
        return self._fuel_ring_materials

    @property
    def fuel_dancoff_factors(self) -> List[float]:
        return self._fuel_dancoff_factors

    @property
    def gap(self) -> Optional[Material]:
        return self._gap

    @property
    def gap_radius(self) -> Optional[float]:
        return self._gap_radius

    @property
    def clad(self) -> Material:
        return self._clad

    @property
    def clad_radius(self) -> float:
        return self._clad_radius

    @property
    def clad_dancoff_factors(self) -> List[float]:
        return self._clad_dancoff_factors

    def load_nuclides(self, ndl: NDLibrary) -> None:
        """
        Loads all the nuclides for all current materials into the data library.

        Parameters
        ----------
        ndl : NDLibrary
            Nuclear data library which should load the nuclides.
        """
        self.fuel.load_nuclides(ndl)
        if self.gap is not None:
            self.gap.load_nuclides(ndl)
        self.clad.load_nuclides(ndl)

    def setup_xs_for_fuel_dancoff_calculation(self) -> None:
        pass
    
    def setup_xs_for_clad_dancoff_calculation(self) -> None:
        pass

    def make_dancoff_moc_cell(
        self, moderator_xs: CrossSection, pitch: float
    ) -> Tuple[SimplePinCell, int, int]:
        """
        Makes a simplified cell suitable for perofrming Dancoff factor
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
        SimplifiedPinCell, int, int
            Pin cell object for MOC calcualtion, ID of the fuel flat source
            region, ID of the clad flat source region.
        """
        if pitch < 2.0 * self.clad_radius:
            raise ValueError(
                "The fuel pin pitch must be > the diameter of the cladding."
            )

        # First we create list of radii and materials
        radii = []
        xs = []

        radii.append(self.fuel_radius)
        xs.append(self._fuel_dancoff_xs)

        if self.gap is not None and self.gap_radius is not None:
            radii.append(self.gap_radius)
            xs.append(self._gap_dancoff_xs)

        radii.append(self.clad_radius)
        xs.append(self._clad_dancoff_xs)

        xs.append(moderator_xs)

        # Make the simple pin cell.
        cell = SimplePinCell(radii, xs, pitch, pitch)

        # Get the FSR IDs for the regions of interest
        cell_fsr_ids = list(cell.get_all_fsr_ids())
        cell_fsr_ids.sort()

        fuel_fsr_id = cell_fsr_ids[0]

        if self.gap is None:
            clad_fsr_id = cell_fsr_ids[1]
        else:
            clad_fsr_id = cell_fsr_ids[2]

        return cell, fuel_fsr_id, clad_fsr_id

    def append_fuel_dancoff_factor(self, D) -> None:
        """
        Saves new Dancoff factor for the fuel that will be used for all
        subsequent cross section updates.

        Parameters
        ----------
        D : float
            New Dancoff factor.
        """
        if D < 0.0 or D > 1.0:
            raise ValueError("Dancoff factor must be in range [0, 1].")
        self._fuel_dancoff_factors.append(D)

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

    def make_moc_cell(self, moderator_xs: CrossSection, pitch: float):
        pass
