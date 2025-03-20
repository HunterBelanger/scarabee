from .._scarabee import (
    NDLibrary,
    Material,
    CrossSection,
    PinCellType,
    SimplePinCell,
    PinCell,
    MOCDriver,
    MixingFraction,
    mix_materials,
)
import numpy as np
from typing import Optional, List
import copy


class FuelPin:
    """
    Represents a generic fuel pin for a PWR.

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
    fuel_dancoff_corrections : list of float
        Dancoff corrections to be used when self-shielding the fuel at each
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
    clad_dancoff_corrections : list of float
        Dancoff corrections to be used when self-shielding the cladding at each
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

        # Get gap related parameters
        if gap is None and gap_radius is not None:
            raise ValueError("Gap material is None but gap radius is defined.")
        elif gap is not None and gap_radius is None:
            raise ValueError("Gap material is defined but gap radius is None.")

        if gap_radius is not None and gap_radius <= fuel_radius:
            raise ValueError("Gap radius must be > fuel radius.")

        self._gap = copy.deepcopy(gap)
        self._gap_radius = gap_radius

        # Get cladding related parameters
        self._clad = copy.deepcopy(clad)

        if clad_radius <= 0.0 or clad_radius <= fuel_radius:
            raise ValueError("Clad radius must be > fuel radius.")
        elif gap_radius is not None and clad_radius <= gap_radius:
            raise ValueError("Clad radius must be > gap radius.")
        self._clad_radius = clad_radius

        # ======================================================================
        # DANCOFF CORRECTION CALCULATION DATA
        # ----------------------------------------------------------------------

        # Initialize empty list of Dancoff corrections for the fuel
        self._fuel_dancoff_corrections: List[float] = []

        # Initialize empty list of Dancoff corrections for the cladding
        self._clad_dancoff_corrections: List[float] = []

        # Initialize empty variables for Dancoff correction calculations.
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

        self._fuel_isolated_dancoff_fsr_ids = []
        self._gap_isolated_dancoff_fsr_ids = []
        self._clad_isolated_dancoff_fsr_ids = []
        self._mod_isolated_dancoff_fsr_ids = []

        self._fuel_full_dancoff_fsr_ids = []
        self._gap_full_dancoff_fsr_ids = []
        self._clad_full_dancoff_fsr_ids = []
        self._mod_full_dancoff_fsr_ids = []

        self._fuel_isolated_dancoff_fsr_inds = []
        self._gap_isolated_dancoff_fsr_inds = []
        self._clad_isolated_dancoff_fsr_inds = []
        self._mod_isolated_dancoff_fsr_inds = []

        self._fuel_full_dancoff_fsr_inds = []
        self._gap_full_dancoff_fsr_inds = []
        self._clad_full_dancoff_fsr_inds = []
        self._mod_full_dancoff_fsr_inds = []

        # ======================================================================
        # TRANSPORT CALCULATION DATA
        # ----------------------------------------------------------------------
        # Lists of the FSR IDs for each fuel ring, used to homogenize flux
        # spectra for depletion. These will be filled by make_moc_cell.
        self._fuel_ring_fsr_ids: List[List[int]] = []
        for r in range(self.num_fuel_rings):
            self._fuel_ring_fsr_ids.append([])
        self._gap_fsr_ids: List[int] = []
        self._clad_fsr_ids: List[int] = []
        self._mod_fsr_ids: List[int] = []

        self._fuel_ring_fsr_inds: List[List[int]] = []
        for r in range(self.num_fuel_rings):
            self._fuel_ring_fsr_inds.append([])
        self._gap_fsr_inds: List[int] = []
        self._clad_fsr_inds: List[int] = []
        self._mod_fsr_inds: List[int] = []

        # Create list of the different radii for fuel pellet
        self._fuel_radii = []
        if self.num_fuel_rings == 1:
            self._fuel_radii.append(self.fuel_radius)
        else:
            V = np.pi * self.fuel_radius * self.fuel_radius
            Vr = V / self.num_fuel_rings
            for ri in range(self.num_fuel_rings):
                Rin = 0.0
                if ri > 0:
                    Rin = self._fuel_radii[-1]
                Rout = np.sqrt((Vr + np.pi * Rin * Rin) / np.pi)
                if Rout > self.fuel_radius:
                    Rout = self.fuel_radius
                self._fuel_radii.append(Rout)

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

        # Holds all the CrossSection objects used for the real transport
        # calculation. These are NOT stored for each depletion step like with
        # the materials.
        self._fuel_ring_xs: List[CrossSection] = []
        self._gap_xs: Optional[CrossSection] = None
        self._clad_xs: Optional[CrossSection] = None

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
    def fuel_dancoff_corrections(self) -> List[float]:
        return self._fuel_dancoff_corrections

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
    def clad_dancoff_corrections(self) -> List[float]:
        return self._clad_dancoff_corrections

    def _check_dx_dy(self, dx, dy, pintype):
        if pintype == PinCellType.Full:
            if dx < 2.0 * self.clad_radius:
                raise ValueError(
                    "The fuel pin cell x width must be > the diameter of the cladding."
                )
            if dy < 2.0 * self.clad_radius:
                raise ValueError(
                    "The fuel pin cell y width must be > the diameter of the cladding."
                )
        elif pintype in [PinCellType.XN, PinCellType.XP]:
            if dx < self.clad_radius:
                raise ValueError(
                    "The fuel pin cell x width must be > the radius of the cladding."
                )
            if dy < 2.0 * self.clad_radius:
                raise ValueError(
                    "The fuel pin cell y width must be > the diameter of the cladding."
                )
        elif pintype in [PinCellType.YN, PinCellType.YP]:
            if dy < self.clad_radius:
                raise ValueError(
                    "The fuel pin cell y width must be > the radius of the cladding."
                )
            if dx < 2.0 * self.clad_radius:
                raise ValueError(
                    "The fuel pin cell x width must be > the diameter of the cladding."
                )
        else:
            if dx < self.clad_radius:
                raise ValueError(
                    "The fuel pin cell x width must be > the radius of the cladding."
                )
            if dy < self.clad_radius:
                raise ValueError(
                    "The fuel pin cell y width must be > the radius of the cladding."
                )

    def load_nuclides(self, ndl: NDLibrary) -> None:
        """
        Loads all the nuclides for all current materials into the data library.

        Parameters
        ----------
        ndl : NDLibrary
            Nuclear data library which should load the nuclides.
        """
        for ring_mats in self.fuel_ring_materials:
            ring_mats[-1].load_nuclides(ndl)

        if self.gap is not None:
            self.gap.load_nuclides(ndl)

        self.clad.load_nuclides(ndl)

    # ==========================================================================
    # Dancoff Correction Related Methods
    def set_xs_for_fuel_dancoff_calculation(self) -> None:
        """
        Sets the 1-group cross sections to calculate the fuel Dancoff correction.
        """
        self._fuel_dancoff_xs.set(
            CrossSection(
                np.array([1.0e5]), np.array([1.0e5]), np.array([[0.0]]), "Fuel"
            )
        )

        if self._gap_dancoff_xs is not None and self.gap is not None:
            self._gap_dancoff_xs.set(
                CrossSection(
                    np.array([self.gap.potential_xs]),
                    np.array([self.gap.potential_xs]),
                    np.array([[0.0]]),
                    "Gap",
                )
            )

        self._clad_dancoff_xs.set(
            CrossSection(
                np.array([self.clad.potential_xs]),
                np.array([self.clad.potential_xs]),
                np.array([[0.0]]),
                "Clad",
            )
        )

    def set_xs_for_clad_dancoff_calculation(self, ndl: NDLibrary) -> None:
        """
        Sets the 1-group cross sections to calculate the clad Dancoff correction.

        Parameters
        ----------
        ndl : NDLibrary
            Nuclear data library for obtaining potential scattering cross
            sections.
        """
        # Create average fuel mixture
        fuel_mats = []
        fuel_vols = []
        for ring in self.fuel_ring_materials:
            fuel_mats.append(ring[-1])
            fuel_vols.append(1.0 / self.num_fuel_rings)
        avg_fuel: Material = mix_materials(
            fuel_mats, fuel_vols, MixingFraction.Volume, ndl
        )

        self._fuel_dancoff_xs.set(
            CrossSection(
                np.array([avg_fuel.potential_xs]),
                np.array([avg_fuel.potential_xs]),
                np.array([[0.0]]),
                "Fuel",
            )
        )

        if self._gap_dancoff_xs is not None and self.gap is not None:
            self._gap_dancoff_xs.set(
                CrossSection(
                    np.array([self.gap.potential_xs]),
                    np.array([self.gap.potential_xs]),
                    np.array([[0.0]]),
                    "Gap",
                )
            )

        self._clad_dancoff_xs.set(
            CrossSection(
                np.array([1.0e5]),
                np.array([1.0e5]),
                np.array([[0.0]]),
                "Clad",
            )
        )

    def make_dancoff_moc_cell(
        self,
        moderator_xs: CrossSection,
        dx: float,
        dy: float,
        pintype: PinCellType,
        isolated: bool,
    ) -> SimplePinCell:
        """
        Makes a simplified cell suitable for performing Dancoff correction
        calculations. The flat source region IDs are stored locally in the
        FuelPin object.

        Parameters
        ----------
        moderator_xs : CrossSection
            One group cross sections for the moderator. Total should equal
            absorption (i.e. no scattering) and should be equal to the
            macroscopic potential cross section.
        dx : float
            Width of the cell along x.
        dy : float
            Width of the cell along y.
        pintype : PinCellType
            How the pin cell should be split (along x, y, or only a quadrant).
        isolated : bool
            If True, the FSR IDs are stored for the isolated pin. Otherwise,
            they are stored for the full pin.

        Returns
        -------
        SimplifiedPinCell
            Pin cell object for MOC Dancoff correction calculation.
        """
        self._check_dx_dy(dx, dy, pintype)

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
        cell = SimplePinCell(radii, xs, dx, dy, pintype)

        # Get the FSR IDs for the regions of interest
        cell_fsr_ids = list(cell.get_all_fsr_ids())
        cell_fsr_ids.sort()

        if isolated:
            self._fuel_isolated_dancoff_fsr_ids.append(cell_fsr_ids[0])
            if self.gap is None:
                self._clad_isolated_dancoff_fsr_ids.append(cell_fsr_ids[1])
                self._mod_isolated_dancoff_fsr_ids.append(cell_fsr_ids[2])
            else:
                self._gap_isolated_dancoff_fsr_ids.append(cell_fsr_ids[1])
                self._clad_isolated_dancoff_fsr_ids.append(cell_fsr_ids[2])
                self._mod_isolated_dancoff_fsr_ids.append(cell_fsr_ids[3])
        else:
            self._fuel_full_dancoff_fsr_ids.append(cell_fsr_ids[0])
            if self.gap is None:
                self._clad_full_dancoff_fsr_ids.append(cell_fsr_ids[1])
                self._mod_full_dancoff_fsr_ids.append(cell_fsr_ids[2])
            else:
                self._gap_full_dancoff_fsr_ids.append(cell_fsr_ids[1])
                self._clad_full_dancoff_fsr_ids.append(cell_fsr_ids[2])
                self._mod_full_dancoff_fsr_ids.append(cell_fsr_ids[3])

        return cell

    def populate_dancoff_fsr_indexes(
        self, isomoc: MOCDriver, fullmoc: MOCDriver
    ) -> None:
        """
        Obtains the flat source region indexes for all of the flat source
        regions used in the Dancoff correction calculations.

        Parameters
        ----------
        isomoc : MOCDriver
            MOC simulation for the isolated pin.
        fullmoc : MOCDriver
            MOC simulation for the full geometry.
        """
        self._fuel_isolated_dancoff_fsr_inds = []
        self._gap_isolated_dancoff_fsr_inds = []
        self._clad_isolated_dancoff_fsr_inds = []
        self._mod_isolated_dancoff_fsr_inds = []

        self._fuel_full_dancoff_fsr_inds = []
        self._gap_full_dancoff_fsr_inds = []
        self._clad_full_dancoff_fsr_inds = []
        self._mod_full_dancoff_fsr_inds = []

        for id in self._fuel_isolated_dancoff_fsr_ids:
            self._fuel_isolated_dancoff_fsr_inds.append(isomoc.get_fsr_indx(id, 0))
        for id in self._gap_isolated_dancoff_fsr_ids:
            self._gap_isolated_dancoff_fsr_inds.append(isomoc.get_fsr_indx(id, 0))
        for id in self._clad_isolated_dancoff_fsr_ids:
            self._clad_isolated_dancoff_fsr_inds.append(isomoc.get_fsr_indx(id, 0))
        for id in self._mod_isolated_dancoff_fsr_ids:
            self._mod_isolated_dancoff_fsr_inds.append(isomoc.get_fsr_indx(id, 0))

        for id in self._fuel_full_dancoff_fsr_ids:
            self._fuel_full_dancoff_fsr_inds.append(fullmoc.get_fsr_indx(id, 0))
        for id in self._gap_full_dancoff_fsr_ids:
            self._gap_full_dancoff_fsr_inds.append(fullmoc.get_fsr_indx(id, 0))
        for id in self._clad_full_dancoff_fsr_ids:
            self._clad_full_dancoff_fsr_inds.append(fullmoc.get_fsr_indx(id, 0))
        for id in self._mod_full_dancoff_fsr_ids:
            self._mod_full_dancoff_fsr_inds.append(fullmoc.get_fsr_indx(id, 0))

    def set_isolated_dancoff_fuel_sources(
        self, isomoc: MOCDriver, moderator: Material
    ) -> None:
        """
        Initializes the fixed sources for the isolated MOC calculation required
        in computing Dancoff corrections. Sources are set for a fuel Dancoff
        corrections calculation.

        Parameters
        ----------
        isomoc : MOCDriver
            MOC simulation for the isolated geometry.
        moderator : Material
            Material definition for the moderator, used to obtain the potential
            scattering cross section.
        """
        # Fuel sources should all be zero !
        for ind in self._fuel_isolated_dancoff_fsr_inds:
            isomoc.set_extern_src(ind, 0, 0.0)

        # Gap sources should all be potential_xs
        if self.gap is not None:
            pot_xs = self.gap.potential_xs
            for ind in self._gap_isolated_dancoff_fsr_inds:
                isomoc.set_extern_src(ind, 0, pot_xs)

        # Clad sources should all be potential_xs
        pot_xs = self.clad.potential_xs
        for ind in self._clad_isolated_dancoff_fsr_inds:
            isomoc.set_extern_src(ind, 0, pot_xs)

        # Moderator sources should all be potential_xs
        pot_xs = moderator.potential_xs
        for ind in self._mod_isolated_dancoff_fsr_inds:
            isomoc.set_extern_src(ind, 0, pot_xs)

    def set_isolated_dancoff_clad_sources(
        self, isomoc: MOCDriver, moderator: Material, ndl: NDLibrary
    ) -> None:
        """
        Initializes the fixed sources for the isolated MOC calculation required
        in computing Dancoff corrections. Sources are set for a clad Dancoff
        correction calculation.

        Parameters
        ----------
        isomoc : MOCDriver
            MOC simulation for the isolated geometry.
        moderator : Material
            Material definition for the moderator, used to obtain the potential
            scattering cross section.
        ndl : NDLibrary
            Nuclear data library for obtaining potential scattering cross
            sections.
        """
        # Create average fuel mixture
        fuel_mats = []
        fuel_vols = []
        for ring in self.fuel_ring_materials:
            fuel_mats.append(ring[-1])
            fuel_vols.append(1.0 / self.num_fuel_rings)
        avg_fuel: Material = mix_materials(
            fuel_mats, fuel_vols, MixingFraction.Volume, ndl
        )

        # Fuel sources should all be potential_xs
        pot_xs = avg_fuel.potential_xs
        for ind in self._fuel_isolated_dancoff_fsr_inds:
            isomoc.set_extern_src(ind, 0, pot_xs)

        # Gap sources should all be potential_xs
        if self.gap is not None:
            pot_xs = self.gap.potential_xs
            for ind in self._gap_isolated_dancoff_fsr_inds:
                isomoc.set_extern_src(ind, 0, pot_xs)

        # Clad sources should all be zero !
        for ind in self._clad_isolated_dancoff_fsr_inds:
            isomoc.set_extern_src(ind, 0, 0.0)

        # Moderator sources should all be potential_xs
        pot_xs = moderator.potential_xs
        for ind in self._mod_isolated_dancoff_fsr_inds:
            isomoc.set_extern_src(ind, 0, pot_xs)

    def set_full_dancoff_fuel_sources(
        self, fullmoc: MOCDriver, moderator: Material
    ) -> None:
        """
        Initializes the fixed sources for the full MOC calculation required
        in computing Dancoff corrections. Sources are set for a fuel Dancoff
        correction calculation.

        Parameters
        ----------
        fullmoc : MOCDriver
            MOC simulation for the full geometry.
        moderator : Material
            Material definition for the moderator, used to obtain the potential
            scattering cross section.
        """
        # Fuel sources should all be zero !
        for ind in self._fuel_full_dancoff_fsr_inds:
            fullmoc.set_extern_src(ind, 0, 0.0)

        # Gap sources should all be potential_xs
        if self.gap is not None:
            pot_xs = self.gap.potential_xs
            for ind in self._gap_full_dancoff_fsr_inds:
                fullmoc.set_extern_src(ind, 0, pot_xs)

        # Clad sources should all be potential_xs
        pot_xs = self.clad.potential_xs
        for ind in self._clad_full_dancoff_fsr_inds:
            fullmoc.set_extern_src(ind, 0, pot_xs)

        # Moderator sources should all be potential_xs
        pot_xs = moderator.potential_xs
        for ind in self._mod_full_dancoff_fsr_inds:
            fullmoc.set_extern_src(ind, 0, pot_xs)

    def set_full_dancoff_clad_sources(
        self, fullmoc: MOCDriver, moderator: Material, ndl: NDLibrary
    ) -> None:
        """
        Initializes the fixed sources for the full MOC calculation required
        in computing Dancoff corrections. Sources are set for a clad Dancoff
        correction calculation.

        Parameters
        ----------
        isomoc : MOCDriver
            MOC simulation for the isolated geometry.
        moderator : Material
            Material definition for the moderator, used to obtain the potential
            scattering cross section.
        ndl : NDLibrary
            Nuclear data library for obtaining potential scattering cross
            sections.
        """
        # Create average fuel mixture
        fuel_mats = []
        fuel_vols = []
        for ring in self.fuel_ring_materials:
            fuel_mats.append(ring[-1])
            fuel_vols.append(1.0 / self.num_fuel_rings)
        avg_fuel: Material = mix_materials(
            fuel_mats, fuel_vols, MixingFraction.Volume, ndl
        )

        # Fuel sources should all be potential_xs
        pot_xs = avg_fuel.potential_xs
        for ind in self._fuel_full_dancoff_fsr_inds:
            fullmoc.set_extern_src(ind, 0, pot_xs)

        # Gap sources should all be potential_xs
        if self.gap is not None:
            pot_xs = self.gap.potential_xs
            for ind in self._gap_full_dancoff_fsr_inds:
                fullmoc.set_extern_src(ind, 0, pot_xs)

        # Clad sources should all be zero !
        for ind in self._clad_full_dancoff_fsr_inds:
            fullmoc.set_extern_src(ind, 0, 0.0)

        # Moderator sources should all be potential_xs
        pot_xs = moderator.potential_xs
        for ind in self._mod_full_dancoff_fsr_inds:
            fullmoc.set_extern_src(ind, 0, pot_xs)

    def compute_fuel_dancoff_correction(
        self, isomoc: MOCDriver, fullmoc: MOCDriver
    ) -> float:
        """
        Computes the Dancoff correction for the fuel region of the fuel pin.

        Parameters
        ----------
        isomoc : MOCDriver
            MOC simulation for the isolated geometry (previously solved).
        fullmoc : MOCDriver
            MOC simulation for the full geometry (previously solved).

        Returns
        -------
        float
            Dancoff correction for the fuel region.
        """
        iso_flux = isomoc.homogenize_flux_spectrum(
            self._fuel_isolated_dancoff_fsr_inds
        )[0]
        full_flux = fullmoc.homogenize_flux_spectrum(self._fuel_full_dancoff_fsr_inds)[
            0
        ]
        return (iso_flux - full_flux) / iso_flux

    def compute_clad_dancoff_correction(
        self, isomoc: MOCDriver, fullmoc: MOCDriver
    ) -> float:
        """
        Computes the Dancoff correction for the cladding region of the fuel pin.

        Parameters
        ----------
        isomoc : MOCDriver
            MOC simulation for the isolated geometry (previously solved).
        fullmoc : MOCDriver
            MOC simulation for the full geometry (previously solved).

        Returns
        -------
        float
            Dancoff correction for the cladding region.
        """
        iso_flux = isomoc.homogenize_flux_spectrum(
            self._clad_isolated_dancoff_fsr_inds
        )[0]
        full_flux = fullmoc.homogenize_flux_spectrum(self._clad_full_dancoff_fsr_inds)[
            0
        ]
        return (iso_flux - full_flux) / iso_flux

    def append_fuel_dancoff_correction(self, C) -> None:
        """
        Saves new Dancoff correction for the fuel that will be used for all
        subsequent cross section updates.

        Parameters
        ----------
        C : float
            New Dancoff correction.
        """
        if C < 0.0 or C > 1.0:
            raise ValueError(
                f"Dancoff correction must be in range [0, 1]. Was provided {C}."
            )
        self._fuel_dancoff_corrections.append(C)

    def append_clad_dancoff_correction(self, C) -> None:
        """
        Saves new Dancoff correction for the cladding that will be used for all
        subsequent cross section updates.

        Parameters
        ----------
        C : float
            New Dancoff correction.
        """
        if C < 0.0 or C > 1.0:
            raise ValueError(
                f"Dancoff correction must be in range [0, 1]. Was provided {C}."
            )
        self._clad_dancoff_corrections.append(C)

    # ==========================================================================
    # Transport Calculation Related Methods
    def set_fuel_xs_for_depletion_step(self, t: int, ndl: NDLibrary) -> None:
        """
        Constructs the CrossSection object for all fuel rings of the pin at the
        specified depletion step.

        Parameters
        ----------
        t : int
            Index for the depletion step.
        ndl : NDLibrary
            Nuclear data library to use for cross sections.
        """
        # Do the fuel cross sections
        if len(self._fuel_ring_xs) == 0:
            # Create initial CrossSection objects
            if self.num_fuel_rings == 1:
                # Compute escape xs
                Ee = 1.0 / (2.0 * self.fuel_radius)
                self._fuel_ring_xs.append(
                    self._fuel_ring_materials[0][t].carlvik_xs(
                        self._fuel_dancoff_corrections[t], Ee, ndl
                    )
                )
                if self._fuel_ring_xs[-1].name == "":
                    self._fuel_ring_xs[-1].set_name("Fuel")
            else:
                # Do each ring
                for ri in range(self.num_fuel_rings):
                    Rin = 0.0
                    if ri > 0:
                        Rin = self._fuel_radii[ri - 1]
                    Rout = self._fuel_radii[ri]
                    self._fuel_ring_xs.append(
                        self._fuel_ring_materials[ri][t].ring_carlvik_xs(
                            self._fuel_dancoff_corrections[t],
                            self.fuel_radius,
                            Rin,
                            Rout,
                            ndl,
                        )
                    )
                    if self._fuel_ring_xs[-1].name == "":
                        self._fuel_ring_xs[-1].set_name("Fuel")

        elif len(self._fuel_ring_xs) == self.num_fuel_rings:
            # Reset XS values. Cannot reassign or pointers will be broken !
            if self.num_fuel_rings == 1:
                # Compute escape xs
                Ee = 1.0 / (2.0 * self.fuel_radius)
                self._fuel_ring_xs[0].set(
                    self._fuel_ring_materials[0][t].carlvik_xs(
                        self._fuel_dancoff_corrections[t], Ee, ndl
                    )
                )
                if self._fuel_ring_xs[0].name == "":
                    self._fuel_ring_xs[0].set_name("Fuel")
            else:
                # Do each ring
                for ri in range(self.num_fuel_rings):
                    Rin = 0.0
                    if ri > 0:
                        Rin = self._fuel_radii[ri - 1]
                    Rout = self._fuel_radii[ri]
                    self._fuel_ring_xs[ri].set(
                        self._fuel_ring_materials[ri][t].ring_carlvik_xs(
                            self._fuel_dancoff_corrections[t],
                            self.fuel_radius,
                            Rin,
                            Rout,
                            ndl,
                        )
                    )
                    if self._fuel_ring_xs[ri].name == "":
                        self._fuel_ring_xs[ri].set_name("Fuel")
        else:
            raise RuntimeError(
                "Number of fuel cross sections does not agree with the number of fuel rings."
            )

    def set_gap_xs(self, ndl: NDLibrary) -> None:
        """
        Constructs the CrossSection object for the gap between the fuel pellet
        and the cladding of the pin.

        Parameters
        ----------
        ndl : NDLibrary
            Nuclear data library to use for cross sections.
        """
        if self.gap is not None:
            if self._gap_xs is None:
                self._gap_xs = self.gap.dilution_xs([1.0e10] * self.gap.size, ndl)
            else:
                self._gap_xs.set(self.gap.dilution_xs([1.0e10] * self.gap.size, ndl))

            if self._gap_xs.name == "":
                self._gap_xs.set_name("Gap")

    def set_clad_xs_for_depletion_step(self, t: int, ndl: NDLibrary) -> None:
        """
        Constructs the CrossSection object for the cladding of the pin at the
        specified depletion step. The depletion step only changes the Dancoff
        correction, not the cladding composition.

        Parameters
        ----------
        t : int
            Index for the depletion step.
        ndl : NDLibrary
            Nuclear data library to use for cross sections.
        """
        # Compute escape xs
        Ee = 0.0
        if self.gap_radius is not None:
            Ee = 1.0 / (2.0 * (self.clad_radius - self.gap_radius))
        else:
            Ee = 1.0 / (2.0 * (self.clad_radius - self.fuel_radius))

        # Get / set the xs
        if self._clad_xs is None:
            self._clad_xs = self.clad.roman_xs(
                self._clad_dancoff_corrections[t], Ee, ndl
            )
        else:
            self._clad_xs.set(
                self.clad.roman_xs(self._clad_dancoff_corrections[t], Ee, ndl)
            )

        if self._clad_xs.name == "":
            self._clad_xs.set_name("Clad")

    def make_moc_cell(
        self,
        moderator_xs: CrossSection,
        dx: float,
        dy: float,
        pintype: PinCellType,
    ) -> PinCell:
        """
        Constructs the pin cell object used in for the global MOC simulation.

        Parameters
        ----------
        moderator_xs : CrossSection
            Cross sections to use for the moderator surrounding the fuel pin.
        dx : float
            Width of the cell along x.
        dy : float
            Width of the cell along y.
        pintype : PinCellType
            How the pin cell should be split (along x, y, or only a quadrant).

        Returns
        -------
        PinCell
            Pin cell suitable for the true MOC calculation.
        """
        if len(self._fuel_ring_xs) != self.num_fuel_rings:
            raise RuntimeError("Fuel cross sections have not yet been built.")
        if self.gap is not None and self._gap_xs is None:
            raise RuntimeError("Gap cross section has not yet been built.")
        if self._clad_xs is None:
            raise RuntimeError("Clad cross section has not yet been built.")
        self._check_dx_dy(dx, dy, pintype)

        # Initialize the radii and cross section lists with the fuel info
        radii = [r for r in self._fuel_radii]
        xss = [xs for xs in self._fuel_ring_xs]

        # Add the gap (if present)
        if self._gap_xs is not None:
            radii.append(self.gap_radius)
            xss.append(self._gap_xs)

        # Add cladding
        radii.append(self.clad_radius)
        xss.append(self._clad_xs)

        # Add another ring of moderator if possible
        if pintype == PinCellType.Full and min(dx, dy) > 2.0 * self.clad_radius:
            radii.append(0.5 * min(dx, dy))
            xss.append(moderator_xs)
        elif (
            pintype in [PinCellType.XN, PinCellType.XP]
            and dx > self.clad_radius
            and dy > 2.0 * self.clad_radius
        ):
            radii.append(min(dx, 0.5 * dy))
            xss.append(moderator_xs)
        elif (
            pintype in [PinCellType.YN, PinCellType.YP]
            and dy > self.clad_radius
            and dx > 2.0 * self.clad_radius
        ):
            radii.append(min(0.5 * dx, dy))
            xss.append(moderator_xs)
        elif dx > self.clad_radius and dy > self.clad_radius:
            radii.append(min(dx, dy))
            xss.append(moderator_xs)

        # Add moderator to the end of materials
        xss.append(moderator_xs)

        return PinCell(radii, xss, dx, dy, pintype)
