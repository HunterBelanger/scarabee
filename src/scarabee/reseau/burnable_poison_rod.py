from .._scarabee import NDLibrary, DepletionChain, Material, CrossSection, MOCDriver
import numpy as np
from typing import Optional, List, Tuple
import copy


class BurnablePoisonRod:
    """
    Represents a burnable poison rod in a PWR which is inserted into an empty
    guide tube. This can be a wet annular burnable absorber (WABA) or a
    borosilicate glass (BSG) burnable absorber.

    Parameters
    ----------
    center : optional Material
        Material at the center of the poison pin. If None, the center will be
        filled with moderator.
    clad : Material
        Material which describes the inner and outer cladding.
    gap : Material
        Material which describes the gap between the cladding and poison.
    poison : Material
        Material which describes the poison.
    center_radius : float
        Outer radius of the center of the poison rod.
    inner_clad_radius : float
        Outer radius of the inner cladding at the center of the rod.
    inner_gap_radius : float
        Outer radius of the inner gap between the cladding and poison.
    poison_radius : float
        Outer radius of the poison.
    outer_gap_radius : float
        Outer radius of the outer gap between the cladding and poison.
    outer_clad_radius : float
        Outer radius of the outer cladding of the poison rod.

    Attributes
    ----------
    center : optional Material
        Material at the center of the poison pin. If None, the center will be
        filled with moderator.
    clad : Material
        Material which describes the inner and outer cladding.
    gap : Material
        Material which describes the gap between the cladding and poison.
    poison_materials : list of Material
        Contains the poison material for each depletion time step.
    center_radius : float
        Outer radius of the center of the poison rod.
    inner_clad_radius : float
        Outer radius of the inner cladding at the center of the rod.
    inner_gap_radius : float
        Outer radius of the inner gap between the cladding and poison.
    poison_radius : float
        Outer radius of the poison.
    outer_gap_radius : float
        Outer radius of the outer gap between the cladding and poison.
    outer_clad_radius : float
        Outer radius of the outer cladding of the poison rod.
    """

    def __init__(
        self,
        center: Optional[Material],
        clad: Material,
        gap: Material,
        poison: Material,
        center_radius: float,
        inner_clad_radius: float,
        inner_gap_radius: float,
        poison_radius: float,
        outer_gap_radius: float,
        outer_clad_radius: float,
    ):
        # Make sure all radii are sorted
        tmp_radii_list = [
            center_radius,
            inner_clad_radius,
            inner_gap_radius,
            poison_radius,
            outer_gap_radius,
            outer_clad_radius,
        ]
        if not sorted(tmp_radii_list):
            raise ValueError("Burnable poison rod radii are not sorted.")

        # Set materials
        self._center = copy.deepcopy(center)
        self._clad = copy.deepcopy(clad)
        self._gap = copy.deepcopy(gap)
        self._poison_materials = [copy.deepcopy(poison)]

        # Set radii
        self._center_radius = center_radius
        self._inner_clad_radius = inner_clad_radius
        self._inner_gap_radius = inner_gap_radius
        self._poison_radius = poison_radius
        self._outer_gap_radius = outer_gap_radius
        self._outer_clad_radius = outer_clad_radius

        # ======================================================================
        # DANCOFF CORRECTION CALCULATION DATA
        # ----------------------------------------------------------------------
        self._center_dancoff_xs: Optional[CrossSection] = None
        if self.center is not None:
            self._center_dancoff_xs = CrossSection(
                np.array([self.center.potential_xs]),
                np.array([self.center.potential_xs]),
                np.array([[0.0]]),
                "BPR Clad",
            )

        self._clad_dancoff_xs: CrossSection = CrossSection(
            np.array([self.clad.potential_xs]),
            np.array([self.clad.potential_xs]),
            np.array([[0.0]]),
            "BPR Clad",
        )

        self._gap_dancoff_xs: CrossSection = CrossSection(
            np.array([self.gap.potential_xs]),
            np.array([self.gap.potential_xs]),
            np.array([[0.0]]),
            "BPR Gap",
        )

        # Make initial poison Dancoff xs with initial material
        self._poison_dancoff_xs: CrossSection = CrossSection(
            np.array([poison.potential_xs]),
            np.array([poison.potential_xs]),
            np.array([[0.0]]),
            "BPR Poison",
        )

        # ======================================================================
        # TRANSPORT CALCULATION DATA
        # ----------------------------------------------------------------------
        self._center_xs: Optional[CrossSection] = None
        self._clad_xs: Optional[CrossSection] = None
        self._gap_xs: Optional[CrossSection] = None
        self._poison_xs: Optional[CrossSection] = None

        # FSR IDs and indices
        self._center_fsr_id: List[int] = []
        self._center_fsr_ind: List[int] = []
        self._clad_fsr_id: List[int] = []
        self._clad_fsr_ind: List[int] = []
        self._gap_fsr_id: List[int] = []
        self._gap_fsr_ind: List[int] = []
        self._poison_fsr_id: List[int] = []
        self._poison_fsr_ind: List[int] = []

        # Flux spectrum of poison needed for depletion
        self._poison_flux_spectrum: np.ndarray = np.array([])

    @property
    def center(self) -> Optional[Material]:
        return self._center

    @property
    def clad(self) -> Material:
        return self._clad

    @property
    def gap(self) -> Material:
        return self._gap

    @property
    def poison_materials(self) -> List[Material]:
        return self._poison_materials

    @property
    def center_radius(self) -> float:
        return self._center_radius

    @property
    def inner_clad_radius(self) -> float:
        return self._inner_clad_radius

    @property
    def inner_gap_radius(self) -> float:
        return self._inner_gap_radius

    @property
    def poison_radius(self) -> float:
        return self._poison_radius

    @property
    def outer_gap_radius(self) -> float:
        return self._outer_gap_radius

    @property
    def outer_clad_radius(self) -> float:
        return self._outer_clad_radius

    def load_nuclides(self, ndl: NDLibrary) -> None:
        """
        Loads all the nuclides for all current materials into the data library.

        Parameters
        ----------
        ndl : NDLibrary
            Nuclear data library which should load the nuclides.
        """
        if self.center is not None:
            self.center.load_nuclides(ndl)
        self.clad.load_nuclides(ndl)
        self.gap.load_nuclides(ndl)
        self.poison.load_nuclides(ndl)

    # ==========================================================================
    # Dancoff Correction Related Methods

    def _make_dancoff_moc_cell(
        self, moderator_xs: CrossSection
    ) -> Tuple[List[float], List[CrossSection]]:
        """
        Returns the list of radii and list of cross sections for the poison pin.
        """
        pass

    def set_xs_for_dancoff_calculation(self) -> None:
        """
        Sets the cross sections for the Dancoff calculations based on the most
        recent poison composition.
        """
        pass

    def set_isolated_dancoff_fuel_sources(
        self, isomoc: MOCDriver, moderator: Material
    ) -> None:
        """
        Initializes the fixed sources for the isolated MOC calculation required
        in computing Dancoff corrections. Sources are set for a fuel Dancoff
        correction calculation.

        Parameters
        ----------
        isomoc : MOCDriver
            MOC simulation for the isolated geometry.
        moderator : Material
            Material definition for the moderator, used to obtain the potential
            scattering cross section.
        """
        pass

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
        pass

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
        pass

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
        pass

    # ==========================================================================
    # Transport Calculation Related Methods

    def set_gap_xs(self, ndl: NDLibrary) -> None:
        """
        Constructs the CrossSection object for the gap between the poison and
        the cladding.

        Parameters
        ----------
        ndl : NDLibrary
            Nuclear data library to use for cross sections.
        """
        pass

    def set_clad_xs(self, ndl: NDLibrary) -> None:
        """
        Constructs the CrossSection object for the cladding of the poison rod.

        Parameters
        ----------
        ndl : NDLibrary
            Nuclear data library to use for cross sections.
        """
        pass

    def set_poison_xs_for_depletion_step(self, t: int, ndl: NDLibrary) -> None:
        """
        Constructs the CrossSection object for the poison at the specified
        depletion step.

        Parameters
        ----------
        t : int
            Index for the depletion step.
        ndl : NDLibrary
            Nuclear data library to use for cross sections.
        """
        pass

    def _make_moc_cell(
        self, moderator_xs: CrossSection
    ) -> Tuple[List[float], List[CrossSection]]:
        """
        Returns the list of radii and list of cross sections for the poison pin.
        """
        pass

    def populate_fsr_indexes(self, moc: MOCDriver) -> None:
        """
        Obtains the flat source region indexes for all of the flat source
        regions used in the full MOC calculations.

        Parameters
        ----------
        moc : MOCDriver
            MOC simulation for the full calculations.
        """
        pass

    def obtain_poison_flux_spectrum(self, moc: MOCDriver) -> None:
        """
        Computes average flux spectrum in the poison from the MOC simulation.

        Parameters
        ----------
        moc : MOCDriver
            MOC simulation for the full calculations.
        """
        pass

    def predict_depletion(
        self, dt: float, chain: DepletionChain, ndl: NDLibrary
    ) -> None:
        """
        Performs the predictor in the integration of the Bateman equation.
        The provided time step should therefore be half of the anticipated full
        time step. The predicted material compositions are appended to the
        materials lists.

        Paramters
        ---------
        dt : float
            Time step for the predictor in seconds.
        chain : DepletionChain
            Depletion chain to use for radioactive decay and transmutation.
        ndl : NDLibrary
            Nuclear data library.
        """
        if dt <= 0:
            raise ValueError("Predictor time step must be > 0.")
        pass

    def correct_depletion(
        self, dt: float, chain: DepletionChain, ndl: NDLibrary
    ) -> None:
        """
        Performs the corrector in the integration of the Bateman equation.
        The provided time step should therefore be the full anticipated time
        step. The corrected material compositions replace the ones where were
        appended in the corrector step.

        Paramters
        ---------
        dt : float
            Time step for the predictor in seconds.
        chain : DepletionChain
            Depletion chain to use for radioactive decay and transmutation.
        ndl : NDLibrary
            Nuclear data library.
        """
        if dt <= 0:
            raise ValueError("Corrector time step must be > 0.")
        pass
