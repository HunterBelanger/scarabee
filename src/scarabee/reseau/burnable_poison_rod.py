from .._scarabee import (
    NDLibrary,
    DepletionChain,
    MaterialComposition,
    Material,
    CrossSection,
    MOCDriver,
    DepletionChain,
    build_depletion_matrix,
)
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

        # Center IDs and indices are only populated if the center material is
        # not None, which indicates that the center should be moderator. If the
        # center is moderator, then those IDs and indices are kept by the
        # parent GuideTube object.

        self._center_isolated_dancoff_fsr_ids = []
        self._center_isolated_dancoff_fsr_inds = []
        self._clad_isolated_dancoff_fsr_ids = []
        self._clad_isolated_dancoff_fsr_inds = []
        self._gap_isolated_dancoff_fsr_ids = []
        self._gap_isolated_dancoff_fsr_inds = []
        self._poison_isolated_dancoff_fsr_ids = []
        self._poison_isolated_dancoff_fsr_inds = []

        self._center_full_dancoff_fsr_ids = []
        self._center_full_dancoff_fsr_inds = []
        self._clad_full_dancoff_fsr_ids = []
        self._clad_full_dancoff_fsr_inds = []
        self._gap_full_dancoff_fsr_ids = []
        self._gap_full_dancoff_fsr_inds = []
        self._poison_full_dancoff_fsr_ids = []
        self._poison_full_dancoff_fsr_inds = []

        # ======================================================================
        # TRANSPORT CALCULATION DATA
        # ----------------------------------------------------------------------
        self._center_xs: Optional[CrossSection] = None
        self._clad_xs: Optional[CrossSection] = None
        self._gap_xs: Optional[CrossSection] = None
        self._poison_xs: Optional[CrossSection] = None

        # FSR IDs and indices. Same rule for center as before.
        self._center_fsr_ids: List[int] = []
        self._center_fsr_inds: List[int] = []
        self._clad_fsr_ids: List[int] = []
        self._clad_fsr_inds: List[int] = []
        self._gap_fsr_ids: List[int] = []
        self._gap_fsr_inds: List[int] = []
        self._poison_fsr_ids: List[int] = []
        self._poison_fsr_inds: List[int] = []

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
        self.poison_materials[-1].load_nuclides(ndl)

    # ==========================================================================
    # Dancoff Correction Related Methods

    def _make_dancoff_moc_cell(
        self, moderator_xs: CrossSection
    ) -> Tuple[List[float], List[CrossSection]]:
        """
        Returns the list of radii and list of cross sections for the poison pin.
        """
        pass

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
        self._center_isolated_dancoff_fsr_inds = []
        self._clad_isolated_dancoff_fsr_inds = []
        self._gap_isolated_dancoff_fsr_inds = []
        self._poison_isolated_dancoff_fsr_inds = []

        self._center_full_dancoff_fsr_inds = []
        self._clad_full_dancoff_fsr_inds = []
        self._gap_full_dancoff_fsr_inds = []
        self._poison_full_dancoff_fsr_inds = []

        for id in self._center_isolated_dancoff_fsr_ids:
            self._center_isolated_dancoff_fsr_inds.append(isomoc.get_fsr_indx(id, 0))
        for id in self._clad_isolated_dancoff_fsr_ids:
            self._clad_isolated_dancoff_fsr_inds.append(isomoc.get_fsr_indx(id, 0))
        for id in self._gap_isolated_dancoff_fsr_ids:
            self._gap_isolated_dancoff_fsr_inds.append(isomoc.get_fsr_indx(id, 0))
        for id in self._poison_isolated_dancoff_fsr_ids:
            self._poison_isolated_dancoff_fsr_inds.append(isomoc.get_fsr_indx(id, 0))

        for id in self._center_full_dancoff_fsr_ids:
            self._center_full_dancoff_fsr_inds.append(fullmoc.get_fsr_indx(id, 0))
        for id in self._clad_full_dancoff_fsr_ids:
            self._clad_full_dancoff_fsr_inds.append(fullmoc.get_fsr_indx(id, 0))
        for id in self._gap_full_dancoff_fsr_ids:
            self._gap_full_dancoff_fsr_inds.append(fullmoc.get_fsr_indx(id, 0))
        for id in self._poison_full_dancoff_fsr_ids:
            self._poison_full_dancoff_fsr_inds.append(fullmoc.get_fsr_indx(id, 0))

    def set_xs_for_dancoff_calculation(self) -> None:
        """
        Sets the cross sections for the Dancoff calculations based on the most
        recent poison composition.
        """
        if self._center_dancoff_xs is not None and self.center is not None:
            self._center_dancoff_xs.set(
                CrossSection(
                    np.array([self.center.potential_xs]),
                    np.array([self.center.potential_xs]),
                    np.array([[0.0]]),
                    "BPR Clad",
                )
            )

        self._clad_dancoff_xs.set(
            CrossSection(
                np.array([self.clad.potential_xs]),
                np.array([self.clad.potential_xs]),
                np.array([[0.0]]),
                "BPR Clad",
            )
        )

        self._gap_dancoff_xs.set(
            CrossSection(
                np.array([self.gap.potential_xs]),
                np.array([self.gap.potential_xs]),
                np.array([[0.0]]),
                "BPR Gap",
            )
        )

        # Make initial poison Dancoff xs with initial material
        self._poison_dancoff_xs.set(
            CrossSection(
                np.array([self.poison_materials[-1].potential_xs]),
                np.array([self.poison_materials[-1].potential_xs]),
                np.array([[0.0]]),
                "BPR Poison",
            )
        )

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
        # We only set sources for the center if it isn't moderator ! Otherwise,
        # the GuideTube class is responsible for taking care of this. This case
        # is indicated by the center material being None. For this case, the
        # list of FSR indices should be empty.
        if self.center is not None:
            pot_xs = self.center.potential_xs
            for ind in self._center_isolated_dancoff_fsr_inds:
                isomoc.set_extern_src(ind, 0, pot_xs)

        pot_xs = self.clad.potential_xs
        for ind in self._clad_isolated_dancoff_fsr_inds:
            isomoc.set_extern_src(ind, 0, pot_xs)

        pot_xs = self.gap.potential_xs
        for ind in self._gap_isolated_dancoff_fsr_inds:
            isomoc.set_extern_src(ind, 0, pot_xs)

        pot_xs = self.poison_materials[-1].potential_xs
        for ind in self._poison_isolated_dancoff_fsr_inds:
            isomoc.set_extern_src(ind, 0, pot_xs)

    def set_isolated_dancoff_clad_sources(
        self, isomoc: MOCDriver, moderator: Material, ndl: NDLibrary
    ) -> None:
        """
        Initializes the fixed sources for the isolated MOC calculation required
        in computing Dancoff corrections. Sources are set for a clad Dancoff
        correction calculation.

        The cladding of a burnable poison pin is no self-shielded. Therefore,
        this method is an alias to set_isolated_dancoff_fuel_sources.

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
        self.set_isolated_dancoff_fuel_sources(isomoc, moderator)

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
        # We only set sources for the center if it isn't moderator ! Otherwise,
        # the GuideTube class is responsible for taking care of this. This case
        # is indicated by the center material being None. For this case, the
        # list of FSR indices should be empty.
        if self.center is not None:
            pot_xs = self.center.potential_xs
            for ind in self._center_full_dancoff_fsr_inds:
                fullmoc.set_extern_src(ind, 0, pot_xs)

        pot_xs = self.clad.potential_xs
        for ind in self._clad_full_dancoff_fsr_inds:
            fullmoc.set_extern_src(ind, 0, pot_xs)

        pot_xs = self.gap.potential_xs
        for ind in self._gap_full_dancoff_fsr_inds:
            fullmoc.set_extern_src(ind, 0, pot_xs)

        pot_xs = self.poison_materials[-1].potential_xs
        for ind in self._poison_full_dancoff_fsr_inds:
            fullmoc.set_extern_src(ind, 0, pot_xs)

    def set_full_dancoff_clad_sources(
        self, fullmoc: MOCDriver, moderator: Material, ndl: NDLibrary
    ) -> None:
        """
        Initializes the fixed sources for the full MOC calculation required
        in computing Dancoff corrections. Sources are set for a clad Dancoff
        correction calculation.

        The cladding of a burnable poison pin is no self-shielded. Therefore,
        this method is an alias to set_isolated_dancoff_fuel_sources.

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
        self.set_full_dancoff_fuel_sources(fullmoc, moderator)

    # ==========================================================================
    # Transport Calculation Related Methods

    def set_center_xs(self, ndl: NDLibrary) -> None:
        """
        Constructs the CrossSection object for the material at the center of
        the poison rod, so long as that material is not moderator.

        Parameters
        ----------
        ndl : NDLibrary
            Nuclear data library to use for cross sections.
        """
        if self.center is None:
            return

        if self._center_xs is None:
            self._center_xs = self.center.dilution_xs([1.0e10] * self.gap.size, ndl)
        else:
            self._center_xs.set(self.center.dilution_xs([1.0e10] * self.gap.size, ndl))

        if self._center_xs.name == "":
            self._center_xs.name = "BPR Center"

    def set_gap_xs(self, ndl: NDLibrary) -> None:
        """
        Constructs the CrossSection object for the gap between the poison and
        the cladding.

        Parameters
        ----------
        ndl : NDLibrary
            Nuclear data library to use for cross sections.
        """
        if self._gap_xs is None:
            self._gap_xs = self.gap.dilution_xs([1.0e10] * self.gap.size, ndl)
        else:
            self._gap_xs.set(self.gap.dilution_xs([1.0e10] * self.gap.size, ndl))

        if self._gap_xs.name == "":
            self._gap_xs.name = "BPR Gap"

    def set_clad_xs(self, ndl: NDLibrary) -> None:
        """
        Constructs the CrossSection object for the cladding of the poison rod.

        The cladding of a poison rod uses infinitely dilute cross sections.

        Parameters
        ----------
        ndl : NDLibrary
            Nuclear data library to use for cross sections.
        """
        if self._clad_xs is None:
            self._clad_xs = self.clad.dilution_xs([1.0e10] * self.clad.size, ndl)
        else:
            self._clad_xs.set(self.clad.dilution_xs([1.0e10] * self.clad.size, ndl))

        if self._clad_xs.name == "":
            self._clad_xs.name = "BPR Clad"

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
        if self._poison_xs is None:
            self._poison_xs = self.poison_materials[t].dilution_xs(
                [1.0e10] * self.poison_materials[t].size, ndl
            )
        else:
            self._poison_xs.set(
                self.poison_materials[t].dilution_xs(
                    [1.0e10] * self.poison_materials[t].size, ndl
                )
            )

        if self._poison_xs.name == "":
            self._poison_xs.name = "BPR Poison"

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
        self._center_fsr_inds = []
        self._clad_fsr_inds = []
        self._gap_fsr_inds = []
        self._poison_fsr_inds = []

        for id in self._center_fsr_ids:
            self._center_fsr_inds.append(moc.get_fsr_indx(id, 0))
        for id in self._clad_fsr_ids:
            self._clad_fsr_inds.append(moc.get_fsr_indx(id, 0))
        for id in self._gap_fsr_ids:
            self._gap_fsr_inds.append(moc.get_fsr_indx(id, 0))
        for id in self._poison_fsr_ids:
            self._poison_fsr_inds.append(moc.get_fsr_indx(id, 0))

    def obtain_poison_flux_spectrum(self, moc: MOCDriver) -> None:
        """
        Computes average flux spectrum in the poison from the MOC simulation.

        Parameters
        ----------
        moc : MOCDriver
            MOC simulation for the full calculations.
        """
        self._poison_flux_spectrum = moc.homogenize_flux_spectrum(self._poison_fsr_inds)

    def normalize_flux_spectrum(self, f) -> None:
        """
        Applies a multiplicative factor to the flux spectra for the poison.
        This permits normalizing the flux to a known assembly power.

        Parameters
        ----------
        f : float
            Normalization factor.
        """
        if f <= 0.0:
            raise ValueError("Normalization factor must be > 0.")

        self._poison_flux_spectrum *= f

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

        # Get the flux and initial material
        flux = self._poison_flux_spectrum
        mat = self._poison_materials[-1]  # Use last available mat !

        # Build depletion matrix and multiply by time step
        dep_matrix = build_depletion_matrix(chain, mat, flux, ndl)
        dep_matrix *= dt

        # At this point, we can clear the xs data from the last material as
        # depletion matrix is now built.
        mat.clear_all_micro_xs_data()

        # Initialize an array with the initial target number densities
        N = np.zeros(dep_matrix.size)
        nuclides = dep_matrix.nuclides
        for i, nuclide in enumerate(nuclides):
            N[i] = mat.atom_density(nuclide)

        # Do the matrix exponential
        dep_matrix.exponential_product(N)

        # Now we can build a new material composition
        new_mat_comp = MaterialComposition()
        for i, nuclide in enumerate(nuclides):
            if N[i] > 0.0:
                new_mat_comp.add_nuclide(nuclide, N[i])

        # Make the new material
        new_mat = Material(new_mat_comp, mat.temperature, ndl)
        self._poison_materials.append(new_mat)

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

        # Get the flux and initial material
        flux = self._poison_flux_spectrum
        mat_pred = self._poison_materials[-1]  # Use last available mat !

        # Build depletion matrix and multiply by time step
        dep_matrix = build_depletion_matrix(chain, mat_pred, flux, ndl)
        dep_matrix *= dt

        # Initialize an array with the initial target number densities
        mat_old = self._poison_materials[-2]  # Go 2 steps back !!
        N = np.zeros(dep_matrix.size)
        nuclides = dep_matrix.nuclides
        for i, nuclide in enumerate(nuclides):
            N[i] = mat_old.atom_density(nuclide)

        # Do the matrix exponential
        dep_matrix.exponential_product(N)

        # Now we can build a new material composition
        new_mat_comp = MaterialComposition()
        for i, nuclide in enumerate(nuclides):
            if N[i] > 0.0:
                new_mat_comp.add_nuclide(nuclide, N[i])

        # Make the new material
        new_mat = Material(new_mat_comp, mat_pred.temperature, ndl)
        self._poison_materials[-1] = new_mat
