from .fuel_pin import FuelPin
from .guide_tube import GuideTube
from .critical_leakage import CriticalLeakage
from ._ensleeve import (
    _ensleeve_quarter,
    _ensleeve_half_top,
    _ensleeve_half_right,
    _ensleeve_full,
)
from .._scarabee import (
    borated_water,
    Material,
    CrossSection,
    NDLibrary,
    DepletionChain,
    PinCellType,
    Vector,
    Direction,
    Cartesian2D,
    MOCDriver,
    BoundaryCondition,
    SimulationMode,
    YamamotoTabuchi6,
    P1CriticalitySpectrum,
    B1CriticalitySpectrum,
    ADF,
    CDF,
    DiffusionCrossSection,
    DiffusionData,
    set_logging_level,
    scarabee_log,
    LogLevel,
)
from enum import Enum
import numpy as np
from typing import Optional, List, Tuple, Union
import copy
from threading import Thread
import matplotlib.pyplot as plt


class Symmetry(Enum):
    """
    Defines the symmetry to use when simulating a PWR fuel assembly.
    """

    Full = 1
    """
    The fuel assembly has no symmetry.
    """

    Half = 2
    """
    The fuel assembly has about either the x or y axis (but not both).
    """

    Quarter = 3
    """
    The fuel assembly has symmetry about the x and y axis.
    """


class PWRAssembly:
    """
    A PWRAssembly instance is responsible for performing all the lattice
    calculations necessary to produce few-group cross sections for a
    single PWR assembly.

    Parameters
    ----------
    shape : tuple of int
        The number of pin cells in the full assembly along x and y.
    pitch : float
        The spacing between fuel pins.
    ndl : NDLibrary
        Nuclear data library used for the calculation.
    cells : list of list of FuelPin or GuideTube
        All of the cells which describe the assembly geometry. Should be
        consistent with the symmetry argument.
    boron_ppm : float
        Moderator boron concentration in parts per million. Default is 800.
    moderator_temp : float
        Moderator temperature in Kelvin. Default is 570.
    moderator_pressure : float
        Moderator pressure in MPa. Default is 15.5.
    symmetry : Symmetry
        Symmetry of the fuel assembly. Default is Symmetry.Full.
    linear_power : float
        Linear power density of the full assembly in kW/cm. This value should
        not be reduced due to symmetry. Default is 42.
    assembly_pitch : optional float
        Spacing between fuel assemblies. If None, assembly pitch is calculated
        from the shape and pitch. Default value is None.
    spacer_grid_width : optional float
        Width of the spacer grid material between pin cells.
        Default value is None.
    spacer_grid : optional Material
        Material defining the spacer grid between pin cells.
    grid_sleeve_width : optional float
        Width of the grid sleeve around the assembly. Default value is None.
    grid_sleeve : optional Material
        Material defining the grid sleeve around the assembly.

    Attributes
    ----------
    shape : tuple of int
        The number of pin cells in the full assembly along x and y.
    pitch : float
        The spacing between fuel pins.
    symmetry : Symmetry
        Symmetry of the fuel assembly.
    assembly_pitch : float
        Spacing between fuel assemblies.
    linear_power : float
        Linear power density of the full assembly in kW/cm.
    initial_heavy_metal_linear_mass : float
        Initial linear density of heavy metal in the full assembly at the
        beginning of life in units of kg/cm.
    boron_ppm : float
        Moderator boron concentration in parts per million.
    moderator_temp : float
        Moderator temperature in Kelvin.
    moderator_pressure : float
        Moderator pressure in MPa.
    moderator : Material
        Material representing the assembly moderator.
    spacer_grid_width : optional float
        Width of the spacer grid material between pin cells.
    spacer_grid : optional Material
        Material defining the spacer grid between pin cells.
    grid_sleeve_width : optional float
        Width of the grid sleeve around the assembly.
    grid_sleeve : optional Material
        Material defining the grid sleeve around the assembly.
    dancoff_moc_track_spacing : float
        Spacing between tracks in the MOC calculations for determining Dancoff
        corrections. Default value is 0.05 cm.
    dancoff_moc_num_angles : int
        Number of azimuthal angles in the MOC calculations for determining
        Dancoff corrections. Default value is 32.
    dancoff_flux_tolerance : float
        Flux convergence tolerance for Dancoff correction calculations. Must be
        in range (0., 1.E-2). Default value is 1.E-5.
    moc_track_spacing : float
        Spacing between tracks in the assembly MOC calculations. Default value
        is 0.05 cm.
    moc_num_angles : int
        Number of azimuthal angles in the assembly MOC calculations. Default
        value is 32.
    flux_tolerance : float
        Flux convergence tolerance for assembly calculations. Must be in range
        (0., 1.E-2). Default value is 1.E-5.
    keff_tolerance : float
        Keff convergence tolerance for assembly calculations. Must be in range
        (0., 1.E-2). Default value is 1.E-5.
    condensation_scheme : list of list of int
        Energy condensation scheme to condense from the group structure of the
        library to the few-groups used in the core solver.
    leakage_model : CriticalLeakage
        Model used to determine the critical leakage flux spectrum, also known
        as the fundamental mode. Default method is homogeneous P1.
    depletion_exposure_steps : optional ndarray
        1D Numpy array of assembly burn-up exposure steps, in units of MWd/kg.
        Default is None.
    depletion_time_steps : optional ndarray
        1D Numpy array of burn-up time steps, in units of days.
        Default is None.
    exposures : ndarray
        1D Numpy array of the total assembly burn-up exposures at which
        material information is available, in units of MWd/kg. Default value
        is an empty array before solve has been called.
    times : ndarray
        1D Numpy array of the total assembly burn-up times at which material
        information is available, in units of days. Default value is an empty
        array before solve has been called.
    keff : float or ndarray
        If depletion was not performed, this is a single float with keff for
        the infinite assembly. If depletion was performed, this is a 1D Numpy
        array for the values of keff at the tabulated burn-up exposures/times.
        Default value is 1 before solve has been called.
    """

    def __init__(
        self,
        shape: Tuple[int, int],
        pitch: float,
        ndl: NDLibrary,
        chain: DepletionChain,
        cells: List[List[Union[FuelPin, GuideTube]]],
        boron_ppm: float = 800.0,
        moderator_temp: float = 570.0,
        moderator_pressure: float = 15.5,
        symmetry: Symmetry = Symmetry.Full,
        linear_power: float = 42.0,
        assembly_pitch: Optional[float] = None,
        spacer_grid_width: Optional[float] = None,
        spacer_grid: Optional[Material] = None,
        grid_sleeve_width: Optional[float] = None,
        grid_sleeve: Optional[Material] = None,
    ):
        self._ndl = ndl
        self._chain = chain
        self._symmetry = symmetry

        if len(shape) != 2:
            raise ValueError("Shape must have 2 entries.")
        if shape[0] <= 0 or shape[1] <= 0:
            raise ValueError("Shape entries must be > 0.")
        if shape[0] != shape[1]:
            raise ValueError("Fuel assemblies must be square.")
        self._shape = shape

        # Compute the number of expected cells along each dimension [y][x]
        # based on the symmetry of the problem being solved.
        expx = self.shape[0]
        expy = self.shape[1]
        if self.symmetry == Symmetry.Half:
            expy = (expy // 2) + (expy % 2)
        elif self.symmetry == Symmetry.Quarter:
            expy = (expy // 2) + (expy % 2)
            expx = (expx // 2) + (expx % 2)
        self._simulated_shape = (expx, expy)

        self._initial_heavy_metal_linear_mass = 0.0

        self._cells_set = False
        self._cells: List[List[Union[FuelPin, GuideTube]]] = [[]]
        self._set_cells(cells)

        if pitch <= 0.5:
            raise ValueError("Pitch must be > 0.5.")
        self._pitch = pitch

        # Compute assembly pitch
        self._assembly_pitch = self.shape[0] * self.pitch
        if assembly_pitch is not None:
            try:
                assembly_pitch = float(assembly_pitch)
            except:
                raise TypeError(
                    "The assembly_pitch argument must be convertible to a float."
                )

            if assembly_pitch <= 0.0:
                raise RuntimeError("Assembly pitch must be > 0.")
            elif assembly_pitch <= self._assembly_pitch:
                raise RuntimeError(
                    "Provided assembly pitch is smaller than that indicated by the pitch and shape."
                )
            self._assembly_pitch = assembly_pitch

        # Check the spacer grid parameters
        self._spacer_grid_width: Optional[float] = None
        self._spacer_grid: Optional[Material] = None
        if spacer_grid is None and spacer_grid_width is not None:
            raise RuntimeError("Spacer grid width is provided but material is not.")
        elif spacer_grid is not None and spacer_grid_width is None:
            raise RuntimeError("Spacer grid material is provided but width is not.")
        else:
            self._spacer_grid_width = spacer_grid_width
            self._spacer_grid = spacer_grid

            if self.spacer_grid_width is not None:
                if self._spacer_grid_width <= 0.0:
                    raise ValueError("Spacer grid width must be > 0.")
                elif self._spacer_grid_width >= self._pitch:
                    raise ValueError("Spacer grid width must be < pitch.")

        # Check the grid sleeve parameters
        self._grid_sleeve_width: Optional[float] = None
        self._grid_sleeve: Optional[Material] = None
        if grid_sleeve is None and grid_sleeve_width is not None:
            raise RuntimeError("Grid sleeve width is provided but material is not.")
        elif grid_sleeve is not None and grid_sleeve_width is None:
            raise RuntimeError("Grid sleeve material is provided but width is not.")
        else:
            self._grid_sleeve_width = grid_sleeve_width
            self._grid_sleeve = grid_sleeve

            if self.grid_sleeve_width is not None:
                if self._grid_sleeve_width <= 0.0:
                    raise ValueError("Grid sleeve width must be > 0.")
                elif self._grid_sleeve_width >= 0.5 * (
                    self._assembly_pitch - self.shape[0] * self.pitch
                ):
                    raise ValueError(
                        "Grid sleeve width must be < the assembly gap width."
                    )

        # Get moderator parameters
        if boron_ppm < 0.0:
            raise ValueError("Boron concentration must be >= 0.")
        self._boron_ppm = boron_ppm

        if moderator_temp <= 0.0:
            raise ValueError("Moderator temperature must be > 0.")
        self._moderator_temp = moderator_temp

        if moderator_pressure <= 0.0:
            raise ValueError("Moderator pressure must be > 0.")
        self._moderator_pressure = moderator_pressure

        if linear_power <= 0.0:
            raise ValueError("Linear power must be > 0.")
        self._linear_power = linear_power

        # Make material for borated water
        self._moderator: Material = borated_water(
            self.boron_ppm, self.moderator_temp, self.moderator_pressure, self._ndl
        )
        self._moderator.name = f"Moderator ({self.boron_ppm} ppm boron)"

        # Set initial boundary conditions
        self._x_min_bc = BoundaryCondition.Periodic
        self._x_max_bc = BoundaryCondition.Periodic
        self._y_min_bc = BoundaryCondition.Periodic
        self._y_max_bc = BoundaryCondition.Periodic

        if self.symmetry != Symmetry.Full:
            self._y_min_bc = BoundaryCondition.Reflective
            self._y_max_bc = BoundaryCondition.Reflective

        if self.symmetry == Symmetry.Quarter:
            self._x_min_bc = BoundaryCondition.Reflective
            self._x_max_bc = BoundaryCondition.Reflective

        # ======================================================================
        # DANCOFF CORRECTION CALCULATION DATA
        # ----------------------------------------------------------------------

        # Make water xs for dancoff calculation
        self._moderator_dancoff_xs: CrossSection = CrossSection(
            np.array([self.moderator.potential_xs]),
            np.array([self.moderator.potential_xs]),
            np.array([[0.0]]),
            "Moderator",
        )

        # Spacer grid and grid sleeve dancoff cross sections
        self._spacer_grid_dancoff_xs: Optional[CrossSection] = None
        self._grid_sleeve_dancoff_xs: Optional[CrossSection] = None
        if self.spacer_grid is not None:
            self._spacer_grid_dancoff_xs: CrossSection = CrossSection(
                np.array([self.spacer_grid.potential_xs]),
                np.array([self.spacer_grid.potential_xs]),
                np.array([[0.0]]),
                "Spacer Grid",
            )
        if self.grid_sleeve is not None:
            self._grid_sleeve_dancoff_xs: CrossSection = CrossSection(
                np.array([self.grid_sleeve.potential_xs]),
                np.array([self.grid_sleeve.potential_xs]),
                np.array([[0.0]]),
                "Grid Sleeve",
            )

        # Isolated cell geometry for Dancoff correction calculations
        self._isolated_dancoff_cells = []
        self._isolated_dancoff_mocs = []

        # Full geometry for Dancoff correction calculations
        self._full_dancoff_cells = []
        self._full_dancoff_geom = None
        self._full_dancoff_moc = None

        # Dancoff correction parameters
        self._dancoff_isolation_scale = 10.0
        self._dancoff_moc_track_spacing = 0.05
        self._dancoff_moc_num_angles = 32
        self._dancoff_flux_tolerance = 1.0e-5

        self._fuel_dancoff_corrections = np.zeros(
            (self._simulated_shape[1], self._simulated_shape[0])
        )
        self._clad_dancoff_corrections = np.zeros(
            (self._simulated_shape[1], self._simulated_shape[0])
        )

        # ======================================================================
        # TRANSPORT CALCULATION DATA
        # ----------------------------------------------------------------------

        # Make water xs for transport calculation
        self._moderator_xs: CrossSection = self.moderator.dilution_xs(
            self.moderator.size * [1.0e10], self._ndl
        )

        # Spacer grid and grid sleeve cross sections
        self._spacer_grid_xs: Optional[CrossSection] = None
        self._grid_sleeve_xs: Optional[CrossSection] = None
        if self.spacer_grid is not None:
            self._spacer_grid_xs: CrossSection = self.spacer_grid.dilution_xs(
                self.spacer_grid.size * [1.0e10], self._ndl
            )
        if self.grid_sleeve is not None:
            self._grid_sleeve_xs: CrossSection = self.grid_sleeve.dilution_xs(
                self.grid_sleeve.size * [1.0e10], self._ndl
            )

        self._moc_track_spacing = 0.05
        self._moc_num_angles = 64
        self._flux_tolerance = 1.0e-5
        self._keff_tolerance = 1.0e-5

        self._asmbly_cells = []
        self._asmbly_geom: Optional[Cartesian2D] = None
        self._asmbly_moc: Optional[MOCDriver] = None

        self._leakage_model: CriticalLeakage = CriticalLeakage.P1
        self._infinite_flux_spectrum = (
            None  # To reset to infinite spectrum in MOC driver
        )

        # Depletion time steps in MWd/kg and depletion time steps in days.
        # Initially starts as None (should be provided by user)
        self._depletion_exposure_steps: Optional[np.ndarray] = None
        self._depletion_time_steps: Optional[np.ndarray] = None

        # Arrays for the depletion exposures (MWd/kg) and times (days)
        self._exposures = np.array([])
        self._times = np.array([])

        # Either a single value or list of values (for each depletion step)
        self._keff: Union[float, List[float]] = 1.0

        # Condensation scheme to make few-group cross sections
        self._condensation_scheme: Optional[List[List[int]]] = None

    @property
    def shape(self):
        return self._shape

    @property
    def pitch(self):
        return self._pitch

    @property
    def symmetry(self):
        return self._symmetry

    @property
    def assembly_pitch(self):
        return self._assembly_pitch

    @property
    def spacer_grid_width(self):
        return self._spacer_grid_width

    @property
    def spacer_grid(self):
        return self._spacer_grid

    @property
    def grid_sleeve_width(self):
        return self._grid_sleeve_width

    @property
    def grid_sleeve(self):
        return self._grid_sleeve

    @property
    def linear_power(self):
        return self._linear_power

    @property
    def initial_heavy_metal_linear_mass(self):
        return self._initial_heavy_metal_linear_mass

    @property
    def boron_ppm(self):
        return self._boron_ppm

    @property
    def moderator_temp(self):
        return self._moderator_temp

    @property
    def moderator_pressure(self):
        return self._moderator_pressure

    @property
    def moderator(self):
        return self._moderator

    @property
    def dancoff_moc_track_spacing(self):
        return self._dancoff_moc_track_spacing

    @dancoff_moc_track_spacing.setter
    def dancoff_moc_track_spacing(self, dts: float):
        if dts <= 0.0 or dts > 0.1:
            raise ValueError("Dancoff track spacing must be in range (0, 0.1].")
        self._dancoff_moc_track_spacing = dts

    @property
    def dancoff_moc_num_angles(self):
        return self._dancoff_moc_num_angles

    @dancoff_moc_num_angles.setter
    def dancoff_moc_num_angles(self, dna: int):
        if dna % 4 != 0:
            raise ValueError(
                "Number of angles for Dancoff correction calculation must be a multiple of 4."
            )
        if dna < 4:
            raise ValueError(
                "Number of angles for Dancoff correction calculation must be > 4."
            )
        self._dancoff_moc_num_angles = int(dna)

    @property
    def dancoff_flux_tolerance(self) -> float:
        return self._dancoff_flux_tolerance

    @dancoff_flux_tolerance.setter
    def dancoff_flux_tolerance(self, tol: float) -> None:
        try:
            tol = float(tol)
        except:
            raise TypeError(
                "Dancoff flux tolerance must be a floating point value in (0,1.E-2)."
            )

        if tol <= 0.0 or tol >= 1.0e-2:
            raise ValueError("Dancoff flux tolerance must be in (0,1.E-2).")

        self._dancoff_flux_tolerance = tol

    @property
    def moc_track_spacing(self):
        return self._moc_track_spacing

    @moc_track_spacing.setter
    def moc_track_spacing(self, dts: float):
        if dts <= 0.0 or dts > 0.1:
            raise ValueError("Track spacing must be in range (0, 0.1].")
        self._moc_track_spacing = dts

    @property
    def moc_num_angles(self):
        return self._moc_num_angles

    @moc_num_angles.setter
    def moc_num_angles(self, dna: int):
        if dna % 4 != 0:
            raise ValueError(
                "Number of angles for MOC calculation must be a multiple of 4."
            )
        if dna < 4:
            raise ValueError("Number of angles for MOC calculation must be > 4.")
        self._moc_num_angles = int(dna)

    @property
    def flux_tolerance(self) -> float:
        return self._flux_tolerance

    @flux_tolerance.setter
    def flux_tolerance(self, tol: float) -> None:
        try:
            tol = float(tol)
        except:
            raise TypeError(
                "Flux tolerance must be a floating point value in (0,1.E-2)."
            )

        if tol <= 0.0 or tol >= 1.0e-2:
            raise ValueError("Flux tolerance must be in (0,1.E-2).")

        self._flux_tolerance = tol

    @property
    def keff_tolerance(self) -> float:
        return self._keff_tolerance

    @keff_tolerance.setter
    def keff_tolerance(self, tol: float) -> None:
        try:
            tol = float(tol)
        except:
            raise TypeError(
                "Keff tolerance must be a floating point value in (0,1.E-2)."
            )

        if tol <= 0.0 or tol >= 1.0e-2:
            raise ValueError("Keff tolerance must be in (0,1.E-2).")

        self._keff_tolerance = tol

    @property
    def cells(self):
        return self._cells

    @property
    def leakage_model(self):
        return self._leakage_model

    @leakage_model.setter
    def leakage_model(self, clm: CriticalLeakage):
        self._leakage_model = clm

    @property
    def depletion_exposure_steps(self) -> Optional[np.ndarray]:
        return self._depletion_exposure_steps

    @depletion_exposure_steps.setter
    def depletion_exposure_steps(self, steps: Optional[np.ndarray]) -> None:
        if steps is None:
            self._depletion_exposure_steps = None
            self._depletion_time_steps = None
            return

        if not isinstance(steps, np.ndarray):
            raise TypeError("Depletion exposure steps must be a 1D Numpy array.")
        if steps.ndim != 1:
            raise ValueError("Depletion exposure steps must be a 1D Numpy array.")
        if steps.size == 0:
            raise ValueError(
                "Depletion exposure steps must have at least one entry. Use None for no depletion."
            )
        for step in steps:
            if step <= 0.0:
                raise ValueError("Depletion exposure steps must be > 0.")
            elif step > 5.0:
                raise ValueError("Depletion exposure steps should be <= 5 MWd/kg.")
        self._depletion_exposure_steps = steps

        self._depletion_time_steps = self._depletion_exposure_steps.copy()
        self._depletion_time_steps *= (
            1.0e3 * (1.0 / self.linear_power) * self.initial_heavy_metal_linear_mass
        )

    @property
    def depletion_time_steps(self) -> Optional[np.ndarray]:
        return self._depletion_time_steps

    @depletion_time_steps.setter
    def depletion_time_steps(self, steps: Optional[np.ndarray]) -> None:
        if steps is None:
            self._depletion_time_steps = None
            self._depletion_exposure_steps = None
            return

        if not isinstance(steps, np.ndarray):
            raise TypeError("Depletion time steps must be a 1D Numpy array.")
        if steps.ndim != 1:
            raise ValueError("Depletion time steps must be a 1D Numpy array.")
        if steps.size == 0:
            raise ValueError(
                "Depletion time steps must have at least one entry. Use None for no depletion."
            )
        for step in steps:
            if step <= 0.0:
                raise ValueError("Depletion time steps must be > 0.")
            elif step > 100:
                raise ValueError("Depletion time steps should be <= 100 days.")
        self._depletion_time_steps = steps

        self._depletion_exposure_steps = self._depletion_time_steps.copy()
        self._depletion_exposure_steps /= (
            1.0e3 * (1.0 / self.linear_power) * self.initial_heavy_metal_linear_mass
        )

    @property
    def exposures(self):
        return self._exposures

    @property
    def times(self):
        return self._times

    @property
    def keff(self):
        return self._keff

    @property
    def condensation_scheme(self):
        return self._condensation_scheme

    @condensation_scheme.setter
    def condensation_scheme(self, cs: List[List[int]]):
        if not isinstance(cs, list):
            raise TypeError("Condensation scheme must be a list of lists of ints.")

        if len(cs) == 0:
            raise TypeError("Condensation scheme must have at least one group.")

        for G in range(len(cs)):
            if not isinstance(cs[G], list):
                raise TypeError("Condensation scheme must be a list of lists of ints.")

            if len(cs[G]) != 2:
                raise TypeError(
                    "Each entry in condensation scheme must have 2 entries."
                )

            try:
                cs[G][0] = int(cs[G][0])
                cs[G][1] = int(cs[G][1])
            except:
                raise TypeError(
                    f"Microgroup indices for macrogroup {G} are not convertible to ints."
                )

            if cs[G][1] < cs[G][0]:
                raise ValueError(
                    f"The microgroup indices in macrogroup {G} are not ordered."
                )

            if G == 0 and cs[G][0] != 0:
                raise ValueError(
                    "The first microgroup index of the 0 macrogroup must be 0."
                )
            elif G == len(cs) - 1 and cs[G][1] != self._ndl.ngroups - 1:
                NG = self._ndl.ngroups - 1
                raise ValueError(
                    f"The last microgroup index of the last macrogroup must be {NG}."
                )

            if G > 0:
                if cs[G][0] != cs[G - 1][1] + 1:
                    raise ValueError("The condensation scheme is not continuous")

        self._condensation_scheme = copy.deepcopy(cs)

    def _set_cells(self, cells: List[List[Union[FuelPin, GuideTube]]]):
        if len(cells) != self._simulated_shape[1]:
            raise ValueError(
                "Shape along y of cells list does not agree with assembly shape and symmetry."
            )
        for j in range(len(cells)):
            if len(cells[j]) != self._simulated_shape[0]:
                raise ValueError(
                    "Shape along x of cells list does not agree with assembly shape and symmetry."
                )

        # Make deep copies of all cells provided by the user. This makes sure
        # that each duplicate FuelPin or GuideTube instance is unique.
        self._cells = []
        self._cells_set = False
        self._initial_heavy_metal_linear_mass = 0.0
        for j in range(len(cells)):
            self._cells.append([])
            for i in range(len(cells[j])):
                self._cells[-1].append(copy.deepcopy(cells[j][i]))

                if isinstance(cells[j][i], FuelPin):
                    lfm = cells[j][i].initial_fissionable_linear_mass

                    # First check for quarter symmetry and being corner pin
                    if (
                        self.symmetry == Symmetry.Quarter
                        and self.shape[0] % 2 == 1
                        and self.shape[1] % 2 == 1
                        and j == self._simulated_shape[1] - 1
                        and i == 0
                    ):
                        lfm *= 0.25
                    # Next, check for being on the side with a half pin in quarter symmetry
                    elif (
                        self.symmetry == Symmetry.Quarter
                        and self.shape[0] % 2 == 1
                        and i == 0
                    ):
                        lfm *= 0.5
                    # Next, check for being on the bottom row with a half pin
                    elif (
                        self.symmetry != Symmetry.Full
                        and self.shape[1] % 2 == 1
                        and j == self._simulated_shape[1] - 1
                    ):
                        lfm *= 0.5

                    self._initial_heavy_metal_linear_mass += lfm

        self._cells_set = True

        if self.symmetry == Symmetry.Half:
            self._initial_heavy_metal_linear_mass *= 2.0
        elif self.symmetry == Symmetry.Quarter:
            self._initial_heavy_metal_linear_mass *= 4.0

        # Convert HM mass from g to kg
        self._initial_heavy_metal_linear_mass *= 1.0e-3

    # ==========================================================================
    # Dancoff Correction Related Methods

    def _init_dancoff_components(self) -> None:
        scarabee_log(
            LogLevel.Info, "Initializing Dancoff correction calculation components."
        )
        set_logging_level(LogLevel.Warning)
        self._init_isolated_dancoff_components()
        self._init_full_dancoff_components()
        self._save_dancoff_fsr_indexes()
        set_logging_level(LogLevel.Info)

    def _init_isolated_dancoff_components(self) -> None:
        # Isolated pitch
        iso_pitch = self._dancoff_isolation_scale * self.pitch

        for j in range(len(self.cells)):
            self._isolated_dancoff_cells.append([])
            self._isolated_dancoff_mocs.append([])
            for i in range(len(self.cells[j])):
                cell = None
                geom = None
                moc = None

                x_min_bc = BoundaryCondition.Vacuum
                x_max_bc = BoundaryCondition.Vacuum
                y_min_bc = BoundaryCondition.Vacuum
                y_max_bc = BoundaryCondition.Vacuum

                # First check for quarter symmetry and being corner pin
                if (
                    self.symmetry == Symmetry.Quarter
                    and self.shape[0] % 2 == 1
                    and self.shape[1] % 2 == 1
                    and j == self._simulated_shape[1] - 1
                    and i == 0
                ):
                    cell = self.cells[j][i].make_dancoff_moc_cell(
                        self._moderator_dancoff_xs,
                        0.5 * iso_pitch,
                        0.5 * iso_pitch,
                        PinCellType.I,
                        True,
                    )
                    geom = Cartesian2D([0.5 * iso_pitch], [0.5 * iso_pitch])
                    geom.set_tiles([cell])
                    x_min_bc = BoundaryCondition.Reflective
                    y_min_bc = BoundaryCondition.Reflective
                # Next, check for being on the side with a half pin in quarter symmetry
                elif (
                    self.symmetry == Symmetry.Quarter
                    and self.shape[0] % 2 == 1
                    and i == 0
                ):
                    cell = self.cells[j][i].make_dancoff_moc_cell(
                        self._moderator_dancoff_xs,
                        0.5 * iso_pitch,
                        iso_pitch,
                        PinCellType.XP,
                        True,
                    )
                    geom = Cartesian2D([0.5 * iso_pitch], [iso_pitch])
                    geom.set_tiles([cell])
                    x_min_bc = BoundaryCondition.Reflective
                # Next, check for being on the bottom row with a half pin
                elif (
                    self.symmetry != Symmetry.Full
                    and self.shape[1] % 2 == 1
                    and j == self._simulated_shape[1] - 1
                ):
                    cell = self.cells[j][i].make_dancoff_moc_cell(
                        self._moderator_dancoff_xs,
                        iso_pitch,
                        0.5 * iso_pitch,
                        PinCellType.YP,
                        True,
                    )
                    geom = Cartesian2D([iso_pitch], [0.5 * iso_pitch])
                    geom.set_tiles([cell])
                    y_min_bc = BoundaryCondition.Reflective
                # Otherwise, we just make the full cell
                else:
                    cell = self.cells[j][i].make_dancoff_moc_cell(
                        self._moderator_dancoff_xs,
                        iso_pitch,
                        iso_pitch,
                        PinCellType.Full,
                        True,
                    )
                    geom = Cartesian2D([iso_pitch], [iso_pitch])
                    geom.set_tiles([cell])

                # Save to cells
                self._isolated_dancoff_cells[-1].append(cell)

                # Make the MOCDriver
                moc = MOCDriver(geom)
                moc.sim_mode = SimulationMode.FixedSource
                moc.x_min_bc = x_min_bc
                moc.x_max_bc = x_max_bc
                moc.y_min_bc = y_min_bc
                moc.y_max_bc = y_max_bc

                # Generate tracks in serial as each call will run with threads
                moc.generate_tracks(
                    self._dancoff_moc_num_angles,
                    self._dancoff_moc_track_spacing,
                    YamamotoTabuchi6(),
                )

                # Save the MOC
                self._isolated_dancoff_mocs[-1].append(moc)

    def _init_full_dancoff_components(self):
        pitch = self.pitch
        if self.spacer_grid_width is not None:
            pitch -= 2.0 * self.spacer_grid_width

        # Get all the cells
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                # Get the geometry bit from the cell
                cell = None

                # First check for quarter symmetry and being corner pin
                if (
                    self.symmetry == Symmetry.Quarter
                    and self.shape[0] % 2 == 1
                    and self.shape[1] % 2 == 1
                    and j == self._simulated_shape[1] - 1
                    and i == 0
                ):
                    cell = self.cells[j][i].make_dancoff_moc_cell(
                        self._moderator_dancoff_xs,
                        0.5 * pitch,
                        0.5 * pitch,
                        PinCellType.I,
                        False,
                    )
                    if self.spacer_grid_width is not None:
                        cell, _ = _ensleeve_quarter(
                            cell,
                            pitch,
                            self.spacer_grid_width,
                            self._spacer_grid_dancoff_xs,
                        )
                # Next, check for being on the side with a half pin in quarter symmetry
                elif (
                    self.symmetry == Symmetry.Quarter
                    and self.shape[0] % 2 == 1
                    and i == 0
                ):
                    cell = self.cells[j][i].make_dancoff_moc_cell(
                        self._moderator_dancoff_xs,
                        0.5 * pitch,
                        pitch,
                        PinCellType.XP,
                        False,
                    )
                    if self.spacer_grid_width is not None:
                        cell, _ = _ensleeve_half_right(
                            cell,
                            pitch,
                            self.spacer_grid_width,
                            self._spacer_grid_dancoff_xs,
                        )
                # Next, check for being on the bottom row with a half pin
                elif (
                    self.symmetry != Symmetry.Full
                    and self.shape[1] % 2 == 1
                    and j == self._simulated_shape[1] - 1
                ):
                    cell = self.cells[j][i].make_dancoff_moc_cell(
                        self._moderator_dancoff_xs,
                        pitch,
                        0.5 * pitch,
                        PinCellType.YP,
                        False,
                    )
                    if self.spacer_grid_width is not None:
                        cell, _ = _ensleeve_half_top(
                            cell,
                            pitch,
                            self.spacer_grid_width,
                            self._spacer_grid_dancoff_xs,
                        )
                # Otherwise, we just make the full cell
                else:
                    cell = self.cells[j][i].make_dancoff_moc_cell(
                        self._moderator_dancoff_xs,
                        pitch,
                        pitch,
                        PinCellType.Full,
                        False,
                    )
                    if self.spacer_grid_width is not None:
                        cell, _ = _ensleeve_full(
                            cell,
                            pitch,
                            self.spacer_grid_width,
                            self._spacer_grid_dancoff_xs,
                        )

                # Save to cells
                self._full_dancoff_cells.append(cell)

        # Construct the Cartesian2D geometry
        dx = self._simulated_shape[0] * [self.pitch]
        dy = self._simulated_shape[1] * [self.pitch]
        if self.symmetry != Symmetry.Full and self.shape[1] % 2 == 1:
            dy[0] *= 0.5
        if self.symmetry == Symmetry.Quarter and self.shape[0] % 2 == 1:
            dx[0] *= 0.5
        self._full_dancoff_geom = Cartesian2D(dx, dy)
        self._full_dancoff_geom.set_tiles(self._full_dancoff_cells)

        # Add gap if needed
        gap_width = 0.5 * (self.assembly_pitch - self.shape[0] * self.pitch)
        if self.grid_sleeve_width is not None:
            gap_width -= self.grid_sleeve_width
        if gap_width > 0.0:
            pins_geom = self._full_dancoff_geom

            if self.symmetry == Symmetry.Quarter:
                if self.grid_sleeve_width is not None:
                    pins_geom, _ = _ensleeve_quarter(
                        pins_geom,
                        self.pitch,
                        self._grid_sleeve_width,
                        self._grid_sleeve_dancoff_xs,
                    )
                self._full_dancoff_geom, _ = _ensleeve_quarter(
                    pins_geom, self.pitch, gap_width, self._moderator_dancoff_xs
                )
            elif self.symmetry == Symmetry.Half:
                if self.grid_sleeve_width is not None:
                    pins_geom, _ = _ensleeve_half_top(
                        pins_geom,
                        self.pitch,
                        self._grid_sleeve_width,
                        self._grid_sleeve_dancoff_xs,
                    )
                self._full_dancoff_geom, _ = _ensleeve_half_top(
                    pins_geom, self.pitch, gap_width, self._moderator_dancoff_xs
                )
            else:
                if self.grid_sleeve_width is not None:
                    pins_geom, _ = _ensleeve_full(
                        pins_geom,
                        self.pitch,
                        self._grid_sleeve_width,
                        self._grid_sleeve_dancoff_xs,
                    )
                self._full_dancoff_geom, _ = _ensleeve_full(
                    pins_geom, self.pitch, gap_width, self._moderator_dancoff_xs
                )

        # Construct the MOC
        self._full_dancoff_moc = MOCDriver(self._full_dancoff_geom)
        self._full_dancoff_moc.sim_mode = SimulationMode.FixedSource
        self._full_dancoff_moc.x_min_bc = self._x_min_bc
        self._full_dancoff_moc.x_max_bc = self._x_max_bc
        self._full_dancoff_moc.y_min_bc = self._y_min_bc
        self._full_dancoff_moc.y_max_bc = self._y_max_bc
        self._full_dancoff_moc.generate_tracks(
            self._dancoff_moc_num_angles,
            self._dancoff_moc_track_spacing,
            YamamotoTabuchi6(),
        )

    def _save_dancoff_fsr_indexes(self):
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i]
                isomoc = self._isolated_dancoff_mocs[j][i]
                cell.populate_dancoff_fsr_indexes(isomoc, self._full_dancoff_moc)

    def _dancoff_components_initialized(self) -> bool:
        if self._full_dancoff_moc is None:
            return False
        return True

    def set_dancoff_moderator_xs(self) -> None:
        """
        Updates the moderator cross section for all Dancoff correction calculations.
        """
        self._moderator_dancoff_xs.set(
            CrossSection(
                np.array([self.moderator.potential_xs]),
                np.array([self.moderator.potential_xs]),
                np.array([[0.0]]),
                "Moderator",
            )
        )

    def set_dancoff_spacer_grid_sleeve_xs(self) -> None:
        """
        Updates the spacer grid and grid sleeve cross sections for all Dancoff
        correction calculations.
        """
        if self.spacer_grid is not None:
            self._spacer_grid_dancoff_xs.set(
                CrossSection(
                    np.array([self.spacer_grid.potential_xs]),
                    np.array([self.spacer_grid.potential_xs]),
                    np.array([[0.0]]),
                    "Spacer Grid",
                )
            )
        if self.grid_sleeve is not None:
            self._grid_sleeve_dancoff_xs.set(
                CrossSection(
                    np.array([self.grid_sleeve.potential_xs]),
                    np.array([self.grid_sleeve.potential_xs]),
                    np.array([[0.0]]),
                    "Grid Sleeve",
                )
            )

    def compute_fuel_dancoff_corrections(self):
        """
        Recomputes all Dancoff corrections for the fuel regions in the problem,
        using the most recent material definitions. All fuel is shelf-shielded
        together, regardless of wether or not is is UO2 or MOX.
        """
        scarabee_log(LogLevel.Info, "Computing Dancoff corrections for the fuel.")
        set_logging_level(LogLevel.Warning)
        if not self._dancoff_components_initialized():
            raise RuntimeError(
                "Dancoff calculation components have not been initialized."
            )

        # Set the xs and sources for all cells
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i]
                isomoc = self._isolated_dancoff_mocs[j][i]

                isomoc.flux_tolerance = self.dancoff_flux_tolerance

                cell.set_xs_for_fuel_dancoff_calculation()

                cell.set_isolated_dancoff_fuel_sources(isomoc, self.moderator)

                cell.set_full_dancoff_fuel_sources(
                    self._full_dancoff_moc, self.moderator
                )
        self._full_dancoff_moc.flux_tolerance = self.dancoff_flux_tolerance

        # Solve all the MOCs in parallel
        threads = []
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i]
                if isinstance(cell, FuelPin):
                    isomoc = self._isolated_dancoff_mocs[j][i]
                    threads.append(Thread(target=isomoc.solve))
                    threads[-1].start()
        threads.append(Thread(target=self._full_dancoff_moc.solve))
        threads[-1].start()
        for t in threads:
            t.join()

        # Go through and let each cell compute the Dancoff correction if it holds
        # a fuel pin.
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i]
                isomoc = self._isolated_dancoff_mocs[j][i]

                if isinstance(cell, FuelPin):
                    C = cell.compute_fuel_dancoff_correction(
                        isomoc, self._full_dancoff_moc
                    )
                    self._fuel_dancoff_corrections[j, i] = C
        set_logging_level(LogLevel.Info)

    def compute_clad_dancoff_corrections(self):
        """
        Recomputes all Dancoff corrections for the fuel pin cladding regions
        in the problem, using the most recent material definitions. All
        cladding is shelf-shielded together.
        """
        scarabee_log(LogLevel.Info, "Computing Dancoff corrections for the cladding.")
        set_logging_level(LogLevel.Warning)
        if self._full_dancoff_moc is None:
            raise RuntimeError(
                "Dancoff calculation components have not been initialized."
            )

        # Set the xs and sources for all cells
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i]
                isomoc = self._isolated_dancoff_mocs[j][i]
                isomoc.flux_tolerance = self.dancoff_flux_tolerance

                cell.set_xs_for_clad_dancoff_calculation(self._ndl)

                cell.set_isolated_dancoff_clad_sources(
                    isomoc, self.moderator, self._ndl
                )

                cell.set_full_dancoff_clad_sources(
                    self._full_dancoff_moc, self.moderator, self._ndl
                )
            self._full_dancoff_moc.flux_tolerance = self.dancoff_flux_tolerance

        # Solve all the MOCs in parallel
        threads = []
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i]
                isomoc = self._isolated_dancoff_mocs[j][i]
                threads.append(Thread(target=isomoc.solve))
                threads[-1].start()
        threads.append(Thread(target=self._full_dancoff_moc.solve))
        threads[-1].start()
        for t in threads:
            t.join()

        # Go through and let each cell compute the Dancoff correction if it holds
        # a fuel pin.
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i]
                isomoc = self._isolated_dancoff_mocs[j][i]

                C = cell.compute_clad_dancoff_correction(isomoc, self._full_dancoff_moc)
                self._clad_dancoff_corrections[j, i] = C
        set_logging_level(LogLevel.Info)

    def apply_dancoff_corrections(self):
        """
        Appends all fuel and cladding Dancoff corrections to the appropriate cell.
        """
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i]
                if isinstance(cell, FuelPin):
                    cell.append_fuel_dancoff_correction(
                        self._fuel_dancoff_corrections[j, i]
                    )
                cell.append_clad_dancoff_correction(
                    self._clad_dancoff_corrections[j, i]
                )

    def self_shield_and_xs_update(self):
        """
        Computes a new set of Dancoff corrections for the fuel and the
        cladding.  After, these are applied to all the cells in the problem.
        """
        if not self._dancoff_components_initialized():
            self._init_dancoff_components()

        # Update the Dancoff cross sections held by the assembly
        self.set_dancoff_moderator_xs()
        self.set_dancoff_spacer_grid_sleeve_xs()

        # Compute Dancoff corrections
        self.compute_fuel_dancoff_corrections()
        self.compute_clad_dancoff_corrections()
        self.apply_dancoff_corrections()

    # ==========================================================================
    # Transport Calculation Related Methods

    def _init_moc(self):
        pitch = self.pitch
        if self.spacer_grid_width is not None:
            pitch -= 2.0 * self.spacer_grid_width

        # Get all the cells
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                # Get the geometry bit from the cell
                cell = None

                # First check for quarter symmetry and being corner pin
                if (
                    self.symmetry == Symmetry.Quarter
                    and self.shape[0] % 2 == 1
                    and self.shape[1] % 2 == 1
                    and j == self._simulated_shape[1] - 1
                    and i == 0
                ):
                    cell = self.cells[j][i].make_moc_cell(
                        self._moderator_xs,
                        0.5 * pitch,
                        0.5 * pitch,
                        PinCellType.I,
                    )
                    if self.spacer_grid_width is not None:
                        cell, _ = _ensleeve_quarter(
                            cell, pitch, self.spacer_grid_width, self._spacer_grid_xs
                        )
                # Next, check for being on the side with a half pin in quarter symmetry
                elif (
                    self.symmetry == Symmetry.Quarter
                    and self.shape[0] % 2 == 1
                    and i == 0
                ):
                    cell = self.cells[j][i].make_moc_cell(
                        self._moderator_xs, 0.5 * pitch, pitch, PinCellType.XP
                    )
                    if self.spacer_grid_width is not None:
                        cell, _ = _ensleeve_half_right(
                            cell, pitch, self.spacer_grid_width, self._spacer_grid_xs
                        )
                # Next, check for being on the bottom row with a half pin
                elif (
                    self.symmetry != Symmetry.Full
                    and self.shape[1] % 2 == 1
                    and j == self._simulated_shape[1] - 1
                ):
                    cell = self.cells[j][i].make_moc_cell(
                        self._moderator_xs, pitch, 0.5 * pitch, PinCellType.YP
                    )
                    if self.spacer_grid_width is not None:
                        cell, _ = _ensleeve_half_top(
                            cell, pitch, self.spacer_grid_width, self._spacer_grid_xs
                        )
                # Otherwise, we just make the full cell
                else:
                    cell = self.cells[j][i].make_moc_cell(
                        self._moderator_xs, pitch, pitch, PinCellType.Full
                    )
                    if self.spacer_grid_width is not None:
                        cell, _ = _ensleeve_full(
                            cell, pitch, self.spacer_grid_width, self._spacer_grid_xs
                        )

                # Save to cells
                self._asmbly_cells.append(cell)

        # Construct the Cartesian2D geometry
        dx = self._simulated_shape[0] * [self.pitch]
        dy = self._simulated_shape[1] * [self.pitch]
        if self.symmetry != Symmetry.Full and self.shape[1] % 2 == 1:
            dy[0] *= 0.5
        if self.symmetry == Symmetry.Quarter and self.shape[0] % 2 == 1:
            dx[0] *= 0.5
        self._asmbly_geom = Cartesian2D(dx, dy)
        self._asmbly_geom.set_tiles(self._asmbly_cells)

        # Add gap if needed
        gap_width = 0.5 * (self.assembly_pitch - self.shape[0] * self.pitch)
        if self.grid_sleeve_width is not None:
            gap_width -= self.grid_sleeve_width
        if gap_width > 0.0:
            pins_geom = self._asmbly_geom

            if self.symmetry == Symmetry.Quarter:
                if self.grid_sleeve_width is not None:
                    pins_geom, _ = _ensleeve_quarter(
                        pins_geom,
                        self.pitch,
                        self._grid_sleeve_width,
                        self._grid_sleeve_xs,
                    )
                self._asmbly_geom, _ = _ensleeve_quarter(
                    pins_geom, self.pitch, gap_width, self._moderator_xs
                )
            elif self.symmetry == Symmetry.Half:
                if self.grid_sleeve_width is not None:
                    pins_geom, _ = _ensleeve_half_top(
                        pins_geom,
                        self.pitch,
                        self._grid_sleeve_width,
                        self._grid_sleeve_xs,
                    )
                self._asmbly_geom, _ = _ensleeve_half_top(
                    pins_geom, self.pitch, gap_width, self._moderator_xs
                )
            else:
                if self.grid_sleeve_width is not None:
                    pins_geom, _ = _ensleeve_full(
                        pins_geom,
                        self.pitch,
                        self._grid_sleeve_width,
                        self._grid_sleeve_xs,
                    )
                self._asmbly_geom, _ = _ensleeve_full(
                    pins_geom, self.pitch, gap_width, self._moderator_xs
                )

        # Construct the MOC
        self._asmbly_moc = MOCDriver(self._asmbly_geom)
        self._asmbly_moc.x_min_bc = self._x_min_bc
        self._asmbly_moc.x_max_bc = self._x_max_bc
        self._asmbly_moc.y_min_bc = self._y_min_bc
        self._asmbly_moc.y_max_bc = self._y_max_bc
        self._asmbly_moc.generate_tracks(
            self._moc_num_angles,
            self._moc_track_spacing,
            YamamotoTabuchi6(),
        )

        self._save_fsr_indexes()

    def _save_fsr_indexes(self):
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i]
                cell.populate_fsr_indexes(self._asmbly_moc)

    def plot(self) -> None:
        """
        Launches the graphical geometry plotter for the assembly calculation.
        """
        if self._asmbly_moc is None:
            raise RuntimeError(
                "Cannot launch plotter. MOCDriver has not been initialized."
            )
        self._asmbly_moc.plot()

    def set_moderator_xs(self) -> None:
        """
        Updates the moderator cross section for transport calculations.
        """
        self._moderator_xs.set(
            self.moderator.dilution_xs(self.moderator.size * [1.0e10], self._ndl)
        )

    def set_spacer_grid_sleeve_xs(self) -> None:
        """
        Updates the spacer grid and grid sleeve cross sections for transport
        calculations.
        """
        if self.spacer_grid is not None:
            self._spacer_grid_xs.set(
                self.spacer_grid.dilution_xs(
                    self.spacer_grid.size * [1.0e10], self._ndl
                )
            )
        if self.grid_sleeve is not None:
            self._grid_sleeve_xs.set(
                self.grid_sleeve.dilution_xs(
                    self.grid_sleeve.size * [1.0e10], self._ndl
                )
            )

    def recompute_all_xs(self) -> None:
        """
        Computes and applies all cross sections using the most recent material
        information and Dancoff corrections.
        """
        self.recompute_all_fuel_xs()
        self.recompute_all_gap_xs()
        self.recompute_all_clad_xs()

        self.recompute_all_guide_tube_fill_xs()

        self.set_moderator_xs()
        self.set_spacer_grid_sleeve_xs()

    def recompute_all_self_shielded_xs(self) -> None:
        """
        Computes and applies all fuel and caldding cross sections using the
        most recent material information and Dancoff corrections.
        """
        self.recompute_all_fuel_xs()
        self.recompute_all_clad_xs()

    def recompute_all_fuel_xs(self) -> None:
        """
        Computes and applies all fuel cross sections using the most recent
        material information and Dancoff corrections.
        """
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i]
                if isinstance(cell, FuelPin):
                    cell.set_fuel_xs_for_depletion_step(-1, self._ndl)

    def recompute_all_clad_xs(self) -> None:
        """
        Computes and applies all cladding cross sections using the most recent
        material information and Dancoff corrections.
        """
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i]
                cell.set_clad_xs_for_depletion_step(-1, self._ndl)

    def recompute_all_gap_xs(self) -> None:
        """
        Computes and applies all gap cross sections using the most recent
        material information.
        """
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i]
                if isinstance(cell, FuelPin):
                    cell.set_gap_xs(self._ndl)

    def recompute_all_guide_tube_fill_xs(self) -> None:
        """
        Computes and applies all cross sections for the fill objects of guide
        tubes. These could be for burnable poison rods or control rods.
        """
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i]
                if isinstance(cell, GuideTube):
                    cell.set_fill_xs_for_depletion_step(-1, self._ndl)

    def apply_leakage_model(self) -> None:
        """
        Applied the critical leakage model to the assembly, modifying the flux
        in the MOC simulation directly.
        """
        # If no leakage, just return
        if self.leakage_model == CriticalLeakage.NoLeakage:
            return

        scarabee_log(LogLevel.Info, "")

        self._infinite_flux_spectrum = self._asmbly_moc.homogenize_flux_spectrum()

        homogenized_moc = self._asmbly_moc.homogenize()

        if self.leakage_model == CriticalLeakage.P1:
            scarabee_log(
                LogLevel.Info, "Performing P1 criticality spectrum calculation"
            )
            critical_spectrum = P1CriticalitySpectrum(homogenized_moc)
        else:
            scarabee_log(
                LogLevel.Info, "Performing B1 criticality spectrum calculation"
            )
            critical_spectrum = B1CriticalitySpectrum(homogenized_moc)

        self._asmbly_moc.apply_criticality_spectrum(critical_spectrum.flux)

        scarabee_log(LogLevel.Info, "Kinf    : {:.5f}".format(critical_spectrum.k_inf))
        scarabee_log(
            LogLevel.Info, "Buckling: {:.5f}".format(critical_spectrum.buckling)
        )

    def apply_infinite_spectrum(self) -> None:
        """
        Undoes the critical flux spectrum adjustment to the MOCDriver.
        This permits subsequent transport calcualtions to converge much faster.
        """
        if self._infinite_flux_spectrum is not None and self._asmbly_moc is not None:
            self._asmbly_moc.apply_criticality_spectrum(self._infinite_flux_spectrum)

    def obtain_flux_spectra(self) -> None:
        """
        Computes the average flux spectrum for material regions whcih are
        depleted. This includes each fuel ring in fuel pins and the poison in
        burnable poison rods.
        """
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                self.cells[j][i].obtain_flux_spectra(self._asmbly_moc)

    def normalize_flux_to_power(self) -> None:
        """
        Normalizes the flux spectra based on the specified linear power density
        for the assembly. It assumes that all power comes from fission.
        """
        assembly_linear_power = 0.0
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i]

                if not isinstance(cell, FuelPin):
                    continue

                pin_linear_power = cell.compute_pin_linear_power(self._ndl)

                # First check for quarter symmetry and being corner pin
                if (
                    self.symmetry == Symmetry.Quarter
                    and self.shape[0] % 2 == 1
                    and self.shape[1] % 2 == 1
                    and j == self._simulated_shape[1] - 1
                    and i == 0
                ):
                    pin_linear_power *= 0.25
                # Next, check for being on the side with a half pin in quarter symmetry
                elif (
                    self.symmetry == Symmetry.Quarter
                    and self.shape[0] % 2 == 1
                    and i == 0
                ):
                    pin_linear_power *= 0.5
                # Next, check for being on the bottom row with a half pin
                elif (
                    self.symmetry != Symmetry.Full
                    and self.shape[1] % 2 == 1
                    and j == self._simulated_shape[1] - 1
                ):
                    pin_linear_power *= 0.5

                assembly_linear_power += pin_linear_power

        if self.symmetry == Symmetry.Half:
            assembly_linear_power *= 2.0
        elif self.symmetry == Symmetry.Quarter:
            assembly_linear_power *= 4.0

        # Compute normalization factor.
        # Multiply by 10^3 to go from kW to W on linear power
        f = 1.0e3 * self.linear_power / assembly_linear_power

        # Normalize flux spectra
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i].normalize_flux_spectrum(f)

    def _compute_few_group_flux(self, r: Vector, u: Direction) -> List[float]:
        """
        Computes the flux in the few-group scheme at the given position
        and direction.

        Parameters
        ----------
        r : Vector
            Position where the flux is evaluated.
        u : Direction
            Direction vector to disambiguate the position.

        Returns
        -------
        list of float
            Values of the few-group flux at the given position.

        Raises
        ------
        RuntimeError
            If the condensation_scheme attribute is not set.
        """
        if self.condensation_scheme is None:
            raise RuntimeError("Energy condensation scheme not set.")

        flux = [0.0 for G in range(len(self.condensation_scheme))]

        for G in range(len(self.condensation_scheme)):
            gmin, gmax = self.condensation_scheme[G][:]

            for g in range(gmin, gmax + 1):
                flux[G] += self._asmbly_moc.flux(r, u, g)

        return flux

    def _compute_average_line_flux(
        self, segments: List[Tuple[int, float]]
    ) -> List[float]:
        """
        Computes the average flux along a set of line segments in the few-group
        structure.

        Paramters
        ---------
        segments : list of tuples of int and float
            List of flat source region index and segment length tuples.

        Returns
        -------
        list of float
            Values of the few-group flux alone the line.

        Raises
        ------
        RuntimeError
            If the condensation_scheme attribute is not set.
        """
        if self.condensation_scheme is None:
            raise RuntimeError("Energy condensation scheme not set.")

        total_length = 0.0
        for s in segments:
            total_length += s[1]
        invs_tot_length = 1.0 / total_length

        flux = [0.0 for G in range(len(self.condensation_scheme))]

        for G in range(len(self.condensation_scheme)):
            gmin, gmax = self.condensation_scheme[G][:]
            for g in range(gmin, gmax + 1):
                for s in segments:
                    flux[G] += s[1] * self._asmbly_moc.flux(s[0], g)

        for G in range(len(self.condensation_scheme)):
            flux[G] *= invs_tot_length

        return flux

    def _compute_adf_cdf(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Computes the assembly and corner discontinuity factors in the few-group
        structure.

        Returns
        -------
        ADF : ndarray
            Assembly Discontinuity Factors
        CDF : ndarray
            Corner Discontinuity Factors

        Raises
        ------
        RuntimeError
            If the condensation_scheme attribute is not set.
        """
        if self.condensation_scheme is None:
            raise RuntimeError("Energy condensation scheme not set.")

        NG = len(self.condensation_scheme)
        moc = self._asmbly_moc

        # First, compute the homogeneous flux
        homog_flux = [0.0 for G in range(NG)]
        total_volume = 0.0
        for i in range(moc.nfsr):
            Vi = moc.volume(i)
            total_volume += Vi

            for G in range(NG):
                gmin, gmax = self.condensation_scheme[G][:]
                for g in range(gmin, gmax + 1):
                    homog_flux[G] += Vi * moc.flux(i, g)
        for G in range(NG):
            homog_flux[G] /= total_volume

        # Create empty arrays for ADFs and CDFs
        adf = np.ones((NG, 6))
        cdf = np.zeros((NG, 4))

        if self.symmetry == Symmetry.Full:
            # Get flux along surfaces
            xn_segments = moc.trace_fsr_segments(
                Vector(moc.x_min + 0.001, moc.y_max), Direction(0.0, -1.0)
            )
            xp_segments = moc.trace_fsr_segments(
                Vector(moc.x_max - 0.001, moc.y_max), Direction(0.0, -1.0)
            )
            yn_segments = moc.trace_fsr_segments(
                Vector(moc.x_min, moc.y_min + 0.001), Direction(1.0, 0.0)
            )
            yp_segments = moc.trace_fsr_segments(
                Vector(moc.x_min, moc.y_max - 0.001), Direction(1.0, 0.0)
            )
            xn_flx = self._compute_average_line_flux(xn_segments)
            xp_flx = self._compute_average_line_flux(xp_segments)
            yn_flx = self._compute_average_line_flux(yn_segments)
            yp_flx = self._compute_average_line_flux(yp_segments)

            # Get flux at corners
            I_flx = self._compute_few_group_flux(
                Vector(moc.x_max - 0.001, moc.y_max - 0.001), Direction(-1.0, -1.0)
            )
            II_flx = self._compute_few_group_flux(
                Vector(moc.x_min + 0.001, moc.y_max - 0.001), Direction(1.0, -1.0)
            )
            III_flx = self._compute_few_group_flux(
                Vector(moc.x_min + 0.001, moc.y_min + 0.001), Direction(1.0, 1.0)
            )
            IV_flx = self._compute_few_group_flux(
                Vector(moc.x_max - 0.001, moc.y_min + 0.001), Direction(-1.0, 1.0)
            )

            for G in range(NG):
                adf[G, ADF.XN] = xn_flx[G] / homog_flux[G]
                adf[G, ADF.XP] = xp_flx[G] / homog_flux[G]
                adf[G, ADF.YN] = yn_flx[G] / homog_flux[G]
                adf[G, ADF.YP] = yp_flx[G] / homog_flux[G]
                # The ADFs on the +/- Z sides will be left at unity

                cdf[G, CDF.I] = I_flx[G] / homog_flux[G]
                cdf[G, CDF.II] = II_flx[G] / homog_flux[G]
                cdf[G, CDF.III] = III_flx[G] / homog_flux[G]
                cdf[G, CDF.IV] = IV_flx[G] / homog_flux[G]

        elif self.symmetry == Symmetry.Half:
            # Get flux along surfaces
            xn_segments = moc.trace_fsr_segments(
                Vector(moc.x_min + 0.001, moc.y_max), Direction(0.0, -1.0)
            )
            xp_segments = moc.trace_fsr_segments(
                Vector(moc.x_max - 0.001, moc.y_max), Direction(0.0, -1.0)
            )
            yp_segments = moc.trace_fsr_segments(
                Vector(moc.x_min, moc.y_max - 0.001), Direction(1.0, 0.0)
            )
            xn_flx = self._compute_average_line_flux(xn_segments)
            xp_flx = self._compute_average_line_flux(xp_segments)
            yp_flx = self._compute_average_line_flux(yp_segments)

            # Get flux at corners
            I_flx = self._compute_few_group_flux(
                Vector(moc.x_max - 0.001, moc.y_max - 0.001), Direction(-1.0, -1.0)
            )
            II_flx = self._compute_few_group_flux(
                Vector(moc.x_min + 0.001, moc.y_max - 0.001), Direction(1.0, -1.0)
            )

            for G in range(NG):
                adf[G, ADF.XN] = xn_flx[G] / homog_flux[G]
                adf[G, ADF.XP] = xp_flx[G] / homog_flux[G]
                adf[G, ADF.YP] = yp_flx[G] / homog_flux[G]
                adf[G, ADF.YN] = adf[G, ADF.YP]
                # The ADFs on the +/- Z sides will be left at unity

                cdf[G, CDF.I] = I_flx[G] / homog_flux[G]
                cdf[G, CDF.II] = II_flx[G] / homog_flux[G]
                cdf[G, CDF.III] = cdf[G, CDF.II]
                cdf[G, CDF.IV] = cdf[G, CDF.I]

        else:  # Quarter symmetry
            # Get flux along surfaces
            xp_segments = moc.trace_fsr_segments(
                Vector(moc.x_max - 0.001, moc.y_max), Direction(0.0, -1.0)
            )
            yp_segments = moc.trace_fsr_segments(
                Vector(moc.x_min, moc.y_max - 0.001), Direction(1.0, 0.0)
            )
            xp_flx = self._compute_average_line_flux(xp_segments)
            yp_flx = self._compute_average_line_flux(yp_segments)

            # Get flux at corners
            I_flx = self._compute_few_group_flux(
                Vector(moc.x_max - 0.001, moc.y_max - 0.001), Direction(-1.0, -1.0)
            )

            for G in range(NG):
                adf[G, ADF.XP] = xp_flx[G] / homog_flux[G]
                adf[G, ADF.XN] = adf[G, ADF.XP]
                adf[G, ADF.YP] = yp_flx[G] / homog_flux[G]
                adf[G, ADF.YN] = adf[G, ADF.YP]
                # The ADFs on the +/- Z sides will be left at unity

                cdf[G, CDF.I] = I_flx[G] / homog_flux[G]
                cdf[G, CDF.II] = cdf[G, CDF.I]
                cdf[G, CDF.III] = cdf[G, CDF.I]
                cdf[G, CDF.IV] = cdf[G, CDF.I]

        return adf, cdf

    def _compute_form_factors(self) -> np.ndarray:
        """
        Computes the one group pin power form factors for the full assembly.

        Returns
        -------
        ndarray
            A 2D Numpy array for the pin power form factors. First index is
            y (from high to low) and the second index is x (from low to high).
        """
        ff = np.zeros((self.shape[1], self.shape[0]))

        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i]

                if not isinstance(cell, FuelPin):
                    continue

                # Pin Power
                pp = cell.compute_pin_linear_power(self._ndl)

                if self.symmetry == Symmetry.Full:
                    ff[j, i] = pp

                elif self.symmetry == Symmetry.Half:
                    ff[j, i] = pp
                    ff[-(j + 1), i] = pp

                elif self.symmetry == Symmetry.Quarter:
                    ox = self.shape[1] // 2

                    ff[j, i + ox] = pp
                    ff[j, -(i + ox + 1)] = pp

                    ff[-(j + 1), i + ox] = pp
                    ff[-(j + 1), -(i + ox + 1)] = pp

        mean_ff = np.mean(ff)
        ff /= mean_ff

        return ff

    def _compute_few_group_xs(self) -> DiffusionCrossSection:
        """
        Computes the few-group diffusion cross sections for the problem.

        Returns
        -------
        DiffusionCrossSection
            Few-group diffusion cross sections for the assembly.

        Raises
        ------
        RuntimeError
            If the condensation_scheme attribute is not set.
        """
        if self.condensation_scheme is None:
            raise RuntimeError("Energy condensation scheme not set.")

        # According to Smith, one should do energy condensation on the
        # diffusion coefficients, and not on the transport cross sections which
        # one could then use to make diffusion coefficients [1]. This is in
        # contradiction to Lattice Physics Computations which states that
        # either method is acceptable [2]. In light of these comments, I have
        # chosen to go with Smith's recommendation of performing energy
        # condensation on the diffusion coefficients.

        homog_xs = self._asmbly_moc.homogenize()
        diff_xs = homog_xs.diffusion_xs()
        flux_spectrum = self._asmbly_moc.homogenize_flux_spectrum()
        return diff_xs.condense(self.condensation_scheme, flux_spectrum)

    def _compute_diffusion_data(self) -> DiffusionData:
        """
        Computes the nodal diffusion data for the assembly.

        Returns
        -------
        DiffusionData
            Few-group diffusion cross sections and discontinuity factors.
        """
        diff_xs = self._compute_few_group_xs()
        adf, cdf = self._compute_adf_cdf()
        ff = self._compute_form_factors()
        return DiffusionData(diff_xs, ff, adf, cdf)

    def _run_assembly_calculation(
        self, self_shield: bool, apply_dancoff_corrections: bool = False
    ) -> None:
        """
        Runs a single MOC calculation, applies critical leakage model, obtains
        the flux spectra for each fuel region, and then normalize the flux
        based on the linear assembly power. Self-shielding can be performed
        if desired. If not, the last computed Dancoff corrections can be
        applied to all fuel pins (to run a time step without explicit
        self-shielding step, using previously known Dancoff corrections).

        Paramters
        ---------
        self_shield : bool
            If True, self-shielding is performed for the fuel and cladding.
        apply_dancoff_corrections : bool, default False
            If self_shield is False and this option is True, the previously
            obtained Dancoff corrections are applied to all cells.
        """
        if self_shield:
            # If we want self-shielding, do that stuff
            self.self_shield_and_xs_update()
        elif apply_dancoff_corrections:
            # Sets dancoff corrections, even if self-shielding wasn't performed
            self.apply_dancoff_corrections()

        self.recompute_all_xs()

        if self._asmbly_moc is None:
            self._init_moc()

        self._asmbly_moc.flux_tolerance = self.flux_tolerance
        self._asmbly_moc.keff_tolerance = self.keff_tolerance

        set_logging_level(LogLevel.Warning)
        self._asmbly_moc.solve()
        set_logging_level(LogLevel.Info)
        scarabee_log(LogLevel.Info, "")
        scarabee_log(LogLevel.Info, "Kinf: {:.5f}".format(self._asmbly_moc.keff))

        self.apply_leakage_model()
        self.obtain_flux_spectra()
        self.normalize_flux_to_power()
        self.apply_infinite_spectrum()

    def _predict_depletion(self, dt: float) -> None:
        # Do all depletions in parallel
        threads = []
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i]
                threads.append(
                    Thread(
                        target=cell.predict_depletion, args=(dt, self._chain, self._ndl)
                    )
                )
                threads[-1].start()
        for t in threads:
            t.join()

    def _correct_depletion(self, dt: float) -> None:
        # Do all depletions in parallel
        threads = []
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i]
                threads.append(
                    Thread(
                        target=cell.correct_depletion, args=(dt, self._chain, self._ndl)
                    )
                )
                threads[-1].start()
        for t in threads:
            t.join()

    def _run_depletion_steps(self) -> None:
        if self._chain is None:
            raise RuntimeError(
                "No depletion chain is present. Cannot run depletion calculation."
            )

        self._keff = np.zeros(self.depletion_exposure_steps.size + 1)
        self._exposures = np.zeros(self.depletion_exposure_steps.size + 1)
        self._times = np.zeros(self.depletion_exposure_steps.size + 1)

        for t, dt in enumerate(self.depletion_time_steps):
            if t > 0:
                self._exposures[t] = (
                    self._exposures[t - 1] + self.depletion_exposure_steps[t - 1]
                )
                self._times[t] = self._times[t - 1] + self.depletion_time_steps[t - 1]

            scarabee_log(LogLevel.Info, "")
            scarabee_log(LogLevel.Info, 60 * "-")
            scarabee_log(LogLevel.Info, "Running Time Step {:}".format(t))
            scarabee_log(
                LogLevel.Info, "Exposure: {:.3E} MWd/kg".format(self._exposures[t])
            )
            scarabee_log(LogLevel.Info, "Time    : {:.3E} days".format(self._times[t]))
            scarabee_log(LogLevel.Info, "")
            # Convert days to seconds
            dt_sec = dt * 60.0 * 60.0 * 24.0

            scarabee_log(LogLevel.Info, "Predictor:")
            # Run initial calcualtion for this time step
            self._run_assembly_calculation(True)
            scarabee_log(LogLevel.Info, "")
            self._keff[t] = self._asmbly_moc.keff

            # Predic isotopes at midpoint of step
            self._predict_depletion(0.5 * dt_sec)

            scarabee_log(LogLevel.Info, "Corrector:")
            # Run the a new transport calcualtion to get rates
            self._run_assembly_calculation(False)

            # Do correction step for isotopes
            self._correct_depletion(dt_sec)

        # Run last step at the end to get keff for our final material compositions
        scarabee_log(LogLevel.Info, "")
        scarabee_log(LogLevel.Info, 60 * "-")
        self._exposures[-1] = self._exposures[-2] + self.depletion_exposure_steps[-1]
        self._times[-1] = self._times[-2] + dt
        scarabee_log(LogLevel.Info, "Running Time Step {:}".format(t))
        scarabee_log(
            LogLevel.Info, "Exposure: {:.3E} MWd/kg".format(self._exposures[-1])
        )
        scarabee_log(LogLevel.Info, "Time    : {:.3E} days".format(self._times[-1]))
        scarabee_log(LogLevel.Info, "")
        self._run_assembly_calculation(True)
        self._keff[-1] = self._asmbly_moc.keff
        scarabee_log(LogLevel.Info, "")

    def solve(self) -> None:
        if self.depletion_exposure_steps is None:
            # Single one-off calulcation
            self._run_assembly_calculation(True)
            self._keff = self._asmbly_moc.keff
        else:
            # Run depletion steps
            self._run_depletion_steps()


# REFERENCES
# [1] K. S. Smith, “Nodal diffusion methods and lattice physics data in LWR
#     analyses: Understanding numerous subtle details,” Prog Nucl Energ,
#     vol. 101, pp. 360–369, 2017, doi: 10.1016/j.pnucene.2017.06.013.
# [2] D. Knott and A. Yamamoto, "Lattice Physics Computations" in
#     Handbook of Nuclear Engineering, 2010, p 1226.
