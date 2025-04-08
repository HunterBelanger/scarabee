from .fuel_pin import FuelPin
from .guide_tube import GuideTube
from .critical_leakage import CriticalLeakage
from .._scarabee import (
    borated_water,
    Material,
    CrossSection,
    NDLibrary,
    PinCellType,
    Cartesian2D,
    MOCDriver,
    BoundaryCondition,
    SimulationMode,
    YamamotoTabuchi6,
    P1CriticalitySpectrum,
    B1CriticalitySpectrum,
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
        Linear power density of the assembly in w/cm. Default is 155.

    Attributes
    ----------
    shape : tuple of int
        The number of pin cells in the full assembly along x and y.
    pitch : float
        The spacing between fuel pins.
    symmetry : Symmetry
        Symmetry of the fuel assembly. Default is Symmetry.Full.
    linear_power : float
        Linear power density of the assembly in kW/cm. Default is 42.
    initial_heavy_metal_linear_mass : float
        Linear density of heavy metal in the assembly at the beginning of life.
        Has units of kg / cm.
    boron_ppm : float
        Moderator boron concentration in parts per million.
    moderator_temp : float
        Moderator temperature in Kelvin.
    moderator_pressure : float
        Moderator pressure in MPa.
    moderator : Material
        Material representing the assembly moderator.
    dancoff_moc_track_spacing : float
        Spacing between tracks in the MOC calculations for determining Dancoff
        corrections. Default value is 0.05 cm.
    dancoff_moc_num_angles : int
        Number of azimuthal angles in the MOC calculations for determining
        Dancoff corrections. Default value is 32.
    moc_track_spacing : float
        Spacing between tracks in the assembly MOC calculations. Default value
        is 0.05 cm.
    moc_num_angles : int
        Number of azimuthal angles in the assembly MOC calculations. Default
        value is 32.
    leakage_model : CriticalLeakage
        Model used to determine the critical leakage flux spectrum, also known
        as the fundamental mode. Default method is homogeneous P1.
    depletion_exposure_steps : ndarray
        1D Numpy array of assembly exposure steps in units of MWd/kg.
    depletion_time_steps : ndarray
        1D Numpy array of time steps in units of days.
    """

    def __init__(
        self,
        shape: Tuple[int, int],
        pitch: float,
        ndl: NDLibrary,
        cells: List[List[Union[FuelPin, GuideTube]]],
        boron_ppm: float = 800.0,
        moderator_temp: float = 570.0,
        moderator_pressure: float = 15.5,
        symmetry: Symmetry = Symmetry.Full,
        linear_power: float = 42.0,
    ):
        self._ndl = ndl
        self._symmetry = symmetry

        if len(shape) != 2:
            raise ValueError("Shape must have 2 entries.")
        if shape[0] <= 0 or shape[1] <= 0:
            raise ValueError("Shape entries must be > 0.")
        if shape[0] != shape[1] and self.symmetry != Symmetry.Full:
            raise ValueError(
                "Can only use full symmetry with rectangular fuel assembly."
            )
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

        self._moc_track_spacing = 0.05
        self._moc_num_angles = 64

        self._asmbly_cells = []
        self._asmbly_geom: Optional[Cartesian2D] = None
        self._asmbly_moc: Optional[MOCDriver] = None

        self._leakage_model: CriticalLeakage = CriticalLeakage.P1
        self._infinite_flux_spectrum = None # To reset to infinite spectrum in MOC driver

        # Depletion time steps in MWd/kg
        self._depletion_exposure_steps = np.array([])
        # Depletion time steps in days
        self._depletion_time_steps = np.array([])

        # Up to 20 MWd/kg, use 0.5 MWd/kg time steps.
        # From 20 MWd/kg up to 40 MWd/kg, use 2 MWd/kg time steps
        depletion_exposures = 40 * [0.5] + 10 * [2.0]
        self.depletion_exposure_steps = np.array(depletion_exposures)

        # Now we add 4 X 0.5 day time steps for accurate Xe equilibrium
        temp_times = 4 * [0.5] + self.depletion_time_steps.tolist()
        self.depletion_time_steps = np.array(temp_times)

        # Arrays for the depletion exposures (MWd/kg) and times (days)
        self._exposures = np.array([])
        self._times = np.array([])

        # Either a single value or list of values (for each depletion step)
        self._keff: Union[float, List[float]] = 1.0

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
    def cells(self):
        return self._cells

    @property
    def leakage_model(self):
        return self._leakage_model

    @leakage_model.setter
    def leakage_model(self, clm: CriticalLeakage):
        self._leakage_model = clm

    @property
    def depletion_exposure_steps(self):
        return self._depletion_exposure_steps

    @depletion_exposure_steps.setter
    def depletion_exposure_steps(self, steps):
        if not isinstance(steps, np.ndarray):
            raise TypeError("Depletion exposure steps must be a 1D Numpy array.")
        if steps.ndim != 1:
            raise ValueError("Depletion exposure steps must be a 1D Numpy array.")
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
    def depletion_time_steps(self):
        return self._depletion_time_steps

    @depletion_time_steps.setter
    def depletion_time_steps(self, steps):
        if not isinstance(steps, np.ndarray):
            raise TypeError("Depletion time steps must be a 1D Numpy array.")
        if steps.ndim != 1:
            raise ValueError("Depletion time steps must be a 1D Numpy array.")
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
                        0.5 * self.pitch,
                        0.5 * self.pitch,
                        PinCellType.I,
                        False,
                    )
                # Next, check for being on the side with a half pin in quarter symmetry
                elif (
                    self.symmetry == Symmetry.Quarter
                    and self.shape[0] % 2 == 1
                    and i == 0
                ):
                    cell = self.cells[j][i].make_dancoff_moc_cell(
                        self._moderator_dancoff_xs,
                        0.5 * self.pitch,
                        self.pitch,
                        PinCellType.XP,
                        False,
                    )
                # Next, check for being on the bottom row with a half pin
                elif (
                    self.symmetry != Symmetry.Full
                    and self.shape[1] % 2 == 1
                    and j == self._simulated_shape[1] - 1
                ):
                    cell = self.cells[j][i].make_dancoff_moc_cell(
                        self._moderator_dancoff_xs,
                        self.pitch,
                        0.5 * self.pitch,
                        PinCellType.YP,
                        False,
                    )
                # Otherwise, we just make the full cell
                else:
                    cell = self.cells[j][i].make_dancoff_moc_cell(
                        self._moderator_dancoff_xs,
                        self.pitch,
                        self.pitch,
                        PinCellType.Full,
                        False,
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

                cell.set_xs_for_fuel_dancoff_calculation()

                cell.set_isolated_dancoff_fuel_sources(isomoc, self.moderator)

                cell.set_full_dancoff_fuel_sources(
                    self._full_dancoff_moc, self.moderator
                )

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

                cell.set_xs_for_clad_dancoff_calculation(self._ndl)

                cell.set_isolated_dancoff_clad_sources(
                    isomoc, self.moderator, self._ndl
                )

                cell.set_full_dancoff_clad_sources(
                    self._full_dancoff_moc, self.moderator, self._ndl
                )

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

        self.compute_fuel_dancoff_corrections()
        self.compute_clad_dancoff_corrections()
        self.apply_dancoff_corrections()

    # ==========================================================================
    # Transport Calculation Related Methods

    def _init_moc(self):
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
                        0.5 * self.pitch,
                        0.5 * self.pitch,
                        PinCellType.I,
                    )
                # Next, check for being on the side with a half pin in quarter symmetry
                elif (
                    self.symmetry == Symmetry.Quarter
                    and self.shape[0] % 2 == 1
                    and i == 0
                ):
                    cell = self.cells[j][i].make_moc_cell(
                        self._moderator_xs, 0.5 * self.pitch, self.pitch, PinCellType.XP
                    )
                # Next, check for being on the bottom row with a half pin
                elif (
                    self.symmetry != Symmetry.Full
                    and self.shape[1] % 2 == 1
                    and j == self._simulated_shape[1] - 1
                ):
                    cell = self.cells[j][i].make_moc_cell(
                        self._moderator_xs, self.pitch, 0.5 * self.pitch, PinCellType.YP
                    )
                # Otherwise, we just make the full cell
                else:
                    cell = self.cells[j][i].make_moc_cell(
                        self._moderator_xs, self.pitch, self.pitch, PinCellType.Full
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

    def recompute_all_xs(self) -> None:
        """
        Computes and applies all cross sections using the most recent material
        information and Dancoff corrections.
        """
        self.recompute_all_fuel_xs()
        self.recompute_all_gap_xs()
        self.recompute_all_clad_xs()

        self.set_moderator_xs()

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

    def obtain_fuel_flux_spectra(self) -> None:
        """
        Computes the average flux spectrum for each fuel ring from the MOC
        simulation, and saves it in the FuelPin instance.
        """
        for j in range(len(self.cells)):
            for i in range(len(self.cells[j])):
                cell = self.cells[j][i]

                if isinstance(cell, FuelPin):
                    cell.obtain_fuel_flux_spectra(self._asmbly_moc)

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
                cell = self.cells[j][i]
                if isinstance(cell, FuelPin):
                    cell.normalize_flux_spectrum(f)

    def _run_assembly_calculation(
        self, self_shield: bool, apply_dancoff_corrections: bool = False
    ) -> None:
        if self_shield:
            # If we want self-shielding, do that stuff
            self.self_shield_and_xs_update()
        elif apply_dancoff_corrections:
            # Sets dancoff corrections, even if self-shielding wasn't performed
            self.apply_dancoff_corrections()

        self.recompute_all_xs()

        if self._asmbly_moc is None:
            self._init_moc()

        self._asmbly_moc.solve()

        self.apply_leakage_model()
        self.obtain_fuel_flux_spectra()
        self.normalize_flux_to_power()
        self.apply_infinite_spectrum()

    def _predict_depletion(self, dt: float) -> None:
        pass

    def _correct_depletion(self, dt: float) -> None:
        pass

    def _run_depletion_steps(self) -> None:
        self._keff = np.zeros(self.depletion_exposure_steps.size + 1)
        self._exposures = np.zeros(self.depletion_exposure_steps.size + 1)
        self._times = np.zeros(self.depletion_exposure_steps.size + 1)

        for t, dt in enumerate(self.depletion_time_steps):
            if t > 0:
                self._exposures[t] = (
                    self._exposures[t - 1] + self.depletion_exposure_steps[t - 1]
                )
                self._times[t] = self._times[t - 1] + dt

            scarabee_log(LogLevel.Info, "")
            scarabee_log(LogLevel.Info, "Running Time Step {:}".format(t))
            scarabee_log(LogLevel.Info, "Exposure: {:.3E} MWd/kg".format(self._exposures[t]))
            scarabee_log(LogLevel.Info, "Time    : {:.3E} days".format(self._times[t]))
            # Convert days to seconds
            dt_sec = dt * 60.0 * 60.0 * 24.0

            # Run initial calcualtion for this time step
            self._run_assembly_calculation(True)
            scarabee_log(LogLevel.Info, "")
            self._keff[t] = self._asmbly_moc.keff

            # Predic isotopes at midpoint of step
            self._predict_depletion(0.5 * dt_sec)

            # Run the a new transport calcualtion to get rates
            self._run_assembly_calculation(False)

            # Do correction step for isotopes
            self._correct_depletion(dt_sec)
            
            scarabee_log(LogLevel.Info, "")

    def solve(self) -> None:
        if self.depletion_exposure_steps.size == 0:
            # Single one-off calulcation
            self._run_assembly_calculation(True)
            self._keff = self._asmbly_moc.keff
        else:
            # Run depletion steps
            self._run_depletion_steps()
