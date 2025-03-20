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
    """ """

    def __init__(
        self,
        shape: Tuple[int, int],
        pitch: float,
        ndl: NDLibrary,
        boron_ppm: float = 800.0,
        moderator_temp: float = 570.0,
        moderator_pressure: float = 15.5,
        symmetry: Symmetry = Symmetry.Full,
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

        self._cells_set = False
        self._cells: List[List[Union[FuelPin, GuideTube]]] = [[]]

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

    @cells.setter
    def cells(self, cells: List[List[Union[FuelPin, GuideTube]]]):
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
        for j in range(len(cells)):
            self._cells.append([])
            for i in range(len(cells[j])):
                self._cells[-1].append(copy.deepcopy(cells[j][i]))
        self._cells_set = True

    @property    
    def leakage_model(self):
        return self._leakage_model

    @leakage_model.setter 
    def leakage_model(self, clm: CriticalLeakage):
        self._leakage_model = clm

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
        # If no leakage, just return
        if self.leakage_model == CriticalLeakage.NoLeakage:
            return

        scarabee_log(LogLevel.Info, "")

        homogenized_moc = self._asmbly_moc.homogenize()

        if self.leakage_model == CriticalLeakage.P1:
            scarabee_log(LogLevel.Info, "Performing P1 criticality spectrum calculation")
            critical_spectrum = P1CriticalitySpectrum(homogenized_moc)
        else:
            scarabee_log(LogLevel.Info, "Performing B1 criticality spectrum calculation")
            critical_spectrum = B1CriticalitySpectrum(homogenized_moc)
        
        self._asmbly_moc.apply_criticality_spectrum(critical_spectrum.flux)

        scarabee_log(LogLevel.Info, "Kinf    : {:.5f}".format(critical_spectrum.k_inf))
        scarabee_log(LogLevel.Info, "Buckling: {:.5f}".format(critical_spectrum.buckling))