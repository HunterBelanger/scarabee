from .fuel_pin import FuelPin
from .guide_tube import GuideTube
from .._scarabee import (
    borated_water,
    Material,
    CrossSection,
    NDLibrary,
    Cartesian2D,
    MOCDriver,
    BoundaryCondition,
    SimulationMode,
    YamamotoTabuchi6,
)
from enum import Enum
import numpy as np
from typing import Optional, List, Tuple, Union
import copy

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

        # Set initial boundary conditions
        self._x_min_bc = BoundaryCondition.Periodic
        self._x_max_bc = BoundaryCondition.Periodic
        self._y_min_bc = BoundaryCondition.Periodic
        self._y_max_bc = BoundaryCondition.Periodic

        if self.symmetry == Symmetry.Half:
            self._y_min_bc = BoundaryCondition.Reflective
            self._y_max_bc = BoundaryCondition.Reflective
        elif self.symmetry == Symmetry.Quarter:
            self._x_min_bc = BoundaryCondition.Reflective
            self._x_max_bc = BoundaryCondition.Reflective
        
        # ======================================================================
        # DANCOFF FACTOR CALCULATION DATA
        # ----------------------------------------------------------------------

        # Make water xs for dancoff calculation
        self._moderator_dancoff_xs: CrossSection = CrossSection(
            np.array([self.moderator.potential_xs]),
            np.array([self.moderator.potential_xs]),
            np.array([[0.0]]),
            "Moderator",
        )

        # Isolated cell geometry for Dancoff factor calculations
        self._isolated_dancoff_cells = []
        self._isolated_dancoff_mocs = []

        # Full geometry for Dancoff factor calculations
        self._full_dancoff_cells = []
        self._full_dancoff_geom = None
        self._full_dancoff_moc = None

        # Dancoff factor parameters
        self._dancoff_moc_track_spacing = 0.05
        self._dancoff_moc_num_angles = 64
        
        # ======================================================================
        # TRANSPORT CALCULATION DATA
        # ----------------------------------------------------------------------

        # Make water xs for transport calculation
        self._moderator_xs: CrossSection = self.moderator.dilution_xs(
            self.moderator.size * [1.0e10], self._ndl
        )

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
        return self._dancoff_track_spacing

    @dancoff_moc_track_spacing.setter
    def dancoff_moc_track_spacing(self, dts: float):
        if dts <= 0.0 or dts > 0.1:
            raise ValueError("Dancoff track spacing must be in range (0, 0.1].")
        self._dancoff_moc_track_spacing = dts

    @property
    def dancoff_moc_num_angles(self):
        return self._dancoff_num_angles

    @dancoff_moc_num_angles.setter
    def dancoff_moc_num_angles(self, dna: int):
        if dna % 4 != 0:
            raise ValueError(
                "Number of angles for Dancoff factor calculation must be a multiple of 4."
            )
        if dna < 4:
            raise ValueError(
                "Number of angles for Dancoff factor calculation must be > 4."
            )
        self._dancoff_moc_num_angles = int(dna)

    @property
    def cells(self):
        return self._cells

    @cells.setter
    def cells(self, cells: List[List[Union[FuelPin, GuideTube]]]):
        # Compute the number of expected cells along each dimension [y][x]
        # based on the symmetry of the problem being solved.
        expx = self.shape[0]
        expy = self.shape[1]
        if self.symmetry == Symmetry.Half:
            expy = (expy // 2) + (expy % 2)
        elif self.symmetry == Symmetry.Quarter:
            expy = (expy // 2) + (expy % 2)
            expx = (expx // 2) + (expx % 2)

        if len(cells) != expy:
            raise ValueError(
                "Shape along y of cells list does not agree with assembly shape and symmetry."
            )
        for j in range(len(cells)):
            if len(cells[j]) != expx:
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

    def set_dancoff_moderator_xs(self) -> None:
        """
        Updates the moderator cross section for all Dancoff factor calculations.
        """
        self._moderator_dancoff_xs.set(
            CrossSection(
                np.array([self.moderator.potential_xs]),
                np.array([self.moderator.potential_xs]),
                np.array([[0.0]]),
                "Moderator",
            )
        )

    def set_moderator_xs(self) -> None:
        """
        Updates the moderator cross section for transport calculations.
        """
        self._moderator_dancoff_xs.set(
            borated_water(
                self.boron_ppm, self.moderator_temp, self.moderator_pressure, self._ndl
            )
        )

    def _init_dancoff_components(self) -> None:
        pass        
