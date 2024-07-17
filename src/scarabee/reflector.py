from _scarabee import *
import numpy as np
from typing import Tuple, List, Optional
from copy import copy

class Reflector:
  def __init__(self, fuel: CrossSection, moderator: CrossSection, assembly_width: float, gap_width: float, baffle_width: float, baffle: Material, ndl: NDLibrary):
    self.fuel = fuel
    self.fuel.name = "Fuel"
    self.moderator = moderator
    self.moderator.name = "Moderator"
    self.assembly_width = assembly_width
    self.gap_width = gap_width
    self.baffle_width = baffle_width
    self.condensation_scheme = []
    self.few_group_condensation_scheme = []

    # No Dancoff correction, as looking at 1D isolated slab for baffle
    Ee = 1.0 / (2.0 * self.baffle_width)
    self.baffle = baffle.roman_xs(0.01, Ee, ndl) 
    self.baffle.name = "Baffle"

    # MOC parameters for assembly calculation
    self.track_spacing = 0.02
    self.num_azimuthal_angles = 32
    self.polar_quadrature = YamamotoTabuchi6()
    self.keff_tolerance = 1.0e-5
    self.flux_tolerance = 1.0e-5


    if self.gap_width + self.baffle_width >= self.assembly_width:
      raise RuntimeError("The assembly width is smaller than the sum of the gap and baffle widths.")

    self._fuel_name = self.fuel.name

  @property
  def track_spacing(self):
      return self._track_spacing

  @track_spacing.setter
  def track_spacing(self, value: float):
      if value <= 0.0:
          raise RuntimeError("Track spacing must be > 0.")

      if value >= 1.0:
          raise RuntimeWarning("Track spacing should be < 1.")

      self._track_spacing = value

  @property
  def num_azimuthal_angles(self):
      return self._num_azimuthal_angles

  @num_azimuthal_angles.setter
  def num_azimuthal_angles(self, value: int):
      if value < 4:
          raise RuntimeError("Number of azimuthal angles must be >= 4.")

      if value % 2 != 0:
          raise RuntimeError("Number of azimuthal angles must be even.")

      self._num_azimuthal_angles = value

  def solve(self):
    self._cylindrical_cell_calc() 
    self._moc_1d_calc()

  def _cylindrical_cell_calc(self):
    # We start by making a cylindrical cell. This is just for condensation.
    radii = []
    mats = []

    # First, we add 5 fuel assemblies worth of rings
    NF = 5 * 17
    dr = self.assembly_width / 17.
    last_rad = 0.
    for i in range(NF):
      radii.append(last_rad + dr)
      last_rad = radii[-1]
      mats.append(self.fuel)

    # We now add one ring for the gap
    gap_regions = [NF]
    radii.append(last_rad + self.gap_width)
    last_rad = radii[-1]
    mats.append(self.moderator)

    # Now we add 6 regions for the baffle
    NB = 6
    baffle_regions = list(range(len(radii), len(radii)+NB))
    dr = self.baffle_width / float(NB)
    for i in range(NB):
      radii.append(last_rad + dr)
      last_rad = radii[-1]
      mats.append(self.baffle)

    # Now we add the outer water reflector regions
    ref_width = self.assembly_width - self.gap_width - self.baffle_width
    NR = int(ref_width / 0.3) + 1
    dr = ref_width / float(NR)
    ref_regions = list(range(len(radii), len(radii)+NR))
    for i in range(NR):
      radii.append(last_rad + dr)
      last_rad = radii[-1]
      mats.append(self.moderator)

    cell = CylindricalCell(radii, mats)
    cell.solve()
    cell_flux = CylindricalFluxSolver(cell)
    cell_flux.albedo = 0.
    cell_flux.solve()

    fuel_xs = cell_flux.homogenize(list(range(NF)))
    fuel_spec = cell_flux.homogenize_flux_spectrum(list(range(NF)))
    self.condensed_fuel_xs = fuel_xs.condense(self.condensation_scheme, fuel_spec)
    self.condensed_fuel_xs.name = "Fuel"

    gap_xs = cell_flux.homogenize(gap_regions)
    gap_spec = cell_flux.homogenize_flux_spectrum(gap_regions)
    self.condensed_gap_xs = gap_xs.condense(self.condensation_scheme, gap_spec)
    self.condensed_gap_xs.name = "Gap Moderator"

    baffle_xs = cell_flux.homogenize(baffle_regions)
    baffle_spec = cell_flux.homogenize_flux_spectrum(baffle_regions)
    self.condensed_baffle_xs = baffle_xs.condense(self.condensation_scheme, baffle_spec)
    self.condensed_baffle_xs.name = "Baffle"

    ref_xs = cell_flux.homogenize(ref_regions)
    ref_spec = cell_flux.homogenize_flux_spectrum(ref_regions)
    self.condensed_reflector_xs = ref_xs.condense(self.condensation_scheme, ref_spec)
    self.condensed_reflector_xs.name = "Reflector Moderator"
  
  def _moc_1d_calc(self):
    # We first determine how wide each cell should be. If we will assume a
    # PWR assembly, and try to divide each "pin cell" into 4 regions.
    dy = 100. * self.assembly_width
    Ncell_assembly = 4 * 17
    dx = self.assembly_width / float(Ncell_assembly) 
    fuel_cell = EmptyCell(self.condensed_fuel_xs, dx, dy)

    # The gap will only have 1 cell
    Ncell_gap = 1
    dx_gap = self.gap_width / Ncell_gap
    gap_cell = EmptyCell(self.condensed_gap_xs, dx_gap, dy)

    # The baffle will be broken up into 3 cells
    Ncell_baffle = 6
    dx_baffle = self.baffle_width / Ncell_baffle
    baffle_cell = EmptyCell(self.condensed_baffle_xs, dx_baffle, dy)

    # Now we divide the rest of the reflector
    refl_width = self.assembly_width - (self.gap_width + self.baffle_width)
    Ncell_refl = int(refl_width / 0.5)
    dx_refl = refl_width / Ncell_refl
    refl_cell = EmptyCell(self.condensed_reflector_xs, dx_refl, dy)

    # According to [1], we should have 5 fuel assemblies, and then the reflector
    cells = (5 * Ncell_assembly)*[fuel_cell] + Ncell_gap*[gap_cell] + Ncell_baffle*[baffle_cell] + Ncell_refl*[refl_cell]
    dxs = (5 * Ncell_assembly)*[dx] + Ncell_gap*[dx_gap] + Ncell_baffle*[dx_baffle] + Ncell_refl*[dx_refl]

    geom = Cartesian2D(dxs, [dy])
    geom.set_tiles(cells)

    moc = MOCDriver(geom)
    moc.x_max_bc = BoundaryCondition.Vacuum
    moc.plot()
    moc.generate_tracks(self.num_azimuthal_angles, self.track_spacing, self.polar_quadrature)
    moc.keff_tolerance = self.keff_tolerance
    moc.flux_tolerance = self.flux_tolerance
    moc.solve()
  

# REFERENCES
# [1] S. Huy, M. Guillo, A. Calloo, C. Brosselard, and D. Couyras,
#     “MULTI-GROUP 1D-REFLECTOR MODELLING FOR EDF PWR,” in
#     PHYSOR 2016, Sun Valley, ID, 2016, pp. 74–83.