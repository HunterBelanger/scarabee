from _scarabee import *
import numpy as np
from typing import Tuple, List, Optional
from copy import copy


class Reflector:
    def __init__(
        self,
        fuel: CrossSection,
        moderator: CrossSection,
        assembly_width: float,
        gap_width: float,
        baffle_width: float,
        baffle: Material,
        ndl: NDLibrary,
    ):
        self.fuel = fuel
        self.fuel.name = "Fuel"
        self.moderator = moderator
        self.moderator.name = "Moderator"
        self.assembly_width = assembly_width
        self.gap_width = gap_width
        self.baffle_width = baffle_width
        self.condensation_scheme = []

        # No Dancoff correction, as looking at 1D isolated slab for baffle
        Ee = 1.0 / (2.0 * self.baffle_width)
        self.baffle = baffle.roman_xs(0., Ee, ndl)
        self.baffle.name = "Baffle"

        # MOC parameters for assembly calculation
        self.track_spacing = 0.02
        self.num_azimuthal_angles = 32
        self.polar_quadrature = YamamotoTabuchi6()
        self.keff_tolerance = 1.0e-5
        self.flux_tolerance = 1.0e-5

        if self.gap_width + self.baffle_width >= self.assembly_width:
            raise RuntimeError(
                "The assembly width is smaller than the sum of the gap and baffle widths."
            )

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
        # We start by making a cylindrical cell. This is just for condensation.
        radii = []
        mats = []

        # According to [1], we should have 5 fuel assemblies, and then the reflector
        # First, we add 5 fuel assemblies worth of rings
        NF = 5 * 17
        dr = self.assembly_width / (17.0)
        last_rad = 0.0
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
        baffle_regions = list(range(len(radii), len(radii) + NB))
        dr = self.baffle_width / float(NB)
        for i in range(NB):
            radii.append(last_rad + dr)
            last_rad = radii[-1]
            mats.append(self.baffle)

        # Now we add the outer water reflector regions
        ref_width = self.assembly_width - self.gap_width - self.baffle_width
        NR = int(ref_width / 0.3) + 1
        dr = ref_width / float(NR)
        ref_regions = list(range(len(radii), len(radii) + NR))
        for i in range(NR):
            radii.append(last_rad + dr)
            last_rad = radii[-1]
            mats.append(self.moderator)

        cell = CylindricalCell(radii, mats)
        cell.solve(parallel=True)
        cell_flux = CylindricalFluxSolver(cell)
        cell_flux.albedo = 0.0
        cell_flux.solve()

        fuel_xs = cell_flux.homogenize(list(range(NF)))
        fuel_spec = cell_flux.homogenize_flux_spectrum(list(range(NF)))
        self.condensed_fuel_xs = fuel_xs.condense(self.condensation_scheme, fuel_spec)
        self.condensed_fuel_xs.name = "Fuel"

        homog_xs = cell_flux.homogenize(list(range(NF, len(radii))))
        homog_spec = cell_flux.homogenize_flux_spectrum(list(range(NF, len(radii))))

        NG = homog_xs.ngroups
        fissile = homog_xs.fissile

        D = np.zeros(NG)
        Ea = np.zeros(NG)
        Es = np.zeros((NG, NG))
        if fissile:
            Ef = np.zeros(NG)
            vEf = np.zeros(NG)
            chi = np.zeros(NG)

        for g in range(NG):
            D[g] = 1.0 / (3.0 * homog_xs.Etr(g))
            Ea[g] = homog_xs.Ea(g)

            if fissile:
                Ef[g] = homog_xs.Ef(g)
                vEf[g] = homog_xs.vEf(g)
                chi[g] = homog_xs.chi(g)

            for gg in range(NG):
                Es[g, gg] = homog_xs.Es_tr(g, gg)

        if fissile:
            diff_xs = DiffusionCrossSection(D, Ea, Es, Ef, vEf, chi)
        else:
            diff_xs = DiffusionCrossSection(D, Ea, Es)

        self.diffusion_xs = diff_xs.condense(
            self.condensation_scheme, homog_spec
        )

        NG = self.diffusion_xs.ngroups

        D = []
        Ea = []
        Ef = []
        vEf = []
        chi = []
        for g in range(NG):
            D.append(self.diffusion_xs.D(g))
            Ea.append(self.diffusion_xs.Ea(g))
            Ef.append(self.diffusion_xs.Ef(g))
            vEf.append(self.diffusion_xs.vEf(g))
            chi.append(self.diffusion_xs.chi(g))

        D_str = "  D: "
        Ea_str = " Ea: "
        Ef_str = " Ef: "
        vEf_str = "vEf: "
        chi_str = "chi: "
        Es_strs = []
        for g in range(NG):
            D_str += "{:.4E}  ".format(D[g])
            Ea_str += "{:.4E}  ".format(Ea[g])

            if self.diffusion_xs.fissile:
                Ef_str += "{:.4E}  ".format(Ef[g])
                vEf_str += "{:.4E}  ".format(vEf[g])
                chi_str += "{:.4E}  ".format(chi[g])

            Es_strs.append("{:} -> g: ".format(g + 1))

            for gg in range(NG):
                Es_strs[-1] += "{:.4E}  ".format(self.diffusion_xs.Es(g, gg))

        scarabee_log(LogLevel.Info, D_str)
        scarabee_log(LogLevel.Info, Ea_str)
        if self.diffusion_xs.fissile:
            scarabee_log(LogLevel.Info, Ef_str)
            scarabee_log(LogLevel.Info, vEf_str)
            scarabee_log(LogLevel.Info, chi_str)
        scarabee_log(LogLevel.Info, "Es:  outgoing group ->")
        for g in range(NG):
            scarabee_log(LogLevel.Info, Es_strs[g])


# REFERENCES
# [1] S. Huy, M. Guillo, A. Calloo, C. Brosselard, and D. Couyras,
#     “MULTI-GROUP 1D-REFLECTOR MODELLING FOR EDF PWR,” in
#     PHYSOR 2016, Sun Valley, ID, 2016, pp. 74–83.
