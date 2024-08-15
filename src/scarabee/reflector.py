from _scarabee import *
import numpy as np

class Reflector:
    """
    A Reflector instance is responsible for performing transport calculations
    necessary to produce few-group cross sections for the reflector of an LWR.
    The core baffle cross sections are self-shielded as an infinite slab,
    using the Roman two-term approximation.

    Parameters
    ----------
    fuel : CrossSection
        Homogenized cross section which is representative of a pin cell.
        This is typically obtained from a previous lattice calcualtion.
    moderator : CrossSection
        Material cross sections for the moderator at desired temperature
        and density.
    assembly_width : float
        Width of a single fuel assembly (and the reflector to be modeled).
    gap_width : float
        Width of the moderator gap between the assembly and the core baffle.
    baffle_width : float
        Width of the core baffle.
    baffle : Material
        Material for the core baffle at desired temperature and density.
    ndl : NDLibrary
        Nuclear data library for constructing the baffle cross sections.

    Attributes
    ----------
    condensation_scheme : list of pairs of ints
        Defines how the energy groups will be condensed from the microgroup
        structure of the nuclear data library, to the few-group structure
        used in the nodal calculation.
    fuel : CrossSection
        Cross sections for a homogenized fuel assembly.
    moderator : CrossSection
        Cross sections for the moderator.
    assembly_width : float
        Width of fuel assembly and reflector.
    gap_width : float
        Width of the moderator gap between a fuel assembly and the core baffle.
    baffle_width : float
        Width of the core baffle.
    baffle : CrossSection
        Self-shielded cross sections for the core baffle.
    keff_tolerance : float
        Convergence criteria for keff. Default is 1.E-5.
    flux_tolerance : float
        Convergence criteria for the flux. Default is 1.E-5.
    diffusion_xs : DiffusionCrossSection
        The few-group diffsuion group constants for the reflector region.
    adf : ndarray
        The assembly discontinuity factors.
    """

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
        self.baffle = baffle.roman_xs(0.0, Ee, ndl)
        self.baffle.name = "Baffle"

        self.keff_tolerance = 1.0e-5
        self.flux_tolerance = 1.0e-5

        if self.gap_width + self.baffle_width >= self.assembly_width:
            raise RuntimeError(
                "The assembly width is smaller than the sum of the gap and baffle widths."
            )

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
        """
        Runs a 1D annular problem to generate few group cross sections for the reflector.
        """
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
        cell_flux.keff_tolerance = self.keff_tolerance
        cell_flux.flux_tolerance = self.flux_tolerance
        cell_flux.solve(parallel=True)
        
        # Here, we compute the ADFs
        homog_flux_spec = cell_flux.homogenize_flux_spectrum(gap_regions+baffle_regions+ref_regions)
        gap_flux_spec = cell_flux.homogenize_flux_spectrum(gap_regions+baffle_regions)

        homog_flux = np.zeros(len(self.condensation_scheme))
        gap_flux = np.zeros(len(self.condensation_scheme))
        self.adf = np.zeros((len(self.condensation_scheme), 4))
        for G in range(len(self.condensation_scheme)):
            g_min = self.condensation_scheme[G][0]
            g_max = self.condensation_scheme[G][1]

            for g in range(g_min, g_max+1):
                homog_flux[G] += homog_flux_spec[g]
                gap_flux[G] += gap_flux_spec[g]
            
            self.adf[G,:] = gap_flux[G] / homog_flux[G]
        
        # Here we compute the cross sections
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

        self.diffusion_xs = diff_xs.condense(self.condensation_scheme, homog_spec)

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
        scarabee_log(LogLevel.Info, "ADF: {:}".format(self.adf[0,:]))
        for g in range(1, NG):
            scarabee_log(LogLevel.Info, "     {:}".format(self.adf[g,:]))

    def save_diffusion_data(self, fname):
        if self.diffusion_xs is None:
            raise RuntimeError("No diffusion cross sections.")

        if self.adf is None:
            raise RuntimeError("No ADFs")

        self.diffusion_data = DiffusionData(self.diffusion_xs)
        self.diffusion_data.adf = self.adf
        self.diffusion_data.save(fname)


# REFERENCES
# [1] S. Huy, M. Guillo, A. Calloo, C. Brosselard, and D. Couyras,
#     “MULTI-GROUP 1D-REFLECTOR MODELLING FOR EDF PWR,” in
#     PHYSOR 2016, Sun Valley, ID, 2016, pp. 74–83.
