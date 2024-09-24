from _scarabee import *
import numpy as np
#import matplotlib.pyplot as plt

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
        structure of the nuclear data library, to the macrogroup structure
        used in the MOC calculation.
    few_group_condensation_scheme : list of pairs of ints
        Defines how the energy groups will be condensed from the macrogroup
        structure of to the few-group structure used in nodal calculations.
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
        self.condensation_scheme = None
        self.few_group_condensation_scheme = None

        # No Dancoff correction, as looking at 1D isolated slab for baffle
        Ee = 1.0 / (2.0 * self.baffle_width)
        self.baffle = baffle.roman_xs(0.0, Ee, ndl)
        self.baffle.name = "Baffle"

        self.keff_tolerance = 1.0e-5
        self.flux_tolerance = 1.0e-5

        self._num_azimuthal_angles = 64
        self._track_spacing = 0.01

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
        if self.condensation_scheme is None:
            raise RuntimeError("Cannot perform reflector calculation without condensation scheme.")
        
        if self.few_group_condensation_scheme is None:
            raise RuntimeError("Cannot perform reflector calculation without few-group condensation scheme.")

        # We start by making a cylindrical cell. This is just for condensation.
        radii = []
        mats = []

        # We add 2 fuel assemblies worth of homogenized core
        NF = 2 * 17
        dr = self.assembly_width / (17.0)
        last_rad = 0.0
        fuel_regions = list(range(NF))
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

        # We now get homogenized xs and spectrums
        fuel_spectrum = cell_flux.homogenize_flux_spectrum(fuel_regions)
        fuel_homog_xs = cell_flux.homogenize(fuel_regions)
        gap_spectrum = cell_flux.homogenize_flux_spectrum(gap_regions)
        gap_homog_xs = cell_flux.homogenize(gap_regions)
        baffle_spectrum = cell_flux.homogenize_flux_spectrum(baffle_regions)
        baffle_homog_xs = cell_flux.homogenize(baffle_regions)
        ref_spectrum = cell_flux.homogenize_flux_spectrum(ref_regions)
        ref_homog_xs = cell_flux.homogenize(ref_regions)

        # Now condense cross sections for the assembly-reflector MOC calculation
        fuel_macro_xs = fuel_homog_xs.condense(self.condensation_scheme, fuel_spectrum)
        gap_macro_xs = gap_homog_xs.condense(self.condensation_scheme, gap_spectrum)
        baffle_macro_xs = baffle_homog_xs.condense(self.condensation_scheme, baffle_spectrum)
        ref_macro_xs = ref_homog_xs.condense(self.condensation_scheme, ref_spectrum)

        # We make a 1D MOC geometry for computing the reflector data
        dy = 2.*self.assembly_width + self.gap_width + self.baffle_width + ref_width
        dx = []
        moc_1d_cells = []
        NF = 2 * 17 * 2 # 17 pins, 2 FSR per pin cell
        dr = (2.*self.assembly_width) / NF
        dx += NF*[dr]
        moc_1d_cells += NF*[EmptyCell(fuel_macro_xs, dr, dy)]
        dr = self.gap_width
        dx += [dr]
        moc_1d_cells += [EmptyCell(gap_macro_xs, dr, dy)]
        dr = self.baffle_width / NB
        dx += NB*[dr]
        moc_1d_cells += NB*[EmptyCell(baffle_macro_xs, dr, dy)]
        dr = ref_width / NR
        dx += NR*[dr]
        moc_1d_cells += NR * [EmptyCell(ref_macro_xs, dr, dy)]
        moc_geom = Cartesian2D(dx, [dy])
        moc_geom.set_tiles(moc_1d_cells)

        moc = MOCDriver(moc_geom)
        moc.x_max_bc = BoundaryCondition.Vacuum
        moc.generate_tracks(self.num_azimuthal_angles, self.track_spacing, YamamotoTabuchi6())
        moc.keff_tolerance = self.keff_tolerance
        moc.flux_tolerance = self.flux_tolerance
        moc.solve()

        few_group_flux = np.zeros((len(self.few_group_condensation_scheme), NF+1+NB+NR))
        for i in range(NF+1+NB+NR):
            for G in range(len(self.few_group_condensation_scheme)):
                g_min = self.few_group_condensation_scheme[G][0]
                g_max = self.few_group_condensation_scheme[G][1]

                for g in range(g_min, g_max+1):
                    few_group_flux[G, i] += moc.flux(i, g)
        dx = np.array(dx)
        x = np.zeros(len(dx))
        for i in range(len(dx)):
            if i == 0:
                x[0] = 0.5*dx[0]
            else:
                x[i] = x[i-1] + 0.5*(dx[i-1] + dx[i])

        # Here we compute the cross sections
        homog_xs = moc.homogenize(list(range(NF, NF+1+NB+NR)))
        homog_spec = moc.homogenize_flux_spectrum(list(range(NF, NF+1+NB+NR)))
        self.diffusion_xs = Reflector._get_diffusion_xs(homog_xs, homog_spec, self.few_group_condensation_scheme)

        homog_xs = moc.homogenize(list(range(0, NF)))
        homog_spec = moc.homogenize_flux_spectrum(list(range(0, NF)))
        fuel_diffusion_xs = Reflector._get_diffusion_xs(homog_xs, homog_spec, self.few_group_condensation_scheme)

        # Do nodal calculation to obtain homogeneous flux
        nodal_tiles = [fuel_diffusion_xs, self.diffusion_xs]
        nodal_geom = DiffusionGeometry(nodal_tiles, [2.*self.assembly_width, self.assembly_width], [8, 4],
                                       [10.], [1], [10.], [1], 1., 0., 1., 1., 1., 1.) 
        nodal_solver = NEMDiffusionDriver(nodal_geom)
        nodal_solver.solve()
        nodal_flux = np.zeros((len(self.few_group_condensation_scheme), len(x)))
        for i in range(len(x)):
            for g in range(len(self.few_group_condensation_scheme)):
                nodal_flux[g, i] = nodal_solver.flux(x[i], 5., 5., g)
        
        # Normalize flux to fast group reflective boundary
        few_group_flux /= few_group_flux[0,0]
        nodal_flux /= nodal_flux[0,0]
        
        # Keep this commented, just in case it's needed later
        #plt.plot(x, few_group_flux[0,:], label="Fast Hetero")
        #plt.plot(x, few_group_flux[1,:], label="Thermal Hetero")
        #plt.plot(x, nodal_flux[0,:], label="Fast Homo")
        #plt.plot(x, nodal_flux[1,:], label="Thermal Homo")
        #plt.legend().set_draggable(True)
        #plt.show()

        # Here, we compute the ADFs
        self.adf = np.zeros((len(self.few_group_condensation_scheme), 4))
        for G in range(len(self.few_group_condensation_scheme)):
            self.adf[G,:] = few_group_flux[G, NF] / nodal_flux[G, NF] 
        
        # Write data to terminal/output file
        self._write_data()
    
    def _get_diffusion_xs(homog_xs, homog_spec, cond_scheme):
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

        return diff_xs.condense(cond_scheme, homog_spec)

    def _write_data(self):
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
