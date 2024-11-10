from _scarabee import *
import numpy as np
#import matplotlib.pyplot as plt

class Reflector:
    """
    A Reflector instance is responsible for performing transport calculations
    necessary to produce few-group cross sections for the reflector of an LWR.
    The core baffle cross sections are self-shielded as an infinite slab,
    using the Roman two-term approximation. The calculation is performed using
    a 1D Sn simulation, in the group structure of the nuclear data library.
    This removes the need to obtain a fine-group spectrum that would be used
    to condense to an intermediate group structure.

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
    few_group_condensation_scheme : list of pairs of ints
        Defines how the energy groups will be condensed from the microgroup
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
        self.few_group_condensation_scheme = None

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

    def solve(self):
        """
        Runs a 1D annular problem to generate few group cross sections for the reflector.
        """
        if self.few_group_condensation_scheme is None:
            raise RuntimeError("Cannot perform reflector calculation without few-group condensation scheme.")

        # We start by making a ReflectorSn to do 1D calculation
        dx = []
        mats = []

        # We add 2 fuel assemblies worth of homogenized core
        NF = 2 * 10 * 17
        dr = 2. * self.assembly_width / (float(NF))
        dx += [dr]*NF
        mats += [self.fuel]*NF

        # We now add one ring for the gap
        NG = 3
        dr = self.gap_width / (float(NG))
        dx += [dr]*NG
        mats += [self.moderator]*NG

        # Now we add 20 regions for the baffle
        NB = 20
        dr = self.baffle_width / float(NB)
        dx += [dr]*NB
        mats += [self.baffle]*NB

        # Now we add the outer water reflector regions
        ref_width = self.assembly_width - self.gap_width - self.baffle_width
        NR = int(ref_width / 0.02) + 1
        dr = ref_width / float(NR)
        dx += [dr]*NR
        mats += [self.moderator]*NR

        ref_sn = ReflectorSN(mats, dx) 
        ref_sn.solve()

        few_group_flux = np.zeros((len(self.few_group_condensation_scheme), NF+NG+NB+NR))
        for i in range(NF+1+NB+NR):
            for G in range(len(self.few_group_condensation_scheme)):
                g_min = self.few_group_condensation_scheme[G][0]
                g_max = self.few_group_condensation_scheme[G][1]
                for g in range(g_min, g_max+1):
                    few_group_flux[G, i] += ref_sn.flux(i, g)
        dx = np.array(dx)
        x = np.zeros(len(dx))
        for i in range(len(dx)):
            if i == 0:
                x[0] = 0.5*dx[0]
            else:
                x[i] = x[i-1] + 0.5*(dx[i-1] + dx[i])
        
        # Normalize flux to fission production
        norm_few_group_flx = 0.
        for i in range(NF):
            v = dx[i]
            for g in range(self.fuel.ngroups):
                norm_few_group_flx += v*ref_sn.flux(i,g)*self.fuel.vEf(g)
        few_group_flux /= norm_few_group_flx

        # Here we compute the cross sections
        ref_homog_xs   = ref_sn.homogenize(list(range(NF, NF+NG+NB+NR)))
        ref_homog_spec = ref_sn.homogenize_flux_spectrum(list(range(NF, NF+NG+NB+NR)))
        ref_homog_diff_xs = ref_homog_xs.diffusion_xs()
        self.diffusion_xs = ref_homog_diff_xs.condense(self.few_group_condensation_scheme, ref_homog_spec)

        fuel_homog_xs   = ref_sn.homogenize(list(range(0, NF)))
        fuel_homog_spec = ref_sn.homogenize_flux_spectrum(list(range(0, NF)))
        fuel_homog_diff_xs = fuel_homog_xs.diffusion_xs()
        fuel_diffusion_xs = fuel_homog_diff_xs.condense(self.few_group_condensation_scheme, fuel_homog_spec)

        # Do nodal calculation to obtain homogeneous flux
        nodal_tiles = [fuel_diffusion_xs, self.diffusion_xs]
        nodal_geom = DiffusionGeometry(nodal_tiles, [2.*self.assembly_width, self.assembly_width], [8, 4],
                                       [1.], [1], [1.], [1], 1., 0., 1., 1., 1., 1.) 
        nodal_solver = NEMDiffusionDriver(nodal_geom)
        nodal_solver.solve()
        nodal_flux = np.zeros((len(self.few_group_condensation_scheme), len(x)))
        for i in range(len(x)):
            for g in range(len(self.few_group_condensation_scheme)):
                nodal_flux[g, i] = nodal_solver.flux(x[i], 0.5, 0.5, g)
        
        # Normalize flux to fission production
        norm_nd_flx = 0.
        for i in range(NF):
            v = dx[i]
            for g in range(len(self.few_group_condensation_scheme)):
                norm_nd_flx += v*nodal_flux[g,i]*fuel_diffusion_xs.vEf(g)
        nodal_flux /= norm_nd_flx
        
        ## Keep this commented, just in case it's needed later
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
