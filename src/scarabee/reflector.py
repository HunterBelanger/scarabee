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
        
        # Here we compute the cross sections
        ref_homog_xs   = ref_sn.homogenize(list(range(NF, NF+NG+NB+NR)))
        ref_homog_spec = ref_sn.homogenize_flux_spectrum(list(range(NF, NF+NG+NB+NR)))
        ref_homog_diff_xs = ref_homog_xs.diffusion_xs()
        self.diffusion_xs = ref_homog_diff_xs.condense(self.few_group_condensation_scheme, ref_homog_spec)

        fuel_homog_xs   = ref_sn.homogenize(list(range(0, NF)))
        fuel_homog_spec = ref_sn.homogenize_flux_spectrum(list(range(0, NF)))
        fuel_homog_diff_xs = fuel_homog_xs.diffusion_xs()
        fuel_diffusion_xs = fuel_homog_diff_xs.condense(self.few_group_condensation_scheme, fuel_homog_spec)

        # Obtain net currents at node boundaries and average flux
        avg_flx_ref = np.zeros(len(self.few_group_condensation_scheme))
        avg_flx_fuel = np.zeros(len(self.few_group_condensation_scheme))
        j_0   = np.zeros(len(self.few_group_condensation_scheme))
        j_mid = np.zeros(len(self.few_group_condensation_scheme))
        j_max = np.zeros(len(self.few_group_condensation_scheme))
        s_mid = NF
        s_max = ref_sn.nsurfaces - 1
        for g in range(len(self.few_group_condensation_scheme)):
            avg_flx_ref[g] = np.mean(few_group_flux[g,NF:])
            avg_flx_fuel[g] = np.mean(few_group_flux[g,:NF])
            
            g_min = self.few_group_condensation_scheme[g][0]
            g_max = self.few_group_condensation_scheme[g][1]
            for gg in range(g_min, g_max+1):
                j_mid[g] += ref_sn.current(s_mid, gg)
                j_max[g] += ref_sn.current(s_max, gg)
        
        # Do nodal calculation to obtain homogeneous flux
        a_fuel = self._nodal_calc(np.sum(dx[:NF]), ref_sn.keff, fuel_diffusion_xs, avg_flx_fuel, j_0, j_mid)
        a_ref = self._nodal_calc(np.sum(dx[NF:]), ref_sn.keff, self.diffusion_xs, avg_flx_ref, j_mid, j_max)
        
        # Compute the ADFs
        self.adf = np.zeros((len(self.few_group_condensation_scheme), 4))
        for G in range(len(self.few_group_condensation_scheme)):
            heter_flx_fuel = few_group_flux[G, NF-1]
            homog_flx_fuel = a_fuel[G, 0] + 0.5*a_fuel[G, 1] + 0.5*a_fuel[G, 2]

            heter_flx_ref = few_group_flux[G, NF]
            homog_flx_ref = a_ref[G,0] - 0.5*a_ref[G, 1] + 0.5*a_ref[G, 2]

            f_fuel = heter_flx_fuel / homog_flx_fuel
            f_ref = heter_flx_ref / homog_flx_ref

            # Normalize to fuel DF
            f_ref = f_ref / f_fuel

            self.adf[G,:] = f_ref
        
        # Write data to terminal/output file
        self._write_data()


    def _nodal_calc(self, dx: float, keff: float, xs: DiffusionCrossSection, avg_flx: np.ndarray, j_neg: np.ndarray, j_pos: np.ndarray):
        """
        Performs the nodal diffusion calculation with reference currents to return the reference nodal flux.

        Paramters
        ---------
        dx : float
             Width of the node.
        keff : float
             Multiplication factor from reference Sn calculation.
        xs : DiffusionCrossSection
             Diffusion cross sections homogenized for the node.
        avg_flx : np.ndarray
             Average flux in each group within the node from reference Sn
             calculation.
        j_neg : np.ndarray
             Reference net current in each group on the negative boundary
             from the reference Sn calculation.
        j_pos : np.ndarray
             Reference net current in each group on the positive boundary
             from the reference Sn calculation.

        Returns
        -------
        np.ndarray
             A 2D numpy array containing the coefficients a_g,i. The first
             index is the group, and the second is the coefficient, which
             ranges from 0 (average flux) to 4.
        """
        NG = xs.ngroups
        Na = NG * 4            # number of a coefficients to solve for
        invs_dx = 1. / dx
        A = np.zeros((Na, Na)) # Matrix to hold all coefficients
        b = np.zeros((Na))     # Results array

        # Load the coefficient matrix and results array
        j = 0 # Matrix row
        for g in range(NG):
            Dg_dx = xs.D(g) * invs_dx
            Erf_g = 0.
            chi_g_keff = xs.chi(g) / keff

            # Each group has 4 equations. See reference [1].
            
            # Eq 2.54
            A[j, g*4 + 2] -= 0.5*Dg_dx*invs_dx
            A[j, g*4 + 0] += (Erf_g / 12.)
            A[j, g*4 + 2] -= 0.1*(Erf_g / 12.)
            for gg in range(NG):
                if gg == g: continue
                A[j, gg*4 + 0] -= (xs.Es(gg, g) / 12.)
                A[j, gg*4 + 2] += 0.1*(xs.Es(gg, g) / 12.)
                A[j, gg*4 + 0] -= chi_g_keff*(xs.vEf(gg) / 12.)
                A[j, gg*4 + 2] += 0.1*chi_g_keff*(xs.vEf(gg) / 12.)
            b[j] = 0.
            j += 1

            # Eq 2.55
            A[j, g*4 + 3] -= 0.2*Dg_dx*invs_dx
            A[j, g*4 + 1] += (Erf_g / 20.)
            A[j, g*4 + 3] -= (Erf_g / (20. * 35.))
            for gg in range(NG):
                if gg == g: continue
                A[j, gg*4 + 1] -= (xs.Es(gg, g) / 20.)
                A[j, gg*4 + 3] += (xs.Es(gg, g) / (20. * 35.))
                A[j, gg*4 + 1] -= chi_g_keff*(xs.vEf(gg) / 20.)
                A[j, gg*4 + 3] += chi_g_keff*(xs.vEf(gg) / (20. * 35.))
            b[j] = 0.
            j += 1

            # Eq 2.57
            A[j, g*4 + 0] -= Dg_dx
            A[j, g*4 + 1] += Dg_dx
            A[j, g*4 + 2] -= 0.5*Dg_dx
            A[j, g*4 + 3] += 0.2*Dg_dx
            b[j] = j_neg[g]
            j += 1

            # Eq 2.58
            A[j, g*4 + 0] -= Dg_dx
            A[j, g*4 + 1] -= Dg_dx
            A[j, g*4 + 2] -= 0.5*Dg_dx
            A[j, g*4 + 3] -= 0.2*Dg_dx
            b[j] = j_pos[g]
            j += 1
        
        a_tmp = np.linalg.solve(A, b)
        a_tmp = np.reshape(a_tmp, (NG, 4))

        a = np.zeros((NG, 5))
        for g in range(NG):
            a[g,0] = avg_flx[g]
            a[g,1:] = a_tmp[g,:]

        return a


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
# [1] S. Machach, “Étude des techniques d’équivalence nodale appliquées aux
#     modèles de réflecteurs dans les réacteurs à eau pressurisée,”
#     Polytechnique Montréal, 2022.