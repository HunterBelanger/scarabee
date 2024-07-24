from _scarabee import *
from .pins import *
import numpy as np
from typing import Tuple
from copy import copy, deepcopy


class PWRAssembly:
    """
    A PWRAssembly instance is responsible for performing all the lattice
    calculations necessary to produce few-group cross sections for a single
    PWR assembly.

    Parameters
    ----------
    pitch : float
        Regular spacing between fuel pins.
    moderator : Material
        Material representing the moderator composition, temperature, and
        density. This is used for all of the moderator in the assembly.
    shape : (int, int)
        A tuple containing the number of pins across and down for the fuel
        assembly.
    ndl : NDLibrary
        The nuclear data library to be used in the calculations.

    Attributes
    ----------
    condensation_scheme : list of pairs of ints
        Defines how the energy groups will be condensed from the microgroup
        structure of the nuclear data library, to the macrogroup structure
        used in the full assembly calculation.
    few_group_condensation_scheme : list of pairs of ints
        Defines how the energy groups will be condensed from the macrogroup
        structure of the assembly calculation to the few group structure for
        the nodal diffusion calculation.
    pins : list of FuelPin or GuideTube
        Description of the pins which make up the assembly geometry.
    dancoff_track_spacing : float
        Spacing between tracks when calculating Dancoff factors.
        Default is 0.05.
    dancoff_num_azimuthal_angles : int
        Number of azimuthal angles when calculating Dancoff factors.
        Default is 64.
    dancoff_polar_quadrature: PolarQuadrature
        The polar quadrature used when calculating Dancoff factors.
        Default is YamamotoTabuchi6.
    dancoff_isolation_factor : float
        The factor used to multiply the pitch when calculating the flux in an
        isolated pin. Default is 20.
    track_spacing : float
        Spacing between tracks in the assembly calculation. Default is 0.02.
    num_azimuthal_angles : int
        Number of azimuthal angles in the assembly calculation. Default is 32.
    polar_quadrature: PolarQuadrature
        The polar quadrature used in the assembly calculation.
        Default is YamamotoTabuchi6.
    keff_tolerance : float
        Convergence criteria for keff. Default is 1.E-5.
    flux_tolerance : float
        Convergence criteria for the flux. Default is 1.E-5.
    plot_assembly : bool
        Indicates wether the GUI plotter for the assembly geometry will be
        activated before performing the calcualtion.
    criticality_spectrum_method : str or None
        The type of leakage approximation used to modify the assembly flux
        spectrum. Acceptable values are "B1", "P1", or None.
    moc_geom : Cartesian2D
        Geometry used in the assembly calculation. Is None until solve has
        been called.
    moc : MOCDriver
        The method of characteristics solver for the assembly calculation.
        Is None until solve has been called.
    """

    def __init__(
        self, pitch: float, moderator: Material, shape: Tuple[int, int], ndl: NDLibrary
    ):
        self.pitch = pitch
        self.ndl = ndl  # Must assign first for calculating moderator xs
        self.moderator = moderator
        self._shape = shape
        self.condensation_scheme = []
        self.few_group_condensation_scheme = []
        self._pins = []

        # MOC parameters for computing dancoff corrections
        self.dancoff_track_spacing = 0.05
        self.dancoff_num_azimuthal_angles = 64
        self.dancoff_isolation_factor = 20.0
        self.dancoff_polar_quadrature = YamamotoTabuchi6()

        # MOC parameters for assembly calculation
        self.track_spacing = 0.02
        self.num_azimuthal_angles = 32
        self.polar_quadrature = YamamotoTabuchi6()
        self.keff_tolerance = 1.0e-5
        self.flux_tolerance = 1.0e-5

        self.plot_assembly = False
        self.moc_geom = None
        self.moc = None

        self.criticality_spectrum_method = "P1"

    @property
    def criticality_spectrum_method(self):
        return self._criticality_spectrum_method

    @criticality_spectrum_method.setter
    def criticality_spectrum_method(self, csm):
        if csm not in ["B1", "b1", "P1", "p1", None]:
            raise RuntimeError("Unknown criticality spectrum method.")

        if csm in ["B1", "b1"]:
            self._criticality_spectrum_method = "B1"
        else:
            self._criticality_spectrum_method = "P1"

    @property
    def pins(self):
        return self._pins

    @pins.setter
    def pins(self, pins):
        if len(pins) != self.shape[0] * self.shape[1]:
            raise RuntimeError(
                "The number of pins does not agree with the assembly shape."
            )

        # First, need a copy of pins, as we will add condensed cross sections
        # and other info to them.
        self._pins = []
        for i in range(len(pins)):
            self._pins.append(copy(pins[i]))
            self._pins[-1].condensed_xs  = deepcopy(pins[i].condensed_xs)
            self._pins[-1].clad_xs_vals  = deepcopy(pins[i].clad_xs_vals)
            self._pins[-1].clad_flux_vals  = deepcopy(pins[i].clad_flux_vals)
            self._pins[-1].clad_a = deepcopy(pins[i].clad_a)
            self._pins[-1].clad_b = deepcopy(pins[i].clad_b)

            if isinstance(pins[i], FuelPin):
                self._pins[-1].fuel_xs_vals  = deepcopy(pins[i].fuel_xs_vals)
                self._pins[-1].fuel_flux_vals  = deepcopy(pins[i].fuel_flux_vals)
                self._pins[-1].fuel_a = deepcopy(pins[i].fuel_a)
                self._pins[-1].fuel_b = deepcopy(pins[i].fuel_b)


    @property
    def shape(self):
        return self._shape

    @property
    def moderator(self):
        return self._moderator

    @moderator.setter
    def moderator(self, value: Material):
        self._moderator = value
        self.moderator_xs = self.moderator.dilution_xs(
            self.moderator.size * [1.0e10], self.ndl
        )
        self.moderator_xs.name = "Moderator"

    @property
    def pitch(self):
        return self._pitch

    @pitch.setter
    def pitch(self, value: float):
        if value <= 0.0:
            raise RuntimeError("Pitch must be > 0.")
        self._pitch = value

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

    @property
    def dancoff_track_spacing(self):
        return self._dancoff_track_spacing

    @dancoff_track_spacing.setter
    def dancoff_track_spacing(self, value: float):
        if value <= 0.0:
            raise RuntimeError("Dancoff factor track spacing must be > 0.")

        if value >= 1.0:
            raise RuntimeWarning("Dancoff factor track spacing should be < 1.")

        self._dancoff_track_spacing = value

    @property
    def dancoff_num_azimuthal_angles(self):
        return self._dancoff_num_azimuthal_angles

    @dancoff_num_azimuthal_angles.setter
    def dancoff_num_azimuthal_angles(self, value: int):
        if value < 4:
            raise RuntimeError(
                "Number of azimuthal angles in Dancoff factor calculation must be >= 4."
            )

        if value % 2 != 0:
            raise RuntimeError(
                "Number of azimuthal angles in Dancoff factor calculation must be even."
            )

        self._dancoff_num_azimuthal_angles = value

    @property
    def dancoff_isolation_factor(self):
        return self._dancoff_isolation_factor

    @dancoff_isolation_factor.setter
    def dancoff_isolation_factor(self, value: float):
        if value <= 1.0:
            raise RuntimeError("Dancoff isolation factor must be > 1.")

        self._dancoff_isolation_factor = value

    @property
    def dancoff_polar_quadrature(self):
        return self._dancoff_polar_quadrature

    @dancoff_polar_quadrature.setter
    def dancoff_polar_quadrature(self, value: PolarQuadrature):
        if value.__class__.__name__ not in [
            "Legendre2",
            "Legendre4",
            "Legendre6",
            "Legendre8",
            "Legendre10",
            "Legendre12",
            "YamamotoTabuchi2",
            "YamamotoTabuchi6",
        ]:
            raise RuntimeError("Unknown polar quadrature type.")

        self._dancoff_polar_quadrature = value

    def solve(self):
        """
        Runs the entire assembly calculation chain, resulting in few group cross sections.
        """
        if len(self.pins) == 0:
            raise RuntimeError("Cannot solve PWR assembly problem. No pins provided.")

        if len(self.condensation_scheme) == 0:
            raise RuntimeError(
                "Cannot solve PWR assembly problem. No energy condensation scheme provided."
            )

        if len(self.few_group_condensation_scheme) == 0:
            raise RuntimeError(
                "Cannot solve PWR assembly problem. No few group energy condensation scheme provided."
            )

        #self._get_fuel_dancoff_corrections()
        self._self_shield_fuel()
        #self._get_clad_dancoff_corrections()
        #self._self_shield_clad()
        #self._pin_cell_calc()
        #self._condense_xs()
        #self._moc()
        #self._criticality_spectrum_calc()
        #self._few_group_xs()
    
    def _self_shield_fuel(self):
        scarabee_log(LogLevel.Info, "")
        scarabee_log(LogLevel.Info, "Performing self-shielding for fuel")
        set_logging_level(LogLevel.Warning)

        Ets = [1.E-5, 1.E-4, 1.E-3, 1.E-2, 1.E-1, 1.E0, 1.E1, 1.E2, 1.E3, 1.E4, 1.E5]

        for Et in Ets:
            set_logging_level(LogLevel.Info)
            scarabee_log(LogLevel.Info, "  Performing XS simulation")
            set_logging_level(LogLevel.Warning)   
            # Now we setup the lattice problem
            fuel_df_pins = []
            for pin in self.pins:
                fuel_df_pins.append(
                    pin.make_fuel_dancoff_cell(pitch=self.pitch, xs=Et, moderator=self.moderator)
                )
            geom = Cartesian2D(self.shape[0] * [self.pitch], self.shape[1] * [self.pitch])
            geom.set_tiles(fuel_df_pins)
            moc = MOCDriver(geom)

            # Set the source
            for i in range(moc.nfsr):
                i_xs = moc.xs(i)
                if i_xs.name != "Fuel":
                    # If we aren't in the fuel, set the source to be the value of
                    # the potential xs (should be Et).
                    moc.set_extern_src(i, 0, i_xs.Et(0))

            # Solve the lattice problem
            moc.generate_tracks(
                self.dancoff_num_azimuthal_angles,
                self.dancoff_track_spacing,
                self.dancoff_polar_quadrature,
            )
            moc.sim_mode = SimulationMode.FixedSource
            moc.flux_tolerance = 1.0e-5
            moc.solve()

            # Now we need to save the flux for each pin
            u = Direction(1.0, 0.0)
            pin_i = 0
            for j in range(self.shape[1]):
                y = moc.y_max - (j + 0.5) * self.pitch
                for i in range(self.shape[0]):
                    x = moc.x_min + (i + 0.5) * self.pitch
                    r = Vector(x, y)
                    xs = moc.xs(r, u)
                    if xs.name == "Fuel":
                        flux = moc.flux(r, u, 0)
                        self.pins[pin_i].fuel_xs_vals.append(Et)
                        self.pins[pin_i].fuel_flux_vals.append(flux)
                    pin_i += 1
            
        for i in range(len(self.pins)):
            if isinstance(self.pins[i], FuelPin):
                self.pins[i].determine_fuel_self_shielding_params()
                print(self.pins[i].fuel_a + self.pins[i].fuel_b)
                self.pins[i].plot_fuel_ss()
 

    def _get_fuel_dancoff_corrections(self):
        scarabee_log(LogLevel.Info, "")
        scarabee_log(LogLevel.Info, "Computing Dancoff factors for fuel")
        set_logging_level(LogLevel.Warning)
        # We first make the system for an isolated fuel pin.
        # We isolate it by multiplying the pitch by 20.
        isolated_fp = None
        for pin in self.pins:
            if isinstance(pin, FuelPin):
                isolated_fp = pin.make_fuel_dancoff_cell(
                    pitch=self.dancoff_isolation_factor * self.pitch,
                    moderator=self.moderator,
                )
        if isolated_fp is None:
            raise RuntimeError("No FuelPin type found in pins.")
        iso_geom = Cartesian2D(
            [self.dancoff_isolation_factor * self.pitch],
            [self.dancoff_isolation_factor * self.pitch],
        )
        iso_geom.set_tiles([isolated_fp])
        iso_moc = MOCDriver(iso_geom)

        # Set the source
        for i in range(iso_moc.nfsr):
            i_xs = iso_moc.xs(i)
            if i_xs.name != "Fuel":
                # If we aren't in the fuel, set the source to be the value of
                # the potential xs (should be Et).
                iso_moc.set_extern_src(i, 0, i_xs.Et(0))

        # Solve the isolated pin problem
        iso_moc.x_min_bc = BoundaryCondition.Vacuum
        iso_moc.x_max_bc = BoundaryCondition.Vacuum
        iso_moc.y_min_bc = BoundaryCondition.Vacuum
        iso_moc.y_max_bc = BoundaryCondition.Vacuum
        iso_moc.generate_tracks(
            self.dancoff_num_azimuthal_angles,
            self.dancoff_track_spacing,
            self.dancoff_polar_quadrature,
        )
        iso_moc.sim_mode = SimulationMode.FixedSource
        iso_moc.flux_tolerance = 1.0e-5
        iso_moc.solve()
        iso_flux = iso_moc.flux(0, 0)

        # Now we setup the lattice problem
        fuel_df_pins = []
        for pin in self.pins:
            fuel_df_pins.append(
                pin.make_fuel_dancoff_cell(pitch=self.pitch, moderator=self.moderator)
            )
        geom = Cartesian2D(self.shape[0] * [self.pitch], self.shape[1] * [self.pitch])
        geom.set_tiles(fuel_df_pins)
        moc = MOCDriver(geom)

        # Set the source
        for i in range(moc.nfsr):
            i_xs = moc.xs(i)
            if i_xs.name != "Fuel":
                # If we aren't in the fuel, set the source to be the value of
                # the potential xs (should be Et).
                moc.set_extern_src(i, 0, i_xs.Et(0))

        # Solve the lattice problem
        moc.generate_tracks(
            self.dancoff_num_azimuthal_angles,
            self.dancoff_track_spacing,
            self.dancoff_polar_quadrature,
        )
        moc.sim_mode = SimulationMode.FixedSource
        moc.flux_tolerance = 1.0e-5
        moc.solve()

        # Now we need to calculate the dancoff correction for each pin
        self.fuel_dancoff_corrections = []
        u = Direction(1.0, 0.0)
        for j in range(self.shape[1]):
            y = moc.y_max - (j + 0.5) * self.pitch
            for i in range(self.shape[0]):
                x = moc.x_min + (i + 0.5) * self.pitch
                r = Vector(x, y)
                xs = moc.xs(r, u)
                if xs.name == "Fuel":
                    flux = moc.flux(r, u, 0)
                    C = (iso_flux - flux) / iso_flux
                    self.fuel_dancoff_corrections.append(C)
                else:
                    self.fuel_dancoff_corrections.append(0.0)

        set_logging_level(LogLevel.Info)
    
    def _self_shield_clad(self):
        scarabee_log(LogLevel.Info, "")
        scarabee_log(LogLevel.Info, "Performing self-shielding for cladding")
        set_logging_level(LogLevel.Warning)

    def _get_clad_dancoff_corrections(self):
        scarabee_log(LogLevel.Info, "")
        scarabee_log(LogLevel.Info, "Computing Dancoff factors for cladding")
        set_logging_level(LogLevel.Warning)
        # ----------------------------------------------------------
        # ISOLATED FUEL PIN
        # We first make the system for an isolated fuel pin.
        # We isolate it by multiplying the pitch by 20.
        isolated_fp = None
        for pin in self.pins:
            if isinstance(pin, FuelPin):
                isolated_fp = pin.make_clad_dancoff_cell(
                    pitch=self.dancoff_isolation_factor * self.pitch,
                    moderator=self.moderator,
                )
                break
        if isolated_fp is None:
            raise RuntimeError("No FuelPin type found in pins.")
        iso_geom_fp = Cartesian2D(
            [self.dancoff_isolation_factor * self.pitch],
            [self.dancoff_isolation_factor * self.pitch],
        )
        iso_geom_fp.set_tiles([isolated_fp])
        iso_moc_fp = MOCDriver(iso_geom_fp)

        # Set the source
        for i in range(iso_moc_fp.nfsr):
            i_xs = iso_moc_fp.xs(i)
            if i_xs.name != "Clad":
                # If we aren't in the fuel, set the source to be the value of
                # the potential xs (should be Et).
                iso_moc_fp.set_extern_src(i, 0, i_xs.Et(0))

        # Solve the isolated pin problem
        iso_moc_fp.x_min_bc = BoundaryCondition.Vacuum
        iso_moc_fp.x_max_bc = BoundaryCondition.Vacuum
        iso_moc_fp.y_min_bc = BoundaryCondition.Vacuum
        iso_moc_fp.y_max_bc = BoundaryCondition.Vacuum
        iso_moc_fp.generate_tracks(
            self.dancoff_num_azimuthal_angles,
            self.dancoff_track_spacing,
            self.dancoff_polar_quadrature,
        )
        iso_moc_fp.sim_mode = SimulationMode.FixedSource
        iso_moc_fp.flux_tolerance = 1.0e-5
        iso_moc_fp.solve()
        for i in range(iso_moc_fp.nfsr):
            i_xs = iso_moc_fp.xs(i)
            if i_xs.name == "Clad":
                iso_flux_fp = iso_moc_fp.flux(i, 0)
                break

        # ----------------------------------------------------------
        # ISOLATED GUIDE TUBE
        # We first make the system for an isolated fuel pin.
        # We isolate it by multiplying the pitch by 20.
        isolated_gt = None
        for pin in self.pins:
            if isinstance(pin, GuideTube):
                isolated_gt = pin.make_clad_dancoff_cell(
                    pitch=self.dancoff_isolation_factor * self.pitch,
                    moderator=self.moderator,
                )
                break
        if isolated_gt is not None:
            iso_geom_gt = Cartesian2D(
                [self.dancoff_isolation_factor * self.pitch],
                [self.dancoff_isolation_factor * self.pitch],
            )
            iso_geom_gt.set_tiles([isolated_gt])
            iso_moc_gt = MOCDriver(iso_geom_gt)

            # Set the source
            for i in range(iso_moc_gt.nfsr):
                i_xs = iso_moc_gt.xs(i)
                if i_xs.name != "Clad":
                    # If we aren't in the fuel, set the source to be the value of
                    # the potential xs (should be Et).
                    iso_moc_gt.set_extern_src(i, 0, i_xs.Et(0))

            # Solve the isolated pin problem
            iso_moc_gt.x_min_bc = BoundaryCondition.Vacuum
            iso_moc_gt.x_max_bc = BoundaryCondition.Vacuum
            iso_moc_gt.y_min_bc = BoundaryCondition.Vacuum
            iso_moc_gt.y_max_bc = BoundaryCondition.Vacuum
            iso_moc_gt.generate_tracks(
                self.dancoff_num_azimuthal_angles,
                self.dancoff_track_spacing,
                self.dancoff_polar_quadrature,
            )
            iso_moc_gt.sim_mode = SimulationMode.FixedSource
            iso_moc_gt.flux_tolerance = 1.0e-5
            iso_moc_gt.solve()
            for i in range(iso_moc_gt.nfsr):
                i_xs = iso_moc_gt.xs(i)
                if i_xs.name == "Clad":
                    iso_flux_gt = iso_moc_gt.flux(i, 0)
                    break

        # ---------------------------------------------------
        # Now we setup the lattice problem
        fuel_df_pins = []
        for pin in self.pins:
            fuel_df_pins.append(
                pin.make_clad_dancoff_cell(pitch=self.pitch, moderator=self.moderator)
            )
        geom = Cartesian2D(self.shape[0] * [self.pitch], self.shape[1] * [self.pitch])
        geom.set_tiles(fuel_df_pins)
        moc = MOCDriver(geom)

        # Set the source
        for i in range(moc.nfsr):
            i_xs = moc.xs(i)
            if i_xs.name != "Clad":
                # If we aren't in the fuel, set the source to be the value of
                # the potential xs (should be Et).
                moc.set_extern_src(i, 0, i_xs.Et(0))

        # Solve the lattice problem
        moc.generate_tracks(
            self.dancoff_num_azimuthal_angles,
            self.dancoff_track_spacing,
            self.dancoff_polar_quadrature,
        )
        moc.sim_mode = SimulationMode.FixedSource
        moc.flux_tolerance = 1.0e-5
        moc.solve()

        # Now we need to calculate the dancoff correction for each pin
        self.clad_dancoff_corrections = []
        u = Direction(1.0, 0.0)
        i_pin = 0
        for j in range(self.shape[1]):
            y = moc.y_max - (j + 0.5) * self.pitch
            for i in range(self.shape[0]):
                x = moc.x_min + (i + 0.5) * self.pitch
                pin = self.pins[i_pin]
                r = Vector(x, y) + pin.clad_offset()
                xs = moc.xs(r, u)

                flux = moc.flux(r, u, 0)

                if xs.name == "Clad":
                    if isinstance(pin, FuelPin):
                        C = (iso_flux_fp - flux) / iso_flux_fp
                    elif isinstance(pin, GuideTube):
                        C = (iso_flux_gt - flux) / iso_flux_gt
                else:
                    C = 0.0
                self.clad_dancoff_corrections.append(C)

                i_pin += 1

        set_logging_level(LogLevel.Info)

    def _pin_cell_calc(self):
        scarabee_log(LogLevel.Info, "")
        scarabee_log(LogLevel.Info, "Performing micro-group pin cell calcuations")
        set_logging_level(LogLevel.Warning)
        self.pin_1d_cells = len(self.pins) * [None]
        self.pin_1d_fluxes = len(self.pins) * [None]

        # First, we do all of the fuel pins
        nfp = 0
        avg_fp = None
        for i in range(len(self.pins)):
            if isinstance(self.pins[i], FuelPin):
                set_logging_level(LogLevel.Info)
                scarabee_log(
                    LogLevel.Info, "  Calculating fuel pin at index {:}".format(i)
                )
                set_logging_level(LogLevel.Warning)

                fuel_dancoff = self.fuel_dancoff_corrections[i]
                clad_dancoff = self.clad_dancoff_corrections[i]
                self.pin_1d_cells[i] = self.pins[i].make_cylindrical_cell(
                    pitch=self.pitch,
                    dancoff_fuel=fuel_dancoff,
                    moderator=self.moderator_xs,
                    ndl=self.ndl,
                    dancoff_clad=clad_dancoff,
                )

                self.pin_1d_cells[i].solve()
                self.pin_1d_fluxes[i] = CylindricalFluxSolver(self.pin_1d_cells[i])
                self.pin_1d_fluxes[i].solve()

                nfp += 1
                if avg_fp is None:
                    avg_fp = self.pin_1d_fluxes[i].homogenize()
                else:
                    avg_fp += self.pin_1d_fluxes[i].homogenize()
        avg_fp *= 1.0 / float(nfp)
        self.average_fuel_pin = avg_fp

        # Now we do all the non fuel pin cells
        buffer_rad = np.sqrt(9.0 * self.pitch * self.pitch / np.pi)
        for i in range(len(self.pins)):
            if not isinstance(self.pins[i], FuelPin):
                set_logging_level(LogLevel.Info)
                scarabee_log(LogLevel.Info, "  Calculating cell at index {:}".format(i))
                set_logging_level(LogLevel.Warning)

                clad_dancoff = self.clad_dancoff_corrections[i]
                self.pin_1d_cells[i] = self.pins[i].make_cylindrical_cell(
                    pitch=self.pitch,
                    moderator=self.moderator_xs,
                    ndl=self.ndl,
                    dancoff_clad=clad_dancoff,
                    buffer=avg_fp,
                    buffer_radius=buffer_rad,
                )

                self.pin_1d_cells[i].solve()
                self.pin_1d_fluxes[i] = CylindricalFluxSolver(self.pin_1d_cells[i])
                self.pin_1d_fluxes[i].solve()

        set_logging_level(LogLevel.Info)

    def _condense_xs(self):
        scarabee_log(LogLevel.Info, "")
        scarabee_log(LogLevel.Info, "Performing pin cell energy condensation")
        set_logging_level(LogLevel.Warning)
        for i in range(len(self.pins)):
            pin = self.pins[i]
            cell_flux = self.pin_1d_fluxes[i]

            NR = cell_flux.nregions
            if not isinstance(pin, FuelPin):
                NR -= 1

            pin.condensed_xs = []
            for r in range(NR):
                xs = cell_flux.xs(r)
                flux_spectrum = cell_flux.homogenize_flux_spectrum([r])
                pin.condensed_xs.append(
                    xs.condense(self.condensation_scheme, flux_spectrum)
                )
                pin.condensed_xs[-1].name = xs.name
        set_logging_level(LogLevel.Info)

    def _moc(self):
        scarabee_log(LogLevel.Info, "")
        scarabee_log(LogLevel.Info, "Performing macrogroup assembly calculation")
        moc_pins = len(self.pins) * [None]
        for i in range(len(self.pins)):
            moc_pins[i] = self.pins[i].make_moc_cell(self.pitch)

        dx = self.shape[0] * [self.pitch]
        dy = self.shape[1] * [self.pitch]
        self.moc_geom = Cartesian2D(dx, dy)
        self.moc_geom.set_tiles(moc_pins)

        self.moc = MOCDriver(self.moc_geom)
        if self.plot_assembly:
            self.moc.plot()
        self.moc.generate_tracks(
            self.num_azimuthal_angles, self.track_spacing, self.polar_quadrature
        )
        self.moc.keff_tolerance = self.keff_tolerance
        self.moc.flux_tolerance = self.flux_tolerance
        self.moc.solve()

    def _criticality_spectrum_calc(self):
        if self.criticality_spectrum_method is None:
            self._criticality_spectrum = None
            return

        scarabee_log(LogLevel.Info, "")
        scarabee_log(
            LogLevel.Info,
            "Performing {:} criticality spectrum calculation".format(
                self.criticality_spectrum_method
            ),
        )
        homogenized_moc = self.moc.homogenize()

        if self.criticality_spectrum_method in ["P1", "p1"]:
            self._criticality_spectrum = P1CriticalitySpectrum(homogenized_moc)
        else:
            self._criticality_spectrum = B1CriticalitySpectrum(homogenized_moc)

        self.moc.apply_criticality_spectrum(self._criticality_spectrum.flux)

        scarabee_log(
            LogLevel.Info, "Kinf    : {:.5f}".format(self._criticality_spectrum.k_inf)
        )
        scarabee_log(
            LogLevel.Info, "Buckling: {:.5f}".format(self._criticality_spectrum.B2)
        )

    def _few_group_xs(self):
        scarabee_log(LogLevel.Info, "")
        scarabee_log(LogLevel.Info, "Generating few group cross sections")

        homog_xs = self.moc.homogenize()
        flux_spectrum = self.moc.homogenize_flux_spectrum()
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

        # According to Smith, one should do energy condensation on the
        # diffusion coefficients, and not on the transport cross sections which
        # one could then use to make diffusion coefficients [1]. This is in
        # contradiction to Lattice Physics Computations which states that
        # either method is acceptable [2]. In light of these comments, I have
        # chosen to go with Smith's recommendation of performing energy
        # condensation on the diffusion coefficients.
        self.diffusion_xs = diff_xs.condense(
            self.few_group_condensation_scheme, flux_spectrum
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
# [1] K. S. Smith, “Nodal diffusion methods and lattice physics data in LWR
#     analyses: Understanding numerous subtle details,” Prog Nucl Energ,
#     vol. 101, pp. 360–369, 2017, doi: 10.1016/j.pnucene.2017.06.013.
# [2] D. Knott and A. Yamamoto, "Lattice Physics Computations" in
#     Handbook of Nuclear Engineering, 2010, p 1226.
