from _scarabee import *
import numpy as np
from typing import Tuple, List, Optional
from enum import Enum

class FuelPin:
    def __init__(self, fuel: Material, fuel_radius: float, clad: Material, clad_width: float, gap: Optional[Material] = None, gap_width: Optional[float] = None):
        self.fuel = fuel
        self.fuel_radius = fuel_radius
        self.clad = clad
        self.clad_width = clad_width
        self.gap = gap
        self.gap_width = gap_width
        self.condensed_xs = []

        if self.gap is not None and self.gap_width is None:
            raise RuntimeError('Fuel gap material is provided, but no gap width.')
        elif self.gap_width is not None and self.gap is None:
            raise RuntimeError('Fuel gap width provided, but no gap material.')

    def clad_offset(self):
        if self.gap_width is not None:
            return Vector(self.fuel_radius + self.gap_width + 0.5*self.clad_width, 0.)
        else:
            return Vector(self.fuel_radius + 0.5*self.clad_width, 0.)

    def make_fuel_dancoff_cell(self, pitch: float, moderator: Material):
        # We first determine all the radii
        radii = []
        radii.append(self.fuel_radius)
        if self.gap is not None:
            radii.append(radii[-1] + self.gap_width)
        radii.append(radii[-1] + self.clad_width)

        mats = []

        # This returns a cell for calculating the fuel pin Dancoff factor.
        # As such, the fuel XS has infinite values.
        Et = np.array([1.E5])
        Ea = np.array([1.E5])
        Es = np.array([[0.]])
        Fuel = CrossSection(Et, Ea, Es, "Fuel")
        mats.append(Fuel)

        if self.gap is not None:
            Et[0] = self.gap.potential_xs
            Ea[0] = self.gap.potential_xs
            gap = CrossSection(Et, Ea, Es, "Gap")
            mats.append(Gap)

        Et[0] = self.clad.potential_xs
        Ea[0] = self.clad.potential_xs
        Clad = CrossSection(Et, Ea, Es, "Clad")
        mats.append(Clad)

        Et[0] = moderator.potential_xs
        Ea[0] = moderator.potential_xs
        Mod = CrossSection(Et, Ea, Es, "Moderator")
        mats.append(Mod)

        return SimplePinCell(radii, mats, pitch, pitch)

    def make_clad_dancoff_cell(self, pitch: float, moderator: Material):
        # We first determine all the radii
        radii = []
        radii.append(self.fuel_radius)
        if self.gap is not None:
            radii.append(radii[-1] + self.gap_width)
        radii.append(radii[-1] + self.clad_width)

        mats = []
 
        Et = np.array([self.fuel.potential_xs])
        Ea = np.array([self.fuel.potential_xs])
        Es = np.array([[0.]])
        Fuel = CrossSection(Et, Ea, Es, "Fuel")
        mats.append(Fuel)

        if self.gap is not None:
            Et[0] = self.gap.potential_xs
            Ea[0] = self.gap.potential_xs
            gap = CrossSection(Et, Ea, Es, "Gap")
            mats.append(Gap)
        
        # This returns a cell for calculating the fuel pin Dancoff factor.
        # As such, the clad XS has infinite values.
        Et[0] = 1.E5
        Ea[0] = 1.E5
        Clad = CrossSection(Et, Ea, Es, "Clad")
        mats.append(Clad)

        Et[0] = moderator.potential_xs
        Ea[0] = moderator.potential_xs
        Mod = CrossSection(Et, Ea, Es, "Moderator")
        mats.append(Mod)

        return SimplePinCell(radii, mats, pitch, pitch)

    def make_cylindrical_cell(self, pitch: float, dancoff_fuel: float, moderator: CrossSection, ndl: NDLibrary, dancoff_clad: Optional[float] = None, clad_dilution = 1.E10):
        # We first determine all the radii
        radii = []
        radii.append(self.fuel_radius)
        if self.gap is not None:
            radii.append(radii[-1] + self.gap_width)
        radii.append(radii[-1] + self.clad_width)
        radii.append(np.sqrt(pitch*pitch / np.pi))
        
        # Next, we determine all the materials.
        # This requires applying self shielding to the fuel and cladding
        mats = []

        # First, treat the fuel
        Ee = 1. / (2. * self.fuel_radius) # Fuel escape xs
        mats.append(self.fuel.carlvik_xs(dancoff_fuel, Ee, ndl))
        mats[-1].name = "Fuel"

        # Next, add the gap (if present)
        if self.gap is not None:
            mats.append(self.gap.dilution_xs(self.gap.size*[1.E10], ndl))
            mats[-1].name = "Gap"

        # Add the cladding
        if dancoff_clad is not None:
            Ee = 1. / (2. * self.clad_width)
            mats.append(self.clad.roman_xs(dancoff_clad, Ee, ndl))
        else:
            mats.append(self.clad.dilution_xs(self.clad.size*[clad_dilution], ndl))
        mats[-1].name = "Clad"

        # Finally, add moderator
        mats.append(moderator)

        return CylindricalCell(radii, mats)

    def make_moc_cell(self, pitch: float):
        radii = []
        radii.append(self.fuel_radius)
        if self.gap is not None:
            radii.append(radii[-1] + self.gap_width)
        radii.append(radii[-1] + self.clad_width)

        mod_width = 0.5*pitch - radii[-1]
        radii.append(radii[-1] + 0.8*mod_width)

        mats = self.condensed_xs.copy()
        mats.append(self.condensed_xs[-1])

        return PinCell(radii, mats, pitch, pitch)

class GuideTube:
    def __init__(self, inner_radius: float, outer_radius: float, clad: Material):
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.clad = clad
        self.condensed_xs = []

        if self.outer_radius <= self.inner_radius:
            raise RuntimeError("Outer radius must be > inner radius.")
    
    def clad_offset(self):
        return Vector(0.5*(self.inner_radius + self.outer_radius), 0.)

    def make_clad_dancoff_cell(self, pitch: float, moderator: Material):
        # We first determine all the radii
        radii = []
        radii.append(self.inner_radius)
        radii.append(self.outer_radius)

        mats = []
 
        Et = np.array([moderator.potential_xs])
        Ea = np.array([moderator.potential_xs])
        Es = np.array([[0.]])
        Mod = CrossSection(Et, Ea, Es, "Moderator")
        mats.append(Mod)

        # This returns a cell for calculating the fuel pin Dancoff factor.
        # As such, the clad XS has infinite values.
        Et[0] = 1.E5
        Ea[0] = 1.E5
        Clad = CrossSection(Et, Ea, Es, "Clad")
        mats.append(Clad)

        mats.append(Mod)

        return SimplePinCell(radii, mats, pitch, pitch)

    def make_fuel_dancoff_cell(self, pitch: float, moderator: Material):
        # We first determine all the radii
        radii = []
        radii.append(self.inner_radius)
        radii.append(self.outer_radius)

        mats = []
 
        Et = np.array([moderator.potential_xs])
        Ea = np.array([moderator.potential_xs])
        Es = np.array([[0.]])
        Mod = CrossSection(Et, Ea, Es, "Moderator")
        mats.append(Mod)

        Et[0] = self.clad.potential_xs
        Ea[0] = self.clad.potential_xs
        Clad = CrossSection(Et, Ea, Es, "Clad")
        mats.append(Clad)

        mats.append(Mod)

        return SimplePinCell(radii, mats, pitch, pitch)

    def make_cylindrical_cell(self, pitch: float, moderator: CrossSection, buffer_radius: float, buffer: CrossSection, ndl: NDLibrary, dancoff_clad: Optional[float] = None, clad_dilution: float = 1.E10):
        # We first determine all the radii
        radii = []
        radii.append(self.inner_radius)
        radii.append(self.outer_radius)
        radii.append(np.sqrt(pitch*pitch / np.pi))
        if radii[-1] >= buffer_radius:
            raise RuntimeError("Buffer radius is smaller than the radius of the cell.")
        radii.append(buffer_radius)

        # Next, we determine all the materials.
        # This requires applying self shielding to the fuel and cladding
        mats = []

        mats.append(moderator)

        # Add the cladding
        if dancoff_clad is not None:
            Ee = 1. / (2. * (self.outer_radius - self.inner_radius))
            mats.append(self.clad.roman_xs(dancoff_clad, Ee, ndl))
        else:
            mats.append(self.clad.dilution_xs(self.clad.size*[clad_dilution], ndl))
        mats[-1].name = "Clad"

        # Add outer moderator
        mats.append(moderator)

        # Add the buffer
        mats.append(buffer)

        return CylindricalCell(radii, mats)

    def make_moc_cell(self, pitch: float):
        r_inner_inner_mod = np.sqrt(0.5*self.inner_radius*self.inner_radius)
        radii = [r_inner_inner_mod, self.inner_radius, self.outer_radius]
        mats = [self.condensed_xs[0]] + self.condensed_xs.copy()
        return PinCell(radii, mats, pitch, pitch)

class PWRAssembly:
    def __init__(self, pitch: float, gap: float, moderator: Material, shape: Tuple[int, int], ndl: NDLibrary):
        self.pitch = pitch
        self.gap = gap
        self.ndl = ndl
        self.moderator = moderator
        self.shape = shape
        self.pins = []
        self.condensation_scheme = []
        
        self.moderator_xs = moderator.dilution_xs(moderator.size*[1.E10], self.ndl)
        self.moderator_xs.name = "Moderator"
        
        # MOC parameters for computing dancoff corrections
        self.dt_dancoff = 0.05
        self.na_dancoff = 64
        self.iso_dancoff = 20.0

        # MOC parameters for assembly calculation
        self.dt = 0.05
        self.na = 32
        self.keff_tolerance = 1.E-5
        self.flux_tolerance = 1.E-5
        self.plot_assembly = True
        self.moc_geom = None
        self.moc = None

    def solve(self):
        self._get_fuel_dancoff_corrections()
        self._get_clad_dancoff_corrections()
        self._pin_cell_calc()
        self._condense_xs()
        self._moc()

        #print()
        #print()
        #print(self.pin_1d_cells[0].radius(self.pin_1d_cells[0].nregions-1))
        #homog_xs = self.pin_1d_fluxes[0].homogenize()
        #rad = np.sqrt(self.pitch*self.pitch/np.pi)
        #print(rad)
        #homo_cell = CylindricalCell([0.25*rad, 0.5*rad, 0.75*rad, rad], [homog_xs, homog_xs, homog_xs, homog_xs])
        #homo_cell.solve()
        #homo_flux = CylindricalFluxSolver(homo_cell)
        #homo_flux.solve()

    def _get_fuel_dancoff_corrections(self):
        # We first make the system for an isolated fuel pin.
        # We isolate it by multiplying the pitch by 20.
        isolated_fp = None
        for pin in self.pins:
            if isinstance(pin, FuelPin):
                isolated_fp = pin.make_fuel_dancoff_cell(pitch=self.iso_dancoff*self.pitch, moderator=self.moderator)
        if isolated_fp is None:
            raise RuntimeError("No FuelPin type found in pins.")
        iso_geom = Cartesian2D([self.iso_dancoff*self.pitch], [self.iso_dancoff*self.pitch])  
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
        iso_moc.generate_tracks(self.na_dancoff, self.dt_dancoff, YamamotoTabuchi6())
        iso_moc.sim_mode = SimulationMode.FixedSource
        iso_moc.flux_tolerance = 1.E-5
        iso_moc.solve()
        iso_flux = iso_moc.flux(0, 0)

        # Now we setup the lattice problem
        fuel_df_pins = []
        for pin in self.pins:
            fuel_df_pins.append(pin.make_fuel_dancoff_cell(pitch=self.pitch, moderator=self.moderator))
        geom = Cartesian2D(self.shape[0]*[self.pitch], self.shape[1]*[self.pitch]) 
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
        moc.generate_tracks(self.na_dancoff, self.dt_dancoff, YamamotoTabuchi6())
        moc.sim_mode = SimulationMode.FixedSource
        moc.flux_tolerance = 1.E-5
        moc.solve()

        # Now we need to calculate the dancoff correction for each pin
        self.fuel_dancoff_corrections = []
        u = Direction(1., 0.)
        for j in range(self.shape[1]):
            y = moc.y_max - self.gap - (j+0.5)*self.pitch
            for i in range(self.shape[0]):
                x = moc.x_min + self.gap + (i+0.5)*self.pitch
                r = Vector(x, y)
                xs = moc.xs(r, u)
                if xs.name == "Fuel":
                    flux = moc.flux(r, u, 0)
                    C = (iso_flux - flux) / iso_flux
                    self.fuel_dancoff_corrections.append(C)
                else:
                    self.fuel_dancoff_corrections.append(0.)

    def _get_clad_dancoff_corrections(self):
        #----------------------------------------------------------
        # ISOLATED FUEL PIN
        # We first make the system for an isolated fuel pin.
        # We isolate it by multiplying the pitch by 20.
        isolated_fp = None
        for pin in self.pins:
            if isinstance(pin, FuelPin):
                isolated_fp = pin.make_clad_dancoff_cell(pitch=self.iso_dancoff*self.pitch, moderator=self.moderator)
                break
        if isolated_fp is None:
            raise RuntimeError("No FuelPin type found in pins.")
        iso_geom_fp = Cartesian2D([self.iso_dancoff*self.pitch], [self.iso_dancoff*self.pitch])  
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
        iso_moc_fp.generate_tracks(self.na_dancoff, self.dt_dancoff, YamamotoTabuchi6())
        iso_moc_fp.sim_mode = SimulationMode.FixedSource
        iso_moc_fp.flux_tolerance = 1.E-5
        iso_moc_fp.solve()
        for i in range(iso_moc_fp.nfsr):
            i_xs = iso_moc_fp.xs(i)
            if i_xs.name == "Clad":
                iso_flux_fp = iso_moc_fp.flux(i, 0)
                break

        #----------------------------------------------------------
        # ISOLATED GUIDE TUBE
        # We first make the system for an isolated fuel pin.
        # We isolate it by multiplying the pitch by 20.
        isolated_gt = None
        for pin in self.pins:
            if isinstance(pin, GuideTube):
                isolated_gt = pin.make_clad_dancoff_cell(pitch=self.iso_dancoff*self.pitch, moderator=self.moderator)
                break
        if isolated_gt is not None:
            iso_geom_gt = Cartesian2D([self.iso_dancoff*self.pitch], [self.iso_dancoff*self.pitch])  
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
            iso_moc_gt.generate_tracks(self.na_dancoff, self.dt_dancoff, YamamotoTabuchi6())
            iso_moc_gt.sim_mode = SimulationMode.FixedSource
            iso_moc_gt.flux_tolerance = 1.E-5
            iso_moc_gt.solve()
            for i in range(iso_moc_gt.nfsr):
                i_xs = iso_moc_gt.xs(i)
                if i_xs.name == "Clad":
                    iso_flux_gt = iso_moc_gt.flux(i, 0)
                    break
        
        #---------------------------------------------------
        # Now we setup the lattice problem
        fuel_df_pins = []
        for pin in self.pins:
            fuel_df_pins.append(pin.make_clad_dancoff_cell(pitch=self.pitch, moderator=self.moderator))
        geom = Cartesian2D(self.shape[0]*[self.pitch], self.shape[1]*[self.pitch]) 
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
        moc.generate_tracks(self.na_dancoff, self.dt_dancoff, YamamotoTabuchi6())
        moc.sim_mode = SimulationMode.FixedSource
        moc.flux_tolerance = 1.E-5
        moc.solve()

        # Now we need to calculate the dancoff correction for each pin
        self.clad_dancoff_corrections = []
        u = Direction(1., 0.)
        i_pin = 0
        for j in range(self.shape[1]):
            y = moc.y_max - self.gap - (j+0.5)*self.pitch
            for i in range(self.shape[0]):
                x = moc.x_min + self.gap + (i+0.5)*self.pitch
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
                    C = 0.
                self.clad_dancoff_corrections.append(C)

                i_pin += 1

    def _pin_cell_calc(self):
        self.pin_1d_cells = len(self.pins) * [None]
        self.pin_1d_fluxes = len(self.pins) * [None]

        # First, we do all of the fuel pins
        nfp = 0
        avg_fp = None
        for i in range(len(self.pins)):
            if isinstance(self.pins[i], FuelPin):
                fuel_dancoff = self.fuel_dancoff_corrections[i]
                clad_dancoff = self.clad_dancoff_corrections[i]
                self.pin_1d_cells[i] = self.pins[i].make_cylindrical_cell(pitch=self.pitch, dancoff_fuel=fuel_dancoff, moderator=self.moderator_xs, ndl=self.ndl, dancoff_clad=clad_dancoff)

                print("Solving pin at index {:}".format(i))
                self.pin_1d_cells[i].solve()
                self.pin_1d_fluxes[i] = CylindricalFluxSolver(self.pin_1d_cells[i])
                self.pin_1d_fluxes[i].solve()
                print()

                nfp += 1
                if avg_fp is None:
                    avg_fp = self.pin_1d_fluxes[i].homogenize()
                else:
                    avg_fp += self.pin_1d_fluxes[i].homogenize()
        avg_fp *= 1. / float(nfp)

        # Now we do all the non fuel pin cells
        buffer_rad = np.sqrt(9.*self.pitch*self.pitch / np.pi)
        for i in range(len(self.pins)):
            if not isinstance(self.pins[i], FuelPin):
                clad_dancoff = self.clad_dancoff_corrections[i]
                self.pin_1d_cells[i] = self.pins[i].make_cylindrical_cell(pitch=self.pitch, moderator=self.moderator_xs, ndl=self.ndl, dancoff_clad=clad_dancoff, buffer=avg_fp, buffer_radius=buffer_rad)
                
                print("Solving pin at index {:}".format(i))
                self.pin_1d_cells[i].solve()
                self.pin_1d_fluxes[i] = CylindricalFluxSolver(self.pin_1d_cells[i])
                self.pin_1d_fluxes[i].solve()
                print()
    
    def _condense_xs(self):
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
                pin.condensed_xs.append(xs.condense(self.condensation_scheme, flux_spectrum))
                pin.condensed_xs[-1].name = xs.name

    def _moc(self):
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
        self.moc.generate_tracks(self.na, self.dt, YamamotoTabuchi6())
        self.moc.keff_tolerance = self.keff_tolerance
        self.moc.flux_tolerance = self.flux_tolerance
        self.moc.solve()
