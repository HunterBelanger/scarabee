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

        if self.gap is not None and self.gap_width is None:
            raise RuntimeError('Fuel gap material is provided, but no gap width.')
        elif self.gap_width is not None and self.gap is None:
            raise RuntimeError('Fuel gap width provided, but no gap material.')

    def make_cylindrical_cell(self, pitch: float, dancoff_fuel: float, moderator: CrossSection, ndl: NDLibrary, dancoff_clad: Optional[float] = None, clad_dillution = 1.E10):
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

        # Next, add the gap (if present)
        if self.gap is not None:
            mats.append(self.gap.dilution_xs(self.gap.size*[1.E10], ndl))

        # Add the cladding
        if dancoff_clad is not None:
            Ee = 1. / (2. * self.clad_width)
            mats.append(self.clad.roman_xs(dancoff_clad, Ee, ndl))
        else:
            mats.append(self.clad.dilution_xs(self.clad.size*[clad_dillution], ndl))

        # Finally, add moderator
        mats.append(moderator.dilution_xs(moderator.size*[1.E10], ndl))

        return CylindricalCell(radii, mats)


class PWRAssembly:
    def __init__(self, pitch: float, gap: float, moderator: Material, shape: Tuple[int, int]):
        self.pitch = pitch
        self.gap = gap
        self.moderator = moderator
        self.shape = shape

