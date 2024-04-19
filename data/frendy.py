import ENDFtk
import shutil
import subprocess
import os
import numpy as np

class FrendyFreeGasMG:
  def __init__(self):
    self.temp = 293.6
    self.dilutions = []
    self.endf_file = ""
    self.label = ""
    self.name = ""

  def _get_endf_info(self):
    # First, get MAT
    tape = ENDFtk.tree.Tape.from_file(self.endf_file)
    self.mat = tape.material_numbers[0]
    mat = tape.MAT(self.mat).parse()

    # Read AWR and ZA from MF1 MT 451
    mf1mt451 = mat.MF(1).MT(451)
    self.awr = mf1mt451.AWR
    self.ZA = mf1mt451.ZA
    self.fissile = mf1mt451.LFI

    # Get potential scattering xs
    mf2mt151 = mat.MF(2).MT(151).isotopes[0]
    rrr = mf2mt151.resonance_ranges[0]
    if rrr.energy_dependent_scattering_radius:
      raise RuntimeError("HELP ENERGY DEPENDENT SCATTERING RADIUS")
    rrr_params = rrr.parameters
    AP = rrr_params.AP
    if AP == 0.:
      # Try getting from the l-values
      AP = rrr_params.l_values[0].APL

    self.pot_xs = 4. * np.pi * AP * AP

  def _frendy_input(self):
    if len(self.dilutions) > 0:
      if self.dilutions[-1] < 1.E10:
        self.dilutions.append(1.E10)

    dil_frmt = len(self.dilutions)*"{:.2E} "
    out = "mg_neutron_mode\n"
    out += "mg_edit_option ( \"1DXS 2 4 16 17 37 22 28 32 33 34 18 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 155 182 191 192 193 197\" \"2DXS 2 4 16 17 37 22 28 32 33 34\" NuChi MGFlux )\n"
    out += "nucl_file_name ({endf})\n".format(endf=self.endf_file)
    out += "mg_file_name {mgfname}\n".format(mgfname=self.name)
    out += "temperature {temp}\n".format(temp=self.temp)
    out += "legendre_order 1\n"
    out += "mg_structure ( shem-cea-281 )\n"
    out += "mg_weighting_spectrum ( fission+1/e+maxwell  )\n"
    if len(self.dilutions) > 0:
      out += "sigma_zero_data ( " + dil_frmt.format(*self.dilutions) + " )\n"
    return out

  def process(self):
    self._get_endf_info()

    frendy_input = self._frendy_input()
    with open("frendy_input", 'w') as fl:
      fl.write(frendy_input)

    subprocess.run(['frendy', 'frendy_input'])

    os.remove('frendy_input')
    os.remove('FMAlternateInputData.txt')
    os.remove(self.endf_file+".ace")
    os.remove(self.endf_file+".ace.dir")

def read_1dxs(fname, nskip=3):
  fl = open(fname, 'r')

  # Skip first lines that have headers / dilutions / temperatures
  for i in range(nskip):
    fl.readline()

  array = []

  for line in fl:
    line = line.strip()
    if len(line) == 0:
      continue

    line = line.split()[4:]
    line.reverse() # Reverse line for dilutions to go from low to high
    for i in range(len(line)):
      line[i] = float(line[i])
    array.append(line)
  fl.close()

  array = np.array(array)
  array = np.swapaxes(array, 0, 1)

  # First index on dilution, second on group
  return array

def read_2dxs(fname, lorder=1):
  fl = open(fname, 'r')

  fl.readline()
  fl.readline()
  fl.readline()
  fl.readline() # Skip P0 header

  array = []
  for l in range(lorder+1):
    lary = []

    line = fl.readline().strip()
    while len(line) > 0:
      line = line.split()[1:]
      for i in range(len(line)):
        line[i] = float(line[i])
      lary.append(line)

      line = fl.readline().strip()
    
    array.append(lary)

    # Skip new Pl header
    fl.readline()
  fl.close()

  array = np.array(array)

  return array


