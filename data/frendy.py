import ENDFtk
import shutil
import subprocess
import os
import numpy as np

_ABS_MT_LIST = [102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 155, 182, 191, 192, 193, 197]
_2D_MT_LIST = [2, 4, 16, 17, 37, 22, 28, 32, 33, 34]
_YIELDS = {2: 1., 4: 1., 16: 2., 17: 3., 37: 4., 22: 1., 28: 1., 32: 1., 33: 1., 34: 1.}

class FrendyMG:
  def __init__(self):
    self.temps = [293.6]
    self.dilutions = [1.E10]
    self.endf_file = None
    self.tsl_file = None
    self.label = ""
    self.name = ""
    self.ngroups = 281
    self.initialized = False
    self.processed = False
    self.resonant = False

  def initialize(self):
    self._get_endf_info()

    if len(self.dilutions) > 1:
      self.resonant = True

    self.Es = np.zeros((len(self.temps), len(self.dilutions), self.ngroups, self.ngroups), dtype=np.float32)
    self.Es1 = np.zeros((len(self.temps), len(self.dilutions), self.ngroups, self.ngroups), dtype=np.float32)
    self.Ea = np.zeros((len(self.temps), len(self.dilutions), self.ngroups), dtype=np.float32)
    if self.fissile:
      self.Ef = np.zeros((len(self.temps), len(self.dilutions), self.ngroups), dtype=np.float32)
      self.nu = np.zeros((len(self.temps), self.ngroups), dtype=np.float32)
      self.chi = np.zeros((len(self.temps), self.ngroups), dtype=np.float32)
    else:
      self.Ef = None
      self.nu = None
      self.chi = None

    self.initialized = True

  def process(self, h5=None):
    if not self.initialized:
      self.initialize()
    
    for i in range(len(self.temps)):
      self._process_temp(i)
    
    self.processed = True

    if h5 is not None:
      self.add_to_hdf5(h5)

  def add_to_hdf5(self, h5):
    grp = h5.create_group(self.name)

    # Save attributes
    grp.attrs['name'] = self.name
    grp.attrs['fissile'] = self.fissile
    grp.attrs['resonant'] = self.resonant
    grp.attrs['awr'] = self.awr
    grp.attrs['ZA'] = self.ZA
    grp.attrs['label'] = self.label
    grp.attrs['potential-xs'] = self.pot_xs
    grp.attrs['temperatures'] = self.temps
    if self.resonant:
      grp.attrs['dilutions'] = self.dilutions

    # Save cross section data
    grp.create_dataset("absorption", data=self.Ea)
    grp.create_dataset("scattering", data=self.Es)
    grp.create_dataset("p1-scattering", data=self.Es1)
    if self.fissile:
      grp.create_dataset("fission", data=self.Ef)
      grp.create_dataset("nu", data=self.nu)
      grp.create_dataset("chi", data=self.chi)

  def _get_endf_info(self):
    # First, get MAT
    tape = ENDFtk.tree.Tape.from_file(self.endf_file)
    self.mat = tape.material_numbers[0]

    # Read AWR and ZA from MF1 MT 451
    mf1mt451 = tape.MAT(self.mat).MF(1).MT(451).parse()
    self.awr = mf1mt451.AWR
    self.ZA = mf1mt451.ZA
    self.fissile = mf1mt451.is_fissile

    # Get potential scattering xs
    mf2mt151 = tape.MAT(self.mat).MF(2).MT(151).parse().isotopes[0]
    rrr = mf2mt151.resonance_ranges[0]
    if rrr.energy_dependent_scattering_radius:
      raise RuntimeError("HELP ENERGY DEPENDENT SCATTERING RADIUS")
    rrr_params = rrr.parameters
    AP = rrr_params.AP
    if AP == 0.:
      # Try getting from the l-values
      AP = rrr_params.l_values[0].APL

    self.pot_xs = 4. * np.pi * AP * AP

  def _frendy_input(self, temp):
    if len(self.dilutions) > 0:
      if self.dilutions[-1] < 1.E10:
        self.dilutions.append(1.E10)

    dil_frmt = len(self.dilutions)*"{:.2E} "
    out = "mg_neutron_mode\n"
    out += "mg_edit_option ( \"1DXS 18 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 155 182 191 192 193 197\" \"2DXS 2 4 16 17 37 22 28 32 33 34\" NuChi MGFlux )\n"
    out += "nucl_file_name ({endf})\n".format(endf=self.endf_file)
    if self.tsl_file is not None:
      out += "nucl_file_name_tsl ({tsl})\n".format(tsl=self.tsl_file)
    out += "mg_file_name {mgfname}\n".format(mgfname=self.name)
    out += "temperature {temp}\n".format(temp=temp)
    out += "legendre_order 1\n"
    out += "mg_structure ( shem-cea-281 )\n"
    out += "mg_weighting_spectrum ( fission+1/e+maxwell  )\n"
    out += "process_gas_xs off\n"
    if len(self.dilutions) > 0:
      out += "sigma_zero_data ( " + dil_frmt.format(*self.dilutions) + " )\n"
    return out

  def _process_temp(self, itemp):
    self._get_endf_info()

    temp = self.temps[itemp]
    frendy_input = self._frendy_input(temp)
    with open("frendy_input", 'w') as fl:
      fl.write(frendy_input)

    subprocess.run(['frendy', 'frendy_input'])

    self._read_temp(itemp)

    try:
      os.remove('frendy_input')
      os.remove('FMAlternateInputData.txt')
      os.remove(os.path.basename(self.endf_file)+".ace")
      os.remove(os.path.basename(self.endf_file)+".ace.dir")

      if self.tsl_file is not None:
        os.remove(os.path.basename(self.tsl_file)+".ace")
        os.remove(os.path.basename(self.tsl_file)+".ace.dir")
    except:
      pass

    for fl in os.listdir():
      if self.name+"_" in fl:
        try:
          os.remove(fl)
        except:
          pass

  def _read_temp(self, itemp):
    # Get contributions to absorption from absorption reaction
    for amt in _ABS_MT_LIST:
      fname = self.name+"_1DXS_" + str(self.ZA) + ".00c_MT" + str(amt) + ".mg"
      if fname in os.listdir():
        xs = read_1dxs(fname)
        self.Ea[itemp,:,:] += xs

    # Check for fission
    if self.fissile:
      fname = self.name+"_NuChi_" + str(self.ZA) + ".00c.mg"
      nu_chi = read_nu_chi(fname)
      self.nu[itemp, :] = nu_chi[0,:]
      self.chi[itemp, :] = nu_chi[3,:]

      fname = self.name+"_1DXS_" + str(self.ZA) + ".00c_MT18.mg"
      xs = read_1dxs(fname)
      self.Ef[itemp,:,:] = xs
      self.Ea[itemp,:,:] += xs

    # Now we get scattering bits
    for d in range(len(self.dilutions)):
      for mt in _2D_MT_LIST:
        # Read xs
        fname = self.name+"_2DXS_" + str(self.ZA) + ".00c_MT" + str(mt) + "_bg" + str(d) + ".mg"

        if fname in os.listdir():
          xs = read_2dxs(fname)
          self.Es[itemp,d,:,:] += xs[0,:,:]
          self.Es1[itemp,d,:,:] += xs[1,:,:]

          v = _YIELDS[mt]
          if v > 1:
            # Subtract from absorption
            for gin in range(self.ngroups):
              self.Ea[itemp,d,gin] -= (v - 1.)*np.sum(xs[0,gin,:])

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

  array = np.array(array, dtype=np.float32)
  array = np.swapaxes(array, 0, 1)

  # First index on dilution, second on group
  return array

def read_nu_chi(fname):
  fl = open(fname, 'r')

  # Skip first lines that have headers / delayed fractions
  for i in range(8):
    fl.readline()

  array = []
  for line in fl:
    line = line.strip()
    if len(line) == 0:
      continue
    line = line.split()[4:]
    for i in range(len(line)):
      line[i] = float(line[i])
    array.append(line)
  fl.close()

  array = np.array(array, dtype=np.float32)
  array = np.swapaxes(array, 0, 1)

  # First index on quantity, second on group
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

  array = np.array(array, dtype=np.float32)

  return array
