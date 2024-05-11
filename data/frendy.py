import ENDFtk
import shutil
import subprocess
import os
import numpy as np

class KRAMXS:
  def __init__(self):
    self.Et = None
    self.Ea = None
    self.Es = None
    self.Es1 = None
    self.Ef = None
    self.nu = None
    self.chi = None

  @property
  def ngroups(self):
    if self.Et is None:
      return 0
    return len(self.Et)

  def __read_line(fl):
    line = fl.readline()
    line = line.strip().split()
    for i in range(len(line)):
      line[i] = float(line[i])
    line = np.array(line)
    return line

  def from_file(fname):
    fl = open(fname, 'r')
    fl.readline() # Skip the XSN 1 header

    # Read scattering matrix first
    Es = []
    ngroups = 100 # This is a guess to start
    line_num = 0
    while line_num < ngroups:
      line_num += 1
      Es.append(KRAMXS.__read_line(fl))
      ngroups = len(Es[-1])
    Es = np.array(Es)
    Es = np.copy(np.swapaxes(Es, 0, 1))

    # Read vEf
    vEf = KRAMXS.__read_line(fl)

    # Read Ea
    Ea = KRAMXS.__read_line(fl)

    # Read Et
    Et = KRAMXS.__read_line(fl)

    # Read Ef
    Ef = KRAMXS.__read_line(fl)

    # Skip FSP 1 line
    fl.readline()

    # Read chi
    chi = KRAMXS.__read_line(fl)

    # Skip ASC 1 line and 1 line
    fl.readline()
    fl.readline()

    # Read P1-scattering matrix
    Es1 = []
    line_num = 0
    while line_num < ngroups:
      line_num += 1
      Es1.append(KRAMXS.__read_line(fl))
    Es1 = np.array(Es1)
    Es1 = np.copy(np.swapaxes(Es1, 0, 1))

    fl.close()

    # Create and return instance
    xs = KRAMXS()
    xs.Et = Et
    xs.Ea = Ea
    xs.Es = Es
    xs.Es1 = Es1
    xs.Ef = Ef
    xs.nu = np.divide(vEf, Ef, out=np.zeros_like(vEf), where=Ef!=0.)
    xs.chi = chi
    return xs

class FrendyMG:
  def __init__(self):
    self.temps = [293.6]
    self.dilutions = None
    self.pot_xs = None
    self.endf_file = None
    self.tsl_file = None
    self.tsl_type = None
    self.label = ""
    self.name = ""
    #self.ngroups = 281
    self.ngroups = 172
    self.initialized = False
    self.processed = False
    self.resonant = False
    self.delete_files = True

  def initialize(self):
    if self.dilutions is not None:
      self.dilutions.sort()
      if len(self.dilutions) == 0 or self.dilutions[-1] < 1.E10:
        self.dilutions.append(1.E10)

    self._get_endf_info()

    self._allocate_arrays()

    self.initialized = True

  def _allocate_arrays(self):
    if self.dilutions is None:
      return

    if len(self.dilutions) > 1:
      self.resonant = True

    self.Es =  np.zeros((len(self.temps), len(self.dilutions), self.ngroups, self.ngroups))
    self.Es1 = np.zeros((len(self.temps), len(self.dilutions), self.ngroups, self.ngroups))
    self.Ea =  np.zeros((len(self.temps), len(self.dilutions), self.ngroups))
    self.flux =  np.zeros((len(self.temps), len(self.dilutions), self.ngroups))
    if self.fissile:
      self.Ef = np.zeros((len(self.temps), len(self.dilutions), self.ngroups))
      self.nu = np.zeros((len(self.temps), self.ngroups))
      self.chi = np.zeros((len(self.temps), self.ngroups))
    else:
      self.Ef = None
      self.nu = None
      self.chi = None

  def process(self, h5=None):
    if not self.initialized:
      self.initialize()

    # Make sure we have all tsl info
    if (self.tsl_file is not None and self.tsl_type is None) or (self.tsl_file is None and self.tsl_type is not None):
      raise RuntimeError("For TSL, must provide both tsl_file and tsl_type.")
    
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
    grp.attrs['dilutions'] = self.dilutions

    # Save cross section data
    grp.create_dataset("absorption", data=self.Ea)
    grp.create_dataset("scatter", data=self.Es)
    grp.create_dataset("p1-scatter", data=self.Es1)
    grp.create_dataset("flux", data=self.flux)
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

    if self.pot_xs is None:
      self.pot_xs = 4. * np.pi * AP * AP

  def _frendy_input(self, temp): 
    out = "mg_neutron_mode\n"
    out += "mg_edit_option ( KRAMXS MGFlux )\n"
    out += "nucl_file_name ({endf})\n".format(endf=self.endf_file)
    if self.tsl_file is not None:
      out += "nucl_file_name_tsl ({tsl})\n".format(tsl=self.tsl_file)
      out += "mg_tsl_data_type {tsl_type}\n".format(tsl_type=self.tsl_type)
    if self.pot_xs is not None:
      out += "potential_scat_xs {pot_xs}\n".format(pot_xs=self.pot_xs)
    out += "mg_file_name {mgfname}\n".format(mgfname=self.name)
    out += "temperature {temp}\n".format(temp=temp)
    out += "legendre_order 1\n"
    #out += "mg_structure ( shem-cea-281 )\n"
    out += "mg_structure ( xmas-nea-lanl-172 )\n"
    out += "mg_weighting_spectrum ( fission+1/e+maxwell  )\n"
    out += "process_gas_xs off\n"
    if self.dilutions is None:
      out += "sigma_zero_data ( auto 0.005 100 1.E-10 factor linear )"
    else:
      dil_frmt = len(self.dilutions)*"{:.2E} "
      out += "sigma_zero_data ( " + dil_frmt.format(*self.dilutions) + " )\n"
    return out

  def _process_temp(self, itemp):
    self._get_endf_info()

    temp = self.temps[itemp]
    frendy_input = self._frendy_input(temp)
    with open("frendy_input", 'w') as fl:
      fl.write(frendy_input)

    subprocess.run(['frendy', 'frendy_input'])

    if itemp == 0 and self.dilutions is None:
      self._get_dilutions()

    self._read_temp(itemp)

    if self.delete_files:
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

  def _get_dilutions(self):
    fname = self.name + "_MGFlux.mg"
    fl = open(fname, 'r')
    fl.readline()
    line = fl.readline()
    fl.close()
    line = line.strip().split()
    line = line[2:]
    for i in range(len(line)):
      line[i] = float(line[i])
    line.sort()
    self.dilutions = line
    self._allocate_arrays()

  def _read_temp(self, itemp):
    # For each dilution, we need to read the KRAMXS file
    for d in range(len(self.dilutions)):
      # Read xs file
      fname = self.name+"_KRAMXS_MACRO_bg" + str(d) + ".mg"
      xs = KRAMXS.from_file(fname)

      # Save values. FRENDY order dilutions from high to low, hence the index
      # shift on d to add them backwards
      self.Ea[itemp,-(d+1),:] = xs.Ea
      self.Es[itemp,-(d+1),:,:] = xs.Es
      self.Es1[itemp,-(d+1),:,:] = xs.Es1
      if self.fissile:
        self.Ef[itemp,-(d+1),:] = xs.Ef
        if d == 0:
          self.nu[itemp,:] = xs.nu
          self.chi[itemp,:] = xs.chi

    # Read flux, and save to file
    fname = self.name + "_MGFlux.mg"
    flux = read_1dxs(fname, 2)
    self.flux[itemp,:,:] = flux

def read_1dxs(fname, nskip):
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
  array = np.copy(np.swapaxes(array, 0, 1))

  # First index on dilution, second on group
  return array