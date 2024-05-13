import ENDFtk
import shutil
import subprocess
import os
import numpy as np

_GROUPS_STRUCTURES = ["WIMS-69", "XMAS-172"]

_GROUP_BOUNDS = {
  "WIMS-69":
    np.array([1.00000E+01, 6.06550E+00, 3.67900E+00, 2.23100E+00, 1.35300E+00,
	            8.21000E-01, 5.00000E-01, 3.02500E-01, 1.83000E-01, 1.11000E-01,
	            6.73400E-02, 4.08500E-02, 2.47800E-02, 1.50300E-02, 9.11800E-03,
	            5.53000E-03, 3.51910E-03, 2.23945E-03, 1.42510E-03, 9.06899E-04,
	            3.67263E-04, 1.48729E-04, 7.55014E-05, 4.80520E-05, 2.77000E-05,
	            1.59680E-05, 9.87700E-06, 4.00000E-06, 3.30000E-06, 2.60000E-06,
	            2.10000E-06, 1.50000E-06, 1.30000E-06, 1.15000E-06, 1.12300E-06,
	            1.09700E-06, 1.07100E-06, 1.04500E-06, 1.02000E-06, 9.96000E-07,
	            9.72000E-07, 9.50000E-07, 9.10000E-07, 8.50000E-07, 7.80000E-07,
	            6.25000E-07, 5.00000E-07, 4.00000E-07, 3.50000E-07, 3.20000E-07,
	            3.00000E-07, 2.80000E-07, 2.50000E-07, 2.20000E-07, 1.80000E-07,
	            1.40000E-07, 1.00000E-07, 8.00000E-08, 6.70000E-08, 5.80000E-08,
	            5.00000E-08, 4.20000E-08, 3.50000E-08, 3.00000E-08, 2.50000E-08,
	            2.00000E-08, 1.50000E-08, 1.00000E-08, 5.00000E-09, 1.00000E-11], dtype=np.float32) * 1.E6,

  "XMAS-172": 
    np.array([1.96403E+01, 1.73325E+01, 1.49182E+01, 1.38403E+01, 1.16183E+01,
	            1.00000E+01, 8.18731E+00, 6.70320E+00, 6.06531E+00, 5.48812E+00,
	            4.49329E+00, 3.67879E+00, 3.01194E+00, 2.46597E+00, 2.23130E+00,
	            2.01897E+00, 1.65299E+00, 1.35335E+00, 1.22456E+00, 1.10803E+00,
	            1.00259E+00, 9.07180E-01, 8.20850E-01, 6.08101E-01, 5.50232E-01,
	            4.97871E-01, 4.50492E-01, 4.07622E-01, 3.01974E-01, 2.73237E-01,
	            2.47235E-01, 1.83156E-01, 1.22773E-01, 1.11090E-01, 8.22975E-02,
	            6.73795E-02, 5.51656E-02, 4.08677E-02, 3.69786E-02, 2.92830E-02,
	            2.73944E-02, 2.47875E-02, 1.66156E-02, 1.50344E-02, 1.11378E-02,
	            9.11882E-03, 7.46586E-03, 5.53084E-03, 5.00451E-03, 3.52662E-03,
	            3.35463E-03, 2.24867E-03, 2.03468E-03, 1.50733E-03, 1.43382E-03,
	            1.23410E-03, 1.01039E-03, 9.14242E-04, 7.48518E-04, 6.77287E-04,
	            4.53999E-04, 3.71703E-04, 3.04325E-04, 2.03995E-04, 1.48625E-04,
	            1.36742E-04, 9.16609E-05, 7.56736E-05, 6.79041E-05, 5.55951E-05,
	            5.15780E-05, 4.82516E-05, 4.55174E-05, 4.01690E-05, 3.72665E-05,
	            3.37201E-05, 3.05113E-05, 2.76077E-05, 2.49805E-05, 2.26033E-05,
	            1.94548E-05, 1.59283E-05, 1.37096E-05, 1.12245E-05, 9.90555E-06,
	            9.18981E-06, 8.31529E-06, 7.52398E-06, 6.16012E-06, 5.34643E-06,
	            5.04348E-06, 4.12925E-06, 4.00000E-06, 3.38075E-06, 3.30000E-06,
	            2.76792E-06, 2.72000E-06, 2.60000E-06, 2.55000E-06, 2.36000E-06,
	            2.13000E-06, 2.10000E-06, 2.02000E-06, 1.93000E-06, 1.84000E-06,
	            1.75500E-06, 1.67000E-06, 1.59000E-06, 1.50000E-06, 1.47500E-06,
	            1.44498E-06, 1.37000E-06, 1.33750E-06, 1.30000E-06, 1.23500E-06,
	            1.17000E-06, 1.15000E-06, 1.12535E-06, 1.11000E-06, 1.09700E-06,
	            1.07100E-06, 1.04500E-06, 1.03500E-06, 1.02000E-06, 9.96000E-07,
	            9.86000E-07, 9.72000E-07, 9.50000E-07, 9.30000E-07, 9.10000E-07,
	            8.60000E-07, 8.50000E-07, 7.90000E-07, 7.80000E-07, 7.05000E-07,
	            6.25000E-07, 5.40000E-07, 5.00000E-07, 4.85000E-07, 4.33000E-07,
	            4.00000E-07, 3.91000E-07, 3.50000E-07, 3.20000E-07, 3.14500E-07,
	            3.00000E-07, 2.80000E-07, 2.48000E-07, 2.20000E-07, 1.89000E-07,
	            1.80000E-07, 1.60000E-07, 1.40000E-07, 1.34000E-07, 1.15000E-07,
	            1.00001E-07, 9.50000E-08, 8.00000E-08, 7.70000E-08, 6.70000E-08,
	            5.80000E-08, 5.00000E-08, 4.20000E-08, 3.50000E-08, 3.00000E-08,
	            2.50000E-08, 2.00000E-08, 1.50000E-08, 1.00000E-08, 6.90000E-09,
	            5.00000E-09, 3.00000E-09, 1.00001E-11], dtype=np.float32) * 1.E6
}

_GROUP_IDS = {"WIMS-69": 9, "XMAS-172": 18}

_DEFAULT_GROUP_STRUCTURE = "XMAS-172"

def set_default_group_structure(name):
  if name not in _GROUPS_STRUCTURES:
    raise RuntimeError("Uknown group structure \"{}\".".format(name))
  _DEFAULT_GROUP_STRUCTURE = name

def get_default_group_structure():
  return _DEFAULT_GROUP_STRUCTURE

def get_default_group_bounds():
  return _GROUP_BOUNDS[_DEFAULT_GROUP_STRUCTURE]

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
  def __init__(self, group_strucutre=None):
    self.temps = [293.6]
    self.dilutions = None
    self.pot_xs = None
    self.endf_file = None
    self.tsl_file = None
    self.tsl_type = None
    self.label = ""
    self.name = ""
    self.group_strucutre = _DEFAULT_GROUP_STRUCTURE
    if group_strucutre is not None:
      if group_strucutre not in _GROUPS_STRUCTURES:
        raise RuntimeError("Unknown group structure \"{}\".".format(group_strucutre))
      self.group_strucutre = group_strucutre
    self.ngroups = len(_GROUP_BOUNDS[self.group_strucutre])-1
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
    out += "mg_structure ( {id} )\n".format(id=_GROUP_IDS[self.group_strucutre])
    out += "mg_weighting_spectrum ( fission+1/e+maxwell  )\n"
    out += "process_gas_xs off\n"
    if self.dilutions is None:
      out += "sigma_zero_data ( auto 0.005 100 1.E-10 rr linear )"
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