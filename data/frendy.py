import ENDFtk
import subprocess
import os
import numpy as np
from scarabee import *

# Estimates of Intermediate-Resonance Lambda values as a function of mass
# number, A, taken from "WIMS-D Library Update", by the IAEA.
_IR_LAMBDA = [1.00, 1.00, 0.99, 0.99, 0.99, 0.99, 0.98, 0.98, 0.98, 0.98, # 010
              0.98, 0.97, 0.97, 0.97, 0.96, 0.96, 0.95, 0.94, 0.93, 0.91, # 020
              0.90, 0.89, 0.87, 0.85, 0.84, 0.82, 0.81, 0.79, 0.78, 0.76, # 030
              0.75, 0.74, 0.72, 0.71, 0.70, 0.69, 0.68, 0.67, 0.66, 0.65, # 040
              0.64, 0.63, 0.62, 0.61, 0.60, 0.59, 0.59, 0.58, 0.57, 0.56, # 050
              0.56, 0.55, 0.54, 0.54, 0.53, 0.52, 0.52, 0.51, 0.51, 0.50, # 060
              0.50, 0.49, 0.48, 0.48, 0.47, 0.47, 0.46, 0.46, 0.46, 0.45, # 070
              0.45, 0.44, 0.44, 0.43, 0.43, 0.42, 0.42, 0.42, 0.41, 0.41, # 080
              0.41, 0.40, 0.40, 0.40, 0.39, 0.39, 0.38, 0.38, 0.38, 0.38, # 090
              0.37, 0.37, 0.37, 0.36, 0.36, 0.36, 0.35, 0.35, 0.35, 0.35, # 100
              0.34, 0.34, 0.34, 0.33, 0.33, 0.33, 0.33, 0.32, 0.32, 0.32, # 110
              0.32, 0.31, 0.31, 0.31, 0.31, 0.31, 0.30, 0.30, 0.30, 0.30, # 120
              0.29, 0.29, 0.29, 0.29, 0.29, 0.28, 0.28, 0.28, 0.28, 0.28, # 130
              0.28, 0.27, 0.27, 0.27, 0.27, 0.27, 0.27, 0.26, 0.26, 0.26, # 140
              0.26, 0.26, 0.26, 0.26, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, # 150
              0.25, 0.25, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, # 160
              0.24, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, # 170
              0.23, 0.23, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, # 180
              0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.21, 0.21, # 190
              0.21, 0.21, 0.21, 0.21, 0.21, 0.21, 0.21, 0.21, 0.21, 0.21, # 200
              0.21, 0.21, 0.21, 0.21, 0.21, 0.21, 0.21, 0.21, 0.21, 0.21, # 210
              0.21, 0.21, 0.21, 0.21, 0.20]                               # 215

def ir_lambda(A):
  if A >= 1 and A <= len(_IR_LAMBDA):
    return _IR_LAMBDA[A-1]
  elif A < 1:
    return _IR_LAMBDA[0]
  else:
    return _IR_LAMBDA[-1]


_GROUPS_STRUCTURES = ["WIMS-69", "APOLLO-99", "XMAS-172", "SHEM-281"]

_GROUP_IDS = {"WIMS-69": "epri-69", "XMAS-172": "xmas-nea-lanl-172", "SHEM-281": "shem-cea-281"}

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
	            2.00000E-08, 1.50000E-08, 1.00000E-08, 5.00000E-09, 1.00000E-11]) * 1.E6,

  "APOLLO-99":
    np.array([1.00000000E+07, 8.18730800E+06, 6.70320100E+06, 5.48811700E+06, 4.49329000E+06,
              3.67879400E+06, 3.01194300E+06, 2.46597100E+06, 2.01896600E+06, 1.65299000E+06,
              1.35335300E+06, 1.10803200E+06, 9.07179875E+05, 6.08101125E+05, 4.07622094E+05,
              2.73237406E+05, 1.83156406E+05, 1.22773500E+05, 8.22975312E+04, 5.51656602E+04,
              3.69786602E+04, 2.47875195E+04, 1.66155801E+04, 1.11377598E+04, 7.46585889E+03,
              5.00451709E+03, 3.35462598E+03, 2.24867407E+03, 1.50733203E+03, 1.01039398E+03,
              6.77287720E+02, 4.53999298E+02, 3.04325012E+02, 2.03995193E+02, 1.36742004E+02,
              9.16609268E+01, 6.79040833E+01, 5.55951691E+01, 4.55174789E+01, 3.72665291E+01,
              3.05112591E+01, 2.49805107E+01, 2.04523201E+01, 1.67449493E+01, 1.37095900E+01,
              1.12244701E+01, 9.18981743E+00, 7.52398682E+00, 6.16012096E+00, 5.04347706E+00,
              4.12925005E+00, 3.38074493E+00, 2.76792002E+00, 2.35999990E+00, 2.13000011E+00,
              2.01999998E+00, 1.92999995E+00, 1.84000003E+00, 1.75500000E+00, 1.66999996E+00,
              1.59000003E+00, 1.50999999E+00, 1.44000006E+00, 1.37000000E+00, 1.30499995E+00,
              1.23500001E+00, 1.16999996E+00, 1.11000001E+00, 1.07000005E+00, 1.03499997E+00,
              9.86000121E-01, 9.30000007E-01, 8.60000014E-01, 7.90000021E-01, 7.04999983E-01,
              6.25000000E-01, 5.40000021E-01, 4.85000014E-01, 4.32999998E-01, 3.91000003E-01,
              3.51999998E-01, 3.14500004E-01, 2.82499999E-01, 2.47999996E-01, 2.19999999E-01,
              1.88999996E-01, 1.59999996E-01, 1.34000003E-01, 1.15000002E-01, 9.49999988E-02,
              7.69999996E-02, 5.90000004E-02, 4.30000015E-02, 2.99999993E-02, 1.99999996E-02,
              1.49999997E-02, 9.99999978E-03, 5.49999997E-03, 3.00000003E-03, 1.10000001E-04]),

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
	            5.00000E-09, 3.00000E-09, 1.00001E-11]) * 1.E6,

  "SHEM-281":
    np.array([1.964030E+07, 1.491823E+07, 1.384029E+07, 1.161833E+07,
              9.999987E+06, 9.048363E+06, 8.187297E+06, 7.408173E+06,
              6.703192E+06, 6.065299E+06, 4.965847E+06, 4.065691E+06,
              3.328707E+06, 2.725314E+06, 2.231299E+06, 1.901387E+06,
              1.636539E+06, 1.405768E+06, 1.336941E+06, 1.286961E+06,
              1.162048E+06, 1.051149E+06, 9.511189E+05, 8.600058E+05,
              7.065112E+05, 5.784425E+05, 4.940018E+05, 4.560211E+05,
              4.125012E+05, 3.838835E+05, 3.206464E+05, 2.678264E+05,
              2.300137E+05, 1.950077E+05, 1.649989E+05, 1.399995E+05,
              1.227732E+05, 1.156235E+05, 9.466450E+04, 8.229736E+04,
              6.737938E+04, 5.516557E+04, 4.991587E+04, 4.086766E+04,
              3.697859E+04, 3.345961E+04, 2.928101E+04, 2.739441E+04,
              2.610010E+04, 2.499908E+04, 2.269941E+04, 1.858471E+04,
              1.620045E+04, 1.489967E+04, 1.360366E+04, 1.113774E+04,
              9.118808E+03, 7.465848E+03, 6.112520E+03, 5.004508E+03,
              4.097345E+03, 3.481068E+03, 2.996183E+03, 2.578838E+03,
              2.219627E+03, 1.910451E+03, 1.614038E+03, 1.345061E+03,
              1.135007E+03, 1.064962E+03, 9.075007E+02, 7.485173E+02,
              6.128342E+02, 5.017462E+02, 4.107950E+02, 3.535746E+02,
              3.199275E+02, 2.837502E+02, 2.417960E+02, 1.979658E+02,
              1.620807E+02, 1.327005E+02, 1.086459E+02, 8.895177E+01,
              7.504548E+01, 6.144204E+01, 5.267255E+01, 4.579131E+01,
              4.399581E+01, 4.016895E+01, 3.372011E+01, 2.760769E+01,
              2.460856E+01, 2.253556E+01, 2.237836E+01, 2.215569E+01,
              2.200114E+01, 2.170178E+01, 2.148585E+01, 2.133597E+01,
              2.122956E+01, 2.114481E+01, 2.106040E+01, 2.097632E+01,
              2.076761E+01, 2.068470E+01, 2.060213E+01, 2.051988E+01,
              2.041754E+01, 2.027512E+01, 2.007338E+01, 1.959735E+01,
              1.939265E+01, 1.919969E+01, 1.908484E+01, 1.795905E+01,
              1.775903E+01, 1.756476E+01, 1.744572E+01, 1.683053E+01,
              1.655014E+01, 1.604977E+01, 1.577923E+01, 1.486626E+01,
              1.473012E+01, 1.459522E+01, 1.447024E+01, 1.425053E+01,
              1.404961E+01, 1.354604E+01, 1.332970E+01, 1.259997E+01,
              1.247210E+01, 1.230855E+01, 1.213015E+01, 1.197947E+01,
              1.181529E+01, 1.170943E+01, 1.158944E+01, 1.126944E+01,
              1.105292E+01, 1.080376E+01, 1.057925E+01, 9.500024E+00,
              9.140311E+00, 8.979950E+00, 8.800375E+00, 8.673690E+00,
              8.524074E+00, 8.300322E+00, 8.130272E+00, 7.970079E+00,
              7.839651E+00, 7.739943E+00, 7.600350E+00, 7.380153E+00,
              7.139869E+00, 6.994292E+00, 6.917776E+00, 6.870208E+00,
              6.835259E+00, 6.810696E+00, 6.791653E+00, 6.776050E+00,
              6.759807E+00, 6.742254E+00, 6.716683E+00, 6.631257E+00,
              6.606106E+00, 6.588293E+00, 6.571843E+00, 6.556090E+00,
              6.539066E+00, 6.514916E+00, 6.481775E+00, 6.432057E+00,
              6.359784E+00, 6.280153E+00, 6.160108E+00, 6.059906E+00,
              5.960142E+00, 5.800211E+00, 5.720146E+00, 5.619790E+00,
              5.530036E+00, 5.488167E+00, 5.410245E+00, 5.380032E+00,
              5.320112E+00, 5.210076E+00, 5.109974E+00, 4.933232E+00,
              4.767845E+00, 4.419800E+00, 4.309812E+00, 4.219828E+00,
              4.000000E+00, 3.882170E+00, 3.712087E+00, 3.543073E+00,
              3.142109E+00, 2.884047E+00, 2.775121E+00, 2.740922E+00,
              2.719898E+00, 2.700115E+00, 2.640041E+00, 2.620053E+00,
              2.590094E+00, 2.550003E+00, 2.469941E+00, 2.330061E+00,
              2.272986E+00, 2.217087E+00, 2.156948E+00, 2.070095E+00,
              1.989920E+00, 1.900077E+00, 1.779966E+00, 1.668949E+00,
              1.588030E+00, 1.519976E+00, 1.443967E+00, 1.410007E+00,
              1.380981E+00, 1.330952E+00, 1.293038E+00, 1.250939E+00,
              1.213968E+00, 1.169989E+00, 1.147969E+00, 1.129974E+00,
              1.116049E+00, 1.103950E+00, 1.091982E+00, 1.077986E+00,
              1.034993E+00, 1.021012E+00, 1.009035E+00, 9.965005E-01,
              9.819591E-01, 9.639598E-01, 9.440222E-01, 9.199779E-01,
              8.800244E-01, 8.200371E-01, 7.199989E-01, 6.249987E-01,
              5.949930E-01, 5.549897E-01, 5.200108E-01, 4.750165E-01,
              4.315786E-01, 3.900011E-01, 3.529935E-01, 3.250079E-01,
              3.050115E-01, 2.799888E-01, 2.549965E-01, 2.311923E-01,
              2.096102E-01, 1.900049E-01, 1.618953E-01, 1.379994E-01,
              1.199949E-01, 1.042977E-01, 8.979683E-02, 7.649686E-02,
              6.519936E-02, 5.549815E-02, 4.730186E-02, 4.029993E-02,
              3.439976E-02, 2.929889E-02, 2.493942E-02, 2.001035E-02,
              1.482996E-02, 1.045050E-02, 7.145263E-03, 4.556021E-03,
              2.499897E-03, 1.100027E-04])
}

_DEFAULT_GROUP_STRUCTURE = "XMAS-172"

_DEFAULT_MAX_LEGENDRE_MOMENT = 1

def set_default_group_structure(name):
  global _DEFAULT_GROUP_STRUCTURE
  if name not in _GROUPS_STRUCTURES:
    raise RuntimeError("Uknown group structure \"{}\".".format(name))
  _DEFAULT_GROUP_STRUCTURE = name

def set_default_max_legendre_moments(l):
  global _DEFAULT_MAX_LEGENDRE_MOMENT
  if l >= 0 and l <= 3:
    _DEFAULT_MAX_LEGENDRE_MOMENT = l
  else:
    raise RuntimeError("Default max legendre moment must be in range [0, 3].")

def get_default_max_legendre_moments():
  return _DEFAULT_MAX_LEGENDRE_MOMENT

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
    self.Es2 = None
    self.Es3 = None
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

  def from_file(fname, max_l):
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
    Es1 = None
    if max_l >= 1:
      Es1 = []
      line_num = 0
      while line_num < ngroups:
        line_num += 1
        Es1.append(KRAMXS.__read_line(fl))
      Es1 = np.array(Es1)
      Es1 = np.copy(np.swapaxes(Es1, 0, 1))

    # Read P2-scattering matrix
    Es2 = None
    if max_l >= 2:
      Es2 = []
      line_num = 0
      while line_num < ngroups:
        line_num += 1
        Es2.append(KRAMXS.__read_line(fl))
      Es2 = np.array(Es2)
      Es2 = np.copy(np.swapaxes(Es2, 0, 1))
    
    # Read P3-scattering matrix
    Es3 = None
    if max_l >= 3:
      Es3 = []
      line_num = 0
      while line_num < ngroups:
        line_num += 1
        Es3.append(KRAMXS.__read_line(fl))
      Es3 = np.array(Es3)
      Es3 = np.copy(np.swapaxes(Es3, 0, 1))

    fl.close()

    # Create and return instance
    xs = KRAMXS()
    xs.Et = Et
    xs.Ea = Ea
    xs.Es = Es

    if Es1 is not None:
      xs.Es1 = Es1
    if Es2 is not None:
      xs.Es2 = Es2
    if Es3 is not None:
      xs.Es3 = Es3

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
    self.max_legendre_moment = _DEFAULT_MAX_LEGENDRE_MOMENT

    if self.max_legendre_moment > 3:
      raise RuntimeError("Only legendre moments up to L=3 are supported.")

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
    
    self.Dtr =  np.zeros((len(self.temps), len(self.dilutions), self.ngroups))
    self.Ea =  np.zeros((len(self.temps), len(self.dilutions), self.ngroups))
    self.Es =  np.zeros((len(self.temps), len(self.dilutions), self.ngroups, self.ngroups))
    if self.fissile:
      self.Ef = np.zeros((len(self.temps), len(self.dilutions), self.ngroups))

      # Nu and chi are only very weakly dependent on temp and dilution.
      # Because of this, we don't tabulate them on temp or dilution.
      self.nu = np.zeros((self.ngroups))
      self.chi = np.zeros((self.ngroups))
    else:
      self.Ef = None
      self.nu = None
      self.chi = None

    if self.max_legendre_moment >= 1:
      self.Es1 = np.zeros((len(self.temps), len(self.dilutions), self.ngroups, self.ngroups))
    else:
      self.Es1 = None

    if self.max_legendre_moment >= 2:
      self.Es2 = np.zeros((len(self.temps), len(self.dilutions), self.ngroups, self.ngroups))
    else:
      self.Es2 = None

    if self.max_legendre_moment >= 3:
      self.Es3 = np.zeros((len(self.temps), len(self.dilutions), self.ngroups, self.ngroups))
    else:
      self.Es3 = None

    # Depletion related reactions
    self.Egamma =  np.zeros((len(self.temps), len(self.dilutions), self.ngroups))
    self.En2n = None
    self.En3n = None
    self.Enp  = None
    self.Ena  = None

  def process(self, h5=None, chi=None):
    if not self.initialized:
      self.initialize()

    # Make sure we have all tsl info
    if (self.tsl_file is not None and self.tsl_type is None) or (self.tsl_file is None and self.tsl_type is not None):
      raise RuntimeError("For TSL, must provide both tsl_file and tsl_type.")
    
    for i in range(len(self.temps)):
      self._process_temp(i)
    
    self.processed = True

    # Truncate threshold reactions to remove zeros
    self._remove_zeros() 

    if chi is not None:
      self.apply_P1_transport_correction(chi)
    else:
      self.apply_inscatter_transport_correction()

    # Apply compression after computing the transport correction ! 
    self._get_compressed_scatter_layout()
    self._compress_scatter_matrices()

    if h5 is not None:
      self.add_to_hdf5(h5)

  def apply_P1_transport_correction(self, chi):
    if not self.processed:
      raise RuntimeError("Cannot apply transport corretion to unprocessed data.")

    for iT in range(len(self.temps)):
      for id in range(len(self.dilutions)):
        # Create a temporary xs set with the provided fission spectrum
        Et = self.Ea[iT, id, :] + np.sum(self.Es[iT, id, :, :], axis=1)
        Es = np.array([self.Es[iT, id, :, :], self.Es1[iT, id, :, :]])
        if self.fissile: 
          TempXS = CrossSection(Et, self.Dtr[iT, id, :], self.Ea[iT, id, :], Es, self.Ef[iT, id, :], self.nu[iT, id, :]*self.Ef[iT, id, :], chi)
        else:
          TempXS = CrossSection(Et, self.Dtr[iT, id, :], self.Ea[iT, id, :], Es, np.zeros(self.ngroups), np.zeros(self.ngroups), chi)

        # We now perform a P1 leakage calculation
        P1_spectrum = P1CriticalitySpectrum(TempXS, 0.0001)

        # We now have diffusion coefficients
        D = P1_spectrum.diff_coeff

        # Compute transport xs
        Etr = 1. / (3. * D)

        # Calculate the delta xs for the transport correction
        self.Dtr[iT, id, :] = Et - Etr
  
  def apply_inscatter_transport_correction(self):
    if not self.processed:
      raise RuntimeError("Cannot apply transport corretion to unprocessed data.")

    if self.Es1 is not None:
      for iT in range(len(self.temps)):
        for id in range(len(self.dilutions)):
          # Calculate the delta xs for the transport correction
          for g in range(self.ngroups):
              self.Dtr[iT, id, g] = np.sum(self.Es1[iT, id, g, :])

  def _get_compressed_scatter_layout(self):
    self.low_grps = []
    self.high_grps = []
    self.data_starts = []

    for g in range(self.ngroups):
      # Find the first outgoing group which isn't 0
      g_low = 0
      for gg in range(self.ngroups):
        not_all_zeros = np.any(self.Es[:,:,g,:gg+1])
        if not_all_zeros:
          g_low = gg
          break

      # Find the last outgoing group which isn't 0
      g_hi = 0
      for gg in range(self.ngroups):
        all_zeros = not np.any(self.Es[:,:,g,gg:])
        if all_zeros:
          g_hi = gg-1
          break
      if g_hi == 0:
        g_hi = self.ngroups - 1
      
      self.low_grps.append(g_low)
      self.high_grps.append(g_hi)

      if g == 0:
        self.data_starts.append(0)
      else:
        self.data_starts.append(self.data_starts[-1] + (self.high_grps[-2] - self.low_grps[-2]) + 1)

    self.len_scatter_matrix_data = self.data_starts[-1] + (self.high_grps[-1] - self.low_grps[-1]) + 1

  def _compress_scatter_matrices(self):
    Es = np.zeros((len(self.temps), len(self.dilutions), self.len_scatter_matrix_data))
    if self.Es1 is not None:
      Es1 = np.zeros((len(self.temps), len(self.dilutions), self.len_scatter_matrix_data))
    if self.Es2 is not None:
      Es2 = np.zeros((len(self.temps), len(self.dilutions), self.len_scatter_matrix_data))
    if self.Es3 is not None:
      Es3 = np.zeros((len(self.temps), len(self.dilutions), self.len_scatter_matrix_data))

    i = 0
    for g in range(self.ngroups):
      g_low = self.low_grps[g]
      g_hi = self.high_grps[g]
      l = g_hi - g_low + 1

      Es[:,:,i:i+l] = self.Es[:,:,g, g_low:g_hi+1]

      if self.Es1 is not None:
        Es1[:,:,i:i+l] = self.Es1[:,:,g, g_low:g_hi+1] 
      if self.Es2 is not None:
        Es2[:,:,i:i+l] = self.Es2[:,:,g, g_low:g_hi+1] 
      if self.Es3 is not None:
        Es3[:,:,i:i+l] = self.Es3[:,:,g, g_low:g_hi+1] 

      i += l
    
    self.Es = Es
    if self.Es1 is not None:
      self.Es1 = Es1
    if self.Es2 is not None:
      self.Es2 = Es2
    if self.Es3 is not None:
      self.Es3 = Es3

  def _remove_zeros(self):
    self.En2n = self._cull_array(self.En2n) 
    self.En3n = self._cull_array(self.En3n)
    self.Enp  = self._cull_array(self.Enp)
    self.Ena  = self._cull_array(self.Ena)
  
  def _cull_array(self, a):
    if a is not None:
      gmax = self.ngroups-1
      for g in range(self.ngroups-1, -1, -1):
        # Check if group is all zeros
        not_all_zeros = np.any(a[:,:,g])
        if not_all_zeros:
          gmax = g
          break
      if gmax != self.ngroups-1:
        return a[:,:,:gmax+1]
      else:
        return a

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
    grp.attrs['ir-lambda'] = self.ir_lambda
    grp.attrs['temperatures'] = self.temps
    grp.attrs['dilutions'] = self.dilutions

    packing = np.zeros((self.ngroups, 3), dtype=np.uint32)
    for g in range(self.ngroups):
       packing[g, 0] = self.data_starts[g]
       packing[g, 1] = self.low_grps[g]
       packing[g, 2] = self.high_grps[g]
    grp.create_dataset('matrix-compression', data=packing)

    # Save cross section data
    grp.create_dataset("transport-correction", data=self.Dtr)
    grp.create_dataset("absorption", data=self.Ea)
    grp.create_dataset("scatter", data=self.Es)
    if self.Es1 is not None:
      grp.create_dataset("p1-scatter", data=self.Es1)
    if self.Es2 is not None:
      grp.create_dataset("p2-scatter", data=self.Es2)
    if self.Es3 is not None:
      grp.create_dataset("p3-scatter", data=self.Es3)
    if self.fissile:
      grp.create_dataset("fission", data=self.Ef)
      grp.create_dataset("nu", data=self.nu)
      grp.create_dataset("chi", data=self.chi)

    # Depletion data 
    grp.create_dataset("(n,gamma)", data=self.Egamma)
    if self.En2n is not None:
      grp.create_dataset("(n,2n)", data=self.En2n)
    if self.En3n is not None:
      grp.create_dataset("(n,3n)", data=self.En3n)
    if self.Enp is not None:
      grp.create_dataset("(n,p)", data=self.Enp)
    if self.Ena is not None:
      grp.create_dataset("(n,a)", data=self.Ena)

  def _get_endf_info(self):
    # First, get MAT
    tape = ENDFtk.tree.Tape.from_file(self.endf_file)
    self.mat = tape.material_numbers[0]

    # Read AWR and ZA from MF1 MT 451
    mf1mt451 = tape.MAT(self.mat).MF(1).MT(451).parse()
    self.awr = mf1mt451.AWR
    self.ZA = mf1mt451.ZA
    self.ir_lambda = ir_lambda(self.ZA % 1000)
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
    # Will write (n,2n), (n,3n), (n,gamma), (n,p), and (n,alpha)
    out = "mg_neutron_mode\n"
    out += "mg_edit_option ( KRAMXS MGFlux 1DXS 16, 17, 102, 103, 107 )\n"
    out += "nucl_file_name ({endf})\n".format(endf=self.endf_file)
    if self.tsl_file is not None:
      out += "nucl_file_name_tsl ({tsl})\n".format(tsl=self.tsl_file)
      out += "mg_tsl_data_type {tsl_type}\n".format(tsl_type=self.tsl_type)
    if self.pot_xs is not None:
      out += "potential_scat_xs {pot_xs}\n".format(pot_xs=self.pot_xs)
    out += "mg_file_name {mgfname}\n".format(mgfname=self.name)
    out += "temperature {temp}\n".format(temp=temp)
    out += "legendre_order {max_l}\n".format(max_l=self.max_legendre_moment)
    if self.group_strucutre in _GROUP_IDS:
      out += "mg_structure ( {id} )\n".format(id=_GROUP_IDS[self.group_strucutre])
    else:
      ebnd_frmt = len(_GROUP_BOUNDS[self.group_strucutre])*"{:.7E}  "
      out += "mg_structure ( " + ebnd_frmt.format(*np.flip(_GROUP_BOUNDS[self.group_strucutre])) + " )\n"
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
      xs = KRAMXS.from_file(fname, self.max_legendre_moment)

      # Save values. FRENDY order dilutions from high to low, hence the index
      # shift on d to add them backwards
      self.Ea[itemp,-(d+1),:] = xs.Ea
      self.Es[itemp,-(d+1),:,:] = xs.Es
      if xs.Es1 is not None:
        self.Es1[itemp,-(d+1),:,:] = xs.Es1
      if xs.Es2 is not None:
        self.Es2[itemp,-(d+1),:,:] = xs.Es2
      if xs.Es3 is not None:
        self.Es3[itemp,-(d+1),:,:] = xs.Es3
      if self.fissile:
        self.Ef[itemp,-(d+1),:] = xs.Ef
        if itemp == 0 and d == 0:
          self.nu = xs.nu
          self.chi = xs.chi
      
    # read in (n,gamma) data
    fname = self.name + "_1DXS_" + str(self.ZA) + ".00c_MT102.mg"
    ngamma = read_1dxs(fname, 3)
    self.Egamma[itemp,:,:] = ngamma

    fls = os.listdir()

    # Check for (n,2n) data
    fname = self.name + "_1DXS_" + str(self.ZA) + ".00c_MT16.mg"
    if fname in fls:
      if itemp == 0:
        self.En2n = np.zeros((len(self.temps), len(self.dilutions), self.ngroups))
      n2n = read_1dxs(fname, 3)
      self.En2n[itemp,:,:] = n2n

    # Check for (n,3n) data
    fname = self.name + "_1DXS_" + str(self.ZA) + ".00c_MT17.mg"
    if fname in fls:
      if itemp == 0:
        self.En3n = np.zeros((len(self.temps), len(self.dilutions), self.ngroups))
      n3n = read_1dxs(fname, 3)
      self.En3n[itemp,:,:] = n3n

    # Check for (n,p) data
    fname = self.name + "_1DXS_" + str(self.ZA) + ".00c_MT103.mg"
    if fname in fls:
      if itemp == 0:
        self.Enp = np.zeros((len(self.temps), len(self.dilutions), self.ngroups))
      nprt = read_1dxs(fname, 3)
      self.Enp[itemp,:,:] = nprt

    # Check for (n,a) data
    fname = self.name + "_1DXS_" + str(self.ZA) + ".00c_MT107.mg"
    if fname in fls:
      if itemp == 0:
        self.Ena = np.zeros((len(self.temps), len(self.dilutions), self.ngroups))
      na = read_1dxs(fname, 3)
      self.Ena[itemp,:,:] = na


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
