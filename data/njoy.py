import ENDFtk
import shutil
import subprocess
import os

class NjoyFreeGasMG:
  def __init__(self):
    self.temps = []
    self.dilutions = []
    self.endf_file = ""
    self.label = ""
    self.mat = None
    self.name = ""

    self._get_endf_info()

  def _get_endf_info(self):
    # First, get MAT
    tape = ENDFtk.tree.Tape.from_file(self.endf_file)
    self.mat = tape.material_numbers[0]
    mat = tape.MAT(self.mat).parse()

    # Read AWR and ZA from MF1 MT 451
    mf1mt451 = mat.MF(1).MT(451)
    self.awr = mf1mt451.AWR
    self.ZA = mf1mt451.ZA

    # Get potential scattering xs
    mf2mt151 = mat.MF(2).MT(451).isotopes[0]
    self.pot_xs = 0.

  def process(self):
    # Copy ENDF file to tape20
    shutil.copy(self.endf_file, "tape20")

    # Process
    njoy_input = self._reconr()
    njoy_input += self._broadr()
    njoy_input += self._purr()
    njoy_input += self._thermr()
    njoy_input += self._groupr()
    njoy_input += "stop\n"

    with open("njoy_input", 'w') as fl:
      fl.write(njoy_input)

    fl = open("njoy_input")
    subprocess.run(['njoy'], stdin=fl)
    fl.close()

    shutil.move("njoy_input", self.name+"_njoy_input")
    shutil.move("output", self.name+"_njoy_output")
    os.remove("tape20")
    os.remove("tape21")
    os.remove("tape22")
    os.remove("tape23")
    os.remove("tape24")

  def _reconr(self):
    out = "reconr\n"
    out += "20 21 /\n"
    out += "'{}'/\n".format(self.label)
    out += "{mat} 0/\n".format(mat=self.mat)
    out += "0.01 /\n"
    out += "0/\n"
    return out

  def _broadr(self):
    tmp_frmt = len(self.temps)*"{:5f} "
    out = "broadr\n"
    out += "20 21 22 /\n"
    out += "{mat} {ntemp} /\n".format(mat=self.mat, ntemp=len(self.temps))
    out += "0.01 /\n"
    out += tmp_frmt.format(*self.temps) + " /\n"
    out += "0 /\n"
    return out

  def _purr(self):
    tmp_frmt = len(self.temps)*"{:5f} "
    dil_frmt = len(self.dilutions)*"{:.2E} "
    out = "purr\n"
    out += "20 22 23 /\n"
    out += "{mat} {ntemp} {ndil} 20 32 /\n".format(mat=self.mat, ntemp=len(self.temps), ndil=len(self.dilutions))
    out += tmp_frmt.format(*self.temps) + " /\n"
    out += dil_frmt.format(*self.dilutions) + " /\n"
    out += "0 /\n"
    return out

  def _thermr(self):
    tmp_frmt = len(self.temps)*"{:5f} "
    out = "thermr\n"
    out += "0 23 24 /\n"
    out += "{mat1} {mat2} 20 {ntemp} 1 0 0 1 221 /\n".format(mat1=self.mat, mat2=self.mat, ntemp=len(self.temps))
    out += tmp_frmt.format(*self.temps) + " /\n"
    out += "0.01 5. /\n"
    return out

  def _groupr(self):
    tmp_frmt = len(self.temps)*"{:5f} "
    dil_frmt = len(self.dilutions)*"{:.2E} "
    out = "groupr\n"
    out += "20 24 /\n"
    out += "{mat} 24 0 4 1 {ntemp} {ndil} /\n".format(mat=self.mat, ntemp=len(self.temps), ndil=len(self.dilutions))
    out += "'{label}' /\n".format(label=self.label)
    out += tmp_frmt.format(*self.temps) + " /\n"
    out += dil_frmt.format(*self.dilutions) + " /\n"
    out += ".10 .025 820.3e3 1.4e6 /\n" # card 8c
    out += "3 /\n"     # All vector XSs
    out += "3 221 /\n" # Thermal vector XS
    out += "3 229 /\n" # Inverse velocity
    out += "3 455 /\n" # Delayed nubar
    out += "5 455 /\n" # Delayed neutron spectrum
    out += "6 /\n"     # All matrix XSs
    out += "6 221 /\n" # Thermal matrix XS
    out += "0 /\n"     # Stop
    return out