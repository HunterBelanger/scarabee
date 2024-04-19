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

  def process(self):
    # First, get MAT
    if self.mat is None:
      tape = ENDFtk.tree.Tape.from_file(self.endf_file)
      self.mat = tape.material_numbers[0]

    # Copy to tape 20
    shutil.copy(self.endf_file, "tape20")

    # Process
    njoy_input = self._reconr()
    njoy_input += self._broadr()
    njoy_input += self._purr()
    njoy_input += self._thermr()
    njoy_input += "stop\n"

    with open("njoy_input", 'w') as fl:
      fl.write(njoy_input)

    fl = open("njoy_input")
    subprocess.run(['njoy'], stdin=fl)
    fl.close()

    shutil.move("tape24", "groupr_out")
    shutil.move("njoy_input", self.name+"_njoy_input")
    shutil.move("output", self.name+"_njoy_output")
    os.remove("tape20")
    os.remove("tape21")
    os.remove("tape22")
    os.remove("tape23")

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
    out += "{mat} {ntemp} {ndil} 20 50 /\n".format(mat=self.mat, ntemp=len(self.temps), ndil=len(self.dilutions))
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
