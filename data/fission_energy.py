import pyPapillonNDL as pndl
import h5py
import sys

def get_isotope(NucName):
  iso = ""

  for c in NucName:
    if c == '_':
      break
    else:
      iso += c
  
  return iso

def get_element(NucName):
  elem = ""

  for c in NucName:
    if c.isalpha():
      elem += c
    else:
      break
  
  return elem

def main():
  if len(sys.argv) != 2:
    raise RuntimeError("Script takes 1 argument: File with MG data library.")
  
  mg_lib = h5py.File(sys.argv[1], 'r+')

  # Go through all nuclides in the library
  for NucName in mg_lib.keys():
    NucGrp = mg_lib[NucName]

    # Skip non-fissile nuclides
    if NucGrp.attrs['fissile'] == False:
      continue

    iso = get_isotope(NucName)
    elem = get_element(iso)

    # Get the name of the ACE file for the nuclide
    ace_fname = "/mnt/c/Users/hunte/Documents/nuclear_data/abeille_ace_files/endf8/data/neutron_dir/" + elem + "/" + iso + "/" + iso + ".293.6.ace.bin"

    # Get the ACE
    ace = pndl.ACE(ace_fname, pndl.ACEType.BINARY)
    nuc = pndl.STNeutron(ace)

    # Get the fission reaction
    fiss = None
    if nuc.has_reaction(18):
      fiss = nuc.reaction(18)
    elif nuc.has_reaction(19):
      fiss = nuc.reaction(19)
    elif nuc.has_reaction(20):
      fiss = nuc.reaction(20)
    elif nuc.has_reaction(21):
      fiss = nuc.reaction(21)
    elif nuc.has_reaction(38):
      fiss = nuc.reaction(38)
    else:
      raise RuntimeError("Could not find fission reaction for {:}".format(NucName))

    # Get energy release (in MeV)
    fiss_q = fiss.q()

    NucGrp.attrs['fission-energy'] = fiss_q


if __name__ == "__main__":
  main()