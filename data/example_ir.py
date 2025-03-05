import pyPapillonNDL as pndl
import matplotlib.pyplot as plt
import numpy as np
from ir_lambda import generate_U238_U235_ir_lambda
import sys
import h5py

def main():
  if len(sys.argv) != 2:
    raise RuntimeError("Script takes 1 argument: File with MG data library.")
  
  mg_lib = h5py.File(sys.argv[1], 'r+')

  mg_energy_bounds = np.array(mg_lib.attrs['group-bounds'])

  ace = pndl.ACE("/mnt/c/Users/BELANH2/Documents/nuclear_data/ENDF-VIII.0/abeille/data/neutron_dir/U/U238/U238.600.0.ace")
  U238 = pndl.STNeutron(ace)
  
  ace = pndl.ACE("/mnt/c/Users/BELANH2/Documents/nuclear_data/ENDF-VIII.0/abeille/data/neutron_dir/U/U235/U235.600.0.ace")
  U235 = pndl.STNeutron(ace)

  for NucName in mg_lib.keys():
    NucGrp = mg_lib[NucName]
    awr = NucGrp.attrs['awr']
    if awr <= 1.:
      continue
    print(">>> Processing {:}".format(NucName))
    sig_pot = NucGrp.attrs['potential-xs']
    ir_lmbda = generate_U238_U235_ir_lambda(U238, U235, awr=awr, sig_pot=sig_pot, mg_energy_bounds=mg_energy_bounds)
    NucGrp.attrs['ir-lambda'] = ir_lmbda

  mg_lib.close() 

if __name__ == "__main__":
  main()
