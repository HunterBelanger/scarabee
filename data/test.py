from frendy import FrendyMG
import h5py

base = "/mnt/c/Users/Hunter/Documents/nuclear_data/ENDF-B-VIII.0_neutrons/"
tslbase = "/mnt/c/Users/Hunter/Documents/nuclear_data/ENDF-B-VIII.0_thermal_scatt/"

h5 = h5py.File('endf8.h5', 'w')

N = FrendyMG()
N.endf_file = base + "n-001_H_001.endf"
N.tsl_file = tslbase + "tsl-HinH2O.endf"
#N.dilutions = [1.E1, 2.E1, 5.E1, 1.E2, 3.E2, 1.E3, 1.E4, 1.E5, 1.E6, 1.E8, 1.E10]
#N.temps = [293.6, 373., 559., 748., 793., 833., 963., 1273., 1773., 2573.]
N.temps = [283.6, 293.6, 300., 323.6, 350., 373.6, 400., 423.6, 450., 473.6, 500., 523.6, 550., 573.6, 600., 623.6, 650., 800.]
N.label = "H1 in H2O from ENDF/B-8.0"
N.name = "H1_H2O"
N.initialize()
N.process(h5)
N.add_to_hdf5(h5)


h5.close()