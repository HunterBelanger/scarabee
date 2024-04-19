from njoy import NjoyFreeGasMG

N = NjoyFreeGasMG()
N.endf_file = "n-092_U_238.endf"
#N.temps = [293.6, 373., 559., 748., 793., 833., 963., 1273., 1773., 2573.]
#N.dilutions = [1.E1, 2.E1, 5.E1, 1.E2, 3.E2, 1.E3, 1.E4, 1.E5, 1.E6, 1.E8]
N.temps = [293.6]
N.label = "U238 from ENDF/B-8.0"
N.name = "U238"

N.process()