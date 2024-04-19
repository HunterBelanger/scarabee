from njoy import NjoyFreeGasMG
from frendy import FrendyFreeGasMG

N = FrendyFreeGasMG()
N.endf_file = "n-092_U_238.endf"
N.dilutions = [1.E1, 2.E1, 5.E1, 1.E2, 3.E2, 1.E3, 1.E4, 1.E5, 1.E6, 1.E8]
N.temp = 293.6
N.label = "U238 from ENDF/B-8.0"
N.name = "U238"
N.process()

print(N.ZA)
print(N.awr)
print(N.pot_xs)

#N.process()