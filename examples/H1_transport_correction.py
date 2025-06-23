from scarabee import *
import numpy as np
import matplotlib.pyplot as plt

ndl = NDLibrary()

# Get U235 so we can steal the fission spectrum
U235_comp = MaterialComposition(Fraction.Atoms)
U235_comp.add_nuclide("U235", 1.)
U235_mat = Material(U235_comp, 293.6, 1., DensityUnits.a_bcm, ndl)
U235 = U235_mat.dilution_xs([1.E10], ndl)
G = U235.ngroups
chi = np.zeros(G)
for g in range(G):
  chi[g] = U235.chi(g)

# We now load H1_H2O
H1_comp = MaterialComposition(Fraction.Atoms)
H1_comp.add_nuclide("H1_H2O", 1.)
H1_mat = Material(H1_comp, 293.6, 1., DensityUnits.a_bcm, ndl)
H1 = H1_mat.dilution_xs([1.E10], ndl)

# We must create a new nuclide with the fission spectrum
Et = np.zeros(G)
Dtr = np.zeros(G)
Ea = np.zeros(G)
Ef = np.zeros(G)
vEf = np.zeros(G)
Es = np.zeros((2, G, G))
for g in range(G):
  Et[g] = H1.Et(g)
  Ea[g] = H1.Ea(g)
  Ef[g] = H1.Ef(g)
  vEf[g] = H1.vEf(g)
  for gg in range(G):
    Es[0, g, gg] = H1.Es(0, g, gg)
    Es[1, g, gg] = H1.Es(1, g, gg)

# Create a new temporary nuclide
TempH1 = CrossSection(Et, Dtr, Ea, Es, Ef, vEf, chi)

# We now perform a P1 leakage calculation
P1_spectrum = P1CriticalitySpectrum(TempH1, 0.0001)

# We now have diffusion coefficients
D = P1_spectrum.diff_coeff

# Compute transport xs
Etr = 1. / (3. * D)

# Now we compute the transport xs ratio
ratio = Etr / Et

# Make the plot
plt.stairs(values=ratio, edges=ndl.group_bounds)
plt.xlabel('Energy [eV]')
plt.ylabel(r'$\Sigma_{tr} / \Sigma_t$ for H$^1$ in H$_2$O')
plt.xscale('log')
plt.show()

# Calculate the delta xs for the transport correction
Delta = Et - Etr

# Correct the xs data for the nuclide
for g in range(G):
  Et[g] -= Delta[g]
  Es[0, g, g] -= Delta[g]

plt.stairs(values=Et, edges=ndl.group_bounds)
plt.xlabel('Energy [eV]')
plt.ylabel('Transport Cross Section [barns]')
plt.xscale('log')
plt.yscale('log')
plt.show()
