from scarabee import *
import numpy as np

Et = np.array([0.32640])
Ea = np.array([0.101184])
Ef = np.array([0.081600])
nu = np.array([3.24])
chi = np.array([1.])
Es = np.array([[0.225216]])
xs = CrossSection(Et, Ea, Es, Ef, nu*Ef, chi)
L = 1.853722

#Et = np.array([0.54628])
#Ea = np.array([0.54628 - 0.464338])
#Ef = np.array([0.054628])
#nu = np.array([1.7])
#chi = np.array([1.])
#Es = np.array([[0.464338]])
#xs = CrossSection(Et, Ea, Es, Ef, nu*Ef, chi)
#L = 10.371065

#Et = np.array([0.33588, 0.54628])
#Ea = np.array([0.33588 - (0.3198 + 0.004555), 0.54628 - 0.42410])
#Ef = np.array([0.002817, 0.097])
#nu = np.array([2.5, 2.5])
#chi = np.array([1., 0.])
#Es = np.array([[0.3198, 0.004555],
#               [0.,     0.42410]])
#xs = CrossSection(Et, Ea, Es, Ef, nu*Ef, chi)
#L = 846.632726

N = 1

dx = N * [L/N]
xsary = N * [xs]

refsn = ReflectorSN(xsary, dx)
refsn.solve()