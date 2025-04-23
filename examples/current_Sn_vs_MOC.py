from scarabee import *
import numpy as np
import matplotlib.pyplot as plt 

Et = np.array([0.54628])
Ea = np.array([0.54628 - 0.464338])
Ef = np.array([0.054628])
nu = np.array([1.7])
chi = np.array([1.])
Es = np.array([[0.464338]])
xs = CrossSection(Et, Ea, Es, Ef, nu*Ef, chi)
#L = 10.371065
L=8.

N = 100

dx = N * [L/N]
xsary = N * [xs]

refsn = ReflectorSN(xsary, dx)
refsn.solve()


distances_sn = [sum(dx[:i]) + dx[i]/2 for i in range(len(dx))]

#Vertical MOC Slices
cell = EmptyCell(xs, L/N, L)
geom = Cartesian2D(dx, [L])
geom.set_tiles([cell]*N)

moc = MOCDriver(geom)
moc.x_min_bc = BoundaryCondition.Reflective
moc.x_max_bc = BoundaryCondition.Vacuum
moc.y_min_bc = BoundaryCondition.Reflective
moc.y_max_bc = BoundaryCondition.Reflective

moc.cmfd = CMFD(dx, [L], [[0,0]])

moc.generate_tracks(64, 0.01, YamamotoTabuchi6())
moc.flux_tolerance = 1.E-5
moc.keff_tolerance = 1.E-5
moc.solve()

sn_fluxes = np.zeros(len(dx))
moc_fluxes = np.zeros(len(dx))
for i in range(len(dx)):
    sn_fluxes[i] = refsn.flux(i,0,0)/dx[i]
    moc_fluxes[i] = moc.flux(i,0,0)/dx[i]

sn_normalize = 1/np.max(sn_fluxes)
moc_normalize = 1/(np.max(moc_fluxes))

sn_fluxes *= sn_normalize
moc_fluxes *= moc_normalize

nsurf_sn = refsn.nsurfaces

exit_current_sn = refsn.current(nsurf_sn-1,0)
#maybe this should be changed in the future to be consistent?
exit_current_moc = moc.cmfd.current(0,nsurf_sn-1)


distances_moc = [sum(dx[:i]) + dx[i]/2 for i in range(len(dx))]

print(f"""
Right side net currents 
Sn: {exit_current_sn}
MOC w/CMFD: {exit_current_moc}

Normalized currents
Sn: {exit_current_sn * sn_normalize}
MOC: {exit_current_moc * moc_normalize}

MOC/Sn: {(exit_current_moc * moc_normalize)/(exit_current_sn * sn_normalize)}
      """)

plt.plot(distances_sn,sn_fluxes,label="Sn fluxes")
plt.plot(distances_moc,moc_fluxes,label="MOC fluxes")
plt.xlabel("Distance [cm]")
plt.ylabel("Normalized flux")
plt.title("Sn and MOC fluxes")
plt.yscale('log')
plt.legend()
plt.show()
plt.close()