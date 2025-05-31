from scarabee import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
'''
This includes a small 3x3 problem surrounded with water to create a fast benchmark for testing changes to the CMFD class
MOC without CMFD:
Iteration  304          keff: 0.95737
'''

Et = np.array([1.77949E-01, 3.29805E-01, 4.80388E-01, 5.54367E-01, 3.11801E-01, 3.95168E-01, 5.64406E-01])
Ea = np.array([8.02480E-03, 3.71740E-03, 2.67690E-02, 9.62360E-02, 3.00200E-02, 1.11260E-01, 2.82780E-01])
Ef = np.array([7.21206E-03, 8.19301E-04, 6.45320E-03, 1.85648E-02, 1.78084E-02, 8.30348E-02, 2.16004E-01])
nu = np.array([2.78145, 2.47443, 2.43383, 2.43380, 2.43380, 2.43380, 2.43380])
chi = np.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0., 0., 0.])
Es = np.array([[1.27537E-01, 4.23780E-02, 9.43740E-06, 5.51630E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00],
               [0.00000E+00, 3.24456E-01, 1.63140E-03, 3.14270E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00],
               [0.00000E+00, 0.00000E+00, 4.50940E-01, 2.67920E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00],
               [0.00000E+00, 0.00000E+00, 0.00000E+00, 4.52565E-01, 5.56640E-03, 0.00000E+00, 0.00000E+00],
               [0.00000E+00, 0.00000E+00, 0.00000E+00, 1.25250E-04, 2.71401E-01, 1.02550E-02, 1.00210E-08],
               [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.29680E-03, 2.65802E-01, 1.68090E-02],
               [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 8.54580E-03, 2.73080E-01]])
UO2xs = CrossSection(Et, Ea, Es, Ef, nu*Ef, chi)


Et = np.array([1.83045E-01, 3.36705E-01, 5.00507E-01, 6.06174E-01, 5.02754E-01, 9.21028E-01, 9.55231E-01])
Ea = np.array([9.48620E-03, 4.65560E-03, 3.62400E-02, 1.32720E-01, 2.08400E-01, 6.58700E-01, 6.90170E-01])
Ef = np.array([8.67209E-03, 1.62426E-03, 1.02716E-02, 3.90447E-02, 1.92576E-02, 3.74888E-01, 4.30599E-01])
nu = np.array([2.90426,     2.91795,     2.86986,     2.87491,     2.87175,     2.86752,     2.87808])
chi = np.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0.,          0.,          0.])
Es = np.array([[1.31504E-01, 4.20460E-02, 8.69720E-06, 5.19380E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00],
               [0.00000E+00, 3.30403E-01, 1.64630E-03, 2.60060E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00],
               [0.00000E+00, 0.00000E+00, 4.61792E-01, 2.47490E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00],
               [0.00000E+00, 0.00000E+00, 0.00000E+00, 4.68021E-01, 5.43300E-03, 0.00000E+00, 0.00000E+00],
               [0.00000E+00, 0.00000E+00, 0.00000E+00, 1.85970E-04, 2.85771E-01, 8.39730E-03, 8.92800E-09],
               [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.39160E-03, 2.47614E-01, 1.23220E-02],
               [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 8.96810E-03, 2.56093E-01]])
M8xs = CrossSection(Et, Ea, Es, Ef, nu*Ef, chi)

Et = np.array([1.59206E-01, 4.12970E-01, 5.90310E-01, 5.84350E-01, 7.18000E-01, 1.25445E+00, 2.65038E+00])
Ea = np.array([6.01050E-04, 1.57930E-05, 3.37160E-04, 1.94060E-03, 5.74160E-03, 1.50010E-02, 3.72390E-02])
Es = np.array([[4.44777E-02, 1.13400E-01, 7.23470E-04, 3.74990E-06, 5.31840E-08, 0.00000E+00, 0.00000E+00],
               [0.00000E+00, 2.82334E-01, 1.29940E-01, 6.23400E-04, 4.80020E-05, 7.44860E-06, 1.04550E-06],
               [0.00000E+00, 0.00000E+00, 3.45256E-01, 2.24570E-01, 1.69990E-02, 2.64430E-03, 5.03440E-04],
               [0.00000E+00, 0.00000E+00, 0.00000E+00, 9.10284E-02, 4.15510E-01, 6.37320E-02, 1.21390E-02],
               [0.00000E+00, 0.00000E+00, 0.00000E+00, 7.14370E-05, 1.39138E-01, 5.11820E-01, 6.12290E-02],
               [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.21570E-03, 6.99913E-01, 5.37320E-01],
               [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.32440E-01, 2.48070E+00]])
H2Oxs = CrossSection(Et, Ea, Es)

# Define Cells
pitch = 1.26

radii = [0.4,   0.54,  0.63]
mats =  [UO2xs, UO2xs, H2Oxs, H2Oxs]
U2 = PinCell(radii, mats, pitch, pitch)

mats =  [M8xs, M8xs, H2Oxs, H2Oxs]
M8 = PinCell(radii, mats, pitch, pitch)
wc = EmptyCell(H2Oxs,1.26,1.26)

dx = [1.26]*5
dy = [1.26]*5
c2d = Cartesian2D(dx, dy)
c2d.set_tiles([M8,U2,M8,wc,wc,
               U2,M8,U2,wc,wc,
               M8,U2,M8,wc,wc,
               wc,wc,wc,wc,wc,
               wc,wc,wc,wc,wc])

#cmfd_groups = [[0,0], [1,1], [2,2], [3,3], [4,4], [5,5], [6,6]]
cmfd_groups = [[0,3],[4,6]]
#cmfd_groups = [[0,6]]
moc_to_cmfd_group_map = np.zeros(len(Et),dtype=np.int32)
for i, pair in enumerate(cmfd_groups):
    for g in range(pair[0],pair[1]+1,1):
        moc_to_cmfd_group_map[g] = int(i)

dx_cmfd = [1.26]*5
dy_cmfd = [1.26]*5

#Note: The CMFD constructor should be updated to include an error message for if dx and dy do not match the size of the problem domain,
# right now it is an unhandled error
moc = MOCDriver(c2d)
moc.x_min_bc = BoundaryCondition.Reflective
moc.x_max_bc = BoundaryCondition.Reflective
moc.y_min_bc = BoundaryCondition.Reflective
moc.y_max_bc = BoundaryCondition.Reflective
moc.cmfd = CMFD(dx_cmfd, dy_cmfd, cmfd_groups)
moc.cmfd.keff_tolerance = 1.e-5
moc.cmfd.flux_tolerance = 1.e-5
#moc.cmfd.set_damping(0.7)
moc.generate_tracks(64, 0.05, YamamotoTabuchi6())
moc.keff_tolerance = 1.e-5
moc.flux_tolerance = 1.e-5
moc.solve()

cmfd_flux = np.zeros((len(dx), len(dy), len(cmfd_groups)))

#Get MOC flux
flux_moc, xm, ym = moc.rasterize_flux(100, 100)

#Get CMFD flux in each cell
for i in range(len(dx)):
    for j in range(len(dy)):
        for g in range(len(cmfd_groups)):
            cmfd_flux[i, j, g] = moc.cmfd.flux(i, j, g)

for g in range(moc_to_cmfd_group_map.size):
    flux_moc_g = flux_moc[g, :, :]
    G = moc_to_cmfd_group_map[g]
    flux_cmfd_g = cmfd_flux[:, :, G].T 

    #Normalize
    flux_moc_norm = flux_moc_g / np.max(flux_moc_g)
    flux_cmfd_norm = flux_cmfd_g / np.max(flux_cmfd_g)

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))

    #Plot MOC flux
    im0 = axs[0].pcolormesh(xm, ym, flux_moc_norm, cmap='jet')
    axs[0].set_title(f"MOC Flux (Group {g+1})")
    cbar0 = fig.colorbar(ScalarMappable(cmap='jet'), ax=axs[0])
    cbar0.set_label("MOC Flux (normalized)")

    #Plot CMFD flux
    im1 = axs[1].imshow(flux_cmfd_norm, origin='lower', cmap='jet')
    axs[1].set_title(f"CMFD Flux (Group {G+1})")
    axs[1].set_xlabel('i')
    axs[1].set_ylabel('j')
    cbar1 = fig.colorbar(ScalarMappable(cmap='jet'), ax=axs[1])
    cbar1.set_label("CMFD Flux (normalized)")

    plt.savefig(f"MOC_CMFD_compare_g{g}.png")
    plt.close()
