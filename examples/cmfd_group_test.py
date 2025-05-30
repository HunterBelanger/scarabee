from scarabee import *
import numpy as np
import matplotlib.pyplot as plt

#For a problem with constant cross sections, CMFD performance should be the same with any group structure
#this takes 63 iterations with any CMFD but 83 with MOC

#constant factor to multiply XS by to simulate different optical thicknesses
tf = 1.0
'''
Et = np.array([0.5]*7)
Ea = np.array([0.2]*7)
Ef = np.array([0.1]*7)
nu = np.array([2.5]*7)
chi = np.array([0.001]*7)
Es = np.array([[0.3/7]*7,
               [0.3/7]*7,
               [0.3/7]*7,
               [0.3/7]*7,
               [0.3/7]*7,
               [0.3/7]*7,
               [0.3/7]*7,])
'''
#There isnt a check in scarabee to see if Ef > Et??
#There is no check in scarabee for sum(chi) != 1, and it only rebalances if XS are added
Et = np.array([0.4, 1.0])
Ea = np.array([0.0, 0.5])
nu = np.array([2.5, 2.5])
Ef = np.array([0.0, 0.01])
chi = np.array([0.01, 0.01])
Es = np.array([[0.4, 0.],
               [0., 0.5]])
UO2 = CrossSection(Et*tf, Ea*tf, Es*tf, Ef*tf, nu*Ef*tf, chi)

Et = np.array([0.2]*2)
Ea = np.array([0.05]*2)
Es = np.array([[0.15/2]*2,
               [0.15/2]*2,])
H2O = CrossSection(Et*tf, Ea*tf, Es*tf)

# Define Cells
radii = [0.54]
mats =  [UO2, H2O]
U = PinCell(radii, mats, 1.26, 1.26)
w = EmptyCell(H2O,1.26,1.26)

dx = [1.26]*5
dy = [1.26]*5
c2d = Cartesian2D(dx, dy)
c2d.set_tiles([U,U,U,w,w,
               U,U,U,w,w,
               U,U,U,w,w,
               w,w,w,w,w,
               w,w,w,w,w])


moc = MOCDriver(c2d)
#moc.cmfd = CMFD(dx, dy, [[0,0], [1,1], [2,2], [3,3], [4,4], [5,5], [6,6]])
#moc.cmfd = CMFD(dx, dy, [[0,4],[5,6]])
#moc.cmfd = CMFD(dx, dy, [[0,1]])
moc.cmfd = CMFD(dx, dy, [[0,0],[1,1]])
moc.generate_tracks(64, 0.01, YamamotoTabuchi6())
moc.keff_tolerance = 1.e-5
moc.flux_tolerance = 1.e-5
moc.solve()
