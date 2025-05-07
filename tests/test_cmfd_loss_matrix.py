import pytest
import pytest_check as check
from scarabee import *
import numpy as np
import scipy.sparse

class TestLossMatrix:

    def test_sparsity_pattern(self):
        #Tests sparsity pattern has correct # of entries and size
        
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
        UO2 = CrossSection(Et, Ea, Es, Ef, nu*Ef, chi)

        Et = np.array([1.59206E-01, 4.12970E-01, 5.90310E-01, 5.84350E-01, 7.18000E-01, 1.25445E+00, 2.65038E+00])
        Ea = np.array([6.01050E-04, 1.57930E-05, 3.37160E-04, 1.94060E-03, 5.74160E-03, 1.50010E-02, 3.72390E-02])
        Es = np.array([[4.44777E-02, 1.13400E-01, 7.23470E-04, 3.74990E-06, 5.31840E-08, 0.00000E+00, 0.00000E+00],
                    [0.00000E+00, 2.82334E-01, 1.29940E-01, 6.23400E-04, 4.80020E-05, 7.44860E-06, 1.04550E-06],
                    [0.00000E+00, 0.00000E+00, 3.45256E-01, 2.24570E-01, 1.69990E-02, 2.64430E-03, 5.03440E-04],
                    [0.00000E+00, 0.00000E+00, 0.00000E+00, 9.10284E-02, 4.15510E-01, 6.37320E-02, 1.21390E-02],
                    [0.00000E+00, 0.00000E+00, 0.00000E+00, 7.14370E-05, 1.39138E-01, 5.11820E-01, 6.12290E-02],
                    [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.21570E-03, 6.99913E-01, 5.37320E-01],
                    [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.32440E-01, 2.48070E+00]])
        H2O = CrossSection(Et, Ea, Es)

        # Define Cells
        radii = [0.54]
        mats =  [UO2, H2O]
        U = PinCell(radii, mats, 1.26, 1.26)
        w = EmptyCell(H2O,1.26,1.26)

        dx = [1.26]*4
        dy = [1.26]*4
        c2d = Cartesian2D(dx, dy)
        c2d.set_tiles([w,U,U,w] + [w,U,U,w] + [w,U,U,w] + [w,U,U,w])

        moc = MOCDriver(c2d,)
        moc.x_min_bc = BoundaryCondition.Reflective
        moc.x_max_bc = BoundaryCondition.Reflective
        moc.y_min_bc = BoundaryCondition.Reflective
        moc.y_max_bc = BoundaryCondition.Reflective
        moc.cmfd = CMFD(dx, dy, [[0,0],[1,1],[2,2],[3,3],[4,4],[5,5],[6,6]])
        moc.cmfd.set_solve(1) #turn off CMFD solve, currently broken
        moc.generate_tracks(64, 0.01, YamamotoTabuchi6())
        #problem doesn't need to converge, just run at least one iteration to make CMFD matrix
        moc.keff_tolerance = 1.e-2
        moc.flux_tolerance = 1.e-2
        moc.solve()

        # after moc.solve()…
        #A is a scipy.sparse.csr_matrix
        A = moc.cmfd.get_loss_matrix()

        #Based on OpenMOC sparsity pattern for this problem, should have  NNZ entries and 7*4*4 by 7*4*4 size 
        check.equal(A.nnz,448+6*4*4*7)
        check.equal(A.shape, (4*4*7,4*4*7))

    #Next test is to check positions of individual components 
