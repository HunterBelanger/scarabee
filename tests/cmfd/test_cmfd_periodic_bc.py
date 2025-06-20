import pytest
import pytest_timeout
import numpy as np
from scarabee import *

"""
Integration test for CMFD with periodic BC 
"""

class TestCMFDPeriodicBC:

    # This problem should never take this long to run
    # if it takes longer, most likely didn't converge
    @pytest.mark.timeout(100)
    def test_UO2_assembly_periodic(self):
        '''
        Test UO2 assembly with Periodic BC
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

        Et = np.array([1.26032E-01, 2.93160E-01, 2.84240E-01, 2.80960E-01, 3.34440E-01, 5.65640E-01, 1.17215E+00])
        Ea = np.array([5.11320E-04, 7.58010E-05, 3.15720E-04, 1.15820E-03, 3.39750E-03, 9.18780E-03, 2.32420E-02])
        Es = np.array([[6.61659E-02, 5.90700E-02, 2.83340E-04, 1.46220E-06, 2.06420E-08, 0.00000E+00, 0.00000E+00],
                    [0.00000E+00, 2.40377E-01, 5.24350E-02, 2.49900E-04, 1.92390E-05, 2.98750E-06, 4.21400E-07],
                    [0.00000E+00, 0.00000E+00, 1.83297E-01, 9.23970E-02, 6.94460E-03, 1.08030E-03, 2.05670E-04],
                    [0.00000E+00, 0.00000E+00, 0.00000E+00, 7.88511E-02, 1.70140E-01, 2.58810E-02, 4.92970E-03],
                    [0.00000E+00, 0.00000E+00, 0.00000E+00, 3.73330E-05, 9.97372E-02, 2.06790E-01, 2.44780E-02],
                    [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 9.17260E-04, 3.16765E-01, 2.38770E-01],
                    [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 4.97920E-02, 1.09912E+00]])
        GT = CrossSection(Et, Ea, Es)

        # Define Cells
        radii = [0.4, 0.54, 0.63]
        mats =  [UO2, UO2,  H2O, H2O]
        U = PinCell(radii, mats, 1.26, 1.26)

        mats =  [GT, GT, H2O, H2O]
        G = PinCell(radii, mats, 1.26, 1.26)

        dx = [1.26]*17
        c2d = Cartesian2D(dx, dx)
        c2d.set_tiles([U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,
                    U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,
                    U,U,U,U,U,G,U,U,G,U,U,G,U,U,U,U,U,
                    U,U,U,G,U,U,U,U,U,U,U,U,U,G,U,U,U,
                    U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,
                    U,U,G,U,U,G,U,U,G,U,U,G,U,U,G,U,U,
                    U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,
                    U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,
                    U,U,G,U,U,G,U,U,G,U,U,G,U,U,G,U,U,
                    U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,
                    U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,
                    U,U,G,U,U,G,U,U,G,U,U,G,U,U,G,U,U,
                    U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,
                    U,U,U,G,U,U,U,U,U,U,U,U,U,G,U,U,U,
                    U,U,U,U,U,G,U,U,G,U,U,G,U,U,U,U,U,
                    U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,
                    U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U,U])

        moc = MOCDriver(c2d)
        moc.x_min_bc = BoundaryCondition.Periodic
        moc.x_max_bc = BoundaryCondition.Periodic
        moc.y_min_bc = BoundaryCondition.Periodic
        moc.y_max_bc = BoundaryCondition.Periodic
        moc.cmfd = CMFD(dx,dx,[[0,1],[2,4],[5,6]])
        moc.generate_tracks(64, 0.05, YamamotoTabuchi6())
        moc.keff_tolerance = 1.E-5
        moc.flux_tolerance = 1.E-5
        moc.solve()

        assert moc.keff == pytest.approx(1.33381, 0.0001)