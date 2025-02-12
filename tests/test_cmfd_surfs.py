import pytest
import pytest_check as check
from scarabee import CMFD, CMFDSurfaceCrossing, CMFDSurfaceCrossingType, Direction, Vector
class TestCMFDSurfaces:
    @pytest.fixture
    
    def make_cmfd(self):
        cmfd = CMFD([1.0,1.0],[1.0,1.0],[(0,1)])
        return cmfd
    
    def test_xp(self, make_cmfd):
        """
        Should return a valid crossing for the positive x surface of cell 0
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0-tol,0.5)
        u = Direction(0.707,0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 0)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.XP)

    def test_xn(self, make_cmfd):
        """
        Should return a valid crossing for the negative x surface of cell 0
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(0.0+tol,0.5)
        u = Direction(-0.707,-0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 0)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.XN)

    def test_yp(self, make_cmfd):
        """
        Should return a valid crossing for the positive y surface of cell 0
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(0.5,1.0-tol)
        u = Direction(0.707,0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 0)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.YP)

    def test_yn(self, make_cmfd):
        """
        Should return a valid crossing for the negative y surface of cell 0
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(0.5,0.0+tol)
        u = Direction(0,-1)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 0)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.YN)


    def test_top_right(self, make_cmfd):
        """
        Should return a valid crossing for the top right surface of cell 0
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0-tol,1.0-tol)
        u = Direction(0.707,0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 0)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.TR)

    def test_TR_to_XP(self, make_cmfd):
        """
        Should return a valid crossing for the XP surface of cell 0
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0-tol,1.0-tol)
        u = Direction(0.707,-0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 0)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.XP)

    def test_TR_to_YP(self, make_cmfd):
        """
        Should return a valid crossing for the YP surface of cell 0
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0-tol,1.0-tol)
        u = Direction(-0.707,0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 0)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.YP)
    
    def test_TR_backwards(self, make_cmfd):
        """
        Should return a valid crossing for the TR surface of cell 0
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0-tol,1.0-tol)
        u = Direction(-0.707,-0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 0)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.TR)

    def test_top_left(self, make_cmfd):
        """
        Should return a valid crossing for the top left surface of cell 1
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0+tol,1.0-tol)
        u = Direction(-0.707,0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 1)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.TL)
    
    def test_TL_to_YP(self, make_cmfd):
        """
        Should return a valid crossing for YP surface of cell 1
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0+tol,1.0-tol)
        u = Direction(0.707,0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 1)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.YP)

    def test_TL_to_XN(self, make_cmfd):
        """
        Should return a valid crossing for XN surface of cell 1
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0+tol,1.0-tol)
        u = Direction(-0.707,-0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 1)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.XN)
    
    def test_TL_backwards(self, make_cmfd):
        """
        Should return a valid crossing for TL surface of cell 1
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0+tol,1.0-tol)
        u = Direction(0.707,-0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 1)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.TL)

    def test_bottom_right(self, make_cmfd):
        """
        Should return a valid crossing for the bottom right surface of cell 2
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0-tol,1.0+tol)
        u = Direction(0.707,-0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 2)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.BR)
    
    def test_BR_to_XP(self, make_cmfd):
        """
        Should return a valid crossing for XP surface of cell 2
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0-tol,1.0+tol)
        u = Direction(0.707,0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 2)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.YN)
    
    def test_BR_to_YN(self, make_cmfd):
        """
        Should return a valid crossing for YN surface of cell 2
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0-tol,1.0+tol)
        u = Direction(-0.707,-0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 2)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.YN)

    def test_BR_backwards(self, make_cmfd):
        """
        Should return a valid crossing for BR surface of cell 2
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0-tol,1.0+tol)
        u = Direction(-0.707,0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 2)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.BR)

    def test_bottom_left(self, make_cmfd):
        """
        Should return a valid crossing for the bottom left surface of cell 3
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0+tol,1.0+tol)
        u = Direction(-0.707,-0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 3)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.BL)
    
    def test_BL_to_XN(self, make_cmfd):
        """
        Should return a valid crossing for XN surface of cell 3
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0+tol,1.0+tol)
        u = Direction(-0.707,0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 3)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.XN)
    
    def test_BL_to_YN(self, make_cmfd):
        """
        Should return a valid crossing for YN surface of cell 3
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0+tol,1.0+tol)
        u = Direction(0.707,-0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 3)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.YN)

    def test_BL_backwards(self, make_cmfd):
        """
        Should return a valid crossing for BL surface of cell 3
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0-tol,1.0+tol)
        u = Direction(0.707,0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 3)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.BL)


