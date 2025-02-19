import pytest
import pytest_check as check
from scarabee import CMFD, CMFDSurfaceCrossing, CMFDSurfaceCrossingType, Direction, Vector

"""
The center of a 4 cell mesh is (1,0)
The y axis runs from [-1.0,1.0]
The x axis runs from [0.0,2.0]
"""
class TestCMFDSurfaces:
    @pytest.fixture
    
    def make_cmfd(self):
        cmfd = CMFD([1.0,1.0],[1.0,1.0],[(0,1)])
        return cmfd
    
    def test_0_to_1(self, make_cmfd):
        """
        Should return a valid crossing for the negative x surface of cell 1
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0-tol,-0.5)
        u = Direction(1.0,0.0)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 1)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.XN)

    def test_0_to_outside_left(self, make_cmfd):
        """
        Should return a valid crossing for the negative x surface of cell 0
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(0.0+tol,-0.5)
        u = Direction(-1.0,0.0)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 0)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.XN)

    def test_0_to_2(self, make_cmfd):
        """
        Should return a valid crossing for the negative y surface of cell 2
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(0.5,0.0-tol)
        u = Direction(0.0,1.0)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 2)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.YN)

    def test_0_to_bottom(self, make_cmfd):
        """
        Should return a valid crossing for the negative y surface of cell 0
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(0.5,-1.0+tol)
        u = Direction(0,-1)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 0)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.YN)


    def test_0_to_3(self, make_cmfd):
        """
        Should return a valid crossing for the top right surface of cell 0
        or BL of 3 
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0-tol,0.0-tol)
        u = Direction(0.707,0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 3)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.BL)

    def test_1_to_2(self, make_cmfd):
        """
        Should return a valid crossing for the top left surface of cell 1
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0+tol,0.0-tol)
        u = Direction(-0.707,0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 2)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.BR)

    def test_2_to_1(self, make_cmfd):
        """
        Should return a valid crossing for the bottom right surface of cell 2
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0-tol,0.0+tol)
        u = Direction(0.707,-0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 1)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.TL)

    def test_3_to_0(self, make_cmfd):
        """
        Should return a valid crossing for the bottom left surface of cell 3
        """
        cmfd = make_cmfd
        tol = 1e-13
        r = Vector(1.0+tol,0.0+tol)
        u = Direction(-0.707,-0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 0)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.TR)




