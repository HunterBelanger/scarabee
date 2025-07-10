import pytest
import pytest_check as check
from scarabee import CMFD, CMFDSurfaceCrossing, CMFDSurfaceCrossingType, Direction, Vector
"""
Unit tests for CMFD::get_surface() method

The test mesh consists of a 2x2 rectangular grid with dx, dy = 1.0,1.0
"""
class TestCMFDSurfaces:
    @pytest.fixture
    
    def make_cmfd(self):
        """
        Makes the test grid. Cell indexes are:
        [2, 3]
        [0, 1]
        """
        cmfd = CMFD([1.0,1.0],[1.0,1.0],[(0,1)])
        return cmfd
    
    def test_cell_XP_to_XN(self, make_cmfd):
        """
        Tests get_surface() for a point, direction pair that goes from 
        the positive x surface of a cell into the negative x surface of 
        an adjacent cell.

        Should return a valid crossing for the negative x surface of cell 1
        """
        cmfd = make_cmfd
        tol = 1e-12
        r = Vector(0.0-tol,-0.5)
        u = Direction(1.0,0.0)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 1)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.XN)

    def test_XN_to_outside_left(self, make_cmfd):
        """
        Tests get_surface() for a point, direction pair that goes 
        from the left of a cell to an exterior boundary
        Should throw an error
        """
        cmfd = make_cmfd
        tol = 1e-12
        r = Vector(-1.0+tol,-0.5)
        u = Direction(-1.0,0.0)
        with pytest.raises(RuntimeError):
            surface = cmfd.get_surface(r,u)

    def test_YP_to_YN(self, make_cmfd):
        """
        Tests get_surface() for a point, direction pair that goes from 
        the positive y surface of a cell into the negative y surface of 
        an adjacent cell.

        Should return a valid crossing for the negative y surface of cell 2
        """
        cmfd = make_cmfd
        tol = 1e-12
        r = Vector(-0.5,0.0-tol)
        u = Direction(0.0,1.0)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 2)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.YN)

    def test_YN_to_outside_bottom(self, make_cmfd):
        """
        Tests get_surface() for a point, direction pair that goes 
        from the bottom of a cell to an exterior boundary
        Should throw an error
        """
        cmfd = make_cmfd
        tol = 1e-12
        r = Vector(-0.5,-1.0+tol)
        u = Direction(0,-1)

        with pytest.raises(RuntimeError):
            surface = cmfd.get_surface(r,u)


    def test_I_to_III(self, make_cmfd):
        """
        Tests get_surface() for a point, direction pair that goes from 
        the top right surface of a cell into the bottom left surface of 
        an adjacent cell.

        Should return a valid crossing for the III surface of cell 3 
        """
        cmfd = make_cmfd
        tol = 1e-12
        r = Vector(0.0-tol,0.0-tol)
        u = Direction(0.707,0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 3)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.III)

    def test_II_to_IV(self, make_cmfd):
        """
        Tests get_surface() for a point, direction pair that goes from 
        the top left surface of a cell into the bottom right surface of 
        an adjacent cell.

        Should return a valid crossing for the IV surface of cell 2
        """
        cmfd = make_cmfd
        tol = 1e-12
        r = Vector(0.0+tol,0.0-tol)
        u = Direction(-0.707,0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 2)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.IV)

    def test_IV_to_II(self, make_cmfd):
        """
        Tests get_surface() for a point, direction pair that goes from 
        the bottom right surface of a cell into the top left surface of 
        an adjacent cell.

        Should return a valid crossing for the top left surface of cell 1
        """
        cmfd = make_cmfd
        tol = 1e-12
        r = Vector(0.0-tol,0.0+tol)
        u = Direction(0.707,-0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 1)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.II)

    def test_III_to_I(self, make_cmfd):
        """
        Tests get_surface() for a point, direction pair that goes from 
        the bottom left surface of a cell into the top right surface of 
        an adjacent cell.

        Should return a valid crossing for the top right surface of cell 0
        """
        cmfd = make_cmfd
        tol = 1e-12
        r = Vector(0.0+tol,0.0+tol)
        u = Direction(-0.707,-0.707)

        surface = cmfd.get_surface(r,u)

        check.equal(surface.is_valid, True)
        check.equal(surface.cell_index, 0)
        check.equal(surface.crossing, CMFDSurfaceCrossingType.I)




