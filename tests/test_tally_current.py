import pytest
import pytest_check as check
from scarabee import CMFD, CMFDSurfaceCrossing, CMFDSurfaceCrossingType, Direction, Vector
"""
Unit tests for CMFD::tally_current() method

The test mesh consists of a 3x3 rectangular grid with dx, dy = 1.0,1.0

Tests for interior sides:
test_interior_xp
test_interior_xn
test_interior_yn
test_interior_yp

Tests for interior corners:
test_interior_TR
test_interior_BR
test_interior_BL
test_interior_TL

Tests for exterior sides:
test_exterior_xp
test_exterior_xn
test_exterior_yp
test_exterior_yn

Tests for exterior corners:
test_exterior_TR_TrueCorner
test_exterior_TR_XPCorner
test_exterior_TR_YPCorner
test_exterior_BR_TrueCorner
test_exterior_BR_XPCorner
test_exterior_BR_YNCorner
test_exterior_BL_TrueCorner
test_exterior_BL_XNCorner
test_exterior_BL_YNCorner
test_exterior_TL_TrueCorner
test_exterior_TL_XNCorner
test_exterior_TL_YPCorner
"""
class TestTallyCurrent:
    @pytest.fixture
    
    def make_cmfd(self):
        """
        Makes the test grid. Cell indexes are:
        [6, 7, 8]
        [3, 4, 5]
        [0, 1, 2]
        """
        cmfd = CMFD([1.0,1.0,1.0],[1.0,1.0,1.0],[(0,1)])
        return cmfd
    
    def test_interior_xp(self, make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(1.0,0.0)
        r = Vector(0.5-tol,0.0)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        test_current = cmfd.current(G,6)

        assert test_current == 1.0

    def test_interior_yp(self, make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(0.0,1.0)
        r = Vector(0.0,0.5-tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        test_current = cmfd.current(G,18)

        assert test_current == 1.0

    def test_interior_xn(self, make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(-1.0,0.0)
        r = Vector(-0.5+tol,0.0)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        test_current = cmfd.current(G,5)

        assert test_current == -1.0
    
    def test_interior_yn(self, make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(0.0,-1.0)
        r = Vector(0.0,-0.5+tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        test_current = cmfd.current(G,17)

        assert test_current == -1.0

    def test_interior_TR(self,make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(0.707,0.707)
        r = Vector(0.5-tol,0.5-tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        tc_1 = cmfd.current(G,6)
        tc_2 = cmfd.current(G,18)
        tc_3 = cmfd.current(G,10)
        tc_4 = cmfd.current(G,22)

        check.equal(tc_1, 0.5)
        check.equal(tc_2, 0.5)
        check.equal(tc_3, 0.5)
        check.equal(tc_4, 0.5)

    def test_interior_BR(self,make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(0.707,-0.707)
        r = Vector(0.5-tol,-0.5+tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        tc_1 = cmfd.current(G,6)
        tc_2 = cmfd.current(G,17)
        tc_3 = cmfd.current(G,21)
        tc_4 = cmfd.current(G,2)

        check.equal(tc_1, 0.5)
        check.equal(tc_2, -0.5)
        check.equal(tc_3, -0.5)
        check.equal(tc_4, 0.5)

    def test_interior_BL(self,make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(-0.707,-0.707)
        r = Vector(-0.5+tol,-0.5+tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        tc_1 = cmfd.current(G,5)
        tc_2 = cmfd.current(G,17)
        tc_3 = cmfd.current(G,1)
        tc_4 = cmfd.current(G,13)

        check.equal(tc_1, -0.5)
        check.equal(tc_2, -0.5)
        check.equal(tc_3, -0.5)
        check.equal(tc_4, -0.5)

    def test_interior_TL(self,make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(-0.707,0.707)
        r = Vector(-0.5+tol,0.5-tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        tc_1 = cmfd.current(G,9)
        tc_2 = cmfd.current(G,18)
        tc_3 = cmfd.current(G,5)
        tc_4 = cmfd.current(G,14)

        check.equal(tc_1, -0.5)
        check.equal(tc_2, 0.5)
        check.equal(tc_3, -0.5)
        check.equal(tc_4, 0.5)
    
    def test_exterior_yp(self,make_cmfd):
        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(0.0,-1.0)
        r = Vector(0.0,1.5-tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        test_current = cmfd.current(G,19)

        check.equal(test_current, -1.0)
        check.equal(surf.is_valid, True)

    def test_exterior_xp(self,make_cmfd):
        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(-1.0,0.0)
        r = Vector(1.5-tol,0.0)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        test_current = cmfd.current(G,7)

        check.equal(test_current, -1.0)
        check.equal(surf.is_valid, True)
    
    def test_exterior_yn(self,make_cmfd):
        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(0.0, 1.0)
        r = Vector(0.0, -1.5+tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        test_current = cmfd.current(G,16)

        check.equal(test_current, 1.0)
        check.equal(surf.is_valid, True)

    def test_exterior_xn(self,make_cmfd):
        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(1.0, 0.0)
        r = Vector(-1.5+tol, 0.0)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        test_current = cmfd.current(G,4)

        check.equal(test_current, 1.0)
        check.equal(surf.is_valid, True)

    def test_exterior_TR_TrueCorner(self,make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(-0.707,-0.707)
        r = Vector(1.5-tol,1.5-tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        tc_1 = cmfd.current(G,11)
        tc_2 = cmfd.current(G,23)

        check.equal(tc_1, -0.5)
        check.equal(tc_2, -0.5)
    
    def test_exterior_TR_XPCorner(self,make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(-0.707,-0.707)
        r = Vector(1.5-tol,0.5-tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        tc_1 = cmfd.current(G,11)
        tc_2 = cmfd.current(G,22)
        tc_3 = cmfd.current(G,7)

        check.equal(tc_1, -0.5)
        check.equal(tc_2, -0.5)
        check.equal(tc_3, -0.5)
    
    def test_exterior_TR_YPCorner(self,make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(-0.707,-0.707)
        r = Vector(0.5-tol,1.5-tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        tc_1 = cmfd.current(G,23)
        tc_2 = cmfd.current(G,10)
        tc_3 = cmfd.current(G,19)

        check.equal(tc_1, -0.5)
        check.equal(tc_2, -0.5)
        check.equal(tc_3, -0.5)
    
    def test_exterior_BR_TrueCorner(self,make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(-0.707,0.707)
        r = Vector(1.5-tol,-1.5+tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        tc_1 = cmfd.current(G,20)
        tc_2 = cmfd.current(G,3)

        check.equal(tc_1, 0.5)
        check.equal(tc_2, -0.5)
    
    def test_exterior_BR_YNCorner(self,make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(-0.707, 0.707)
        r = Vector(0.5-tol,-1.5+tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        tc_1 = cmfd.current(G,2)
        tc_2 = cmfd.current(G,20)
        tc_3 = cmfd.current(G,16)

        check.equal(tc_1, -0.5)
        check.equal(tc_2, 0.5)
        check.equal(tc_3, 0.5)
    
    def test_exterior_BR_XPCorner(self,make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(-0.707, 0.707)
        r = Vector(1.5-tol,-0.5+tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        tc_1 = cmfd.current(G,7)
        tc_2 = cmfd.current(G,21)
        tc_3 = cmfd.current(G,3)

        check.equal(tc_1, -0.5)
        check.equal(tc_2, 0.5)
        check.equal(tc_3, -0.5)
    
    def test_exterior_BL_TrueCorner(self,make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(0.707,0.707)
        r = Vector(-1.5+tol,-1.5+tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        tc_1 = cmfd.current(G,0)
        tc_2 = cmfd.current(G,12)

        check.equal(tc_1, 0.5)
        check.equal(tc_2, 0.5)

    def test_exterior_BL_YNCorner(self,make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(0.707,0.707)
        r = Vector(-0.5+tol,-1.5+tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        tc_1 = cmfd.current(G,1)
        tc_2 = cmfd.current(G,16)
        tc_3 = cmfd.current(G,12)

        check.equal(tc_1, 0.5)
        check.equal(tc_2, 0.5)
        check.equal(tc_3, 0.5)
    
    def test_exterior_BL_XNCorner(self,make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(0.707,0.707)
        r = Vector(-1.5+tol,-0.5+tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        tc_1 = cmfd.current(G,4)
        tc_2 = cmfd.current(G,13)
        tc_3 = cmfd.current(G,0)

        check.equal(tc_1, 0.5)
        check.equal(tc_2, 0.5)
        check.equal(tc_3, 0.5)
    
    def test_exterior_TL_TrueCorner(self,make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(0.707,-0.707)
        r = Vector(-1.5+tol,1.5-tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        tc_1 = cmfd.current(G,8)
        tc_2 = cmfd.current(G,15)

        check.equal(tc_1, 0.5)
        check.equal(tc_2, -0.5)
    
    def test_exterior_TL_YPCorner(self,make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(0.707,-0.707)
        r = Vector(-0.5+tol,1.5-tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        tc_1 = cmfd.current(G,19)
        tc_2 = cmfd.current(G,9)
        tc_3 = cmfd.current(G,15)

        check.equal(tc_1, -0.5)
        check.equal(tc_2, 0.5)
        check.equal(tc_3, -0.5)
    
    def test_exterior_TL_XNCorner(self,make_cmfd):

        cmfd = make_cmfd
        aflx = 1.0
        tol = 1e-12
        u = Direction(0.707,-0.707)
        r = Vector(-1.5+tol,0.5-tol)
        G = 0

        surf = cmfd.get_surface(r,u)
        cmfd.tally_current(aflx,u,G,surf)

        tc_1 = cmfd.current(G,14)
        tc_2 = cmfd.current(G,4)
        tc_3 = cmfd.current(G,8)

        check.equal(tc_1, -0.5)
        check.equal(tc_2, 0.5)
        check.equal(tc_3, 0.5)

    
