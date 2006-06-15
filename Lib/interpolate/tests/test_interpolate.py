from numpy.testing import *
from numpy import mgrid, pi, sin, ogrid
set_package_path()
from interpolate import interp1d, interp2d
restore_path()

class test_interp2d(ScipyTestCase):
    def test_interp2d(self):
        y, x = mgrid[0:pi:20j, 0:pi:21j]
        z = sin(x+y)
        I = interp2d(x, y, z)
        assert_almost_equal(I(1.0, 1.0), sin(2.0), decimal=2)

	v,u = ogrid[0:pi:24j, 0:pi:25j]
        assert_almost_equal(I(u.ravel(), v.ravel()), sin(v+u), decimal=2)

if __name__ == "__main__":
    ScipyTest().run()
