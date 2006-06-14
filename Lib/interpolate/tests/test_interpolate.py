from numpy.testing import *
from numpy import mgrid, pi, sin
set_package_path()
from interpolate import interp1d, interp2d
restore_path()

class test_interp2d(ScipyTestCase):
    def test_interp2d(self):
        x, y = mgrid[0:pi:20j, 0:pi:21j]
        z = sin(x*y)
        I = interp2d(x, y, z)
        assert_almost_equal(I(1.0, 1.0), sin(1.0), decimal=2)
# These don't work yet; you'll also have to figure out the accuracy
#        assert_almost_equal(I(x.ravel(), y.ravel()), z)
#        assert_almost_equal(I(x, y), z)

if __name__ == "__main__":
    ScipyTest().run()
