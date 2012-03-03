#!/usr/bin/env python
# Created by Pearu Peterson, June 2003
""" Test functions for interpolate.fitpack2 module
"""
__usage__ = """
Build interpolate:
  python setup_interpolate.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.interpolate.test(<level>)'
Run tests if interpolate is not installed:
  python tests/test_fitpack.py [<level>]
"""
#import libwadpy

from numpy.testing import assert_equal, assert_almost_equal, assert_array_equal, \
        assert_array_almost_equal, assert_allclose, TestCase, run_module_suite
from numpy import array, diff, shape
from scipy.interpolate.fitpack2 import UnivariateSpline, LSQBivariateSpline, \
    SmoothBivariateSpline, RectBivariateSpline

class TestUnivariateSpline(TestCase):
    def test_linear_constant(self):
        x = [1,2,3]
        y = [3,3,3]
        lut = UnivariateSpline(x,y,k=1)
        assert_array_almost_equal(lut.get_knots(),[1,3])
        assert_array_almost_equal(lut.get_coeffs(),[3,3])
        assert_almost_equal(lut.get_residual(),0.0)
        assert_array_almost_equal(lut([1,1.5,2]),[3,3,3])

    def test_preserve_shape(self):
        x = [1, 2, 3]
        y = [0, 2, 4]
        lut = UnivariateSpline(x, y, k=1)
        arg = 2
        assert_equal(shape(arg), shape(lut(arg)))
        assert_equal(shape(arg), shape(lut(arg, nu=1)))
        arg = [1.5, 2, 2.5]
        assert_equal(shape(arg), shape(lut(arg)))
        assert_equal(shape(arg), shape(lut(arg, nu=1)))

    def test_linear_1d(self):
        x = [1,2,3]
        y = [0,2,4]
        lut = UnivariateSpline(x,y,k=1)
        assert_array_almost_equal(lut.get_knots(),[1,3])
        assert_array_almost_equal(lut.get_coeffs(),[0,4])
        assert_almost_equal(lut.get_residual(),0.0)
        assert_array_almost_equal(lut([1,1.5,2]),[0,1,2])

    def test_subclassing(self):
        # See #731

        class ZeroSpline(UnivariateSpline):
            def __call__(self, x):
                return 0*array(x)

        sp = ZeroSpline([1,2,3,4,5], [3,2,3,2,3], k=2)
        assert_array_equal(sp([1.5, 2.5]), [0., 0.])

    def test_empty_input(self):
        """Test whether empty input returns an empty output. Ticket 1014"""
        x = [1,3,5,7,9]
        y = [0,4,9,12,21]
        spl = UnivariateSpline(x, y, k=3)
        assert_array_equal(spl([]), array([]))

    def test_resize_regression(self):
        """Regression test for #1375."""
        x = [-1., -0.65016502, -0.58856235, -0.26903553, -0.17370892,
             -0.10011001, 0., 0.10011001, 0.17370892, 0.26903553, 0.58856235,
             0.65016502, 1.]
        y = [1.,0.62928599, 0.5797223, 0.39965815, 0.36322694, 0.3508061,
             0.35214793, 0.3508061, 0.36322694, 0.39965815, 0.5797223,
             0.62928599, 1.]
        w = [1.00000000e+12, 6.88875973e+02, 4.89314737e+02, 4.26864807e+02,
             6.07746770e+02, 4.51341444e+02, 3.17480210e+02, 4.51341444e+02,
             6.07746770e+02, 4.26864807e+02, 4.89314737e+02, 6.88875973e+02,
             1.00000000e+12]
        spl = UnivariateSpline(x=x, y=y, w=w, s=None)
        desired = array([ 0.35100374,  0.51715855,  0.87789547,  0.98719344])
        assert_allclose(spl([0.1, 0.5, 0.9, 0.99]), desired, atol=5e-4)


class TestLSQBivariateSpline(TestCase):
    def test_linear_constant(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [3,3,3,3,3,3,3,3,3]
        s = 0.1
        tx = [1+s,3-s]
        ty = [1+s,3-s]
        lut = LSQBivariateSpline(x,y,z,tx,ty,kx=1,ky=1)

        assert_almost_equal(lut(2,2), 3.)

    def test_bilinearity(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [0,7,8,3,4,7,1,3,4]
        s = 0.1
        tx = [1+s,3-s]
        ty = [1+s,3-s]
        lut = LSQBivariateSpline(x,y,z,tx,ty,kx=1,ky=1)

        tx, ty = lut.get_knots()

        for xa, xb in zip(tx[:-1], tx[1:]):
            for ya, yb in zip(ty[:-1], ty[1:]):
                for t in [0.1, 0.5, 0.9]:
                    for s in [0.3, 0.4, 0.7]:
                        xp = xa*(1-t) + xb*t
                        yp = ya*(1-s) + yb*s
                        zp = (+ lut(xa, ya)*(1-t)*(1-s)
                              + lut(xb, ya)*t*(1-s)
                              + lut(xa, yb)*(1-t)*s
                              + lut(xb, yb)*t*s)
                        assert_almost_equal(lut(xp,yp), zp)

    def test_integral(self):
        x = [1,1,1,2,2,2,8,8,8]
        y = [1,2,3,1,2,3,1,2,3]
        z = array([0,7,8,3,4,7,1,3,4])

        s = 0.1
        tx = [1+s,3-s]
        ty = [1+s,3-s]
        lut = LSQBivariateSpline(x,y,z,tx,ty,kx=1,ky=1)
        tx, ty = lut.get_knots()

        tz = lut(tx, ty)
        trpz = .25*(diff(tx)[:,None]*diff(ty)[None,:]
                    *(tz[:-1,:-1]+tz[1:,:-1]+tz[:-1,1:]+tz[1:,1:])).sum()

        assert_almost_equal(lut.integral(tx[0], tx[-1], ty[0], ty[-1]), trpz)

    def test_empty_input(self):
        """Test whether empty inputs returns an empty output. Ticket 1014"""
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [3,3,3,3,3,3,3,3,3]
        s = 0.1
        tx = [1+s,3-s]
        ty = [1+s,3-s]
        lut = LSQBivariateSpline(x,y,z,tx,ty,kx=1,ky=1)

        assert_array_equal(lut([], []), array([]))

class TestSmoothBivariateSpline(TestCase):
    def test_linear_constant(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [3,3,3,3,3,3,3,3,3]
        lut = SmoothBivariateSpline(x,y,z,kx=1,ky=1)
        assert_array_almost_equal(lut.get_knots(),([1,1,3,3],[1,1,3,3]))
        assert_array_almost_equal(lut.get_coeffs(),[3,3,3,3])
        assert_almost_equal(lut.get_residual(),0.0)
        assert_array_almost_equal(lut([1,1.5,2],[1,1.5]),[[3,3],[3,3],[3,3]])

    def test_linear_1d(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [0,0,0,2,2,2,4,4,4]
        lut = SmoothBivariateSpline(x,y,z,kx=1,ky=1)
        assert_array_almost_equal(lut.get_knots(),([1,1,3,3],[1,1,3,3]))
        assert_array_almost_equal(lut.get_coeffs(),[0,0,4,4])
        assert_almost_equal(lut.get_residual(),0.0)
        assert_array_almost_equal(lut([1,1.5,2],[1,1.5]),[[0,0],[1,1],[2,2]])

    def test_integral(self):
        x = [1,1,1,2,2,2,4,4,4]
        y = [1,2,3,1,2,3,1,2,3]
        z = array([0,7,8,3,4,7,1,3,4])

        lut = SmoothBivariateSpline(x,y,z,kx=1,ky=1,s=0)
        tx = [1,2,4]
        ty = [1,2,3]

        tz = lut(tx, ty)
        trpz = .25*(diff(tx)[:,None]*diff(ty)[None,:]
                    *(tz[:-1,:-1]+tz[1:,:-1]+tz[:-1,1:]+tz[1:,1:])).sum()
        assert_almost_equal(lut.integral(tx[0], tx[-1], ty[0], ty[-1]), trpz)

        lut2 = SmoothBivariateSpline(x,y,z,kx=2,ky=2,s=0)
        assert_almost_equal(lut2.integral(tx[0], tx[-1], ty[0], ty[-1]), trpz,
                            decimal=0) # the quadratures give 23.75 and 23.85

        tz = lut(tx[:-1], ty[:-1])
        trpz = .25*(diff(tx[:-1])[:,None]*diff(ty[:-1])[None,:]
                    *(tz[:-1,:-1]+tz[1:,:-1]+tz[:-1,1:]+tz[1:,1:])).sum()
        assert_almost_equal(lut.integral(tx[0], tx[-2], ty[0], ty[-2]), trpz)

class TestRectBivariateSpline(TestCase):
    def test_defaults(self):
        x = array([1,2,3,4,5])
        y = array([1,2,3,4,5])
        z = array([[1,2,1,2,1],[1,2,1,2,1],[1,2,3,2,1],[1,2,2,2,1],[1,2,1,2,1]])
        lut = RectBivariateSpline(x,y,z)
        assert_array_almost_equal(lut(x,y),z)

    def test_evaluate(self):
        x = array([1,2,3,4,5])
        y = array([1,2,3,4,5])
        z = array([[1,2,1,2,1],[1,2,1,2,1],[1,2,3,2,1],[1,2,2,2,1],[1,2,1,2,1]])
        lut = RectBivariateSpline(x,y,z)

        xi = [1, 2.3, 5.3, 0.5, 3.3, 1.2, 3]
        yi = [1, 3.3, 1.2, 4.0, 5.0, 1.0, 3]
        zi = lut.ev(xi, yi)
        zi2 = array([lut(xp, yp)[0,0] for xp, yp in zip(xi, yi)])

        assert_almost_equal(zi, zi2)

if __name__ == "__main__":
    run_module_suite()
