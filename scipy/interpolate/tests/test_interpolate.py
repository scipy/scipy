from numpy.testing import *
from numpy import mgrid, pi, sin, ogrid, poly1d
import numpy as np

from scipy.interpolate import interp1d, interp2d, lagrange


class TestInterp2D(TestCase):
    def test_interp2d(self):
        y, x = mgrid[0:pi:20j, 0:pi:21j]
        z = sin(x+y)
        I = interp2d(x, y, z)
        assert_almost_equal(I(1.0, 1.0), sin(2.0), decimal=2)

        v,u = ogrid[0:pi:24j, 0:pi:25j]
        assert_almost_equal(I(u.ravel(), v.ravel()), sin(v+u), decimal=2)


class TestInterp1D(TestCase):

    def setUp(self):
        self.x10 = np.arange(10.)
        self.y10 = np.arange(10.)
        self.x25 = self.x10.reshape((2,5))
        self.x2 = np.arange(2.)
        self.y2 = np.arange(2.)
        self.x1 = np.array([0.])
        self.y1 = np.array([0.])

        self.y210 = np.arange(20.).reshape((2, 10))
        self.y102 = np.arange(20.).reshape((10, 2))

        self.fill_value = -100.0

    def test_validation(self):
        """ Make sure that appropriate exceptions are raised when invalid values
        are given to the constructor.
        """

        # These should all work.
        interp1d(self.x10, self.y10, kind='linear')
        interp1d(self.x10, self.y10, kind='cubic')
        interp1d(self.x10, self.y10, kind='slinear')
        interp1d(self.x10, self.y10, kind='quadratic')
        interp1d(self.x10, self.y10, kind='zero')
        interp1d(self.x10, self.y10, kind=0)
        interp1d(self.x10, self.y10, kind=1)
        interp1d(self.x10, self.y10, kind=2)
        interp1d(self.x10, self.y10, kind=3)

        # x array must be 1D.
        self.assertRaises(ValueError, interp1d, self.x25, self.y10)

        # y array cannot be a scalar.
        self.assertRaises(ValueError, interp1d, self.x10, np.array(0))

        # Check for x and y arrays having the same length.
        self.assertRaises(ValueError, interp1d, self.x10, self.y2)
        self.assertRaises(ValueError, interp1d, self.x2, self.y10)
        self.assertRaises(ValueError, interp1d, self.x10, self.y102)
        interp1d(self.x10, self.y210)
        interp1d(self.x10, self.y102, axis=0)

        # Check for x and y having at least 1 element.
        self.assertRaises(ValueError, interp1d, self.x1, self.y10)
        self.assertRaises(ValueError, interp1d, self.x10, self.y1)
        self.assertRaises(ValueError, interp1d, self.x1, self.y1)


    def test_init(self):
        """ Check that the attributes are initialized appropriately by the
        constructor.
        """

        self.assert_(interp1d(self.x10, self.y10).copy)
        self.assert_(not interp1d(self.x10, self.y10, copy=False).copy)
        self.assert_(interp1d(self.x10, self.y10).bounds_error)
        self.assert_(not interp1d(self.x10, self.y10, bounds_error=False).bounds_error)
        self.assert_(np.isnan(interp1d(self.x10, self.y10).fill_value))
        self.assertEqual(
            interp1d(self.x10, self.y10, fill_value=3.0).fill_value,
            3.0,
        )
        self.assertEqual(
            interp1d(self.x10, self.y10).axis,
            0,
        )
        self.assertEqual(
            interp1d(self.x10, self.y210).axis,
            1,
        )
        self.assertEqual(
            interp1d(self.x10, self.y102, axis=0).axis,
            0,
        )
        assert_array_equal(
            interp1d(self.x10, self.y10).x,
            self.x10,
        )
        assert_array_equal(
            interp1d(self.x10, self.y10).y,
            self.y10,
        )
        assert_array_equal(
            interp1d(self.x10, self.y210).y,
            self.y210,
        )


    def test_linear(self):
        """ Check the actual implementation of linear interpolation.
        """

        interp10 = interp1d(self.x10, self.y10)
        assert_array_almost_equal(
            interp10(self.x10),
            self.y10,
        )
        assert_array_almost_equal(
            interp10(1.2),
            np.array([1.2]),
        )
        assert_array_almost_equal(
            interp10([2.4, 5.6, 6.0]),
            np.array([2.4, 5.6, 6.0]),
        )

    def test_cubic(self):
        """ Check the actual implementation of spline interpolation.
        """

        interp10 = interp1d(self.x10, self.y10, kind='cubic')
        assert_array_almost_equal(
            interp10(self.x10),
            self.y10,
        )
        assert_array_almost_equal(
            interp10(1.2),
            np.array([1.2]),
        )
        assert_array_almost_equal(
            interp10([2.4, 5.6, 6.0]),
            np.array([2.4, 5.6, 6.0]),
        )

    def _bounds_check(self, kind='linear'):
        """ Test that our handling of out-of-bounds input is correct.
        """

        extrap10 = interp1d(self.x10, self.y10, fill_value=self.fill_value,
            bounds_error=False, kind=kind)
        assert_array_equal(
            extrap10(11.2),
            np.array([self.fill_value]),
        )
        assert_array_equal(
            extrap10(-3.4),
            np.array([self.fill_value]),
        )
        assert_array_equal(
            extrap10._check_bounds(np.array([-1.0, 0.0, 5.0, 9.0, 11.0])),
            np.array([True, False, False, False, True]),
        )

        raises_bounds_error = interp1d(self.x10, self.y10, bounds_error=True,
                                       kind=kind)
        self.assertRaises(ValueError, raises_bounds_error, -1.0)
        self.assertRaises(ValueError, raises_bounds_error, 11.0)
        raises_bounds_error([0.0, 5.0, 9.0])

    def test_bounds(self):
        for kind in ('linear', 'cubic'):
            self._bounds_check(kind=kind)

    def _nd_check(self, kind='linear'):
        """ Check the behavior when the inputs and outputs are multidimensional.
        """
        # Multidimensional input.
        interp10 = interp1d(self.x10, self.y10, kind=kind)
        assert_array_almost_equal(
            interp10(np.array([[3.4, 5.6], [2.4, 7.8]])),
            np.array([[3.4, 5.6], [2.4, 7.8]]),
        )

        # Multidimensional outputs.
        interp210 = interp1d(self.x10, self.y210, kind=kind)
        assert_array_almost_equal(
            interp210(1.5),
            np.array([[1.5], [11.5]]),
        )
        assert_array_almost_equal(
            interp210(np.array([1.5, 2.4])),
            np.array([[1.5, 2.4],
                      [11.5, 12.4]]),
        )

        interp102 = interp1d(self.x10, self.y102, axis=0, kind=kind)
        assert_array_almost_equal(
            interp102(1.5),
            np.array([[3.0, 4.0]]),
        )
        assert_array_almost_equal(
            interp102(np.array([1.5, 2.4])),
            np.array([[3.0, 4.0],
                      [4.8, 5.8]]),
        )

        # Both at the same time!
        x_new = np.array([[3.4, 5.6], [2.4, 7.8]])
        assert_array_almost_equal(
            interp210(x_new),
            np.array([[[3.4, 5.6], [2.4, 7.8]],
                      [[13.4, 15.6], [12.4, 17.8]]]),
        )
        assert_array_almost_equal(
            interp102(x_new),
            np.array([[[6.8, 7.8], [11.2, 12.2]],
                      [[4.8, 5.8], [15.6, 16.6]]]),
        )

        # Check large ndim output
        a = [4, 5, 6, 7]
        y = np.arange(np.prod(a)).reshape(*a)
        for n, s in enumerate(a):
            x = np.arange(s)
            z = interp1d(x, y, axis=n, kind=kind)
            assert_array_almost_equal(z(x), y)

            x2 = np.arange(2*3*1).reshape((2,3,1)) / 12.
            b = list(a)
            b[n:n+1] = [2,3,1]
            assert_array_almost_equal(z(x2).shape, b)

    def test_nd(self):
        for kind in ('linear', 'cubic'):
            self._nd_check(kind=kind)

class TestLagrange(TestCase):

    def test_lagrange(self):
        p = poly1d([5,2,1,4,3])
        xs = np.arange(len(p.coeffs))
        ys = p(xs)
        pl = lagrange(xs,ys)
        assert_array_almost_equal(p.coeffs,pl.coeffs)

if __name__ == "__main__":
    run_module_suite()
