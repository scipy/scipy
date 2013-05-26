from __future__ import division, print_function, absolute_import

from numpy.testing import assert_, assert_equal, assert_almost_equal, \
        assert_array_almost_equal, assert_raises, assert_array_equal, \
        dec, TestCase, run_module_suite
from numpy import mgrid, pi, sin, ogrid, poly1d, linspace
import numpy as np

from scipy.interpolate import interp1d, interp2d, lagrange

from scipy.lib._gcutils import assert_deallocated


class TestInterp2D(TestCase):
    def test_interp2d(self):
        y, x = mgrid[0:2:20j, 0:pi:21j]
        z = sin(x+0.5*y)
        I = interp2d(x, y, z)
        assert_almost_equal(I(1.0, 2.0), sin(2.0), decimal=2)

        v,u = ogrid[0:2:24j, 0:pi:25j]
        assert_almost_equal(I(u.ravel(), v.ravel()), sin(u+0.5*v), decimal=2)

    def test_interp2d_meshgrid_input(self):
        # Ticket #703
        x = linspace(0, 2, 16)
        y = linspace(0, pi, 21)
        z = sin(x[None,:] + y[:,None]/2.)
        I = interp2d(x, y, z)
        assert_almost_equal(I(1.0, 2.0), sin(2.0), decimal=2)

    def test_interp2d_meshgrid_input_unsorted(self):
        np.random.seed(1234)
        x = linspace(0, 2, 16)
        y = linspace(0, pi, 21)

        z = sin(x[None,:] + y[:,None]/2.)
        ip1 = interp2d(x.copy(), y.copy(), z, kind='cubic')

        np.random.shuffle(x)
        z = sin(x[None,:] + y[:,None]/2.)
        ip2 = interp2d(x.copy(), y.copy(), z, kind='cubic')

        np.random.shuffle(x)
        np.random.shuffle(y)
        z = sin(x[None,:] + y[:,None]/2.)
        ip3 = interp2d(x, y, z, kind='cubic')

        x = linspace(0, 2, 31)
        y = linspace(0, pi, 30)

        assert_equal(ip1(x, y), ip2(x, y))
        assert_equal(ip1(x, y), ip3(x, y))

    def test_interp2d_linear(self):
        # Ticket #898
        a = np.zeros([5, 5])
        a[2, 2] = 1.0
        x = y = np.arange(5)
        b = interp2d(x, y, a, 'linear')
        assert_almost_equal(b(2.0, 1.5), np.array([0.5]), decimal=2)
        assert_almost_equal(b(2.0, 2.5), np.array([0.5]), decimal=2)

    def test_interp2d_bounds(self):
        x = np.linspace(0, 1, 5)
        y = np.linspace(0, 2, 7)
        z = x[:,None]**2 + y[None,:]

        ix = np.linspace(-1, 3, 31)
        iy = np.linspace(-1, 3, 33)

        b = interp2d(x, y, z, bounds_error=True)
        assert_raises(ValueError, b, ix, iy)

        b = interp2d(x, y, z, fill_value=np.nan)
        iz = b(ix, iy)
        mx = (ix < 0) | (ix > 1)
        my = (iy < 0) | (iy > 2)
        assert_(np.isnan(iz[my,:]).all())
        assert_(np.isnan(iz[:,mx]).all())
        assert_(np.isfinite(iz[~my,:][:,~mx]).all())


class TestInterp1D(object):

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
        interp1d(self.x10, self.y10, kind='nearest')
        interp1d(self.x10, self.y10, kind=0)
        interp1d(self.x10, self.y10, kind=1)
        interp1d(self.x10, self.y10, kind=2)
        interp1d(self.x10, self.y10, kind=3)

        # x array must be 1D.
        assert_raises(ValueError, interp1d, self.x25, self.y10)

        # y array cannot be a scalar.
        assert_raises(ValueError, interp1d, self.x10, np.array(0))

        # Check for x and y arrays having the same length.
        assert_raises(ValueError, interp1d, self.x10, self.y2)
        assert_raises(ValueError, interp1d, self.x2, self.y10)
        assert_raises(ValueError, interp1d, self.x10, self.y102)
        interp1d(self.x10, self.y210)
        interp1d(self.x10, self.y102, axis=0)

        # Check for x and y having at least 1 element.
        assert_raises(ValueError, interp1d, self.x1, self.y10)
        assert_raises(ValueError, interp1d, self.x10, self.y1)
        assert_raises(ValueError, interp1d, self.x1, self.y1)

    def test_init(self):
        """ Check that the attributes are initialized appropriately by the
        constructor.
        """

        assert_(interp1d(self.x10, self.y10).copy)
        assert_(not interp1d(self.x10, self.y10, copy=False).copy)
        assert_(interp1d(self.x10, self.y10).bounds_error)
        assert_(not interp1d(self.x10, self.y10, bounds_error=False).bounds_error)
        assert_(np.isnan(interp1d(self.x10, self.y10).fill_value))
        assert_equal(
            interp1d(self.x10, self.y10, fill_value=3.0).fill_value,
            3.0,
        )
        assert_equal(
            interp1d(self.x10, self.y10).axis,
            0,
        )
        assert_equal(
            interp1d(self.x10, self.y210).axis,
            1,
        )
        assert_equal(
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

    def test_nearest(self):
        """Check the actual implementation of nearest-neighbour interpolation.
        """

        interp10 = interp1d(self.x10, self.y10, kind='nearest')
        assert_array_almost_equal(
            interp10(self.x10),
            self.y10,
        )
        assert_array_almost_equal(
            interp10(1.2),
            np.array(1.),
        )
        assert_array_almost_equal(
            interp10([2.4, 5.6, 6.0]),
            np.array([2., 6., 6.]),
        )

    @dec.knownfailureif(True, "zero-order splines fail for the last point")
    def test_zero(self):
        """Check the actual implementation of zero-order spline interpolation.
        """
        interp10 = interp1d(self.x10, self.y10, kind='zero')
        assert_array_almost_equal(interp10(self.x10), self.y10)
        assert_array_almost_equal(interp10(1.2), np.array(1.))
        assert_array_almost_equal(interp10([2.4, 5.6, 6.0]),
                                  np.array([2., 6., 6.]))

    def _bounds_check(self, kind='linear'):
        """ Test that our handling of out-of-bounds input is correct.
        """

        extrap10 = interp1d(self.x10, self.y10, fill_value=self.fill_value,
            bounds_error=False, kind=kind)
        assert_array_equal(
            extrap10(11.2),
            np.array(self.fill_value),
        )
        assert_array_equal(
            extrap10(-3.4),
            np.array(self.fill_value),
        )
        assert_array_equal(
            extrap10([[[11.2], [-3.4], [12.6], [19.3]]]),
            np.array(self.fill_value),
        )
        assert_array_equal(
            extrap10._check_bounds(np.array([-1.0, 0.0, 5.0, 9.0, 11.0])),
            np.array([True, False, False, False, True]),
        )

        raises_bounds_error = interp1d(self.x10, self.y10, bounds_error=True,
                                       kind=kind)
        assert_raises(ValueError, raises_bounds_error, -1.0)
        assert_raises(ValueError, raises_bounds_error, 11.0)
        raises_bounds_error([0.0, 5.0, 9.0])

    def _bounds_check_int_nan_fill(self, kind='linear'):
        x = np.arange(10).astype(np.int_)
        y = np.arange(10).astype(np.int_)
        c = interp1d(x, y, kind=kind, fill_value=np.nan, bounds_error=False)
        yi = c(x - 1)
        assert_(np.isnan(yi[0]))
        assert_array_almost_equal(yi, np.r_[np.nan, y[:-1]])

    def test_bounds(self):
        for kind in ('linear', 'cubic', 'nearest',
                     'slinear', 'zero', 'quadratic'):
            self._bounds_check(kind)
            self._bounds_check_int_nan_fill(kind)

    def _nd_check_interp(self, kind='linear'):
        """Check the behavior when the inputs and outputs are multidimensional.
        """

        # Multidimensional input.
        interp10 = interp1d(self.x10, self.y10, kind=kind)
        assert_array_almost_equal(
            interp10(np.array([[3., 5.], [2., 7.]])),
            np.array([[3., 5.], [2., 7.]]),
        )

        # Scalar input -> 0-dim scalar array output
        assert_(isinstance(interp10(1.2), np.ndarray))
        assert_equal(interp10(1.2).shape, ())

        # Multidimensional outputs.
        interp210 = interp1d(self.x10, self.y210, kind=kind)
        assert_array_almost_equal(
            interp210(1.),
            np.array([1., 11.]),
        )
        assert_array_almost_equal(
            interp210(np.array([1., 2.])),
            np.array([[1., 2.],
                      [11., 12.]]),
        )

        interp102 = interp1d(self.x10, self.y102, axis=0, kind=kind)
        assert_array_almost_equal(
            interp102(1.),
            np.array([2.0, 3.0]),
        )
        assert_array_almost_equal(
            interp102(np.array([1., 3.])),
            np.array([[2., 3.],
                      [6., 7.]]),
        )

        # Both at the same time!
        x_new = np.array([[3., 5.], [2., 7.]])
        assert_array_almost_equal(
            interp210(x_new),
            np.array([[[3., 5.], [2., 7.]],
                      [[13., 15.], [12., 17.]]]),
        )
        assert_array_almost_equal(
            interp102(x_new),
            np.array([[[6., 7.], [10., 11.]],
                      [[4., 5.], [14., 15.]]]),
        )

    def _nd_check_shape(self, kind='linear'):
        # Check large ndim output shape
        a = [4, 5, 6, 7]
        y = np.arange(np.prod(a)).reshape(*a)
        for n, s in enumerate(a):
            x = np.arange(s)
            z = interp1d(x, y, axis=n, kind=kind)
            assert_array_almost_equal(z(x), y, err_msg=kind)

            x2 = np.arange(2*3*1).reshape((2,3,1)) / 12.
            b = list(a)
            b[n:n+1] = [2,3,1]
            assert_array_almost_equal(z(x2).shape, b, err_msg=kind)

    def test_nd(self):
        for kind in ('linear', 'cubic', 'slinear', 'quadratic', 'nearest'):
            self._nd_check_interp(kind)
            self._nd_check_shape(kind)

    def _check_complex(self, dtype=np.complex_, kind='linear'):
        x = np.array([1, 2.5, 3, 3.1, 4, 6.4, 7.9, 8.0, 9.5, 10])
        y = x * x ** (1 + 2j)
        y = y.astype(dtype)

        # simple test
        c = interp1d(x, y, kind=kind)
        assert_array_almost_equal(y[:-1], c(x)[:-1])

        # check against interpolating real+imag separately
        xi = np.linspace(1, 10, 31)
        cr = interp1d(x, y.real, kind=kind)
        ci = interp1d(x, y.imag, kind=kind)
        assert_array_almost_equal(c(xi).real, cr(xi))
        assert_array_almost_equal(c(xi).imag, ci(xi))

    def test_complex(self):
        for kind in ('linear', 'nearest', 'cubic', 'slinear', 'quadratic',
                     'zero'):
            self._check_complex(np.complex64, kind)
            self._check_complex(np.complex128, kind)

    @dec.knownfailureif(True, "zero-order splines fail for the last point")
    def test_nd_zero_spline(self):
        # zero-order splines don't get the last point right,
        # see test_zero above
        #yield self._nd_check_interp, 'zero'
        #yield self._nd_check_interp, 'zero'
        pass

    def test_circular_refs(self):
        # Test interp1d can be automatically garbage collected
        x = np.linspace(0, 1)
        y = np.linspace(0, 1)
        # Confirm interp can be released from memory after use
        with assert_deallocated(interp1d, x, y) as interp:
            new_y = interp([0.1, 0.2])
            del interp


class TestLagrange(TestCase):

    def test_lagrange(self):
        p = poly1d([5,2,1,4,3])
        xs = np.arange(len(p.coeffs))
        ys = p(xs)
        pl = lagrange(xs,ys)
        assert_array_almost_equal(p.coeffs,pl.coeffs)

if __name__ == "__main__":
    run_module_suite()
