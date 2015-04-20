from __future__ import division, print_function, absolute_import

import warnings

from numpy.testing import (assert_almost_equal, assert_array_equal,
        TestCase, run_module_suite, assert_allclose, assert_equal, assert_)
from scipy.interpolate import (KroghInterpolator, krogh_interpolate,
        BarycentricInterpolator, barycentric_interpolate,
        PiecewisePolynomial, piecewise_polynomial_interpolate,
        approximate_taylor_polynomial, pchip, PchipInterpolator)
from scipy._lib.six import xrange
import scipy
import numpy as np
from scipy.interpolate import splrep, splev


def check_shape(interpolator_cls, x_shape, y_shape, deriv_shape=None, axis=0):
    np.random.seed(1234)

    x = [-1, 0, 1]
    s = list(range(1, len(y_shape)+1))
    s.insert(axis % (len(y_shape)+1), 0)
    y = np.random.rand(*((3,) + y_shape)).transpose(s)

    # Cython code chokes on y.shape = (0, 3) etc, skip them
    if y.size == 0:
        return

    xi = np.zeros(x_shape)
    yi = interpolator_cls(x, y, axis=axis)(xi)

    target_shape = ((deriv_shape or ()) + y.shape[:axis]
                    + x_shape + y.shape[axis:][1:])
    assert_equal(yi.shape, target_shape)

    # check it works also with lists
    if x_shape and y.size > 0:
        interpolator_cls(list(x), list(y), axis=axis)(list(xi))

    # check also values
    if xi.size > 0 and deriv_shape is None:
        bs_shape = (y.shape[:axis] + ((1,)*len(x_shape)) + y.shape[axis:][1:])
        yv = y[((slice(None,None,None),)*(axis % y.ndim))+(1,)].reshape(bs_shape)

        yi, y = np.broadcast_arrays(yi, yv)
        assert_allclose(yi, y)

SHAPES = [(), (0,), (1,), (3,2,5)]


def test_shapes():
    for ip in [KroghInterpolator, BarycentricInterpolator, pchip]:
        for s1 in SHAPES:
            for s2 in SHAPES:
                for axis in range(-len(s2), len(s2)):
                    yield check_shape, ip, s1, s2, None, axis


def test_derivs_shapes():
    def krogh_derivs(x, y, axis=0):
        return KroghInterpolator(x, y, axis).derivatives

    for s1 in SHAPES:
        for s2 in SHAPES:
            for axis in range(-len(s2), len(s2)):
                yield check_shape, krogh_derivs, s1, s2, (3,), axis


def test_deriv_shapes():
    def krogh_deriv(x, y, axis=0):
        return KroghInterpolator(x, y, axis).derivative

    def pchip_deriv(x, y, axis=0):
        return pchip(x, y, axis).derivative()

    def pchip_deriv2(x, y, axis=0):
        return pchip(x, y, axis).derivative(2)

    def pchip_deriv_inplace(x, y, axis=0):
        class P(PchipInterpolator):
            def __call__(self, x):
                return PchipInterpolator.__call__(self, x, 1)
            pass
        return P(x, y, axis)

    for ip in [krogh_deriv, pchip_deriv, pchip_deriv2, pchip_deriv_inplace]:
        for s1 in SHAPES:
            for s2 in SHAPES:
                for axis in range(-len(s2), len(s2)):
                    yield check_shape, ip, s1, s2, (), axis


def _check_complex(ip):
    x = [1, 2, 3, 4]
    y = [1, 2, 1j, 3]
    p = ip(x, y)
    assert_allclose(y, p(x))


def test_complex():
    for ip in [KroghInterpolator, BarycentricInterpolator, pchip]:
        yield _check_complex, ip


class CheckKrogh(TestCase):
    def setUp(self):
        self.true_poly = scipy.poly1d([-2,3,1,5,-4])
        self.test_xs = np.linspace(-1,1,100)
        self.xs = np.linspace(-1,1,5)
        self.ys = self.true_poly(self.xs)

    def test_lagrange(self):
        P = KroghInterpolator(self.xs,self.ys)
        assert_almost_equal(self.true_poly(self.test_xs),P(self.test_xs))

    def test_scalar(self):
        P = KroghInterpolator(self.xs,self.ys)
        assert_almost_equal(self.true_poly(7),P(7))
        assert_almost_equal(self.true_poly(np.array(7)), P(np.array(7)))

    def test_derivatives(self):
        P = KroghInterpolator(self.xs,self.ys)
        D = P.derivatives(self.test_xs)
        for i in xrange(D.shape[0]):
            assert_almost_equal(self.true_poly.deriv(i)(self.test_xs),
                                D[i])

    def test_low_derivatives(self):
        P = KroghInterpolator(self.xs,self.ys)
        D = P.derivatives(self.test_xs,len(self.xs)+2)
        for i in xrange(D.shape[0]):
            assert_almost_equal(self.true_poly.deriv(i)(self.test_xs),
                                D[i])

    def test_derivative(self):
        P = KroghInterpolator(self.xs,self.ys)
        m = 10
        r = P.derivatives(self.test_xs,m)
        for i in xrange(m):
            assert_almost_equal(P.derivative(self.test_xs,i),r[i])

    def test_high_derivative(self):
        P = KroghInterpolator(self.xs,self.ys)
        for i in xrange(len(self.xs),2*len(self.xs)):
            assert_almost_equal(P.derivative(self.test_xs,i),
                                np.zeros(len(self.test_xs)))

    def test_hermite(self):
        xs = [0,0,0,1,1,1,2]
        ys = [self.true_poly(0),
              self.true_poly.deriv(1)(0),
              self.true_poly.deriv(2)(0),
              self.true_poly(1),
              self.true_poly.deriv(1)(1),
              self.true_poly.deriv(2)(1),
              self.true_poly(2)]
        P = KroghInterpolator(self.xs,self.ys)
        assert_almost_equal(self.true_poly(self.test_xs),P(self.test_xs))

    def test_vector(self):
        xs = [0, 1, 2]
        ys = np.array([[0,1],[1,0],[2,1]])
        P = KroghInterpolator(xs,ys)
        Pi = [KroghInterpolator(xs,ys[:,i]) for i in xrange(ys.shape[1])]
        test_xs = np.linspace(-1,3,100)
        assert_almost_equal(P(test_xs),
                np.rollaxis(np.asarray([p(test_xs) for p in Pi]),-1))
        assert_almost_equal(P.derivatives(test_xs),
                np.transpose(np.asarray([p.derivatives(test_xs) for p in Pi]),
                    (1,2,0)))

    def test_empty(self):
        P = KroghInterpolator(self.xs,self.ys)
        assert_array_equal(P([]), [])

    def test_shapes_scalarvalue(self):
        P = KroghInterpolator(self.xs,self.ys)
        assert_array_equal(np.shape(P(0)), ())
        assert_array_equal(np.shape(P(np.array(0))), ())
        assert_array_equal(np.shape(P([0])), (1,))
        assert_array_equal(np.shape(P([0,1])), (2,))

    def test_shapes_scalarvalue_derivative(self):
        P = KroghInterpolator(self.xs,self.ys)
        n = P.n
        assert_array_equal(np.shape(P.derivatives(0)), (n,))
        assert_array_equal(np.shape(P.derivatives(np.array(0))), (n,))
        assert_array_equal(np.shape(P.derivatives([0])), (n,1))
        assert_array_equal(np.shape(P.derivatives([0,1])), (n,2))

    def test_shapes_vectorvalue(self):
        P = KroghInterpolator(self.xs,np.outer(self.ys,np.arange(3)))
        assert_array_equal(np.shape(P(0)), (3,))
        assert_array_equal(np.shape(P([0])), (1,3))
        assert_array_equal(np.shape(P([0,1])), (2,3))

    def test_shapes_1d_vectorvalue(self):
        P = KroghInterpolator(self.xs,np.outer(self.ys,[1]))
        assert_array_equal(np.shape(P(0)), (1,))
        assert_array_equal(np.shape(P([0])), (1,1))
        assert_array_equal(np.shape(P([0,1])), (2,1))

    def test_shapes_vectorvalue_derivative(self):
        P = KroghInterpolator(self.xs,np.outer(self.ys,np.arange(3)))
        n = P.n
        assert_array_equal(np.shape(P.derivatives(0)), (n,3))
        assert_array_equal(np.shape(P.derivatives([0])), (n,1,3))
        assert_array_equal(np.shape(P.derivatives([0,1])), (n,2,3))

    def test_wrapper(self):
        P = KroghInterpolator(self.xs,self.ys)
        assert_almost_equal(P(self.test_xs),krogh_interpolate(self.xs,self.ys,self.test_xs))
        assert_almost_equal(P.derivative(self.test_xs,2),krogh_interpolate(self.xs,self.ys,self.test_xs,der=2))
        assert_almost_equal(P.derivatives(self.test_xs,2),krogh_interpolate(self.xs,self.ys,self.test_xs,der=[0,1]))

    def test_int_inputs(self):
        # Check input args are cast correctly to floats, gh-3669
        x = [0, 234,468,702,936,1170,1404,2340,3744,6084,8424,13104,60000]
        offset_cdf = np.array([-0.95, -0.86114777, -0.8147762, -0.64072425, -0.48002351,
                               -0.34925329, -0.26503107, -0.13148093, -0.12988833, -0.12979296,
                               -0.12973574, -0.08582937, 0.05])
        f = KroghInterpolator(x, offset_cdf)

        assert_allclose(abs((f(x) - offset_cdf) / f.derivative(x, 1)), 0, atol=1e-10)


class CheckTaylor(TestCase):
    def test_exponential(self):
        degree = 5
        p = approximate_taylor_polynomial(np.exp, 0, degree, 1, 15)
        for i in xrange(degree+1):
            assert_almost_equal(p(0),1)
            p = p.deriv()
        assert_almost_equal(p(0),0)


class CheckBarycentric(TestCase):
    def setUp(self):
        self.true_poly = scipy.poly1d([-2,3,1,5,-4])
        self.test_xs = np.linspace(-1,1,100)
        self.xs = np.linspace(-1,1,5)
        self.ys = self.true_poly(self.xs)

    def test_lagrange(self):
        P = BarycentricInterpolator(self.xs,self.ys)
        assert_almost_equal(self.true_poly(self.test_xs),P(self.test_xs))

    def test_scalar(self):
        P = BarycentricInterpolator(self.xs,self.ys)
        assert_almost_equal(self.true_poly(7),P(7))
        assert_almost_equal(self.true_poly(np.array(7)),P(np.array(7)))

    def test_delayed(self):
        P = BarycentricInterpolator(self.xs)
        P.set_yi(self.ys)
        assert_almost_equal(self.true_poly(self.test_xs),P(self.test_xs))

    def test_append(self):
        P = BarycentricInterpolator(self.xs[:3],self.ys[:3])
        P.add_xi(self.xs[3:],self.ys[3:])
        assert_almost_equal(self.true_poly(self.test_xs),P(self.test_xs))

    def test_vector(self):
        xs = [0, 1, 2]
        ys = np.array([[0,1],[1,0],[2,1]])
        P = BarycentricInterpolator(xs,ys)
        Pi = [BarycentricInterpolator(xs,ys[:,i]) for i in xrange(ys.shape[1])]
        test_xs = np.linspace(-1,3,100)
        assert_almost_equal(P(test_xs),
                np.rollaxis(np.asarray([p(test_xs) for p in Pi]),-1))

    def test_shapes_scalarvalue(self):
        P = BarycentricInterpolator(self.xs,self.ys)
        assert_array_equal(np.shape(P(0)), ())
        assert_array_equal(np.shape(P(np.array(0))), ())
        assert_array_equal(np.shape(P([0])), (1,))
        assert_array_equal(np.shape(P([0,1])), (2,))

    def test_shapes_vectorvalue(self):
        P = BarycentricInterpolator(self.xs,np.outer(self.ys,np.arange(3)))
        assert_array_equal(np.shape(P(0)), (3,))
        assert_array_equal(np.shape(P([0])), (1,3))
        assert_array_equal(np.shape(P([0,1])), (2,3))

    def test_shapes_1d_vectorvalue(self):
        P = BarycentricInterpolator(self.xs,np.outer(self.ys,[1]))
        assert_array_equal(np.shape(P(0)), (1,))
        assert_array_equal(np.shape(P([0])), (1,1))
        assert_array_equal(np.shape(P([0,1])), (2,1))

    def test_wrapper(self):
        P = BarycentricInterpolator(self.xs,self.ys)
        assert_almost_equal(P(self.test_xs),barycentric_interpolate(self.xs,self.ys,self.test_xs))


class CheckPiecewise(TestCase):
    def setUp(self):
        self.tck = splrep([0,1,2,3,4,5], [0,10,-1,3,7,2], s=0)
        self.test_xs = np.linspace(-1,6,100)
        self.spline_ys = splev(self.test_xs, self.tck)
        self.spline_yps = splev(self.test_xs, self.tck, der=1)
        self.xi = np.unique(self.tck[0])
        self.yi = [[splev(x, self.tck, der=j) for j in xrange(3)] for x in self.xi]

    def test_construction(self):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            P = PiecewisePolynomial(self.xi, self.yi, 3)

        assert_almost_equal(P(self.test_xs), self.spline_ys)

    def test_scalar(self):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            P = PiecewisePolynomial(self.xi,self.yi,3)

        assert_almost_equal(P(self.test_xs[0]),self.spline_ys[0])
        assert_almost_equal(P.derivative(self.test_xs[0],1),self.spline_yps[0])
        assert_almost_equal(P(np.array(self.test_xs[0])),self.spline_ys[0])
        assert_almost_equal(P.derivative(np.array(self.test_xs[0]),1),
                            self.spline_yps[0])

    def test_derivative(self):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            P = PiecewisePolynomial(self.xi,self.yi,3)

        assert_almost_equal(P.derivative(self.test_xs,1),self.spline_yps)

    def test_derivatives(self):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            P = PiecewisePolynomial(self.xi,self.yi,3)

        m = 4
        r = P.derivatives(self.test_xs,m)
        #print r.shape, r
        for i in xrange(m):
            assert_almost_equal(P.derivative(self.test_xs,i),r[i])

    def test_vector(self):
        xs = [0, 1, 2]
        ys = [[[0,1]],[[1,0],[-1,-1]],[[2,1]]]
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            P = PiecewisePolynomial(xs,ys)
            Pi = [PiecewisePolynomial(xs,[[yd[i] for yd in y] for y in ys])
                for i in xrange(len(ys[0][0]))]

        test_xs = np.linspace(-1,3,100)
        assert_almost_equal(P(test_xs),
                np.rollaxis(np.asarray([p(test_xs) for p in Pi]),-1))
        assert_almost_equal(P.derivative(test_xs,1),
                np.transpose(np.asarray([p.derivative(test_xs,1) for p in Pi]),
                    (1,0)))

    def test_incremental(self):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            P = PiecewisePolynomial([self.xi[0]], [self.yi[0]], 3)

        for i in xrange(1,len(self.xi)):
            P.append(self.xi[i],self.yi[i],3)
        assert_almost_equal(P(self.test_xs),self.spline_ys)

    def test_shapes_scalarvalue(self):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            P = PiecewisePolynomial(self.xi,self.yi,4)

        assert_array_equal(np.shape(P(0)), ())
        assert_array_equal(np.shape(P(np.array(0))), ())
        assert_array_equal(np.shape(P([0])), (1,))
        assert_array_equal(np.shape(P([0,1])), (2,))

    def test_shapes_scalarvalue_derivative(self):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            P = PiecewisePolynomial(self.xi,self.yi,4)

        n = 4
        assert_array_equal(np.shape(P.derivative(0,1)), ())
        assert_array_equal(np.shape(P.derivative(np.array(0),1)), ())
        assert_array_equal(np.shape(P.derivative([0],1)), (1,))
        assert_array_equal(np.shape(P.derivative([0,1],1)), (2,))

    def test_shapes_vectorvalue(self):
        yi = np.multiply.outer(np.asarray(self.yi),np.arange(3))
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            P = PiecewisePolynomial(self.xi,yi,4)

        assert_array_equal(np.shape(P(0)), (3,))
        assert_array_equal(np.shape(P([0])), (1,3))
        assert_array_equal(np.shape(P([0,1])), (2,3))

    def test_shapes_vectorvalue_1d(self):
        yi = np.multiply.outer(np.asarray(self.yi),np.arange(1))
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            P = PiecewisePolynomial(self.xi,yi,4)

        assert_array_equal(np.shape(P(0)), (1,))
        assert_array_equal(np.shape(P([0])), (1,1))
        assert_array_equal(np.shape(P([0,1])), (2,1))

    def test_shapes_vectorvalue_derivative(self):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            P = PiecewisePolynomial(self.xi, np.multiply.outer(self.yi,
                                                               np.arange(3)),4)

        n = 4
        assert_array_equal(np.shape(P.derivative(0,1)), (3,))
        assert_array_equal(np.shape(P.derivative([0],1)), (1,3))
        assert_array_equal(np.shape(P.derivative([0,1],1)), (2,3))

    def test_wrapper(self):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            P = PiecewisePolynomial(self.xi,self.yi)

        assert_almost_equal(P(self.test_xs),
                            piecewise_polynomial_interpolate(self.xi, self.yi,
                                                             self.test_xs))
        assert_almost_equal(P.derivative(self.test_xs,2),
                            piecewise_polynomial_interpolate(self.xi,
                                                             self.yi,
                                                             self.test_xs,
                                                             der=2))
        assert_almost_equal(P.derivatives(self.test_xs,2),
                            piecewise_polynomial_interpolate(self.xi,
                                                             self.yi,
                                                             self.test_xs,
                                                             der=[0,1]))


class TestPCHIP(TestCase):
    def _make_random(self, npts=20):
        np.random.seed(1234)
        xi = np.sort(np.random.random(npts))
        yi = np.random.random(npts)
        return pchip(xi, yi), xi, yi

    def test_overshoot(self):
        # PCHIP should not overshoot
        p, xi, yi = self._make_random()
        for i in range(len(xi)-1):
            x1, x2 = xi[i], xi[i+1]
            y1, y2 = yi[i], yi[i+1]
            if y1 > y2:
                y1, y2 = y2, y1
            xp = np.linspace(x1, x2, 10)
            yp = p(xp)
            assert_(((y1 <= yp) & (yp <= y2)).all())

    def test_monotone(self):
        # PCHIP should preserve monotonicty
        p, xi, yi = self._make_random()
        for i in range(len(xi)-1):
            x1, x2 = xi[i], xi[i+1]
            y1, y2 = yi[i], yi[i+1]
            xp = np.linspace(x1, x2, 10)
            yp = p(xp)
            assert_(((y2-y1) * (yp[1:] - yp[:1]) > 0).all())

    def test_cast(self):
        # regression test for integer input data, see gh-3453
        data = np.array([[0, 4, 12, 27, 47, 60, 79, 87, 99, 100],
                         [-33, -33, -19, -2, 12, 26, 38, 45, 53, 55]])
        xx = np.arange(100)
        curve = pchip(data[0], data[1])(xx)

        data1 = data * 1.0
        curve1 = pchip(data1[0], data1[1])(xx)

        assert_allclose(curve, curve1, atol=1e-14, rtol=1e-14)

if __name__ == '__main__':
    run_module_suite()
