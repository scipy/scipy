from __future__ import division, print_function, absolute_import

import warnings

import numpy as np

from numpy.testing import (
    assert_almost_equal, assert_array_equal, TestCase, run_module_suite,
    assert_allclose, assert_equal, assert_, assert_raises)

from scipy.interpolate import (
    KroghInterpolator, krogh_interpolate,
    BarycentricInterpolator, barycentric_interpolate,
    approximate_taylor_polynomial, pchip, PchipInterpolator,
    Akima1DInterpolator, CubicSpline)
from scipy._lib.six import xrange


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
        bs_shape = y.shape[:axis] + (1,)*len(x_shape) + y.shape[axis:][1:]
        yv = y[((slice(None,),)*(axis % y.ndim)) + (1,)]
        yv = yv.reshape(bs_shape)

        yi, y = np.broadcast_arrays(yi, yv)
        assert_allclose(yi, y)

SHAPES = [(), (0,), (1,), (3, 2, 5)]


def test_shapes():
    for ip in [KroghInterpolator, BarycentricInterpolator, pchip,
               Akima1DInterpolator, CubicSpline]:
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

    def pchip_antideriv(x, y, axis=0):
        return pchip(x, y, axis).derivative()

    def pchip_antideriv2(x, y, axis=0):
        return pchip(x, y, axis).derivative(2)

    def pchip_deriv_inplace(x, y, axis=0):
        class P(PchipInterpolator):
            def __call__(self, x):
                return PchipInterpolator.__call__(self, x, 1)
            pass
        return P(x, y, axis)

    def akima_deriv(x, y, axis=0):
        return Akima1DInterpolator(x, y, axis).derivative()

    def akima_antideriv(x, y, axis=0):
        return Akima1DInterpolator(x, y, axis).antiderivative()

    def cspline_deriv(x, y, axis=0):
        return CubicSpline(x, y, axis).derivative()

    def cspline_antideriv(x, y, axis=0):
        return CubicSpline(x, y, axis).antiderivative()

    for ip in [krogh_deriv, pchip_deriv, pchip_deriv2, pchip_deriv_inplace,
               pchip_antideriv, pchip_antideriv2, akima_deriv, akima_antideriv,
               cspline_deriv, cspline_antideriv]:
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
    for ip in [KroghInterpolator, BarycentricInterpolator, pchip, CubicSpline]:
        yield _check_complex, ip


class CheckKrogh(TestCase):
    def setUp(self):
        self.true_poly = np.poly1d([-2,3,1,5,-4])
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
        self.true_poly = np.poly1d([-2,3,1,5,-4])
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

    def test_nag(self):
        # Example from NAG C implementation,
        # http://nag.com/numeric/cl/nagdoc_cl25/html/e01/e01bec.html
        # suggested in gh-5326 as a smoke test for the way the derivatives
        # are computed (see also gh-3453)
        from scipy._lib.six import StringIO
        dataStr = '''
          7.99   0.00000E+0
          8.09   0.27643E-4
          8.19   0.43750E-1
          8.70   0.16918E+0
          9.20   0.46943E+0
         10.00   0.94374E+0
         12.00   0.99864E+0
         15.00   0.99992E+0
         20.00   0.99999E+0
        '''
        data = np.loadtxt(StringIO(dataStr))
        pch = pchip(data[:,0], data[:,1])

        resultStr = '''
           7.9900       0.0000
           9.1910       0.4640
          10.3920       0.9645
          11.5930       0.9965
          12.7940       0.9992
          13.9950       0.9998
          15.1960       0.9999
          16.3970       1.0000
          17.5980       1.0000
          18.7990       1.0000
          20.0000       1.0000
        '''
        result = np.loadtxt(StringIO(resultStr))
        assert_allclose(result[:,1], pch(result[:,0]), rtol=0., atol=5e-5)

    def test_endslopes(self):
        # this is a smoke test for gh-3453: PCHIP interpolator should not
        # set edge slopes to zero if the data do not suggest zero edge derivatives
        x = np.array([0.0, 0.1, 0.25, 0.35])
        y1 = np.array([279.35, 0.5e3, 1.0e3, 2.5e3])
        y2 = np.array([279.35, 2.5e3, 1.50e3, 1.0e3])
        for pp in (pchip(x, y1), pchip(x, y2)):
            for t in (x[0], x[-1]):
                assert_(pp(t, 1) != 0)

    def test_all_zeros(self):
        x = np.arange(10)
        y = np.zeros_like(x)

        # this should work and not generate any warnings
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            pch = pchip(x, y)

        xx = np.linspace(0, 9, 101)
        assert_equal(pch(xx), 0.)


class TestCubicSpline(object):
    def test_general(self):
        # This test can be reproduced in Octave/MATLAB with the following code:
        # x = [-1, 0, 0.5, 2, 4, 4.5, 5.5, 9];
        # y = [0, -0.5, 2, 3, 2.5, 1, 1, 0.5];
        # x_eval = [-0.5, 0.1, 1, 3, 4.1, 5, 7];
        # y_true = spline(x, y, x_eval);

        x = np.array([-1, 0, 0.5, 2, 4, 4.5, 5.5, 9])
        y = np.array([0, -0.5, 2, 3, 2.5, 1, 1, 0.5])

        x_test = np.array([-0.5, 0.1, 1, 3, 4.1, 5, 7])
        y_true = np.array([
            -2.111801336869092, 0.019677035288960, 3.164826064645384,
            3.384966086626085, 2.210195906837064, 0.464510284120531,
            4.956649059572205])

        cs = CubicSpline(x, y)
        assert_allclose(cs(x_test), y_true, rtol=1e-12)

        Y = np.vstack((y, y)).T
        cs = CubicSpline(x, Y)
        Y_true = np.vstack((y_true, y_true)).T
        assert_allclose(cs(x_test), Y_true, rtol=1e-12)

        Y = np.vstack((y, y))
        cs = CubicSpline(x, Y, axis=1)
        Y_true = np.vstack((y_true, y_true))
        assert_allclose(cs(x_test), Y_true, rtol=1e-12)

        Y = np.empty((2, 2, 8))
        Y[0, 0, :] = y
        Y[0, 1, :] = y
        Y[1, 0, :] = y
        Y[1, 1, :] = y
        cs = CubicSpline(x, Y, axis=-1)
        Y_true = np.empty((2, 2, 7))
        Y_true[0, 0, :] = y_true
        Y_true[0, 1, :] = y_true
        Y_true[1, 0, :] = y_true
        Y_true[1, 1, :] = y_true
        assert_allclose(cs(x_test), Y_true, rtol=1e-12)

    def test_corner_cases(self):
        x = np.array([-1, 0, 1])
        y = np.array([2, 0, 1])
        x_test = np.array([-0.5, 0.5])
        y_true = np.array([5/8, 1/8])
        cs = CubicSpline(x, y)
        assert_allclose(cs(x_test), y_true, rtol=1e-12)

        x = np.array([-1, 1])
        y = np.array([-1, 3])
        y_true = np.array([0, 2])
        cs = CubicSpline(x, y)
        assert_allclose(cs(x_test), y_true, rtol=1e-12)

    def _check_continuity(self, S, alpha=1e-8):
        x = S.x
        dx = -alpha * np.diff(x)
        for i in range(3):
            D = S.derivative()
            t = S(x[1:]) + D(x[1:]) * dx
            v = S(x[1:] + dx)
            assert_allclose(t, v, rtol=alpha*1e3)
            S = D

    def test_dtypes(self):
        x = np.array([0, 1, 2, 3], dtype=int)
        y = np.array([-5, 2, 3, 1], dtype=int)
        cs = CubicSpline(x, y)
        assert_allclose(y, cs(x), rtol=1e-12)
        self._check_continuity(cs)

        y = np.array([-1+1j, 0.0, 1-1j, 0.5-1.5j])
        cs = CubicSpline(x, y)
        assert_allclose(y, cs(x), rtol=1e-12)
        self._check_continuity(cs)

    def test_small_dx(self):
        rng = np.random.RandomState(0)
        x = np.sort(rng.uniform(size=100))
        y = 1e4 + rng.uniform(size=100)
        cs = CubicSpline(x, y)
        assert_allclose(cs(x), y, rtol=1e-12)
        self._check_continuity(cs)

    def test_incorrect_inputs(self):
        x = np.array([1, 2, 3])
        y = np.array([1, 2, 3])
        xc = np.array([1 + 1j, 2, 3])
        xn = np.array([np.nan, 2, 3])
        xo = np.array([2, 1, 3])
        yn = np.array([np.nan, 2, 3])

        assert_raises(ValueError, CubicSpline, xc, y)
        assert_raises(ValueError, CubicSpline, xn, y)
        assert_raises(ValueError, CubicSpline, x, yn)
        assert_raises(ValueError, CubicSpline, xo, y)


if __name__ == '__main__':
    run_module_suite()
