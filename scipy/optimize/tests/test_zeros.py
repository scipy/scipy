from __future__ import division, print_function, absolute_import

from math import sqrt, exp, sin, cos

from numpy.testing import (assert_warns, assert_, 
                           assert_allclose,
                           assert_equal)
import numpy as np

from scipy.optimize import zeros as cc
from scipy.optimize import zeros

# Import testing parameters
from scipy.optimize._tstutils import functions, fstrings


class TestBasic(object):
    def run_check(self, method, name):
        a = .5
        b = sqrt(3)
        xtol = 4*np.finfo(float).eps
        rtol = 4*np.finfo(float).eps
        for function, fname in zip(functions, fstrings):
            zero, r = method(function, a, b, xtol=xtol, rtol=rtol,
                             full_output=True)
            assert_(r.converged)
            assert_allclose(zero, 1.0, atol=xtol, rtol=rtol,
                err_msg='method %s, function %s' % (name, fname))

    def test_bisect(self):
        self.run_check(cc.bisect, 'bisect')

    def test_ridder(self):
        self.run_check(cc.ridder, 'ridder')

    def test_brentq(self):
        self.run_check(cc.brentq, 'brentq')

    def test_brenth(self):
        self.run_check(cc.brenth, 'brenth')

    def test_newton(self):
        f1 = lambda x: x**2 - 2*x - 1
        f1_1 = lambda x: 2*x - 2
        f1_2 = lambda x: 2.0 + 0*x

        f2 = lambda x: exp(x) - cos(x)
        f2_1 = lambda x: exp(x) + sin(x)
        f2_2 = lambda x: exp(x) + cos(x)

        for f, f_1, f_2 in [(f1, f1_1, f1_2), (f2, f2_1, f2_2)]:
            x = zeros.newton(f, 3, tol=1e-6)
            assert_allclose(f(x), 0, atol=1e-6)
            x = zeros.newton(f, 3, fprime=f_1, tol=1e-6)
            assert_allclose(f(x), 0, atol=1e-6)
            x = zeros.newton(f, 3, fprime=f_1, fprime2=f_2, tol=1e-6)
            assert_allclose(f(x), 0, atol=1e-6)

    def test_newton_array(self):
        """test newton with array"""

        def f1(x, *a):
            b = a[0] + x * a[3]
            return a[1] - a[2] * (np.exp(b / a[5]) - 1.0) - b / a[4] - x

        def f1_1(x, *a):
            b = a[3] / a[5]
            return -a[2] * np.exp(a[0] / a[5] + x * b) * b - a[3] / a[4] - 1

        def f1_2(x, *a):
            b = a[3] / a[5]
            return -a[2] * np.exp(a[0] / a[5] + x * b) * b ** 2

        a0 = np.array([
            5.32725221, 5.48673747, 5.49539973,
            5.36387202, 4.80237316, 1.43764452,
            5.23063958, 5.46094772, 5.50512718,
            5.42046290
        ])
        a1 = (np.sin(range(10)) + 1.0) * 7.0
        args = (a0, a1, 1e-09, 0.004, 10, 0.27456)
        x0 = 7.0
        x = zeros.newton(f1, x0, f1_1, args)
        x_expected = (
            6.17264965, 11.7702805, 12.2219954,
            7.11017681, 1.18151293, 0.143707955,
            4.31928228, 10.5419107, 12.7552490,
            8.91225749
        )
        assert_allclose(x, x_expected)
        # test halley's
        x = zeros.newton(f1, x0, f1_1, args, fprime2=f1_2)
        assert_allclose(x, x_expected)
        # test secant
        x = zeros.newton(lambda y, z: z - y ** 2, 4.0, args=([15.0, 17.0],))
        assert_allclose(x, (3.872983346207417, 4.123105625617661))
        # test derivative zero warning
        assert_warns(RuntimeWarning, zeros.newton,
                     lambda y: y ** 2, [0., 0.], lambda y: 2 * y)

    def test_deriv_zero_warning(self):
        func = lambda x: x**2
        dfunc = lambda x: 2*x
        assert_warns(RuntimeWarning, cc.newton, func, 0.0, dfunc)


def test_gh_5555():
    root = 0.1

    def f(x):
        return x - root

    methods = [cc.bisect, cc.ridder]
    xtol = 4*np.finfo(float).eps
    rtol = 4*np.finfo(float).eps
    for method in methods:
        res = method(f, -1e8, 1e7, xtol=xtol, rtol=rtol)
        assert_allclose(root, res, atol=xtol, rtol=rtol,
                        err_msg='method %s' % method.__name__)


def test_gh_5557():
    # Show that without the changes in 5557 brentq and brenth might
    # only achieve a tolerance of 2*(xtol + rtol*|res|).

    # f linearly interpolates (0, -0.1), (0.5, -0.1), and (1,
    # 0.4). The important parts are that |f(0)| < |f(1)| (so that
    # brent takes 0 as the initial guess), |f(0)| < atol (so that
    # brent accepts 0 as the root), and that the exact root of f lies
    # more than atol away from 0 (so that brent doesn't achieve the
    # desired tolerance).
    def f(x):
        if x < 0.5:
            return -0.1
        else:
            return x - 0.6

    atol = 0.51
    rtol = 4*np.finfo(float).eps
    methods = [cc.brentq, cc.brenth]
    for method in methods:
        res = method(f, 0, 1, xtol=atol, rtol=rtol)
        assert_allclose(0.6, res, atol=atol, rtol=rtol)


class TestRootResults:
    def test_repr(self):
        r = zeros.RootResults(root=1.0,
                              iterations=44,
                              function_calls=46,
                              flag=0)
        expected_repr = ("      converged: True\n           flag: 'converged'"
                         "\n function_calls: 46\n     iterations: 44\n"
                         "           root: 1.0")
        assert_equal(repr(r), expected_repr)

