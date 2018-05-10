from __future__ import division, print_function, absolute_import

from math import sqrt, exp, sin, cos

from numpy.testing import (assert_warns, assert_,
                           assert_allclose,
                           assert_equal)
from numpy import finfo

from scipy.optimize import zeros as cc
from scipy.optimize import zeros

# Import testing parameters
from scipy.optimize._tstutils import functions, fstrings

try:
    from numpy.testing import assert_raises_regex as _assert_raises_regex
    # assert_raises_regex first appeared in NumPy 1.9
except ImportError:
    from numpy.testing import assert_raises

    def _assert_raises_regex(exception_class, _, *args, **kwargs):
        assert_raises(exception_class, *args, **kwargs)


class TestBasic(object):
    def run_check(self, method, name):
        a = .5
        b = sqrt(3)
        xtol = 4*finfo(float).eps
        rtol = 4*finfo(float).eps
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

    def test_newton_full_output(self):
        # Test the full_output capability, both when converging and not.
        # Use the f2 family.
        f2 = lambda x: exp(x) - cos(x)
        f2_1 = lambda x: exp(x) + sin(x)
        f2_2 = lambda x: exp(x) + cos(x)

        x0 = 3
        expected_counts = [(11, 12), (8, 16), (5, 15)]

        x, r = zeros.newton(f2, x0, tol=1e-6, full_output=True, disp=False)
        assert_(r.converged)
        assert_equal(x, r.root)
        assert_equal((r.iterations, r.function_calls), expected_counts[0])
        assert(r.function_calls <= r.iterations + 1)

        # Now repeat, allowing one fewer iteration
        iters = r.iterations - 1
        x, r = zeros.newton(f2, x0, tol=1e-6, maxiter=iters,
                            full_output=True, disp=False)
        assert_(not r.converged)
        assert_equal(x, r.root)
        assert_equal(r.iterations, iters)
        # Check that the correct Exception is raised and
        # validate the start of the message.
        _assert_raises_regex(
            RuntimeError,
            'Failed to converge after %d iterations, value is .*' % (iters),
            zeros.newton, f2, x0, tol=1e-6, maxiter=iters,
            full_output=True, disp=True)

        x, r = zeros.newton(f2, x0, fprime=f2_1, tol=1e-6,
                            full_output=True, disp=False)
        assert_(r.converged)
        assert_equal(x, r.root)
        assert_equal((r.iterations, r.function_calls), expected_counts[1])
        assert_equal(r.function_calls, 2*r.iterations)
        iters = r.iterations - 1
        x, r = zeros.newton(f2, x0, fprime=f2_1, tol=1e-6,
                            maxiter=iters, full_output=True, disp=False)
        assert_(not r.converged)
        assert_equal(x, r.root)
        assert_equal(r.iterations, iters)

        x, r = zeros.newton(f2, x0, fprime=f2_1, fprime2=f2_2, tol=1e-6,
                            full_output=True, disp=False)
        assert_(r.converged)
        assert_equal(x, r.root)
        assert_equal((r.iterations, r.function_calls), expected_counts[2])
        assert_equal(r.function_calls, 3*r.iterations)
        iters = r.iterations - 1
        x, r = zeros.newton(f2, x0, fprime=f2_1, fprime2=f2_2, tol=1e-6,
                            maxiter=iters, full_output=True, disp=False)
        assert_(not r.converged)
        assert_equal(x, r.root)
        assert_equal(r.iterations, iters)

    def test_deriv_zero_warning(self):
        func = lambda x: x**2
        dfunc = lambda x: 2*x
        assert_warns(RuntimeWarning, cc.newton, func, 0.0, dfunc)


def test_gh_5555():
    root = 0.1

    def f(x):
        return x - root

    methods = [cc.bisect, cc.ridder]
    xtol = 4*finfo(float).eps
    rtol = 4*finfo(float).eps
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
    rtol = 4*finfo(float).eps
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


def test_complex_halley():
    """Test Halley's works with complex roots"""
    def f(x, *a):
        return a[0] * x**2 + a[1] * x + a[2]

    def f_1(x, *a):
        return 2 * a[0] * x + a[1]

    def f_2(x, *a):
        return 2 * a[0]

    z = complex(1.0, 2.0)
    coeffs = (2.0, 3.0, 4.0)
    y = zeros.newton(f, z, args=coeffs, fprime=f_1, fprime2=f_2, tol=1e-6)
    # (-0.75000000000000078+1.1989578808281789j)
    assert_allclose(f(y, *coeffs), 0, atol=1e-6)
