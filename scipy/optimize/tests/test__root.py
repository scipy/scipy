"""
Unit tests for optimization routines from _root.py.
"""
from __future__ import division, print_function, absolute_import

from numpy.testing import assert_, assert_allclose, dec
import numpy as np

from scipy.optimize import root, fixpoint
from scipy.special import lambertw


class TestRoot(object):
    def test_tol_parameter(self):
        # Check that the minimize() tol= argument does something
        def func(z):
            x, y = z
            return np.array([x**3 - 1, y**3 - 1])

        def dfunc(z):
            x, y = z
            return np.array([[3*x**2, 0], [0, 3*y**2]])

        for method in ['hybr', 'lm', 'broyden1', 'broyden2', 'anderson',
                       'diagbroyden', 'krylov']:
            if method in ('linearmixing', 'excitingmixing'):
                # doesn't converge
                continue

            if method in ('hybr', 'lm'):
                jac = dfunc
            else:
                jac = None

            sol1 = root(func, [1.1,1.1], jac=jac, tol=1e-4, method=method)
            sol2 = root(func, [1.1,1.1], jac=jac, tol=0.5, method=method)
            msg = "%s: %s vs. %s" % (method, func(sol1.x), func(sol2.x))
            assert_(sol1.success, msg)
            assert_(sol2.success, msg)
            assert_(abs(func(sol1.x)).max() < abs(func(sol2.x)).max(),
                    msg)

    def test_minimize_scalar_coerce_args_param(self):
        # github issue #3503
        def func(z, f=1):
            x, y = z
            return np.array([x**3 - 1, y**3 - f])
        root(func, [1.1, 1.1], args=1.5)


class TestFixPoint(object):

    def _check_scalar_trivial(self, method):
        # f(x) = 2x; fixed point should be x=0
        def func(x):
            return 2.0*x
        x0 = 1.0
        sol = fixpoint(func, x0, method=method)
        assert_(sol.success)
        assert_allclose(sol.x, 0.0, atol=1e-15)

    def _check_scalar_basic1(self, method):
        # f(x) = x**2; x0=1.05; fixed point should be x=1
        def func(x):
            return x**2
        x0 = 1.05
        sol = fixpoint(func, x0, method=method)
        assert_(sol.success)
        assert_allclose(sol.x, 1.0)

    def _check_scalar_basic2(self, method):
        # f(x) = x**0.5; x0=1.05; fixed point should be x=1 or x=0
        def func(x):
            return x**0.5
        x0 = 1.05
        sol = fixpoint(func, x0, method=method)
        assert_(sol.success)
        assert_allclose(sol.x, 1.0)

    def _check_array_trivial(self, method):
        def func(x):
            return 2.0*x
        x0 = [0.3, 0.15]
        olderr = np.seterr(all='ignore')
        try:
            sol = fixpoint(func, x0, method=method)
        finally:
            np.seterr(**olderr)
        assert_(sol.success)
        assert_allclose(sol.x, [0.0, 0.0], atol=1e-15)

    def _check_array_basic1(self, method):
        # f(x) = c * x**2; fixed point should be x=1/c or x=0
        opts = dict()

        def func(x, c):
            return c * x**2

        c = np.array([0.75, 1.0, 1.25])
        x0 = [1.1, 1.15, 0.9]
        olderr = np.seterr(all='ignore')
        try:
            sol = fixpoint(func, x0, args=(c,), method=method, options=opts)
        finally:
            np.seterr(**olderr)
        if not sol.success:
            print(sol)
        assert_(sol.success)
        assert_(np.all((abs(sol.x - 1.0/c) < 1e-6) | (abs(sol.x) < 1e-6)))

    def _check_array_basic2(self, method):
        # f(x) = c * x**0.5; fixed point should be x=c**2 or x=0
        def func(x, c):
            return c * x**0.5
        c = np.array([0.75, 1.0, 1.25])
        x0 = [0.8, 1.1, 1.1]
        sol = fixpoint(func, x0, args=(c,), method=method)
        assert_(sol.success)
        assert_(np.all((abs(sol.x - c**2) < 1e-6) | (abs(sol.x) < 1e-6)))

    def _check_lambertw(self, method):
        # python-list/2010-December/594592.html
        if method == 'steffensen':
            opts = dict(maxiter=500)
        elif method == 'squarem':
            opts = dict(maxfev=500)
        sol = fixpoint(lambda xx: np.exp(-2.0*xx)/2.0, 1.0, tol=1e-12,
                       method=method, args=(), options=opts)
        assert_(sol.success)
        assert_allclose(sol.x, np.exp(-2.0*sol.x)/2.0)
        assert_allclose(sol.x, lambertw(1)/2)

    def test_methods(self):
        for method in ['steffensen', 'squarem']:
            if method != 'squarem':
                # No convergence for squarem --- these functions are
                # not contractions
                yield self._check_scalar_trivial, method
                yield self._check_scalar_basic1, method
                yield self._check_array_trivial, method
                yield self._check_array_basic1, method

            yield self._check_scalar_basic2, method
            yield self._check_array_basic2, method
            yield self._check_lambertw, method
