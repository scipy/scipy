"""
Unit tests for optimization routines from _root.py.
"""
from __future__ import division, print_function, absolute_import

from numpy.testing import assert_, assert_almost_equal, assert_raises
import numpy as np

from scipy.optimize import root, root_scalar


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

    def test_root_scalar(self):
        f = lambda x : x**3
        a, b = -1, 1
        # These should work
        for method in ['brentq', 'brenth', 'ridder', 'bisect']:
            x, r = root_scalar(f, a, b, method=method, full_output=True)
            msg = "Failed to find a root for %s" % method
            assert_(r.converged, msg)
            assert_(np.allclose(x, 0), msg)
        # Make sure it catches a bogus method
        assert_raises(ValueError, root_scalar, f, a, b, method=None)
