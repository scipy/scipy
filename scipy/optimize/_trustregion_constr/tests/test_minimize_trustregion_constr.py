from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.linalg import block_diag
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import LinearOperator
from numpy.testing import (TestCase, assert_array_almost_equal,
                           assert_array_equal, assert_array_less,
                           assert_raises, assert_equal, assert_,
                           run_module_suite, assert_allclose, assert_warns,
                           dec)
from scipy.optimize._trustregion_constr.minimize_trustregion_constr import ScalarFunction


class ExScalarFunction:

    def __init__(self):
        self.nfev = 0
        self.ngev = 0
        self.nhev = 0

    def fun(self, x):
        self.nfev += 1
        return 2*(x[0]**2 + x[1]**2 - 1) - x[0]

    def grad(self, x):
        self.ngev += 1
        return np.array([4*x[0]-1, 4*x[1]])

    def hess(self, x):
        self.nhev += 1
        return 4*np.eye(2)


class TestScalarFunction(TestCase):

    def test_finite_difference_grad(self):
        ex = ExScalarFunction()
        nfev = 0
        ngev = 0

        x0 = [1.0, 0.0]
        analit = ScalarFunction(ex.fun, x0, (), ex.grad, None)
        nfev += 1
        ngev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        approx = ScalarFunction(ex.fun, x0, (), '2-point', None)
        nfev += 3
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(analit.f, approx.f)
        assert_array_almost_equal(analit.g, approx.g)

        x = [2.0, 1.0]
        g_analit = analit.grad(x)
        ngev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        g_approx = approx.grad(x)
        nfev += 3
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_almost_equal(g_analit, g_approx)

        x = [2.5, 0.3]
        f_analit = analit.fun(x)
        g_analit = analit.grad(x)
        nfev += 1
        ngev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        f_approx = approx.fun(x)
        g_approx = approx.grad(x)
        nfev += 3
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_almost_equal(f_analit, f_approx)
        assert_array_almost_equal(g_analit, g_approx)

    def test_finite_difference_hess(self):
        ex = ExScalarFunction()
        nfev = 0
        ngev = 0
        nhev = 0

        x0 = [1.0, 0.0]
        analit = ScalarFunction(ex.fun, x0, (), ex.grad, ex.hess)
        nfev += 1
        ngev += 1
        nhev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        approx = ScalarFunction(ex.fun, x0, (), ex.grad, '2-point')
        nfev += 1
        ngev += 3
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        assert_array_equal(analit.f, approx.f)
        assert_array_almost_equal(analit.g, approx.g)
        assert_array_almost_equal(analit.H, approx.H)

        x = [2.0, 1.0]
        H_analit = analit.hess(x)
        nhev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        H_approx = approx.hess(x)
        ngev += 3
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        assert_array_almost_equal(H_analit, H_approx)

        x = [2.5, 0.3]
        g_analit = analit.grad(x)
        H_analit = analit.hess(x)
        ngev += 1
        nhev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        g_approx = approx.grad(x)
        H_approx = approx.hess(x)
        ngev += 3
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        assert_array_almost_equal(g_analit, g_approx)
        assert_array_almost_equal(H_analit, H_approx)

    def test_finite_difference_hess(self):
        ex = ExScalarFunction()
        nfev = 0
        ngev = 0
        nhev = 0

        x0 = [1.0, 0.0]
        analit = ScalarFunction(ex.fun, x0, (), ex.grad, ex.hess)
        nfev += 1
        ngev += 1
        nhev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        approx = ScalarFunction(ex.fun, x0, (), ex.grad, '2-point')
        nfev += 1
        ngev += 3
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        assert_array_equal(analit.f, approx.f)
        assert_array_almost_equal(analit.g, approx.g)
        assert_array_almost_equal(analit.H, approx.H)

        x = [2.0, 1.0]
        H_analit = analit.hess(x)
        nhev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        H_approx = approx.hess(x)
        ngev += 3
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        assert_array_almost_equal(H_analit, H_approx)

        x = [2.5, 0.3]
        g_analit = analit.grad(x)
        H_analit = analit.hess(x)
        ngev += 1
        nhev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        g_approx = approx.grad(x)
        H_approx = approx.hess(x)
        ngev += 3
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        assert_array_almost_equal(g_analit, g_approx)
        assert_array_almost_equal(H_analit, H_approx)

    def test_finite_difference_hess_linear_operator(self):
        ex = ExScalarFunction()
        nfev = 0
        ngev = 0
        nhev = 0
        finite_diff_options = {"as_linear_operator": True}

        x0 = [1.0, 0.0]
        analit = ScalarFunction(ex.fun, x0, (), ex.grad, ex.hess)
        nfev += 1
        ngev += 1
        nhev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        approx = ScalarFunction(ex.fun, x0, (), ex.grad, '2-point', finite_diff_options)
        assert_(isinstance(approx.H, LinearOperator))
        for v in ([1.0, 2.0], [3.0, 4.0], [5.0, 2.0]):
            assert_array_equal(analit.f, approx.f)
            assert_array_almost_equal(analit.g, approx.g)
            assert_array_almost_equal(analit.H.dot(v), approx.H.dot(v))
        nfev += 1
        ngev += 4
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)

        x = [2.0, 1.0]
        H_analit = analit.hess(x)
        nhev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        H_approx = approx.hess(x)
        assert_(isinstance(H_approx, LinearOperator))
        for v in ([1.0, 2.0], [3.0, 4.0], [5.0, 2.0]):
            assert_array_almost_equal(H_analit.dot(v), H_approx.dot(v))
        ngev += 4
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)

        x = [2.5, 0.3]
        g_analit = analit.grad(x)
        H_analit = analit.hess(x)
        ngev += 1
        nhev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        g_approx = approx.grad(x)
        H_approx = approx.hess(x)
        assert_(isinstance(H_approx, LinearOperator))
        for v in ([1.0, 2.0], [3.0, 4.0], [5.0, 2.0]):
            assert_array_almost_equal(H_analit.dot(v), H_approx.dot(v))
        ngev += 4
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
