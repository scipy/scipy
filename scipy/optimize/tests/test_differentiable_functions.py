from __future__ import division, print_function, absolute_import
import numpy as np
from numpy.testing import (TestCase, assert_array_almost_equal,
                           assert_array_equal, assert_)
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator
from scipy.optimize._differentiable_functions import (ScalarFunction,
                                                      VectorFunction,
                                                      LinearVectorFunction,
                                                      IdentityVectorFunction)


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
        analit = ScalarFunction(ex.fun, x0, (), ex.grad, None, None)
        nfev += 1
        ngev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        approx = ScalarFunction(ex.fun, x0, (), '2-point', None, None)
        nfev += 3
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(analit.f, approx.f)
        assert_array_almost_equal(analit.g, approx.g)

        x = [10, 0.3]
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

        x = [2.0, 1.0]
        analit.update_x(x)
        g_analit = analit.evaluate_grad()
        ngev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        approx.update_x(x)
        g_approx = approx.evaluate_grad()
        nfev += 3
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_almost_equal(g_analit, g_approx)

        x = [2.5, 0.3]
        analit.update_x(x)
        f_analit = analit.evaluate_fun()
        g_analit = analit.evaluate_grad()
        nfev += 1
        ngev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        approx.update_x(x)
        f_approx = approx.evaluate_fun()
        g_approx = approx.evaluate_grad()
        nfev += 3
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_almost_equal(f_analit, f_approx)
        assert_array_almost_equal(g_analit, g_approx)

        x = [2, 0.3]
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

    def test_finite_difference_hess_linear_operator(self):
        ex = ExScalarFunction()
        nfev = 0
        ngev = 0
        nhev = 0

        x0 = [1.0, 0.0]
        analit = ScalarFunction(ex.fun, x0, (), ex.grad, ex.hess, None)
        nfev += 1
        ngev += 1
        nhev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        approx = ScalarFunction(ex.fun, x0, (), ex.grad, '2-point', None)
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
        analit.update_x(x)
        H_analit = analit.evaluate_hess()
        nhev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        approx.update_x(x)
        H_approx = approx.evaluate_hess()
        assert_(isinstance(H_approx, LinearOperator))
        for v in ([1.0, 2.0], [3.0, 4.0], [5.0, 2.0]):
            assert_array_almost_equal(H_analit.dot(v), H_approx.dot(v))
        ngev += 4
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)

        x = [2.1, 1.2]
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
        analit.update_x(x)
        _ = analit.evaluate_grad()
        H_analit = analit.evaluate_hess()
        ngev += 1
        nhev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        approx.update_x(x)
        _ = approx.evaluate_grad()
        H_approx = approx.evaluate_hess()
        assert_(isinstance(H_approx, LinearOperator))
        for v in ([1.0, 2.0], [3.0, 4.0], [5.0, 2.0]):
            assert_array_almost_equal(H_analit.dot(v), H_approx.dot(v))
        ngev += 4
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)

        x = [5.2, 2.3]
        _ = analit.grad(x)
        H_analit = analit.hess(x)
        ngev += 1
        nhev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)
        approx.update_x(x)
        _ = approx.grad(x)
        H_approx = approx.hess(x)
        assert_(isinstance(H_approx, LinearOperator))
        for v in ([1.0, 2.0], [3.0, 4.0], [5.0, 2.0]):
            assert_array_almost_equal(H_analit.dot(v), H_approx.dot(v))
        ngev += 4
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.ngev, ngev)
        assert_array_equal(ex.nhev, nhev)


class ExVectorialFunction:

    def __init__(self):
        self.nfev = 0
        self.njev = 0
        self.nhev = 0

    def fun(self, x):
        self.nfev += 1
        return np.array([2*(x[0]**2 + x[1]**2 - 1) - x[0],
                         4*(x[0]**3 + x[1]**2 - 4) - 3*x[0]])

    def jac(self, x):
        self.njev += 1
        return np.array([[4*x[0]-1, 4*x[1]],
                         [12*x[0]**2-3, 8*x[1]]])

    def hess(self, x, v):
        self.nhev += 1
        return v[0]*4*np.eye(2) + v[1]*np.array([[24*x[0], 0],
                                                 [0, 8]])

class TestVectorialFunction(TestCase):

    def test_finite_difference_jac(self):
        ex = ExVectorialFunction()
        nfev = 0
        njev = 0

        x0 = [1.0, 0.0]
        v0 = [0.0, 1.0]
        analit = VectorFunction(ex.fun, x0, ex.jac, None, None, None,
                                (-np.inf, np.inf), None)
        nfev += 1
        njev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        approx = VectorFunction(ex.fun, x0, '2-point', None, None, None,
                                (-np.inf, np.inf), None)
        nfev += 3
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        assert_array_equal(analit.f, approx.f)
        assert_array_almost_equal(analit.J, approx.J)

        x = [10, 0.3]
        f_analit = analit.fun(x)
        J_analit = analit.jac(x)
        nfev += 1
        njev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        f_approx = approx.fun(x)
        J_approx = approx.jac(x)
        nfev += 3
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        assert_array_almost_equal(f_analit, f_approx)
        assert_array_almost_equal(J_analit, J_approx, decimal=4)

        x = [2.0, 1.0]
        analit.update_x(x)
        J_analit = analit.evaluate_jac()
        njev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        approx.update_x(x)
        J_approx = approx.evaluate_jac()
        nfev += 3
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        assert_array_almost_equal(J_analit, J_approx)

        x = [2.5, 0.3]
        analit.update_x(x)
        f_analit = analit.evaluate_fun()
        J_analit = analit.evaluate_jac()
        nfev += 1
        njev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        approx.update_x(x)
        f_approx = approx.evaluate_fun()
        J_approx = approx.evaluate_jac()
        nfev += 3
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        assert_array_almost_equal(f_analit, f_approx)
        assert_array_almost_equal(J_analit, J_approx)

        x = [2, 0.3]
        f_analit = analit.fun(x)
        J_analit = analit.jac(x)
        nfev += 1
        njev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        f_approx = approx.fun(x)
        J_approx = approx.jac(x)
        nfev += 3
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        assert_array_almost_equal(f_analit, f_approx)
        assert_array_almost_equal(J_analit, J_approx)

    def test_finite_difference_hess_linear_operator(self):
        ex = ExVectorialFunction()
        nfev = 0
        njev = 0
        nhev = 0

        x0 = [1.0, 0.0]
        v0 = [1.0, 2.0]
        analit = VectorFunction(ex.fun, x0, ex.jac, ex.hess, None, None,
                                (-np.inf, np.inf), None)
        nfev += 1
        njev += 1
        nhev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        assert_array_equal(ex.nhev, nhev)
        approx = VectorFunction(ex.fun, x0, ex.jac, '2-point', None, None,
                                (-np.inf, np.inf), None)
        assert_(isinstance(approx.H, LinearOperator))
        for p in ([1.0, 2.0], [3.0, 4.0], [5.0, 2.0]):
            assert_array_equal(analit.f, approx.f)
            assert_array_almost_equal(analit.J, approx.J)
            assert_array_almost_equal(analit.H.dot(p), approx.H.dot(p))
        nfev += 1
        njev += 4
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        assert_array_equal(ex.nhev, nhev)

        x = [2.0, 1.0]
        analit.update_x(x)
        H_analit = analit.evaluate_hess()
        nhev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        assert_array_equal(ex.nhev, nhev)
        approx.update_x(x)
        H_approx = approx.evaluate_hess()
        assert_(isinstance(H_approx, LinearOperator))
        for p in ([1.0, 2.0], [3.0, 4.0], [5.0, 2.0]):
            assert_array_almost_equal(H_analit.dot(p), H_approx.dot(p))
        njev += 4
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        assert_array_equal(ex.nhev, nhev)

        x = [2.1, 1.2]
        v = [1.0, 1.0]
        H_analit = analit.hess(x, v)
        nhev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        assert_array_equal(ex.nhev, nhev)
        H_approx = approx.hess(x, v)
        assert_(isinstance(H_approx, LinearOperator))
        for v in ([1.0, 2.0], [3.0, 4.0], [5.0, 2.0]):
            assert_array_almost_equal(H_analit.dot(v), H_approx.dot(v))
        njev += 4
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        assert_array_equal(ex.nhev, nhev)

        x = [2.5, 0.3]
        analit.update_x(x)
        _ = analit.evaluate_jac()
        H_analit = analit.evaluate_hess()
        njev += 1
        nhev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        assert_array_equal(ex.nhev, nhev)
        approx.update_x(x)
        _ = approx.evaluate_jac()
        H_approx = approx.evaluate_hess()
        assert_(isinstance(H_approx, LinearOperator))
        for v in ([1.0, 2.0], [3.0, 4.0], [5.0, 2.0]):
            assert_array_almost_equal(H_analit.dot(v), H_approx.dot(v), decimal=4)
        njev += 4
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        assert_array_equal(ex.nhev, nhev)

        x = [5.2, 2.3]
        v = [2.3, 5.2]
        _ = analit.jac(x)
        H_analit = analit.hess(x, v)
        njev += 1
        nhev += 1
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        assert_array_equal(ex.nhev, nhev)
        approx.update_x(x)
        _ = approx.jac(x)
        H_approx = approx.hess(x, v)
        assert_(isinstance(H_approx, LinearOperator))
        for v in ([1.0, 2.0], [3.0, 4.0], [5.0, 2.0]):
            assert_array_almost_equal(H_analit.dot(v), H_approx.dot(v), decimal=4)
        njev += 4
        assert_array_equal(ex.nfev, nfev)
        assert_array_equal(ex.njev, njev)
        assert_array_equal(ex.nhev, nhev)


def test_LinearVectorFunction():
    A_dense = np.array([
        [-1, 2, 0],
        [0, 4, 2]
    ])
    A_sparse = csr_matrix(A_dense)
    x = np.array([1, -1, 0])
    v = np.array([-1, 1])
    Ax = np.array([-3, -4])

    f1 = LinearVectorFunction(A_dense, None)
    assert_(not f1.sparse_jacobian)

    f2 = LinearVectorFunction(A_dense, True)
    assert_(f2.sparse_jacobian)

    f3 = LinearVectorFunction(A_dense, False)
    assert_(not f3.sparse_jacobian)

    f4 = LinearVectorFunction(A_sparse, None)
    assert_(f4.sparse_jacobian)

    f5 = LinearVectorFunction(A_sparse, True)
    assert_(f5.sparse_jacobian)

    f6 = LinearVectorFunction(A_sparse, False)
    assert_(not f6.sparse_jacobian)

    assert_array_equal(f1.fun(x), Ax)
    assert_array_equal(f2.fun(x), Ax)
    assert_array_equal(f1.jac(x), A_dense)
    assert_array_equal(f2.jac(x).toarray(), A_sparse.toarray())
    assert_array_equal(f1.hess(x, v).toarray(), np.zeros((3, 3)))


def test_IdentityVectorFunction():
    f1 = IdentityVectorFunction(3, None)
    f2 = IdentityVectorFunction(3, False)
    f3 = IdentityVectorFunction(3, True)

    assert_(f1.sparse_jacobian)
    assert_(not f2.sparse_jacobian)
    assert_(f3.sparse_jacobian)

    x = np.array([-1, 2, 1])
    v = np.array([-2, 3, 0])

    assert_array_equal(f1.fun(x), x)
    assert_array_equal(f2.fun(x), x)

    assert_array_equal(f1.jac(x).toarray(), np.eye(3))
    assert_array_equal(f2.jac(x), np.eye(3))

    assert_array_equal(f1.hess(x, v).toarray(), np.zeros((3, 3)))
