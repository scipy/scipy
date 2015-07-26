from itertools import product
import warnings

import numpy as np
from numpy.testing import (run_module_suite, assert_, assert_allclose,
                           assert_raises, assert_equal)

from scipy.optimize import least_squares
from scipy.sparse import issparse, lil_matrix
from scipy.sparse.linalg import aslinearoperator
from scipy.optimize._lsq_common import (
    step_size_to_bound, find_active_constraints, make_strictly_feasible,
    intersect_trust_region, build_quadratic_1d, minimize_quadratic_1d,
    evaluate_quadratic)
from scipy.optimize._lsq_trf import scaling_vector


# Tests of functions from ._lsq_common.


class TestBounds(object):
    def test_step_size_to_bounds(self):
        lb = np.array([-1.0, 2.5, 10.0])
        ub = np.array([1.0, 5.0, 100.0])
        x = np.array([0.0, 2.5, 12.0])

        s = np.array([0.1, 0.0, 0.0])
        step, hits = step_size_to_bound(x, s, lb, ub)
        assert_equal(step, 10)
        assert_equal(hits, [1, 0, 0])

        s = np.array([0.01, 0.05, -1.0])
        step, hits = step_size_to_bound(x, s, lb, ub)
        assert_equal(step, 2)
        assert_equal(hits, [0, 0, -1])

        s = np.array([10.0, -0.0001, 100.0])
        step, hits = step_size_to_bound(x, s, lb, ub)
        assert_equal(step, np.array(-0))
        assert_equal(hits, [0, -1, 0])

        s = np.array([1.0, 0.5, -2.0])
        step, hits = step_size_to_bound(x, s, lb, ub)
        assert_equal(step, 1.0)
        assert_equal(hits, [1, 0, -1])

        s = np.zeros(3)
        step, hits = step_size_to_bound(x, s, lb, ub)
        assert_equal(step, np.inf)
        assert_equal(hits, [0, 0, 0])

    def test_find_active_constraints(self):
        lb = np.array([0.0, -10.0, 1.0])
        ub = np.array([1.0, 0.0, 100.0])

        x = np.array([0.5, -5.0, 2.0])
        active = find_active_constraints(x, lb, ub)
        assert_equal(active, [0, 0, 0])

        x = np.array([0.0, 0.0, 10.0])
        active = find_active_constraints(x, lb, ub)
        assert_equal(active, [-1, 1, 0])

        active = find_active_constraints(x, lb, ub, rtol=0)
        assert_equal(active, [-1, 1, 0])

        x = np.array([1e-9, -1e-8, 100 - 1e-9])
        active = find_active_constraints(x, lb, ub)
        assert_equal(active, [0, 0, 1])

        active = find_active_constraints(x, lb, ub, rtol=1.5e-9)
        assert_equal(active, [-1, 0, 1])

        lb = np.array([1.0, -np.inf, -np.inf])
        ub = np.array([np.inf, 10.0, np.inf])

        x = np.ones(3)
        active = find_active_constraints(x, lb, ub)
        assert_equal(active, [-1, 0, 0])

        # Handles out-of-bound cases.
        x = np.array([0.0, 11.0, 0.0])
        active = find_active_constraints(x, lb, ub)
        assert_equal(active, [-1, 1, 0])

        active = find_active_constraints(x, lb, ub, rtol=0)
        assert_equal(active, [-1, 1, 0])

    def test_make_strictly_feasible(self):
        lb = np.array([-0.5, -0.8, 2.0])
        ub = np.array([0.8, 1.0, 3.0])

        x = np.array([-0.5, 0.0, 2 + 1e-10])

        x_new = make_strictly_feasible(x, lb, ub, rstep=0)
        assert_(x_new[0] > -0.5)
        assert_equal(x_new[1:], x[1:])

        x_new = make_strictly_feasible(x, lb, ub, rstep=1e-4)
        assert_equal(x_new, [-0.5 + 1e-4, 0.0, 2 * (1 + 1e-4)])

        x = np.array([-0.5, -1, 3.1])
        x_new = make_strictly_feasible(x, lb, ub)
        assert_(np.all((x_new >= lb) & (x_new <= ub)))

        x_new = make_strictly_feasible(x, lb, ub, rstep=0)
        assert_(np.all((x_new >= lb) & (x_new <= ub)))

        lb = np.array([-1, 100.0])
        ub = np.array([1, 100.0 + 1e-10])
        x = np.array([0, 100.0])
        x_new = make_strictly_feasible(x, lb, ub, rstep=1e-8)
        assert_equal(x_new, [0, 100.0 + 0.5e-10])

    def test_scaling_vector(self):
        lb = np.array([-np.inf, -5.0, 1.0, -np.inf])
        ub = np.array([1.0, np.inf, 10.0, np.inf])
        x = np.array([0.5, 2.0, 5.0, 0.0])
        g = np.array([1.0, 0.1, -10.0, 0.0])
        v, dv = scaling_vector(x, g, lb, ub)
        assert_equal(v, [1.0, 7.0, 5.0, 1.0])
        assert_equal(dv, [0.0, 1.0, -1.0, 0.0])


class TestQuadraticFunction(object):
    def __init__(self):
        self.J = np.array([
            [0.1, 0.2],
            [-1.0, 1.0],
            [0.5, 0.2]])
        self.g = np.array([0.8, -2.0])
        self.diag = np.array([1.0, 2.0])

    def test_build_quadratic_1d(self):
        s = np.zeros(2)
        a, b = build_quadratic_1d(self.J, self.g, s)
        assert_equal(a, 0)
        assert_equal(b, 0)

        a, b = build_quadratic_1d(self.J, self.g, s, diag=self.diag)
        assert_equal(a, 0)
        assert_equal(b, 0)

        s = np.array([1.0, -1.0])
        a, b = build_quadratic_1d(self.J, self.g, s)
        assert_equal(a, 2.05)
        assert_equal(b, 2.8)

        a, b = build_quadratic_1d(self.J, self.g, s, diag=self.diag)
        assert_equal(a, 3.55)
        assert_equal(b, 2.8)

        s0 = np.array([0.5, 0.5])
        a, b = build_quadratic_1d(self.J, self.g, s, diag=self.diag, s0=s0)
        assert_equal(a, 3.55)
        assert_allclose(b, 2.39)

    def test_minimize_quadratic_1d(self):
        a = 5
        b = -1

        t, y = minimize_quadratic_1d(a, b, 1, 2)
        assert_equal(t, 1)
        assert_equal(y, a * t**2 + b * t)

        t, y = minimize_quadratic_1d(a, b, -2, -1)
        assert_equal(t, -1)
        assert_equal(y, a * t**2 + b * t)

        t, y = minimize_quadratic_1d(a, b, -1, 1)
        assert_equal(t, 0.1)
        assert_equal(y, a * t**2 + b * t)

    def test_evaluate_quadratic(self):
        s = np.array([1.0, -1.0])

        value = evaluate_quadratic(self.J, self.g, s)
        assert_equal(value, 4.85)

        value = evaluate_quadratic(self.J, self.g, s, diag=self.diag)
        assert_equal(value, 6.35)

        s = np.array([[1.0, -1.0],
                     [1.0, 1.0],
                     [0.0, 0.0]])

        values = evaluate_quadratic(self.J, self.g, s)
        assert_allclose(values, [4.85, -0.91, 0.0])

        values = evaluate_quadratic(self.J, self.g, s, diag=self.diag)
        assert_allclose(values, [6.35, 0.59, 0.0])


class TestTrustRegion(object):
    def test_intersect(self):
        Delta = 1.0

        x = np.zeros(3)
        s = np.array([1.0, 0.0, 0.0])
        t_neg, t_pos = intersect_trust_region(x, s, Delta)
        assert_equal(t_neg, -1)
        assert_equal(t_pos, 1)

        s = np.array([-1.0, 1.0, -1.0])
        t_neg, t_pos = intersect_trust_region(x, s, Delta)
        assert_allclose(t_neg, -3**-0.5)
        assert_allclose(t_pos, 3**-0.5)

        x = np.array([0.5, -0.5, 0])
        s = np.array([0, 0, 1.0])
        t_neg, t_pos = intersect_trust_region(x, s, Delta)
        assert_allclose(t_neg, -2**-0.5)
        assert_allclose(t_pos, 2**-0.5)

        x = np.ones(3)
        assert_raises(ValueError, intersect_trust_region, x, s, Delta)

        x = np.zeros(3)
        s = np.zeros(3)
        assert_raises(ValueError, intersect_trust_region, x, s, Delta)


# Tests of least_squares.


def fun_trivial(x, a=0):
    return (x - a)**2 + 5.0


def jac_trivial(x, a=0.0):
    return 2 * (x - a)


def fun_2d_trivial(x):
    return np.array([x[0], x[1]])


def jac_2d_trivial(x):
    return np.identity(2)


def fun_rosenbrock(x):
    return np.array([10 * (x[1] - x[0]**2), (1 - x[0])])


def jac_rosenbrock(x):
    return np.array([
        [-20 * x[0], 10],
        [-1, 0]
    ])


def jac_rosenbrock_bad_dim(x):
    return np.array([
        [-20 * x[0], 10],
        [-1, 0],
        [0.0, 0.0]
    ])


def fun_rosenbrock_cropped(x):
    return fun_rosenbrock(x)[0]


def jac_rosenbrock_cropped(x):
    return jac_rosenbrock(x)[0]


# When x is 1-d array, return is 2-d array.
def fun_wrong_dimensions(x):
    return np.array([x, x**2, x**3])


def jac_wrong_dimensions(x, a=0.0):
    return np.atleast_3d(jac_trivial(x, a=a))


class BroydenTridiagonal(object):
    def __init__(self, n=100, mode='sparse'):
        np.random.seed(0)

        self.n = n

        self.x0 = -np.ones(n)
        self.lb = np.linspace(-2, -1.5, n)
        self.ub = np.linspace(-0.8, 0.0, n)

        self.lb += 0.1 * np.random.randn(n)
        self.ub += 0.1 * np.random.randn(n)

        self.x0 += 0.1 * np.random.randn(n)
        self.x0 = make_strictly_feasible(self.x0, self.lb, self.ub)

        if mode == 'sparse':
            self.sparsity = lil_matrix((n, n), dtype=int)
            i = np.arange(n)
            self.sparsity[i, i] = 1
            i = np.arange(1, n)
            self.sparsity[i, i - 1] = 1
            i = np.arange(n - 1)
            self.sparsity[i, i + 1] = 1

            self.jac = self._jac
        elif mode == 'operator':
            self.jac = lambda x: aslinearoperator(self._jac(x))
        elif mode == 'dense':
            self.sparsity = None
            self.jac = lambda x: self._jac(x).toarray()
        else:
            assert_(False)

    def fun(self, x):
        f = (3 - x) * x + 1
        f[1:] -= x[:-1]
        f[:-1] -= 2 * x[1:]
        return f

    def _jac(self, x):
        J = lil_matrix((self.n, self.n))
        i = np.arange(self.n)
        J[i, i] = 3 - 2 * x
        i = np.arange(1, self.n)
        J[i, i - 1] = -1
        i = np.arange(self.n - 1)
        J[i, i + 1] = -2
        return J


class BaseMixin(object):
    def test_basic(self):
        # Test that the basic calling sequence works.
        res = least_squares(fun_trivial, 2., method=self.method)
        assert_allclose(res.x, 0, atol=1e-4)
        assert_allclose(res.fun, fun_trivial(res.x))

    def test_args_kwargs(self):
        # Test that args and kwargs are passed correctly to the functions.
        a = 3.0
        for jac in ['2-point', '3-point', 'cs', jac_trivial]:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                res = least_squares(fun_trivial, 2.0, jac, args=(a,),
                                    method=self.method)
                assert_allclose(res.x, a, rtol=1e-4)

                assert_raises(TypeError, least_squares, fun_trivial, 2.0,
                              args=(3, 4,), method=self.method)

                res = least_squares(fun_trivial, 2.0, jac, kwargs={'a': a},
                                    method=self.method)
                assert_allclose(res.x, a, rtol=1e-4)
                assert_raises(TypeError, least_squares, fun_trivial, 2.0,
                              kwargs={'kaboom': 3}, method=self.method)

    def test_jac_options(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            for jac in ['2-point', '3-point', 'cs', jac_trivial]:
                res = least_squares(fun_trivial, 2.0, jac, method=self.method)
                assert_allclose(res.x, 0, atol=1e-4)

        assert_raises(ValueError, least_squares, fun_trivial, 2.0, jac='oops',
                      method=self.method)

    def test_nfev_options(self):
        for max_nfev in [None, 20]:
            res = least_squares(fun_trivial, 2.0, max_nfev=max_nfev,
                                method=self.method)
            assert_allclose(res.x, 0, atol=1e-4)

    def test_scaling_options(self):
        for scaling in [1.0, np.array([2.0]), 'jac']:
            res = least_squares(fun_trivial, 2.0, scaling=scaling)
            assert_allclose(res.x, 0)
        assert_raises(ValueError, least_squares, fun_trivial,
                      2.0, scaling='auto', method=self.method)
        assert_raises(ValueError, least_squares, fun_trivial,
                      2.0, scaling=-1.0, method=self.method)
        assert_raises(ValueError, least_squares, fun_trivial,
                      2.0, scaling=None, method=self.method)
        assert_raises(ValueError, least_squares, fun_trivial,
                      2.0, scaling=1.0+2.0j, method=self.method)

    def test_diff_step(self):
        # res1 and res2 should be equivalent.
        # res2 and res3 should be different.
        res1 = least_squares(fun_trivial, 2.0, diff_step=1e-1,
                             method=self.method)
        res2 = least_squares(fun_trivial, 2.0, diff_step=-1e-1,
                             method=self.method)
        res3 = least_squares(fun_trivial, 2.0,
                             diff_step=None, method=self.method)
        assert_allclose(res1.x, 0, atol=1e-4)
        assert_allclose(res2.x, 0, atol=1e-4)
        assert_allclose(res3.x, 0, atol=1e-4)
        assert_equal(res1.x, res2.x)
        assert_equal(res1.nfev, res2.nfev)
        assert_(res2.nfev != res3.nfev)

    def test_incorrect_options_usage(self):
        assert_raises(TypeError, least_squares, fun_trivial, 2.0,
                      method=self.method, options={'no_such_option': 100})
        assert_raises(TypeError, least_squares, fun_trivial, 2.0,
                      method=self.method, options={'max_nfev': 100})

    def test_full_result(self):
        res = least_squares(fun_trivial, 2.0, method=self.method)
        assert_equal(res.x, np.array([0]))
        assert_equal(res.cost, 12.5)
        assert_equal(res.fun, np.array([5]))
        assert_equal(res.jac, np.array([[0]]))
        assert_equal(res.optimality, 0)
        assert_equal(res.active_mask, np.array([0]))
        assert_(res.nfev < 10)
        if self.method == 'lm':
            assert_(res.njev is None)
        else:
            assert_(res.njev < 10)
        assert_(res.status > 0)
        assert_(res.success)

    def test_full_result_single_fev(self):
        # MINPACK checks the number of nfev after the iteration,
        # so it's hard to tell what he is going to compute.
        if self.method == 'lm':
            return

        res = least_squares(fun_trivial, 2.0, method=self.method,
                            max_nfev=1)
        assert_equal(res.x, np.array([2]))
        assert_equal(res.cost, 40.5)
        assert_equal(res.fun, np.array([9]))
        assert_equal(res.jac, np.array([[4]]))
        assert_equal(res.optimality, 36)
        assert_equal(res.active_mask, np.array([0]))
        assert_equal(res.nfev, 1)
        assert_equal(res.njev, 1)
        assert_equal(res.status, 0)
        assert_equal(res.success, 0)

    def test_rosenbrock(self):
        x0 = [-2, 1]
        x_opt = [1, 1]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            for jac, scaling, tr_solver in product(
                    ['2-point', '3-point', 'cs', jac_rosenbrock],
                    [1.0, np.array([1.0, 5.0]), 'jac'],
                    ['exact', 'lsmr']):
                res = least_squares(fun_rosenbrock, x0, jac, scaling=scaling,
                                    tr_solver=tr_solver, method=self.method)
                assert_allclose(res.x, x_opt)

    def test_rosenbrock_cropped(self):
        x0 = [-2, 1]
        if self.method == 'lm':
            assert_raises(ValueError, least_squares, fun_rosenbrock_cropped,
                          x0, method='lm')
        else:
            for jac, scaling, tr_solver in product(
                    ['2-point', '3-point', 'cs', jac_rosenbrock_cropped],
                    [1.0, np.array([1.0, 5.0]), 'jac'],
                    ['exact', 'lsmr']):
                res = least_squares(
                    fun_rosenbrock_cropped, x0, jac, scaling=scaling,
                    tr_solver=tr_solver, method=self.method)
                assert_allclose(res.cost, 0, atol=1e-14)

    def test_fun_wrong_dimensions(self):
        assert_raises(RuntimeError, least_squares, fun_wrong_dimensions,
                      2.0, method=self.method)

    def test_jac_wrong_dimensions(self):
        assert_raises(ValueError, least_squares, fun_trivial,
                      2.0, jac_wrong_dimensions, method=self.method)

    def test_fun_and_jac_inconsistent_dimensions(self):
        x0 = [1, 2]
        assert_raises(ValueError, least_squares, fun_rosenbrock, x0,
                      jac_rosenbrock_bad_dim, method=self.method)

    def test_x0_multidimensional(self):
        x0 = np.ones(4).reshape(2, 2)
        assert_raises(ValueError, least_squares, fun_trivial, x0,
                      method=self.method)


class BoundsMixin(object):
    def test_inconsistent(self):
        assert_raises(ValueError, least_squares, fun_trivial, 2.0,
                      bounds=(10.0, 0.0), method=self.method)

    def test_infeasible(self):
        assert_raises(ValueError, least_squares, fun_trivial, 2.0,
                      bounds=(3., 4), method=self.method)

    def test_wrong_number(self):
        assert_raises(ValueError, least_squares, fun_trivial, 2.,
                      bounds=(1., 2, 3), method=self.method)

    def test_inconsistent_shape(self):
        assert_raises(ValueError, least_squares, fun_trivial, 2.0,
                      bounds=(1.0, [2.0, 3.0]), method=self.method)
        # 1-D array wont't be broadcasted
        assert_raises(ValueError, least_squares, fun_rosenbrock, [1.0, 2.0],
                      bounds=([0.0], [3.0, 4.0]), method=self.method)

    def test_in_bounds(self):
        for jac in ['2-point', '3-point', 'cs', jac_trivial]:
            res = least_squares(fun_trivial, 2.0, jac=jac,
                                bounds=(-1.0, 3.0), method=self.method)
            assert_allclose(res.x, 0.0, atol=1e-4)
            assert_equal(res.active_mask, [0])
            assert_(-1 <= res.x <= 3)
            res = least_squares(fun_trivial, 2.0, jac=jac,
                                bounds=(0.5, 3.0), method=self.method)
            assert_allclose(res.x, 0.5, atol=1e-4)
            assert_equal(res.active_mask, [-1])
            assert_(0.5 <= res.x <= 3)

    def test_bounds_shape(self):
        for jac in ['2-point', '3-point', 'cs', jac_2d_trivial]:
            x0 = [1.0, 1.0]
            res = least_squares(fun_2d_trivial, x0, jac=jac)
            assert_allclose(res.x, [0.0, 0.0])
            res = least_squares(fun_2d_trivial, x0, jac=jac,
                                bounds=(0.5, [2.0, 2.0]), method=self.method)
            assert_allclose(res.x, [0.5, 0.5])
            res = least_squares(fun_2d_trivial, x0, jac=jac,
                                bounds=([0.3, 0.2], 3.0), method=self.method)
            assert_allclose(res.x, [0.3, 0.2])
            res = least_squares(
                fun_2d_trivial, x0, jac=jac, bounds=([-1, 0.5], [1.0, 3.0]),
                method=self.method)
            assert_allclose(res.x, [0.0, 0.5], atol=1e-5)

    def test_rosenbrock_bounds(self):
        x0_1 = np.array([-2.0, 1.0])
        x0_2 = np.array([2.0, 2.0])
        x0_3 = np.array([-2.0, 2.0])
        x0_4 = np.array([0.0, 2.0])
        x0_5 = np.array([-1.2, 1.0])
        problems = [
            (x0_1, ([-np.inf, -1.5], np.inf)),
            (x0_2, ([-np.inf, 1.5], np.inf)),
            (x0_3, ([-np.inf, 1.5], np.inf)),
            (x0_4, ([-np.inf, 1.5], [1.0, np.inf])),
            (x0_2, ([1.0, 1.5], [3.0, 3.0])),
            (x0_5, ([-50.0, 0.0], [0.5, 100]))
        ]
        for x0, bounds in problems:
            for jac, scaling, tr_solver in product(
                    ['2-point', '3-point', 'cs', jac_rosenbrock],
                    [1.0, [1.0, 2.0], 'jac'],
                    ['exact', 'lsmr']):
                res = least_squares(fun_rosenbrock, x0, jac, bounds,
                                    scaling=scaling, tr_solver=tr_solver,
                                    method=self.method)
                assert_allclose(res.optimality, 0.0, atol=1e-5)


class SparseMixin(object):
    def test_exact_tr_solver(self):
        p = BroydenTridiagonal()
        assert_raises(ValueError, least_squares, p.fun, p.x0, p.jac,
                      tr_solver='exact', method=self.method)
        assert_raises(ValueError, least_squares, p.fun, p.x0,
                      tr_solver='exact', jac_sparsity=p.sparsity,
                      method=self.method)

    def test_equivalence(self):
        sparse = BroydenTridiagonal(mode='sparse')
        dense = BroydenTridiagonal(mode='dense')
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            res_sparse = least_squares(
                sparse.fun, sparse.x0, jac=sparse.jac,
                method=self.method)
            res_dense = least_squares(
                dense.fun, dense.x0, jac=sparse.jac,
                method=self.method)
            assert_equal(res_sparse.nfev, res_dense.nfev)
            assert_allclose(res_sparse.x, res_dense.x, atol=1e-20)
            assert_allclose(res_sparse.cost, 0, atol=1e-20)
            assert_allclose(res_dense.cost, 0, atol=1e-20)

    def test_tr_options(self):
        p = BroydenTridiagonal()
        res = least_squares(p.fun, p.x0, p.jac, method=self.method,
                            tr_options={'btol': 1e-10})
        assert_allclose(res.cost, 0, atol=1e-20)

    def test_wrong_parameters(self):
        p = BroydenTridiagonal()
        assert_raises(ValueError, least_squares, p.fun, p.x0, p.jac,
                      tr_solver='best', method=self.method)
        assert_raises(TypeError, least_squares, p.fun, p.x0, p.jac,
                      tr_solver='lsmr', tr_options={'tol': 1e-10})

    def test_solver_selection(self):
        sparse = BroydenTridiagonal(mode='sparse')
        dense = BroydenTridiagonal(mode='dense')
        res_sparse = least_squares(sparse.fun, sparse.x0, jac=sparse.jac,
                                   method=self.method)
        res_dense = least_squares(dense.fun, dense.x0, jac=dense.jac,
                                  method=self.method)
        assert_allclose(res_sparse.cost, 0, atol=1e-20)
        assert_allclose(res_dense.cost, 0, atol=1e-20)
        assert_(issparse(res_sparse.jac))
        assert_(isinstance(res_dense.jac, np.ndarray))

    def test_numerical_jac(self):
        p = BroydenTridiagonal()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            for jac in ['2-point', '3-point', 'cs']:
                res_dense = least_squares(p.fun, p.x0, jac, method=self.method)
                res_sparse = least_squares(
                    p.fun, p.x0, jac,method=self.method,
                    jac_sparsity=p.sparsity)
                assert_equal(res_dense.nfev, res_sparse.nfev)
                assert_allclose(res_dense.x, res_sparse.x, atol=1e-20)
                assert_allclose(res_dense.cost, 0, atol=1e-20)
                assert_allclose(res_sparse.cost, 0, atol=1e-20)

    def test_with_bounds(self):
        p = BroydenTridiagonal()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            for jac, jac_sparsity in product(
                    [p.jac, '2-point', '3-point', 'cs'], [None, p.sparsity]):
                res_1 = least_squares(
                    p.fun, p.x0, jac, bounds=(p.lb, np.inf),
                    method=self.method,jac_sparsity=jac_sparsity)
                res_2 = least_squares(
                    p.fun, p.x0, jac, bounds=(-np.inf, p.ub),
                    method=self.method, jac_sparsity=jac_sparsity)
                res_3 = least_squares(
                    p.fun, p.x0, jac, bounds=(p.lb, p.ub),
                    method=self.method, jac_sparsity=jac_sparsity)
                assert_allclose(res_1.optimality, 0, atol=1e-10)
                assert_allclose(res_2.optimality, 0, atol=1e-10)
                assert_allclose(res_3.optimality, 0, atol=1e-10)

    def test_wrong_jac_sparsity(self):
        p = BroydenTridiagonal()
        sparsity = p.sparsity[:-1]
        assert_raises(ValueError, least_squares, p.fun, p.x0,
                      jac_sparsity=sparsity, method=self.method)

    def test_linear_operator(self):
        p = BroydenTridiagonal(mode='operator')
        res = least_squares(p.fun, p.x0, p.jac, method=self.method)
        assert_allclose(res.cost, 0.0, atol=1e-20)
        assert_raises(ValueError, least_squares, p.fun, p.x0, p.jac,
                      method=self.method, tr_solver='exact')

    def test_scaling(self):
        p = BroydenTridiagonal()
        res = least_squares(p.fun, p.x0, p.jac, method=self.method,
                            scaling='jac')
        assert_allclose(res.cost, 0.0, atol=1e-20)

        p = BroydenTridiagonal(mode='operator')
        assert_raises(ValueError, least_squares, p.fun, p.x0, p.jac,
                      method=self.method, scaling='jac')


class TestDogbox(BaseMixin, BoundsMixin, SparseMixin):
    method = 'dogbox'


class TestTRF(BaseMixin, BoundsMixin, SparseMixin):
    method = 'trf'

    def test_lsmr_regularization(self):
        p = BroydenTridiagonal()
        for regularize in [True, False]:
            res = least_squares(p.fun, p.x0, p.jac, method='trf',
                                tr_options={'regularize': regularize})
            assert_allclose(res.cost, 0, atol=1e-20)


class TestLM(BaseMixin):
    method = 'lm'

    def test_bounds_not_supported(self):
        assert_raises(ValueError, least_squares, fun_trivial,
                      2.0, bounds=(-3.0, 3.0), method='lm')

    def test_m_less_n_not_supported(self):
        x0 = [-2, 1]
        assert_raises(ValueError, least_squares, fun_rosenbrock_cropped, x0,
                      method='lm')

    def test_sparse_not_supported(self):
        p = BroydenTridiagonal()
        assert_raises(ValueError, least_squares, p.fun, p.x0, p.jac,
                      method='lm')

    def test_jac_sparsity_not_supported(self):
        assert_raises(ValueError, least_squares, fun_trivial, 2.0,
                      jac_sparsity=[1], method='lm')

    def test_LinearOperator_not_supported(self):
        p = BroydenTridiagonal(mode="operator")
        assert_raises(ValueError, least_squares, p.fun, p.x0, p.jac,
                      method='lm')


def test_basic():
    # test that 'method' arg is really optional
    res = least_squares(fun_trivial, 2.0)
    assert_allclose(res.x, 0, atol=1e-10)


if __name__ == "__main__":
    run_module_suite()
