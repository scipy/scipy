import math
from itertools import product

import numpy as np
from numpy.testing import assert_allclose, assert_equal, assert_
import pytest
from pytest import raises as assert_raises

from scipy._lib._testutils import IS_WASM

from scipy._lib._util import (MapWrapper, _ScalarFunctionWrapper,
                              _item_for_scalar_function)
from scipy.sparse import csr_array, csc_array, lil_array
from scipy._lib._array_api import (xp_result_type, make_xp_test_case, array_namespace,
                                   xp_assert_close)
from scipy._external import array_api_extra as xpx

from scipy.optimize._numdiff import (
    _adjust_scheme_to_bounds, approx_derivative, check_derivative,
    group_columns, _eps_for_method, _compute_absolute_step)
from scipy.optimize import rosen


def test_group_columns():
    structure = [
        [1, 1, 0, 0, 0, 0],
        [1, 1, 1, 0, 0, 0],
        [0, 1, 1, 1, 0, 0],
        [0, 0, 1, 1, 1, 0],
        [0, 0, 0, 1, 1, 1],
        [0, 0, 0, 0, 1, 1],
        [0, 0, 0, 0, 0, 0]
    ]
    for transform in [np.asarray, csr_array, csc_array, lil_array]:
        A = transform(structure)
        order = np.arange(6)
        groups_true = np.array([0, 1, 2, 0, 1, 2])
        groups = group_columns(A, order)
        assert_equal(groups, groups_true)

        order = [1, 2, 4, 3, 5, 0]
        groups_true = np.array([2, 0, 1, 2, 0, 1])
        groups = group_columns(A, order)
        assert_equal(groups, groups_true)

    # Test repeatability.
    groups_1 = group_columns(A)
    groups_2 = group_columns(A)
    assert_equal(groups_1, groups_2)


@make_xp_test_case(
    _eps_for_method
)
def test_correct_fp_eps(xp):
    # check that relative step size is correct for FP size
    EPS = xp.finfo(xp.float64).eps
    relative_step = {"2-point": EPS**0.5,
                     "3-point": EPS**(1/3),
                     "cs": EPS**0.5}
    for method in ['2-point', '3-point', 'cs']:
        x = _eps_for_method(xp.float64, xp.float64, method, xp=xp)
        y = xp.asarray(relative_step[method], dtype=xp.float64)
        xp_assert_close(
            x,
            y,
        )
        assert x.dtype == xp.float64

    # check another FP size
    EPS = xp.finfo(xp.float32).eps
    relative_step = {"2-point": EPS**0.5,
                    "3-point": EPS**(1/3),
                     "cs": EPS**0.5}

    for method in ['2-point', '3-point', 'cs']:
        x = _eps_for_method(xp.float64, xp.float32, method, xp=xp)
        y = xp.asarray(relative_step[method], dtype=xp.float32)
        xp_assert_close(
            x,
            y
        )
        assert x.dtype == xp.float32
        x = _eps_for_method(xp.float32, xp.float64, method, xp=xp)
        xp_assert_close(
            x,
            y
        )
        assert x.dtype == xp.float32
        x = _eps_for_method(xp.float32, xp.float32, method, xp=xp)
        xp_assert_close(
            x,
            y
        )
        assert x.dtype == xp.float32


class TestAdjustSchemeToBounds:
    def test_no_bounds(self):
        x0 = np.zeros(3)
        h = np.full(3, 1e-2)
        inf_lower = np.empty_like(x0)
        inf_upper = np.empty_like(x0)
        inf_lower.fill(-np.inf)
        inf_upper.fill(np.inf)

        h_adjusted, one_sided = _adjust_scheme_to_bounds(
            x0, h, 1, '1-sided', inf_lower, inf_upper)
        assert_allclose(h_adjusted, h)
        assert_(np.all(one_sided))

        h_adjusted, one_sided = _adjust_scheme_to_bounds(
            x0, h, 2, '1-sided', inf_lower, inf_upper)
        assert_allclose(h_adjusted, h)
        assert_(np.all(one_sided))

        h_adjusted, one_sided = _adjust_scheme_to_bounds(
            x0, h, 1, '2-sided', inf_lower, inf_upper)
        assert_allclose(h_adjusted, h)
        assert_(np.all(~one_sided))

        h_adjusted, one_sided = _adjust_scheme_to_bounds(
            x0, h, 2, '2-sided', inf_lower, inf_upper)
        assert_allclose(h_adjusted, h)
        assert_(np.all(~one_sided))

    def test_with_bound(self):
        x0 = np.array([0.0, 0.85, -0.85])
        lb = -np.ones(3)
        ub = np.ones(3)
        h = np.array([1, 1, -1]) * 1e-1

        h_adjusted, _ = _adjust_scheme_to_bounds(x0, h, 1, '1-sided', lb, ub)
        assert_allclose(h_adjusted, h)

        h_adjusted, _ = _adjust_scheme_to_bounds(x0, h, 2, '1-sided', lb, ub)
        assert_allclose(h_adjusted, np.array([1, -1, 1]) * 1e-1)

        h_adjusted, one_sided = _adjust_scheme_to_bounds(
            x0, h, 1, '2-sided', lb, ub)
        assert_allclose(h_adjusted, np.abs(h))
        assert_(np.all(~one_sided))

        h_adjusted, one_sided = _adjust_scheme_to_bounds(
            x0, h, 2, '2-sided', lb, ub)
        assert_allclose(h_adjusted, np.array([1, -1, 1]) * 1e-1)
        assert_equal(one_sided, np.array([False, True, True]))

    def test_tight_bounds(self):
        lb = np.array([-0.03, -0.03])
        ub = np.array([0.05, 0.05])
        x0 = np.array([0.0, 0.03])
        h = np.array([-0.1, -0.1])

        h_adjusted, _ = _adjust_scheme_to_bounds(x0, h, 1, '1-sided', lb, ub)
        assert_allclose(h_adjusted, np.array([0.05, -0.06]))

        h_adjusted, _ = _adjust_scheme_to_bounds(x0, h, 2, '1-sided', lb, ub)
        assert_allclose(h_adjusted, np.array([0.025, -0.03]))

        h_adjusted, one_sided = _adjust_scheme_to_bounds(
            x0, h, 1, '2-sided', lb, ub)
        assert_allclose(h_adjusted, np.array([0.03, -0.03]))
        assert_equal(one_sided, np.array([False, True]))

        h_adjusted, one_sided = _adjust_scheme_to_bounds(
            x0, h, 2, '2-sided', lb, ub)
        assert_allclose(h_adjusted, np.array([0.015, -0.015]))
        assert_equal(one_sided, np.array([False, True]))


@make_xp_test_case(
    approx_derivative
)
def test_xp_smoke(xp):
    # a quick smoke test to check approx_derivative with different array-api
    dtypes = [xp.float32, xp.float64]

    def fun(x):
        return xp.asarray(2*(x[0]**2 + x[1]**2 - 1) - x[0], dtype=x.dtype)

    def grad(x):
        return xp.stack([4*x[0]-1, 4*x[1]])

    def hess(x):
        return 4 * xp.eye(2, dtype=x.dtype)

    for dtype, method, abs_step in product(
        dtypes,
        ['2-point', '3-point'],
        [False, True]
    ):
        x0 = xp.asarray([3.0, 4.0], dtype=dtype)
        f0 = fun(x0)

        # we don't need to check specified rel_steps, as this is tested by
        # test__compute_absolute_step
        if abs_step:
            _h = _compute_absolute_step(None, x0, f0, method, xp=xp)
        else:
            _h = None

        g = approx_derivative(fun, x0, method=method, abs_step=_h)
        xp_assert_close(g, grad(x0))

        h = approx_derivative(grad, x0, method=method, abs_step=_h)
        xp_assert_close(h, hess(x0))

        # test bounds setting. If we set a lb/ub equal to x0, then
        # _adjust_scheme_to_bounds is exercised.
        _b = xp.asarray([3.0, 4.0], dtype=dtype)
        g = approx_derivative(fun, x0, bounds=(-xp.inf, _b), method=method, abs_step=_h)
        xp_assert_close(g, grad(x0))

        h = approx_derivative(grad, x0, bounds=(_b, xp.inf), method=method, abs_step=_h)
        xp_assert_close(h, hess(x0))


class TestFunctions:
    def fun_scalar_scalar(self, x):
        xp = array_namespace(x)
        return xp.sinh(x)

    def jac_scalar_scalar(self, x):
        xp = array_namespace(x)
        return xp.cosh(x)

    def fun_scalar_vector(self, x):
        xp = array_namespace(x)
        return xp.stack([x[0] ** 2, xp.tan(x[0]), xp.exp(x[0])])

    def jac_scalar_vector(self, x):
        xp = array_namespace(x)
        return xp.reshape(
            xp.stack([2 * x[0], xp.cos(x[0]) ** -2, xp.exp(x[0])]),
            (-1, 1)
        )

    def fun_vector_scalar(self, x):
        xp = array_namespace(x)
        return xp.sin(x[0] * x[1]) * xp.log(x[0])

    def jac_vector_scalar(self, x):
        xp = array_namespace(x)
        return xp.stack([
            x[1] * xp.cos(x[0] * x[1]) * xp.log(x[0]) +
            xp.sin(x[0] * x[1]) / x[0],
            x[0] * xp.cos(x[0] * x[1]) * xp.log(x[0])],
        )

    def fun_vector_vector(self, x):
        xp = array_namespace(x)
        return xp.stack([
            x[0] * xp.sin(x[1]),
            x[1] * xp.cos(x[0]),
            x[0] ** 3 * x[1] ** -0.5]
        )

    def jac_vector_vector(self, x):
        xp = array_namespace(x)
        return xp.stack([
            [xp.sin(x[1]), x[0] * xp.cos(x[1])],
            [-x[1] * xp.sin(x[0]), xp.cos(x[0])],
            [3 * x[0] ** 2 * x[1] ** -0.5, -0.5 * x[0] ** 3 * x[1] ** -1.5]],
        )


class TestApproxDerivativesDense(TestFunctions):

    def wrong_dimensions_fun(self, x):
        return np.array([x**2, np.tan(x), np.exp(x)])

    def fun_vector_vector_with_arg(self, x, arg):
        """Used to test passing custom arguments with check_derivative()"""
        assert arg == 42
        return np.array([
            x[0] * np.sin(x[1]),
            x[1] * np.cos(x[0]),
            x[0] ** 3 * x[1] ** -0.5
        ])

    def jac_vector_vector_with_arg(self, x, arg):
        """Used to test passing custom arguments with check_derivative()"""
        assert arg == 42
        return np.array([
            [np.sin(x[1]), x[0] * np.cos(x[1])],
            [-x[1] * np.sin(x[0]), np.cos(x[0])],
            [3 * x[0] ** 2 * x[1] ** -0.5, -0.5 * x[0] ** 3 * x[1] ** -1.5]
        ])

    def fun_parametrized(self, x, c0, c1=1.0):
        return np.array([np.exp(c0 * x[0]), np.exp(c1 * x[1])])

    def jac_parametrized(self, x, c0, c1=0.1):
        return np.array([
            [c0 * np.exp(c0 * x[0]), 0],
            [0, c1 * np.exp(c1 * x[1])]
        ])

    def fun_with_nan(self, x):
        return x if np.abs(x) <= 1e-8 else np.nan

    def jac_with_nan(self, x):
        return 1.0 if np.abs(x) <= 1e-8 else np.nan

    def fun_zero_jacobian(self, x):
        return np.array([x[0] * x[1], np.cos(x[0] * x[1])])

    def jac_zero_jacobian(self, x):
        return np.array([
            [x[1], x[0]],
            [-x[1] * np.sin(x[0] * x[1]), -x[0] * np.sin(x[0] * x[1])]
        ])

    def jac_non_numpy(self, x):
        # x can be a scalar or an array [val].
        # Cast to true scalar before handing over to math.exp
        xp = np.asarray(x).item()
        return math.exp(xp)

    @make_xp_test_case(
        approx_derivative
    )
    def test_scalar_scalar(self, xp):
        x0 = _item_for_scalar_function(xp.asarray(1.0, dtype=xp.float64))
        jac_diff_2 = approx_derivative(self.fun_scalar_scalar, x0,
                                       method='2-point')
        jac_diff_3 = approx_derivative(self.fun_scalar_scalar, x0)
        jac_diff_4 = approx_derivative(self.fun_scalar_scalar, x0,
                                       method='cs')
        jac_true = xpx.atleast_nd(self.jac_scalar_scalar(x0), ndim=1, xp=xp)
        xp_assert_close(jac_diff_2, jac_true, rtol=1e-6)
        xp_assert_close(jac_diff_3, jac_true, rtol=1e-9)
        xp_assert_close(jac_diff_4, jac_true, rtol=1e-12)

    @make_xp_test_case(
        approx_derivative
    )
    def test_scalar_scalar_abs_step(self, xp):
        # can approx_derivative use abs_step?
        x0 = _item_for_scalar_function(xp.asarray(1.0, dtype=xp.float64))
        jac_diff_2 = approx_derivative(self.fun_scalar_scalar, x0,
                                       method='2-point', abs_step=1.49e-8)
        jac_diff_3 = approx_derivative(self.fun_scalar_scalar, x0,
                                       abs_step=1.49e-8)
        jac_diff_4 = approx_derivative(self.fun_scalar_scalar, x0,
                                       method='cs', abs_step=1.49e-8)
        jac_true = xpx.atleast_nd(self.jac_scalar_scalar(x0), ndim=1, xp=xp)
        xp_assert_close(jac_diff_2, jac_true, rtol=1e-6)
        xp_assert_close(jac_diff_3, jac_true, rtol=1e-9)
        xp_assert_close(jac_diff_4, jac_true, rtol=1e-12)

    @make_xp_test_case(
        approx_derivative
    )
    def test_scalar_vector(self, xp):
        x0 = _item_for_scalar_function(xp.asarray(0.5, dtype=xp.float64))
        with MapWrapper(1 if IS_WASM else 2) as mapper:
            jac_diff_2 = approx_derivative(self.fun_scalar_vector, x0,
                                           method='2-point', workers=mapper)
        jac_diff_3 = approx_derivative(self.fun_scalar_vector, x0, workers=map)
        jac_diff_4 = approx_derivative(self.fun_scalar_vector, x0,
                                       method='cs', workers=None)
        jac_true = self.jac_scalar_vector(
            xpx.atleast_nd(x0, ndim=1, xp=xp)
        )
        xp_assert_close(jac_diff_2, jac_true, rtol=1e-6)
        xp_assert_close(jac_diff_3, jac_true, rtol=1e-9)
        xp_assert_close(jac_diff_4, jac_true, rtol=1e-12)

    @pytest.mark.fail_slow(5.0)
    def test_workers_evaluations_and_nfev(self):
        # check that nfev consumed by approx_derivative is tracked properly
        # and that parallel evaluation is same as series
        n_workers = 1 if IS_WASM else 2
        x0 = [0.5, 1.5, 2.0]
        with MapWrapper(n_workers) as mapper:
            md2, mdct2 = approx_derivative(rosen, x0,
                                           method='2-point', workers=mapper,
                                           full_output=True)
            md3, mdct3 = approx_derivative(rosen, x0,
                                           workers=mapper, full_output=True)
        # supply a number for workers. This is not normally recommended
        # for upstream workers as setting up processes incurs a large overhead
        md4, mdct4 = approx_derivative(rosen, x0,
                                       method='cs', workers=n_workers,
                                       full_output=True)

        sfr = _ScalarFunctionWrapper(rosen)
        d2, dct2 = approx_derivative(sfr, x0, method='2-point', full_output=True)
        assert_equal(dct2['nfev'], sfr.nfev)
        sfr.nfev = 0
        d3, dct3 = approx_derivative(sfr, x0, full_output=True)
        assert_equal(dct3['nfev'], sfr.nfev)
        sfr.nfev = 0
        d4, dct4 = approx_derivative(sfr, x0, method='cs', full_output=True)
        assert_equal(dct4['nfev'], sfr.nfev)

        assert_equal(mdct2['nfev'], dct2['nfev'])
        assert_equal(mdct3['nfev'], dct3['nfev'])
        assert_equal(mdct4['nfev'], dct4['nfev'])
        # also check that gradients are equivalent
        assert_equal(md2, d2)
        assert_equal(md3, d3)
        assert_equal(md4, d4)

    @make_xp_test_case(
        approx_derivative
    )
    def test_vector_scalar(self, xp):
        x0 = xp.asarray([100.0, -0.5], dtype=xp.float64)
        jac_diff_2 = approx_derivative(self.fun_vector_scalar, x0,
                                       method='2-point')
        jac_diff_3 = approx_derivative(self.fun_vector_scalar, x0)
        jac_diff_4 = approx_derivative(self.fun_vector_scalar, x0,
                                       method='cs')
        jac_true = self.jac_vector_scalar(x0)
        xp_assert_close(jac_diff_2, jac_true, rtol=1e-6)
        xp_assert_close(jac_diff_3, jac_true, rtol=1e-7)
        xp_assert_close(jac_diff_4, jac_true, rtol=1e-12)

    @make_xp_test_case(
        approx_derivative
    )
    def test_vector_scalar_abs_step(self, xp):
        # can approx_derivative use abs_step?
        x0 = xp.asarray([100.0, -0.5], dtype=xp.float64)
        jac_diff_2 = approx_derivative(self.fun_vector_scalar, x0,
                                       method='2-point', abs_step=1.49e-8)
        jac_diff_3 = approx_derivative(self.fun_vector_scalar, x0,
                                       abs_step=1.49e-8, rel_step=np.inf)
        jac_diff_4 = approx_derivative(self.fun_vector_scalar, x0,
                                       method='cs', abs_step=1.49e-8)
        jac_true = xp.asarray(self.jac_vector_scalar([100.0, -0.5]))
        xp_assert_close(jac_diff_2, jac_true, rtol=1e-6)
        xp_assert_close(jac_diff_3, jac_true, rtol=3e-9)
        xp_assert_close(jac_diff_4, jac_true, rtol=1e-12)

    @make_xp_test_case(
        approx_derivative
    )
    def test_vector_vector(self, xp):
        x0 = xp.asarray([-100.0, 0.2], dtype=xp.float64)
        jac_diff_2 = approx_derivative(self.fun_vector_vector, x0,
                                       method='2-point')
        jac_diff_3 = approx_derivative(self.fun_vector_vector, x0)
        with MapWrapper(1 if IS_WASM else 2) as mapper:
            jac_diff_4 = approx_derivative(self.fun_vector_vector, x0,
                                           method='cs', workers=mapper)
        output = self.jac_vector_vector(np.array([-100.0, 0.2]))
        jac_true = xp.asarray(output, dtype=xp.float64)
        xp_assert_close(jac_diff_2, jac_true, rtol=1e-5)
        xp_assert_close(jac_diff_3, jac_true, rtol=1e-6)
        xp_assert_close(jac_diff_4, jac_true, rtol=2e-8)

    def test_wrong_dimensions(self):
        x0 = 1.0
        assert_raises(RuntimeError, approx_derivative,
                      self.wrong_dimensions_fun, x0)
        f0 = self.wrong_dimensions_fun(np.atleast_1d(x0))
        assert_raises(ValueError, approx_derivative,
                      self.wrong_dimensions_fun, x0, f0=f0)

    @make_xp_test_case(
        approx_derivative
    )
    def test_custom_rel_step(self, xp):
        x0 = xp.asarray([-0.1, 0.1], dtype=xp.float64)
        jac_diff_2 = approx_derivative(self.fun_vector_vector, x0,
                                       method='2-point', rel_step=1e-4)
        jac_diff_3 = approx_derivative(self.fun_vector_vector, x0,
                                       rel_step=1e-4)
        jac_true = xp.asarray(self.jac_vector_vector([-0.1, 0.1]))
        xp_assert_close(jac_diff_2, jac_true, rtol=1e-2)
        xp_assert_close(jac_diff_3, jac_true, rtol=1e-4)

    def test_options(self):
        x0 = np.array([1.0, 1.0])
        c0 = -1.0
        c1 = 1.0
        lb = 0.0
        ub = 2.0
        f0 = self.fun_parametrized(x0, c0, c1=c1)
        rel_step = np.array([-1e-6, 1e-7])
        jac_true = self.jac_parametrized(x0, c0, c1)
        jac_diff_2 = approx_derivative(
            self.fun_parametrized, x0, method='2-point', rel_step=rel_step,
            f0=f0, args=(c0,), kwargs=dict(c1=c1), bounds=(lb, ub))
        jac_diff_3 = approx_derivative(
            self.fun_parametrized, x0, rel_step=rel_step,
            f0=f0, args=(c0,), kwargs=dict(c1=c1), bounds=(lb, ub))
        assert_allclose(jac_diff_2, jac_true, rtol=1e-6)
        assert_allclose(jac_diff_3, jac_true, rtol=1e-9)

    @make_xp_test_case(
        approx_derivative
    )
    def test_with_bounds_2_point(self, xp):
        lb = -xp.ones(2)
        ub = xp.ones(2)

        x0 = xp.asarray([-2.0, 0.2], dtype=xp.float64)
        assert_raises(ValueError, approx_derivative,
                      self.fun_vector_vector, x0, bounds=(lb, ub))

        x0 = xp.asarray([-1.0, 1.0], dtype=xp.float64)
        jac_diff = approx_derivative(self.fun_vector_vector, x0,
                                     method='2-point', bounds=(lb, ub))
        jac_true = xp.asarray(self.jac_vector_vector([-1.0, 1.0]))
        xp_assert_close(jac_diff, jac_true, rtol=1e-6)

    @make_xp_test_case(
        approx_derivative
    )
    def test_with_bounds_3_point(self, xp):
        lb = xp.asarray([1.0, 1.0], dtype=xp.float64)
        ub = xp.asarray([2.0, 2.0], dtype=xp.float64)

        x0 = xp.asarray([1.0, 2.0], dtype=xp.float64)
        jac_true = xp.asarray(
            self.jac_vector_vector(np.array([1.0, 2.0])), dtype=xp.float64
        )

        jac_diff = approx_derivative(self.fun_vector_vector, x0)
        xp_assert_close(jac_diff, jac_true, rtol=1e-9)

        jac_diff = approx_derivative(self.fun_vector_vector, x0,
                                     bounds=(lb, xp.inf))
        xp_assert_close(jac_diff, jac_true, rtol=1e-9)

        jac_diff = approx_derivative(self.fun_vector_vector, x0,
                                     bounds=(-xp.inf, ub))
        xp_assert_close(jac_diff, jac_true, rtol=1e-9)

        jac_diff = approx_derivative(self.fun_vector_vector, x0,
                                     bounds=(lb, ub))
        xp_assert_close(jac_diff, jac_true, rtol=1e-9)

    def test_tight_bounds(self):
        x0 = np.array([10.0, 10.0])
        lb = x0 - 3e-9
        ub = x0 + 2e-9
        jac_true = self.jac_vector_vector(x0)
        jac_diff = approx_derivative(
            self.fun_vector_vector, x0, method='2-point', bounds=(lb, ub))
        assert_allclose(jac_diff, jac_true, rtol=1e-6)
        jac_diff = approx_derivative(
            self.fun_vector_vector, x0, method='2-point',
            rel_step=1e-6, bounds=(lb, ub))
        assert_allclose(jac_diff, jac_true, rtol=1e-6)

        jac_diff = approx_derivative(
            self.fun_vector_vector, x0, bounds=(lb, ub))
        assert_allclose(jac_diff, jac_true, rtol=1e-6)
        jac_diff = approx_derivative(
            self.fun_vector_vector, x0, rel_step=1e-6, bounds=(lb, ub))
        assert_allclose(jac_true, jac_diff, rtol=1e-6)

    def test_bound_switches(self):
        lb = -1e-8
        ub = 1e-8
        x0 = 0.0
        jac_true = self.jac_with_nan(x0)
        jac_diff_2 = approx_derivative(
            self.fun_with_nan, x0, method='2-point', rel_step=1e-6,
            bounds=(lb, ub))
        jac_diff_3 = approx_derivative(
            self.fun_with_nan, x0, rel_step=1e-6, bounds=(lb, ub))
        assert_allclose(jac_diff_2, jac_true, rtol=1e-6)
        assert_allclose(jac_diff_3, jac_true, rtol=1e-9)

        x0 = 1e-8
        jac_true = self.jac_with_nan(x0)
        jac_diff_2 = approx_derivative(
            self.fun_with_nan, x0, method='2-point', rel_step=1e-6,
            bounds=(lb, ub))
        jac_diff_3 = approx_derivative(
            self.fun_with_nan, x0, rel_step=1e-6, bounds=(lb, ub))
        assert_allclose(jac_diff_2, jac_true, rtol=1e-6)
        assert_allclose(jac_diff_3, jac_true, rtol=1e-9)

    def test_non_numpy(self):
        x0 = 1.0
        jac_true = self.jac_non_numpy(x0)
        jac_diff_2 = approx_derivative(self.jac_non_numpy, x0,
                                       method='2-point')
        jac_diff_3 = approx_derivative(self.jac_non_numpy, x0)
        assert_allclose(jac_diff_2, jac_true, rtol=1e-6)
        assert_allclose(jac_diff_3, jac_true, rtol=1e-8)

        # math.exp cannot handle complex arguments, hence this raises
        assert_raises(TypeError, approx_derivative, self.jac_non_numpy, x0,
                      **dict(method='cs'))

    def test_fp(self):
        # checks that approx_derivative works for FP size other than 64.
        # Example is derived from the minimal working example in gh12991.
        rng = np.random.default_rng(1)

        def func(p, x):
            return p[0] + p[1] * x

        def err(p, x, y):
            return func(p, x) - y

        x = np.linspace(0, 1, 100, dtype=np.float64)
        y = rng.random(size=100, dtype=np.float64)
        p0 = np.array([-1.0, -1.0])

        jac_fp64 = approx_derivative(err, p0, method='2-point', args=(x, y))

        # parameter vector is float32, func output is float64
        jac_fp = approx_derivative(err, p0.astype(np.float32),
                                   method='2-point', args=(x, y))
        assert err(p0, x, y).dtype == np.float64
        assert_allclose(jac_fp, jac_fp64, atol=1e-3)

        # parameter vector is float64, func output is float32
        def err_fp32(p):
            assert p.dtype == np.float32
            return err(p, x, y).astype(np.float32)

        jac_fp = approx_derivative(err_fp32, p0.astype(np.float32),
                                   method='2-point')
        assert_allclose(jac_fp, jac_fp64, atol=1e-3)
        assert jac_fp.dtype == np.float32

        # check upper bound of error on the derivative for 2-point
        def f(x):
            return np.sin(x)
        def g(x):
            return np.cos(x)
        def hess(x):
            return -np.sin(x)

        def calc_atol(h, x0, f, hess, EPS):
            # truncation error
            t0 = h / 2 * max(np.abs(hess(x0)), np.abs(hess(x0 + h)))
            # roundoff error. There may be a divisor (>1) missing from
            # the following line, so this contribution is possibly
            # overestimated
            t1 = EPS / h * max(np.abs(f(x0)), np.abs(f(x0 + h)))
            return t0 + t1

        for dtype in [np.float16, np.float32, np.float64]:
            EPS = np.finfo(dtype).eps
            x0 = np.array(1.0).astype(dtype)
            h = _compute_absolute_step(None, x0, f(x0), '2-point')
            atol = calc_atol(h, x0, f, hess, EPS)
            err = approx_derivative(f, x0, method='2-point',
                                    abs_step=h) - g(x0)
            assert abs(err) < atol


    @pytest.mark.parametrize("x0_dtype", (np.float16, np.float32, np.float64))
    @pytest.mark.parametrize("f0_dtype", (np.float16, np.float32, np.float64))
    @pytest.mark.parametrize("method", ['2-point', '3-point'])
    def test_check_dtype(self, x0_dtype, f0_dtype, method):
        # the output of approx_derivative should be the promoted
        # type of x0, f0.

        # both are of the same dtype
        x = np.array([2.0, 3.0, 4.0], dtype=x0_dtype)

        def f(x):
            return f0_dtype(rosen(x))

        promoted_type = xp_result_type(
        x,
            f(x),
            force_floating=True,
            xp=np
        )
        g = approx_derivative(f, x, method=method)
        assert g.dtype == promoted_type

        # setting abs_step or rel_step shouldn't change output dtype
        g = approx_derivative(f, x, rel_step=np.float16(0.1), method=method)
        assert g.dtype == promoted_type

        g = approx_derivative(f, x, abs_step=np.float16(0.1), method=method)
        assert g.dtype == promoted_type

    def test_check_derivative(self):
        x0 = np.array([-10.0, 10])
        accuracy = check_derivative(self.fun_vector_vector,
                                    self.jac_vector_vector, x0)
        assert_(accuracy < 1e-9)
        accuracy = check_derivative(self.fun_vector_vector,
                                    self.jac_vector_vector, x0)
        assert_(accuracy < 1e-6)

        x0 = np.array([0.0, 0.0])
        accuracy = check_derivative(self.fun_zero_jacobian,
                                    self.jac_zero_jacobian, x0)
        assert_(accuracy == 0)
        accuracy = check_derivative(self.fun_zero_jacobian,
                                    self.jac_zero_jacobian, x0)
        assert_(accuracy == 0)

    def test_check_derivative_with_kwargs(self):
        x0 = np.array([-10.0, 10])
        accuracy = check_derivative(self.fun_vector_vector_with_arg,
                                    self.jac_vector_vector_with_arg,
                                    x0,
                                    kwargs={'arg': 42})
        assert_(accuracy < 1e-9)


class TestApproxDerivativeSparse:
    # Example from Numerical Optimization 2nd edition, p. 198.
    def setup_method(self):
        self.rng = np.random.default_rng(121091202)
        self.n = 50
        self.lb = -0.1 * (1 + np.arange(self.n))
        self.ub = 0.1 * (1 + np.arange(self.n))
        self.x0 = np.empty(self.n)
        self.x0[::2] = (1 - 1e-7) * self.lb[::2]
        self.x0[1::2] = (1 - 1e-7) * self.ub[1::2]

        self.J_true = self.jac(self.x0)

    def fun(self, x):
        e = x[1:]**3 - x[:-1]**2
        return np.hstack((0, 3 * e)) + np.hstack((2 * e, 0))

    def jac(self, x):
        n = x.size
        J = np.zeros((n, n))
        J[0, 0] = -4 * x[0]
        J[0, 1] = 6 * x[1]**2
        for i in range(1, n - 1):
            J[i, i - 1] = -6 * x[i-1]
            J[i, i] = 9 * x[i]**2 - 4 * x[i]
            J[i, i + 1] = 6 * x[i+1]**2
        J[-1, -1] = 9 * x[-1]**2
        J[-1, -2] = -6 * x[-2]

        return J

    def structure(self, n):
        A = np.zeros((n, n), dtype=int)
        A[0, 0] = 1
        A[0, 1] = 1
        for i in range(1, n - 1):
            A[i, i - 1: i + 2] = 1
        A[-1, -1] = 1
        A[-1, -2] = 1

        return A

    @pytest.mark.fail_slow(5)
    def test_all(self):
        A = self.structure(self.n)
        order = np.arange(self.n)
        groups_1 = group_columns(A, order)
        self.rng.shuffle(order)
        groups_2 = group_columns(A, order)

        with MapWrapper(1 if IS_WASM else 2) as mapper:
            for method, groups, l, u, mf in product(
                    ['2-point', '3-point', 'cs'], [groups_1, groups_2],
                    [-np.inf, self.lb], [np.inf, self.ub], [map, mapper]):
                J = approx_derivative(self.fun, self.x0, method=method,
                                      bounds=(l, u), sparsity=(A, groups),
                                      workers=mf)
                assert_(isinstance(J, csr_array))
                assert_allclose(J.toarray(), self.J_true, rtol=1e-6)

                rel_step = np.full_like(self.x0, 1e-8)
                rel_step[::2] *= -1
                J = approx_derivative(self.fun, self.x0, method=method,
                                      rel_step=rel_step, sparsity=(A, groups),
                                      workers=mf)
                assert_allclose(J.toarray(), self.J_true, rtol=1e-5)

    def test_no_precomputed_groups(self):
        A = self.structure(self.n)
        J = approx_derivative(self.fun, self.x0, sparsity=A)
        assert_allclose(J.toarray(), self.J_true, rtol=1e-6)

    def test_equivalence(self):
        structure = np.ones((self.n, self.n), dtype=int)
        groups = np.arange(self.n)
        for method in ['2-point', '3-point', 'cs']:
            J_dense = approx_derivative(self.fun, self.x0, method=method)
            J_sparse = approx_derivative(
                self.fun, self.x0, sparsity=(structure, groups), method=method)
            assert_allclose(J_dense, J_sparse.toarray(),
                            rtol=5e-16, atol=7e-15)

    def test_check_derivative(self):
        def jac(x):
            return csr_array(self.jac(x))

        accuracy = check_derivative(self.fun, jac, self.x0,
                                    bounds=(self.lb, self.ub))
        assert_(accuracy < 1e-9)

        accuracy = check_derivative(self.fun, jac, self.x0,
                                    bounds=(self.lb, self.ub))
        assert_(accuracy < 1e-9)


class TestApproxDerivativeLinearOperator(TestFunctions):

    @make_xp_test_case(
        approx_derivative
    )
    def test_scalar_scalar(self, xp):
        x0 = _item_for_scalar_function(xp.asarray(1.0, dtype=xp.float64))
        jac_diff_2 = approx_derivative(self.fun_scalar_scalar, x0,
                                       method='2-point',
                                       as_linear_operator=True)
        jac_diff_3 = approx_derivative(self.fun_scalar_scalar, x0,
                                       as_linear_operator=True)
        jac_diff_4 = approx_derivative(self.fun_scalar_scalar, x0,
                                       method='cs',
                                       as_linear_operator=True)
        jac_true = self.jac_scalar_scalar(x0)
        rng = np.random.default_rng(11290049580398)
        for i in range(10):
            p = rng.uniform(-10, 10, size=(1,))
            p = xp.asarray(p, dtype=x0.dtype)

            xp_assert_close(jac_diff_2.dot(p), jac_true*p,
                            rtol=1e-5)
            xp_assert_close(jac_diff_3.dot(p), jac_true*p,
                            rtol=5e-6)
            xp_assert_close(jac_diff_4.dot(p), jac_true*p,
                            rtol=5e-6)

    @make_xp_test_case(
        approx_derivative
    )
    def test_scalar_vector(self, xp):
        x0 = _item_for_scalar_function(xp.asarray(0.5, dtype=xp.float64))

        jac_diff_2 = approx_derivative(self.fun_scalar_vector, x0,
                                       method='2-point',
                                       as_linear_operator=True)
        jac_diff_3 = approx_derivative(self.fun_scalar_vector, x0,
                                       as_linear_operator=True)
        jac_diff_4 = approx_derivative(self.fun_scalar_vector, x0,
                                       method='cs',
                                       as_linear_operator=True)

        jac_true = self.jac_scalar_vector(
            xpx.atleast_nd(xp.asarray(x0, dtype=xp.float64), ndim=1, xp=xp)
        )
        rng = np.random.default_rng(11290049580398)
        for i in range(10):
            p = rng.uniform(-10, 10, size=(1,))
            p = xp.asarray(p, dtype=x0.dtype)
            xp_assert_close(jac_diff_2.dot(p), xp.tensordot(jac_true, p, axes=1),
                            rtol=1e-5)
            xp_assert_close(jac_diff_3.dot(p), xp.tensordot(jac_true, p, axes=1),
                            rtol=5e-6)
            xp_assert_close(jac_diff_4.dot(p), xp.tensordot(jac_true, p, axes=1),
                            rtol=5e-6)

    @make_xp_test_case(
        approx_derivative
    )
    def test_vector_scalar(self, xp):
        x0 = xp.asarray([100.0, -0.5], dtype=xp.float64)
        jac_diff_2 = approx_derivative(self.fun_vector_scalar, x0,
                                       method='2-point',
                                       as_linear_operator=True)
        jac_diff_3 = approx_derivative(self.fun_vector_scalar, x0,
                                       as_linear_operator=True)
        jac_diff_4 = approx_derivative(self.fun_vector_scalar, x0,
                                       method='cs',
                                       as_linear_operator=True)
        jac_true = xp.asarray(self.jac_vector_scalar([100.0, -0.5]))
        rng = np.random.default_rng(11290049580398)
        for i in range(10):
            p = rng.uniform(-10, 10, size=x0.shape)
            p = xp.asarray(p, dtype=x0.dtype)
            expected = xpx.atleast_nd(
                xp.tensordot(jac_true, p, axes=1), ndim=1, xp=xp
            )
            xp_assert_close(jac_diff_2.dot(p), expected,
                            rtol=1e-5)
            xp_assert_close(jac_diff_3.dot(p), expected,
                            rtol=5e-6)
            xp_assert_close(jac_diff_4.dot(p), expected,
                            rtol=1e-7)

    @make_xp_test_case(
        approx_derivative
    )
    def test_vector_vector(self, xp):
        x0 = xp.asarray([-100.0, 0.2], dtype=xp.float64)
        jac_diff_2 = approx_derivative(self.fun_vector_vector, x0,
                                       method='2-point',
                                       as_linear_operator=True)
        jac_diff_3 = approx_derivative(self.fun_vector_vector, x0,
                                       as_linear_operator=True)
        jac_diff_4 = approx_derivative(self.fun_vector_vector, x0,
                                       method='cs',
                                       as_linear_operator=True)
        jac_true = xp.asarray(self.jac_vector_vector([-100.0, 0.2]))
        rng = np.random.default_rng(11290049580398)
        for i in range(10):
            p = rng.uniform(-10, 10, size=x0.shape)
            p = xp.asarray(p, dtype=x0.dtype)
            expected = xpx.atleast_nd(
                xp.tensordot(jac_true, p, axes=1), ndim=1, xp=xp
            )
            x = jac_diff_2.dot(p)
            xp_assert_close(x, expected, rtol=1e-5)
            x = jac_diff_3.dot(p)
            xp_assert_close(x, expected, rtol=1e-6)
            x = jac_diff_4.dot(p)
            xp_assert_close(x, expected, rtol=1.5e-7)

    def test_exception(self):
        x0 = np.array([-100.0, 0.2])
        assert_raises(ValueError, approx_derivative,
                      self.fun_vector_vector, x0,
                      method='2-point', bounds=(1, np.inf))


def test_absolute_step_sign():
    # test for gh12487
    # if an absolute step is specified for 2-point differences make sure that
    # the side corresponds to the step. i.e. if step is positive then forward
    # differences should be used, if step is negative then backwards
    # differences should be used.

    # function has double discontinuity at x = [-1, -1]
    # first component is \/, second component is /\
    def f(x):
        return -np.abs(x[0] + 1) + np.abs(x[1] + 1)

    # check that the forward difference is used
    grad = approx_derivative(f, [-1, -1], method='2-point', abs_step=1e-8)
    assert_allclose(grad, [-1.0, 1.0])

    # check that the backwards difference is used
    grad = approx_derivative(f, [-1, -1], method='2-point', abs_step=-1e-8)
    assert_allclose(grad, [1.0, -1.0])

    # check that the forwards difference is used with a step for both
    # parameters
    grad = approx_derivative(
        f, [-1, -1], method='2-point', abs_step=[1e-8, 1e-8]
    )
    assert_allclose(grad, [-1.0, 1.0])

    # check that we can mix forward/backwards steps.
    grad = approx_derivative(
        f, [-1, -1], method='2-point', abs_step=[1e-8, -1e-8]
     )
    assert_allclose(grad, [-1.0, -1.0])
    grad = approx_derivative(
        f, [-1, -1], method='2-point', abs_step=[-1e-8, 1e-8]
    )
    assert_allclose(grad, [1.0, 1.0])

    # the forward step should reverse to a backwards step if it runs into a
    # bound
    # This is kind of tested in TestAdjustSchemeToBounds, but only for a lower level
    # function.
    grad = approx_derivative(
        f, [-1, -1], method='2-point', abs_step=1e-8,
        bounds=(-np.inf, -1)
    )
    assert_allclose(grad, [1.0, -1.0])

    grad = approx_derivative(
        f, [-1, -1], method='2-point', abs_step=-1e-8, bounds=(-1, np.inf)
    )
    assert_allclose(grad, [-1.0, 1.0])


@make_xp_test_case(
    _compute_absolute_step
)
def test__compute_absolute_step(xp):
    # tests calculation of absolute step from rel_step
    methods = ['2-point', '3-point', 'cs']

    x0 = xp.asarray([1e-5, 0, 1, 1e5], dtype=xp.float64)

    EPS = xp.finfo(xp.float64).eps
    relative_step = {
        "2-point": EPS**0.5,
        "3-point": EPS**(1/3),
        "cs": EPS**0.5
    }
    f0 = xp.asarray(1.0, dtype=xp.float64)

    for method in methods:
        rel_step = relative_step[method]
        correct_step = xp.stack(
            [xp.asarray(rel_step),
             xp.asarray(rel_step * 1.),
             xp.asarray(rel_step * 1.),
             xp.asarray(rel_step * xp.abs(x0[3]))],
        )

        abs_step = _compute_absolute_step(None, x0, f0, method, xp=xp)
        xp_assert_close(abs_step, correct_step)

        sign_x0 = xp.astype(-x0 >= 0, xp.float64) * 2 - 1
        abs_step = _compute_absolute_step(None, -x0, f0, method, xp=xp)
        xp_assert_close(abs_step, sign_x0 * correct_step)

    # if a relative step is provided it should be used
    rel_step = xp.asarray([0.1, 1, 10, 100], dtype=xp.float64)
    correct_step = xp.stack(
        [xp.asarray(rel_step[0] * x0[0]),
         xp.asarray(relative_step['2-point']),
         xp.asarray(rel_step[2] * 1.),
         xp.asarray(rel_step[3] * xp.abs(x0[3]))],
    )

    abs_step = _compute_absolute_step(rel_step, x0, f0, '2-point')
    xp_assert_close(abs_step, correct_step)

    sign_x0 = xp.astype(-x0 >= 0, xp.float64) * 2 - 1
    abs_step = _compute_absolute_step(rel_step, -x0, f0, '2-point')
    xp_assert_close(abs_step, sign_x0 * correct_step)

    # the dtype of absolute step should be the same as x0
    #def _compute_absolute_step(rel_step, x0, f0, method):
    x0 = xp.asarray([1e-5, 0, 1, 1e5], dtype=xp.float32)
    abs_step = _compute_absolute_step(None, x0, f0, '3-point')
    assert abs_step.dtype == xp.float32
    abs_step = _compute_absolute_step(None, x0, f0, '2-point')
    assert abs_step.dtype == xp.float32

    # check for 16-bit x0 and f0
    if hasattr(xp, "float16"):
        x0 = xp.asarray([0.1, 0, 1, 50], dtype=xp.float16)
        abs_step = _compute_absolute_step(None, x0, f0, '3-point')
        assert abs_step.dtype == xp.float16
        abs_step = _compute_absolute_step(None, x0, f0, '2-point')
        assert abs_step.dtype == xp.float16

        x0 = xp.asarray([1e-3, 0, 1, 10], dtype=xp.float64)
        f0 = xp.asarray(1.0, dtype=xp.float16)
        abs_step = _compute_absolute_step(rel_step, x0, f0, '2-point')
        assert abs_step.dtype == xp.float64

        x0 = xp.asarray([1e-5, 0, 1, 1e5], dtype=xp.float32)
        abs_step = _compute_absolute_step(rel_step, x0, f0, '2-point')
        assert abs_step.dtype == xp.float32

        x0 = xp.asarray([1e-3, 0, 1, 10], dtype=xp.float16)
        abs_step = _compute_absolute_step(rel_step, x0, f0, '2-point')
        assert abs_step.dtype == xp.float16
