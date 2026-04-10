import numpy as np
import pytest

from scipy._lib._array_api import (
    is_numpy, make_xp_test_case, xp_assert_close, xp_assert_equal,
    lazy_xp_function
)
from scipy._external import array_api_extra as xpx
from scipy.optimize._pava import pava
from scipy.optimize import isotonic_regression

lazy_xp_function(isotonic_regression)
lazy_xp_function(pava, jax_jit=False)

@make_xp_test_case(isotonic_regression)
class TestIsotonicRegression:
    @pytest.mark.parametrize(
        ("y", "w", "msg"),
        [
            ([[0, 1]], None,
             "array has incorrect number of dimensions: 2; expected 1"),
            ([0, 1], [[1, 2]],
             "Input arrays y and w must have one dimension of equal length"),
            ([0, 1], [1],
             "Input arrays y and w must have one dimension of equal length"),
            (1, [1, 2],
             "Input arrays y and w must have one dimension of equal length"),
            ([1, 2], 1,
             "Input arrays y and w must have one dimension of equal length"),
            ([0, 1], [0, 1],
             "Weights w must be strictly positive"),
        ]
    )
    def test_raise_error(self, y, w, msg, xp):
        y, w = map(lambda x: xp.asarray(x) if x is not None else x, [y, w])
        with pytest.raises(ValueError, match=msg):
            isotonic_regression(y=y, weights=w)

    
    def test_simple_pava(self, xp):
        # Test case of Busing 2020
        # https://doi.org/10.18637/jss.v102.c01
        y = xp.asarray([8, 4, 8, 2, 2, 0, 8], dtype=xp.float64)
        w = xp.ones_like(y)
        dtype = np.intp if is_numpy(xp) else xp.int64
        r = xp.full(shape=y.shape[0] + 1, fill_value=-1, dtype=dtype)
        pava(y, w, r)
        xp_assert_close(y, xp.asarray([4, 4, 4, 4, 4, 4, 8], dtype=xp.float64))
        # Only first 2 elements of w are changed.
        xp_assert_close(w, xp.asarray([6, 1, 1, 1, 1, 1, 1], dtype=xp.float64))
        # Only first 3 elements of r are changed.
        xp_assert_close(r, xp.asarray([0, 6, 7, -1, -1, -1, -1, -1], dtype=dtype))

    @pytest.mark.parametrize("y_dtype", [np.float64, np.float32, np.int64, np.int32])
    @pytest.mark.parametrize("w_dtype", [np.float64, np.float32, np.int64, np.int32])
    @pytest.mark.parametrize("w", [None, "ones"])
    def test_simple_isotonic_regression(self, w, w_dtype, y_dtype, xp):
        # Test case of Busing 2020
        # https://doi.org/10.18637/jss.v102.c01
        y = xp.asarray(np.array([8, 4, 8, 2, 2, 0, 8], dtype=y_dtype))
        if w is not None:
            w = xp.asarray(np.ones_like(y, dtype=w_dtype))
        res = isotonic_regression(y, weights=w)
        xp_assert_close(res.x, xp.asarray([4, 4, 4, 4, 4, 4, 8], dtype=xp.float64))
        xp_assert_close(res.weights, xp.asarray([6, 1], dtype=xp.float64))
        blocks_dtype = xpx.default_dtype(xp, kind="integral")
        xp_assert_close(res.blocks, xp.asarray([0, 6, 7], dtype=blocks_dtype))
        # Assert that y was not overwritten
        xp_assert_equal(y, xp.asarray(np.asarray([8, 4, 8, 2, 2, 0, 8], y_dtype)))

    @pytest.mark.parametrize("increasing", [True, False])
    def test_linspace(self, increasing, xp):
        n = 10
        y = xp.linspace(0, 1, n) if increasing else xp.linspace(1, 0, n)
        res = isotonic_regression(y, increasing=increasing)
        xp_assert_close(res.x, xp.astype(y, xp.float64))
        xp_assert_close(res.blocks, xp.arange(n + 1))

    def test_weights(self, xp):
        w = xp.asarray([1, 2, 5, 0.5, 0.5, 0.5, 1, 3])
        y = xp.asarray([3, 2, 1, 10, 9, 8, 20, 10])
        res = isotonic_regression(y, weights=w)
        xp_assert_close(
            res.x,
            xp.asarray([12/8, 12/8, 12/8, 9, 9, 9, 50/4, 50/4], dtype=xp.float64)
        )
        xp_assert_close(res.weights, xp.asarray([8, 1.5, 4], dtype=xp.float64))
        blocks_dtype = xpx.default_dtype(xp, kind="integral")
        xp_assert_close(res.blocks, xp.asarray([0, 3, 6, 8], dtype=blocks_dtype))

        # weights are like repeated observations, we repeat the 3rd element 5
        # times.
        w2 = xp.asarray([1, 2, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 1, 3])
        y2 = xp.asarray([3, 2, 1, 1, 1, 1, 1, 10, 9, 8, 20, 10])
        res2 = isotonic_regression(y2, weights=w2)
        xp_assert_close(
            xp.diff(res2.x[0:7]), xp.asarray(0),
            check_dtype=False, check_shape=False
        )
        xp_assert_close(res2.x[4:], res.x)
        xp_assert_close(res2.weights, res.weights)
        xp_assert_close(res2.blocks[1:] - 4, res.blocks[1:])

    def test_against_R_monotone(self, xp):
        y = xp.asarray([0, 6, 8, 3, 5, 2, 1, 7, 9, 4])
        res = isotonic_regression(y)
        # R code
        # library(monotone)
        # options(digits=8)
        # monotone(c(0, 6, 8, 3, 5, 2, 1, 7, 9, 4))
        x_R = xp.asarray([
            0, 4.1666667, 4.1666667, 4.1666667, 4.1666667, 4.1666667,
            4.1666667, 6.6666667, 6.6666667, 6.6666667,
        ], dtype=xp.float64)
        xp_assert_close(res.x, x_R)
        xp_assert_equal(res.blocks, xp.asarray([0, 1, 7, 10]))

        n = 100
        y = xp.linspace(0, 1, num=n, endpoint=False)
        y = 5 * y + xp.sin(10 * y)
        res = isotonic_regression(y)
        # R code
        # library(monotone)
        # n <- 100
        # y <- 5 * ((1:n)-1)/n + sin(10 * ((1:n)-1)/n)
        # options(digits=8)
        # monotone(y)
        x_R = xp.asarray([
            0.00000000, 0.14983342, 0.29866933, 0.44552021, 0.58941834, 0.72942554,
            0.86464247, 0.99421769, 1.11735609, 1.23332691, 1.34147098, 1.44120736,
            1.53203909, 1.57081100, 1.57081100, 1.57081100, 1.57081100, 1.57081100,
            1.57081100, 1.57081100, 1.57081100, 1.57081100, 1.57081100, 1.57081100,
            1.57081100, 1.57081100, 1.57081100, 1.57081100, 1.57081100, 1.57081100,
            1.57081100, 1.57081100, 1.57081100, 1.57081100, 1.57081100, 1.57081100,
            1.57081100, 1.57081100, 1.57081100, 1.57081100, 1.57081100, 1.57081100,
            1.57081100, 1.57081100, 1.57081100, 1.57081100, 1.57081100, 1.57081100,
            1.57081100, 1.57081100, 1.57081100, 1.62418532, 1.71654534, 1.81773256,
            1.92723551, 2.04445967, 2.16873336, 2.29931446, 2.43539782, 2.57612334,
            2.72058450, 2.86783750, 3.01691060, 3.16681390, 3.31654920, 3.46511999,
            3.61154136, 3.75484992, 3.89411335, 4.02843976, 4.15698660, 4.27896904,
            4.39366786, 4.50043662, 4.59870810, 4.68799998, 4.76791967, 4.83816823,
            4.86564130, 4.86564130, 4.86564130, 4.86564130, 4.86564130, 4.86564130,
            4.86564130, 4.86564130, 4.86564130, 4.86564130, 4.86564130, 4.86564130,
            4.86564130, 4.86564130, 4.86564130, 4.86564130, 4.86564130, 4.86564130,
            4.86564130, 4.86564130, 4.86564130, 4.86564130,
        ], dtype=xp.float64)
        atol = 10 * xp.finfo(y.dtype).eps
        xp_assert_close(res.x, x_R, atol=atol)

        # Test increasing
        assert xp.all(xp.diff(res.x) >= 0)

        # Test balance property: sum(y) == sum(x)
        xp_assert_close(xp.sum(res.x), xp.sum(y, dtype=xp.float64))

        # Reverse order
        res_inv = isotonic_regression(-y, increasing=False)
        xp_assert_close(-res_inv.x, res.x)
        xp_assert_equal(res_inv.blocks, res.blocks)

    def test_readonly(self, xp):
        x = xp.arange(3, dtype=xp.float64)
        w = xp.ones(3, dtype=xp.float64)

        if is_numpy(xp):
            x.flags.writeable = False
            w.flags.writeable = False

        res = isotonic_regression(x, weights=w)
        assert xp.all(xp.isfinite(res.x))
        assert xp.all(xp.isfinite(res.weights))
        assert xp.all(xp.isfinite(res.blocks))

    def test_non_contiguous_arrays(self, xp):
        x = xp.arange(10, dtype=xp.float64)[::3]
        w = xp.ones(10, dtype=xp.float64)[::3]
        if is_numpy(xp):
            assert not x.flags.c_contiguous
            assert not x.flags.f_contiguous
            assert not w.flags.c_contiguous
            assert not w.flags.f_contiguous

        res = isotonic_regression(x, weights=w)
        assert xp.all(xp.isfinite(res.x))
        assert xp.all(xp.isfinite(res.weights))
        assert xp.all(xp.isfinite(res.blocks))
