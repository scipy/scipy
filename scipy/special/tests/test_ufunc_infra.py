import numpy as np
import pytest

from numpy.testing import assert_equal
from numpy._core._exceptions import UFuncTypeError
# ufunc using special._ufuncs_tools._with_cache_optimization
from scipy.special import mathieu_sem
# raw ufunc without cache optimization
from scipy.special._ufuncs import _mathieu_sem


# Tests that ufunc kwargs still work when _with_cache_optimization is applied
class TestWithCacheOptimization:
    @pytest.mark.parametrize("m_shape,q_shape,x_shape,where_shape", [
        ((5, 1), (1, 1), (1, 10), ()),
        ((2, 3), (2, 3), (2, 3), ()),
        ((1, 5), (1, 5), (10, 5), ()),
        ((5, 1), (1, 1), (1, 10), (1, 10)),
        ((2, 3), (2, 3), (2, 3), (2, 3)),
        ((1, 5), (1, 5), (10, 5), (10, 5)),
        ((5, 1), (1, 1), (1, 10), (5, 1)),
        ((2, 3), (2, 3), (2, 3), (2, 1)),
        ((1, 5), (1, 5), (10, 5), (1, 5)),
        ((5, 1), (5, 1), (5, 1), (3, 1, 1)),
        ((1, 10), (1, 10), (1, 10), (5, 1)),
        ((5,), (5,), (5,), (2, 3, 1, 1)),
        ((5, 1), (5, 1), (1, 1), (1, 10)),
        ((5, 1, 1), (1, 4, 1), (1, 1, 3), (2, 1, 1, 1)),
    ])
    @pytest.mark.parametrize("out_order", ["C", "F"])
    def test_out(self, m_shape, q_shape, x_shape, where_shape, out_order):
        # the call to mathieu_sem in this test raised a segfault on scipy/main prior to
        # https://github.com/scipy/scipy/pull/25127 -- https://github.com/scipy/xsf/pull/146.
        rng = np.random.default_rng(1234)

        m = rng.integers(1, 20, m_shape)
        q = rng.uniform(0, 10, q_shape)
        x = rng.uniform(0, 90, x_shape)

        if where_shape == ():
            batch_shape = np.broadcast_shapes(m_shape, q_shape, x_shape)
            where = True
        else:
            where = rng.choice([True, False], size=where_shape)
            batch_shape = np.broadcast_shapes(m_shape, q_shape, x_shape, where_shape)

        out0 = np.full(batch_shape, np.nan, order=out_order)
        out1 = np.full(batch_shape, np.nan, order=out_order)
        res0, res1 = mathieu_sem(m, q, x, out=(out0, out1), where=where)
        assert res0 is out0 and res1 is out1

        expected0 = np.full(batch_shape, np.nan, order=out_order)
        expected1 = np.full(batch_shape, np.nan, order=out_order)
        _mathieu_sem(m, q, x, out=(expected0, expected1), where=where)
        assert_equal((res0, res1), (expected0, expected1))
        assert res0.flags == expected0.flags
        assert res1.flags == expected1.flags

    def test_out_broadcasts(self):
        m = [[1, 3, 4]]
        q = [[2, 2, 2]]
        x = [[30, 60, 90]]
        out0 = np.full((3, 3), np.nan)
        out1 = np.full((3, 3), np.nan)
        mathieu_sem(m, q, x, out=(out0, out1))
        expected0 = np.full((3, 3), np.nan)
        expected1 = np.full((3, 3), np.nan)
        mathieu_sem(m, q, x, out=(expected0, expected1))
        assert_equal((out0, out1), (expected0, expected1))

    @pytest.mark.parametrize("subok", [True, False])
    def test_subok(self, subok):
        m = np.ma.masked_array([1, 2, 3])
        q = np.ma.masked_array([2.1, 3.2, 4.3])
        x = np.ma.masked_array([10, 20, 30])
        res0, res1 = mathieu_sem(m, q, x, subok=subok)
        expected0, expected1 = _mathieu_sem(m, q, x, subok=subok)
        assert type(res0) is type(expected0) and type(res1) is type(expected1)

    @pytest.mark.parametrize("order_m", ["C", "F"])
    @pytest.mark.parametrize("order_q", ["C", "F"])
    @pytest.mark.parametrize("order_x", ["C", "F"])
    # order="K" works differently with _with_cache_optimization
    @pytest.mark.parametrize("order", ["C", "F", "A"])
    @pytest.mark.parametrize("m_shape,q_shape,x_shape", [
        ((5, 1), (1, 1), (1, 10)),
        ((2, 3), (2, 3), (2, 3)),
        ((1, 5), (1, 5), (10, 5)),
        ((5, 1), (1, 1), (1, 10)),
        ((2, 3), (2, 3), (2, 3)),
        ((1, 5), (1, 5), (10, 5)),
        ((5, 1), (1, 1), (1, 10)),
        ((2, 3), (2, 3), (2, 3)),
        ((1, 5), (1, 5), (10, 5)),
        ((5, 1), (5, 1), (5, 1)),
        ((1, 10), (1, 10), (1, 10)),
        ((5,), (5,), (5,)),
        ((5, 1), (5, 1), (1, 1)),
        ((5, 1, 1), (1, 4, 1), (1, 1, 3)),
    ])
    def test_order(self,order_m, order_q, order_x, order, m_shape, q_shape, x_shape):
        rng = np.random.default_rng(1234)
        m = np.asarray(rng.integers(1, 20, m_shape), copy=True, order=order_m)
        q = np.asarray(rng.uniform(0, 10, q_shape), copy=True, order=order_q)
        x = np.asarray(rng.uniform(0, 90, x_shape), copy=True, order=order_x)
        res0, res1 = mathieu_sem(m, q, x, order=order)
        expected0, expected1 = _mathieu_sem(m, q, x, order=order)
        assert_equal((res0, res1), (expected0, expected1))

        assert res0.flags == expected0.flags
        assert res1.flags == expected1.flags

    def test_casting(self):
        m = np.float16(2.0)
        q = np.float16(1.0)
        x = np.float16(45.0)
        with pytest.raises(UFuncTypeError):
            res0, res1 = mathieu_sem(m, q, x, casting="no")

        res0, res1 = mathieu_sem(m, q, x, casting="safe")
        expected0, expected1 = _mathieu_sem(m, q, x, casting="safe")
        assert_equal((res0, res1), (expected0, expected1))

    @pytest.mark.parametrize("dtype_kwarg", [np.float32, np.float64])
    @pytest.mark.parametrize("input_dtype", [np.float16, np.float32, np.float64])
    def test_dtype(self, dtype_kwarg, input_dtype):
        m = np.asarray([1, 2, 3], dtype=input_dtype)
        q = np.asarray([2, 3, 4], dtype=input_dtype)
        x = np.asarray([30, 60, 90], dtype=input_dtype)

        res0, res1 = mathieu_sem(m, q, x, dtype=dtype_kwarg)
        expected0, expected1 = _mathieu_sem(m, q, x, dtype=dtype_kwarg)
        assert_equal((res0, res1), (expected0, expected1))
        assert res0.dtype == res1.dtype == expected0.dtype == expected1.dtype

    @pytest.mark.parametrize(
        "signature", [
            (np.float32, np.float32, np.float32, np.float32, np.float32),
            (np.float64, np.float64, np.float64, np.float64, np.float64),
        ],
    )
    def test_signature(self, signature):
        m = [1, 2, 3]
        q = [2, 3, 4]
        x = [30, 60, 90]
        res0, res1 = mathieu_sem(m, q, x, signature=signature)
        expected0, expected1 = _mathieu_sem(m, q, x, signature=signature)
        assert_equal((res0, res1), (expected0, expected1))
        assert res0.dtype == res1.dtype == expected0.dtype == expected1.dtype
