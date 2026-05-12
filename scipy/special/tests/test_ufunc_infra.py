import numpy as np
import pytest

from numpy.testing import assert_equal
# ufunc using special._ufuncs_tools._with_cache_optimization
from scipy.special import mathieu_sem
# raw ufunc without cache optimization
from scipy.special._ufuncs import _mathieu_sem


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
def test_out(m_shape, q_shape, x_shape, where_shape, out_order):
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


@pytest.mark.parametrize("subok", [True, False])
def test_subok(subok):
    m = np.ma.masked_array([1, 2, 3])
    q = np.ma.masked_array([2.1, 3.2, 4.3])
    x = np.ma.masked_array([10, 20, 30])
    res0, res1 = mathieu_sem(m, q, x, subok=subok)
    expected0, expected1 = _mathieu_sem(m, q, x, subok=subok)
    assert type(res0) == type(res1) == type(expected0) == type(expected1)


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
def test_order(order_m, order_q, order_x, order, m_shape, q_shape, x_shape):
    rng = np.random.default_rng(1234)
    m = np.asarray(rng.integers(1, 20, m_shape), copy=True, order=order_m)
    q = np.asarray(rng.uniform(0, 10, q_shape), copy=True, order=order_q)
    x = np.asarray(rng.uniform(0, 90, x_shape), copy=True, order=order_x)
    res0, res1 = mathieu_sem(m, q, x, order=order)
    expected0, expected1 = _mathieu_sem(m, q, x, order=order)
    assert_equal((res0, res1), (expected0, expected1))
    
    assert res0.flags == expected0.flags
    assert res1.flags == expected1.flags    
