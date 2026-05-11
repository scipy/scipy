import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_equal

from scipy.special import mathieu_sem
from scipy.special._ufuncs import _mathieu_sem


class TestMathieuCacheOptimization:
    def _make_buffer(self, shape, order, contig):
        if contig:
            return np.full(shape, np.nan, order=order)
        return np.full((shape[0]*2, *shape[1:]), np.nan, order=order)[::2]


    @pytest.mark.parametrize("order_out0", ["C", "F"])
    @pytest.mark.parametrize("order_out1", ["C", "F"])
    @pytest.mark.parametrize("contig_out0", [True, False])
    @pytest.mark.parametrize("contig_out1", [True, False])
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
        # Cases where broadcasting with where changes the out shape.
        ((5, 1), (5, 1), (5, 1), (3, 1, 1)),
        ((1, 10), (1, 10), (1, 10), (5, 1)),
        ((5,), (5,), (5,), (2, 3, 1, 1)),
        ((5, 1), (5, 1), (1, 1), (1, 10)),
        ((5, 1, 1), (1, 4, 1), (1, 1, 3), (2, 1, 1, 1)),
    ])
    # order should not impact results when out is passed.
    @pytest.mark.parametrize("order", ["A", "C", "F", "K"])
    def test_out(
            self,
            order_out0,
            order_out1,
            contig_out0,
            contig_out1,
            m_shape,
            q_shape,
            x_shape,
            where_shape,
            order,
    ):
        rng = np.random.default_rng(1234)

        m = rng.uniform(1, 5, m_shape)
        q = rng.uniform(0, 10, q_shape)
        x = rng.uniform(0, 90, x_shape)

        if where_shape == ():
            batch_shape = np.broadcast_shapes(m_shape, q_shape, x_shape)
            where = True
        else:
            where = rng.choice([True, False], size=where_shape)
            batch_shape = np.broadcast_shapes(m_shape, q_shape, x_shape, where_shape)

        out0 = self._make_buffer(batch_shape, order_out0, contig_out0)
        out1 = self._make_buffer(batch_shape, order_out1, contig_out1)
        res0, res1 = mathieu_sem(m, q, x, out=(out0, out1), order=order)
        assert res0 is out0 and res1 is out1

        expected0 = self._make_buffer(batch_shape, order_out0, contig_out0)
        expected1 = self._make_buffer(batch_shape, order_out1, contig_out1)
        _mathieu_sem(m, q, x, out=(expected0, expected1), order=order)
        assert_equal((out0, out1), (expected0, expected1))
        assert res0.flags == expected0.flags
        assert res1.flags == expected1.flags
