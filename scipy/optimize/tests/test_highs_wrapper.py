"""
Unit test for the highspy wrapper
"""

import numpy as np
import scipy.optimize._highspy._core as _h # type: ignore[import-not-found]
from scipy.optimize._highspy import _highs_wrapper
from numpy.testing import assert_allclose


def test_get_marg_bnd():
    col_dual = [1, 2, 3, 4, 5]
    # Make one of each possible value
    col_status = [
        _h.HighsBasisStatus.kLower,
        _h.HighsBasisStatus.kUpper,
        _h.HighsBasisStatus.kBasic,
        _h.HighsBasisStatus.kZero,
        _h.HighsBasisStatus.kNonbasic,
    ]
    desired = np.array([
        [1, 0, 0, 0, 0],
        [0, 2, 0, 0, 0],
    ])
    actual = _highs_wrapper._get_marg_bnds(col_dual, col_status)
    assert_allclose(actual, desired)

    # In actual code, this function gets called with lists.
    # But let's try arrays too.
    col_dual_np = np.array(col_dual)
    col_status_np = np.array(col_status)
    actual = _highs_wrapper._get_marg_bnds(col_dual_np, col_status_np)
    assert_allclose(actual, desired)
