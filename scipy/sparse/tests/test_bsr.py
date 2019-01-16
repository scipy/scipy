from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_equal, assert_array_equal
from scipy.sparse import bsr_matrix

def test_bsr_eliminate_zeros():
    N = 12
    np.random.seed(0)
    m = bsr_matrix(np.random.random((N, N)), blocksize=(2, 3))

    # eliminate some blocks, but not all
    m.data[m.data <= 0.9] = 0
    m.eliminate_zeros()
    assert_equal(m.nnz, 66)
    assert_array_equal(m.data.shape, (11, 2, 3))

    # eliminate all remaining blocks (github issue #9687)
    m.data[m.data <= 1.0] = 0
    m.eliminate_zeros()
    assert_equal(m.nnz, 0)
    assert_array_equal(m.data.shape, (0, 2, 3))
    assert_array_equal(m.todense(), np.zeros((12,12)))

    # test fast path
    m.eliminate_zeros()
    assert_equal(m.nnz, 0)
    assert_array_equal(m.data.shape, (0, 2, 3))
    assert_array_equal(m.todense(), np.zeros((12,12)))

