from __future__ import division, print_function, absolute_import

import scipy.special as sc
import numpy as np
from numpy.testing import assert_almost_equal


def test_softmax_fixtures():
    assert_almost_equal(sc.softmax([1000, 0, 0, 0]), np.array([1, 0, 0, 0]))
    assert_almost_equal(sc.softmax([1, 1]), np.array([.5, .5]))


def test_softmax_multi_axes():
    assert_almost_equal(sc.softmax([[1000, 0], [1000, 0]], axis=0), np.array([[.5, .5], [.5, .5]]))
    assert_almost_equal(sc.softmax([[1000, 0], [1000, 0]], axis=1), np.array([[1, 0], [1, 0]]))


def test_softmax_dims():
    assert sc.softmax([1, 0]).ndim == 1
    assert sc.softmax([[1, 0], [2, 3]]).ndim == 2
