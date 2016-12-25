from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_almost_equal, run_module_suite

from scipy.sparse.linalg import minres


def test_minres():
    a = np.random.randn(25).reshape(5,5)
    b = np.random.randn(5)
    c = np.random.randn(5)
    x1 = minres(a, b)[0] + c
    x2 = minres(a, b, x0=c)[0]
    assert_almost_equal(x1, x2)

if __name__ == "__main__":
    # Comment out the next line to run unit tests only
    run_module_suite()
