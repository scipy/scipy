"""Test functions for the sparse.linalg.norm module
"""

from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (assert_raises, assert_allclose, assert_equal, assert_,
        decorators, TestCase, run_module_suite)

from numpy import matrix, abs
from scipy.sparse import *
from scipy.sparse.linalg import norm

class TestNorm(TestCase):
    def test_norm(self):
        a = np.arange(9) - 4
        b = a.reshape((3, 3))
        b = csr_matrix(b)
        
        #Frobenius norm is the default
        assert_equal(norm(b), 7.745966692414834)        
        assert_equal(norm(b, 'fro'), 7.745966692414834)
        
        assert_equal(norm(b, np.inf), 9)
        assert_equal(norm(b, -np.inf), 2)
        assert_equal(norm(b, 1), 7)
        assert_equal(norm(b, -1), 6)
        
        #_multi_svd_norm is not implemented for sparse matrix
        assert_raises(NotImplementedError, norm, b, 2)
        assert_raises(NotImplementedError, norm, b, -2)
