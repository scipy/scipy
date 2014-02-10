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
                
    def test_norm_axis(self):
        a = np.array([[ 1, 2, 3],
                      [-1, 1, 4]])
                
        c = csr_matrix(a)        
        #Frobenius norm
        assert_equal(norm(c, axis=0), np.sqrt(np.power(np.asmatrix(a), 2).sum(axis=0)))
        assert_equal(norm(c, axis=1), np.sqrt(np.power(np.asmatrix(a), 2).sum(axis=1)))
        
        assert_equal(norm(c, np.inf, axis=0), max(abs(np.asmatrix(a)).sum(axis=0)))
        assert_equal(norm(c, np.inf, axis=1), max(abs(np.asmatrix(a)).sum(axis=1)))
                
        assert_equal(norm(c, -np.inf, axis=0), min(abs(np.asmatrix(a)).sum(axis=0)))
        assert_equal(norm(c, -np.inf, axis=1), min(abs(np.asmatrix(a)).sum(axis=1)))
                        
        assert_equal(norm(c, 1, axis=0), abs(np.asmatrix(a)).sum(axis=0))
        assert_equal(norm(c, 1, axis=1), abs(np.asmatrix(a)).sum(axis=1))
        
        assert_equal(norm(c, -1, axis=0), min(abs(np.asmatrix(a)).sum(axis=0))  )
        assert_equal(norm(c, -1, axis=1), min(abs(np.asmatrix(a)).sum(axis=1))  )
        
         #_multi_svd_norm is not implemented for sparse matrix
        assert_raises(NotImplementedError, norm, c, 2, 0)
        #assert_raises(NotImplementedError, norm, c, -2, 0)
