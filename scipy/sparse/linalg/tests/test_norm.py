"""Test functions for the sparse.linalg.norm module
"""

from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.linalg import norm as npnorm
from numpy.testing import (assert_raises, assert_equal, assert_allclose,
        TestCase)

import scipy.sparse
from scipy.sparse.linalg import norm as spnorm


class TestNorm(TestCase):
    def test_norm(self):
        a = np.arange(9) - 4
        b = a.reshape((3, 3))
        b = scipy.sparse.csr_matrix(b)

        #Frobenius norm is the default
        assert_equal(spnorm(b), 7.745966692414834)        
        assert_equal(spnorm(b, 'fro'), 7.745966692414834)

        assert_equal(spnorm(b, np.inf), 9)
        assert_equal(spnorm(b, -np.inf), 2)
        assert_equal(spnorm(b, 1), 7)
        assert_equal(spnorm(b, -1), 6)

        #_multi_svd_norm is not implemented for sparse matrix
        assert_raises(NotImplementedError, spnorm, b, 2)
        assert_raises(NotImplementedError, spnorm, b, -2)


class TestVsNumpyNorm(TestCase):
    _sparse_types = (
            scipy.sparse.bsr_matrix,
            scipy.sparse.coo_matrix,
            scipy.sparse.csc_matrix,
            scipy.sparse.csr_matrix,
            scipy.sparse.dia_matrix,
            scipy.sparse.dok_matrix,
            scipy.sparse.lil_matrix,
            )
    _test_matrices = (
            (np.arange(9) - 4).reshape((3, 3)),
            [
                [ 1, 2, 3],
                [-1, 1, 4]],
            [
                [ 1, 0, 3],
                [-1, 1, 4j]],
            )
    def test_sparse_matrix_norms(self):
        for sparse_type in self._sparse_types:
            for M in self._test_matrices:
                S = sparse_type(M)
                assert_allclose(spnorm(S), npnorm(M))
                assert_allclose(spnorm(S, 'fro'), npnorm(M, 'fro'))
                assert_allclose(spnorm(S, np.inf), npnorm(M, np.inf))
                assert_allclose(spnorm(S, -np.inf), npnorm(M, -np.inf))
                assert_allclose(spnorm(S, 1), npnorm(M, 1))
                assert_allclose(spnorm(S, -1), npnorm(M, -1))
