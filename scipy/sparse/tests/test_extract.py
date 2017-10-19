"""test sparse matrix construction functions"""

from __future__ import division, print_function, absolute_import

from itertools import product

from numpy.testing import assert_equal
from scipy.sparse import csr_matrix
from scipy._lib._version import NumpyVersion

import numpy as np
from scipy.sparse import extract, spmatrix


class TestExtract(object):
    def setup_method(self):
        self.cases = [
            csr_matrix([[1,2]]),
            csr_matrix([[1,0]]),
            csr_matrix([[0,0]]),
            csr_matrix([[1],[2]]),
            csr_matrix([[1],[0]]),
            csr_matrix([[0],[0]]),
            csr_matrix([[1,2],[3,4]]),
            csr_matrix([[0,1],[0,0]]),
            csr_matrix([[0,0],[1,0]]),
            csr_matrix([[0,0],[0,0]]),
            csr_matrix([[1,2,0,0,3],[4,5,0,6,7],[0,0,8,9,0]]),
            csr_matrix([[1,2,0,0,3],[4,5,0,6,7],[0,0,8,9,0]]).T,
            csr_matrix([[0,2,2,0,3],[4,4,0,2,7],[0,0,4,2,0]]),
        ]

    def find(self):
        for A in self.cases:
            I,J,V = extract.find(A)
            assert_equal(A.toarray(), csr_matrix(((I,J),V), shape=A.shape))

    def test_tril(self):
        for A in self.cases:
            B = A.toarray()
            for k in [-3,-2,-1,0,1,2,3]:
                assert_equal(extract.tril(A,k=k).toarray(), np.tril(B,k=k))

    def test_triu(self):
        for A in self.cases:
            B = A.toarray()
            for k in [-3,-2,-1,0,1,2,3]:
                assert_equal(extract.triu(A,k=k).toarray(), np.triu(B,k=k))

    def test_unique(self):
        for A in self.cases:
            B = A.toarray()
            if NumpyVersion(np.__version__) <= '1.8.2':
                # return_counts is not compatible with these versions of Numpy.
                return_counts_options = [False]
            else:
                return_counts_options = [True, False]

            for return_indices, return_inverse, return_counts in (
                    product(*([[True, False]] * 2 + return_counts_options))):
                sparse_result = extract.unique(
                    A, return_indices, return_inverse, return_counts)
                np_result = np.unique(
                    B, return_indices, return_inverse, return_counts)
                if not isinstance(sparse_result, tuple):  # Just the uniques
                    assert_equal(sparse_result, np_result)
                    continue

                sparse_result = list(reversed(sparse_result))
                np_result = list(reversed(np_result))

                sparse_uniques = sparse_result.pop()
                np_uniques = np_result.pop()
                assert_equal(sparse_uniques, np_uniques)

                if return_indices:  # Compare indices of first unique values
                    assert_equal(sparse_result.pop(), np_result.pop())

                if return_inverse:
                    # Special check for inverse indices, since they have
                    # different structures. Instead of comparing return values,
                    # compare the reconstructed matrices.
                    sparse_inv = sparse_result.pop()
                    sparse_recon = np.zeros(np.prod(A.shape), dtype=B.dtype)
                    sparse_recon[sparse_inv[0]] = (
                        sparse_uniques[sparse_inv[1]])  # don't bother w reshape

                    np_inv = np_result.pop()
                    np_recon = np_uniques[np_inv]

                    assert_equal(sparse_recon, np_recon)

                if return_counts:  # Compare counts
                    assert_equal(sparse_result.pop(), np_result.pop())
