"""test sparse matrix construction functions"""

from __future__ import division, print_function, absolute_import

from itertools import product

from numpy.testing import assert_equal
from scipy.sparse import csr_matrix

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
            for return_index, return_inverse, return_counts in (
                    product(*([[True, False]] * 3))):
                sparse_result = list(reversed(extract.unique(
                    A, return_index, return_inverse, return_counts)))
                np_result = list(reversed(np.unique(
                    B, return_index, return_inverse, return_counts)))

                sparse_uniques = sparse_result.pop()
                np_uniques = np_result.pop()
                assert_equal(sparse_uniques, np_uniques)

                if return_index: # Compare indices of first unique values
                    assert_equal(sparse_result.pop(), np_result.pop())

                if return_inverse:
                    # Special check for inverse indices, since they have
                    # different structure. Instead of comparing return values,
                    # compare the reconstructed matrices.
                    sparse_inverse = sparse_result.pop()
                    sparse_reconstructed = np.zeros(np.prod(A.shape))
                    sparse_reconstructed[sparse_inverse[0]] = (
                        sparse_uniques[sparse_inverse[1]])
                    sparse_reconstructed.reshape(A.shape)

                    np_inverse = np_result.pop()
                    np_reconstructed = np_uniques[np_inverse]

                    assert_equal(sparse_reconstructed, np_reconstructed)

                if return_counts: # Compare counts
                    assert_equal(sparse_result.pop(), np_result.pop())
