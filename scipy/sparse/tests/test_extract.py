"""test sparse matrix construction functions"""

from __future__ import division, print_function, absolute_import

from numpy.testing import TestCase, assert_equal
from scipy.sparse import csr_matrix

import numpy as np
from scipy.sparse import extract


class TestExtract(TestCase):
    def setUp(self):
        cases = []

        cases.append(csr_matrix([[1,2]]))
        cases.append(csr_matrix([[1,0]]))
        cases.append(csr_matrix([[0,0]]))
        cases.append(csr_matrix([[1],[2]]))
        cases.append(csr_matrix([[1],[0]]))
        cases.append(csr_matrix([[0],[0]]))
        cases.append(csr_matrix([[1,2],[3,4]]))
        cases.append(csr_matrix([[0,1],[0,0]]))
        cases.append(csr_matrix([[0,0],[1,0]]))
        cases.append(csr_matrix([[0,0],[0,0]]))
        cases.append(csr_matrix([[1,2,0,0,3],[4,5,0,6,7],[0,0,8,9,0]]))
        cases.append(csr_matrix([[1,2,0,0,3],[4,5,0,6,7],[0,0,8,9,0]]).T)

        self.cases = cases

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
