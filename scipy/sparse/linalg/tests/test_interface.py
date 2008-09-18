#!/usr/bin/env python
""" Test functions for the sparse.linalg.interface module
"""

from numpy.testing import *

import numpy
from numpy import array, matrix, dtype
from scipy.sparse import csr_matrix

from scipy.sparse.linalg.interface import *


class TestInterface(TestCase):
    def test_aslinearoperator(self):
        cases = []

        cases.append( matrix([[1,2,3],[4,5,6]]) )
        cases.append( array([[1,2,3],[4,5,6]]) )
        cases.append( csr_matrix([[1,2,3],[4,5,6]]) )

        class matlike:
            def __init__(self):
                self.dtype = dtype('int')
                self.shape = (2,3)
            def matvec(self,x):
                y = array([ 1*x[0] + 2*x[1] + 3*x[2],
                            4*x[0] + 5*x[1] + 6*x[2]])
                if len(x.shape) == 2:
                    y = y.reshape(-1,1)
                return y

            def rmatvec(self,x):
                return array([ 1*x[0] + 4*x[1],
                               2*x[0] + 5*x[1],
                               3*x[0] + 6*x[1]])

        cases.append( matlike() )


        for M in cases:
            A = aslinearoperator(M)
            M,N = A.shape

            assert_equal(A.matvec(array([1,2,3])),       [14,32])
            assert_equal(A.matvec(array([[1],[2],[3]])), [[14],[32]])

            assert_equal(A.rmatvec(array([1,2])),     [9,12,15])
            assert_equal(A.rmatvec(array([[1],[2]])), [[9],[12],[15]])

            assert_equal(A.matmat(array([[1,4],[2,5],[3,6]])), \
                    [[14,32],[32,77]] )

            if hasattr(M,'dtype'):
                assert_equal(A.dtype, M.dtype)

if __name__ == "__main__":
    nose.run(argv=['', __file__])
