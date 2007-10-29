from numpy.testing import *

import numpy
import scipy
from scipy import matrix,array,diag,zeros
from scipy.sparse import csr_matrix


set_package_path()
from scipy.sandbox.multigrid.utils import infinity_norm, diag_sparse, \
                                          expand_into_blocks
restore_path()


class TestUtils(NumpyTestCase):
    def check_infinity_norm(self):
        A = matrix([[-4]])
        assert_equal(infinity_norm(csr_matrix(A)),4)

        A = matrix([[1,0,-5],[-2,5,0]])
        assert_equal(infinity_norm(csr_matrix(A)),7)

        A = matrix([[0,1],[0,-5]])
        assert_equal(infinity_norm(csr_matrix(A)),5)

        A = matrix([[1.3,-4.7,0],[-2.23,5.5,0],[9,0,-2]])
        assert_equal(infinity_norm(csr_matrix(A)),11)

    def check_diag_sparse(self):
        #check sparse -> array
        A = matrix([[-4]])
        assert_equal(diag_sparse(csr_matrix(A)),[-4])

        A = matrix([[1,0,-5],[-2,5,0]])
        assert_equal(diag_sparse(csr_matrix(A)),[1,5])

        A = matrix([[0,1],[0,-5]])
        assert_equal(diag_sparse(csr_matrix(A)),[0,-5])

        A = matrix([[1.3,-4.7,0],[-2.23,5.5,0],[9,0,-2]])
        assert_equal(diag_sparse(csr_matrix(A)),[1.3,5.5,-2])

        #check array -> sparse
        A = matrix([[-4]])
        assert_equal(diag_sparse(array([-4])).todense(),csr_matrix(A).todense())

        A = matrix([[1,0],[0,5]])
        assert_equal(diag_sparse(array([1,5])).todense(),csr_matrix(A).todense())

        A = matrix([[0,0],[0,-5]])
        assert_equal(diag_sparse(array([0,-5])).todense(),csr_matrix(A).todense())

        A = matrix([[1.3,0,0],[0,5.5,0],[0,0,-2]])
        assert_equal(diag_sparse(array([1.3,5.5,-2])).todense(),csr_matrix(A).todense())

    def check_expand_into_blocks(self):
        cases = []
        cases.append( ( matrix([[1]]), (1,2) ) )
        cases.append( ( matrix([[1]]), (2,1) ) )
        cases.append( ( matrix([[1]]), (2,2) ) )
        cases.append( ( matrix([[1,2]]), (1,2) ) )
        cases.append( ( matrix([[1,2],[3,4]]), (2,2) ) )
        cases.append( ( matrix([[0,0],[0,0]]), (3,1) ) )
        cases.append( ( matrix([[0,1,0],[0,2,3]]), (3,2) ) )
        cases.append( ( matrix([[1,0,0],[2,0,3]]), (2,5) ) )

        for A,dims in cases:
            m,n = dims
            result = expand_into_blocks(csr_matrix(A),m,n).todense()

            expected = zeros((m*A.shape[0],n*A.shape[1]))
            for i in range(m):
                for j in range(n):
                    expected[i::m,j::n] = A

            assert_equal(expected,result)


if __name__ == '__main__':
    NumpyTest().run()
