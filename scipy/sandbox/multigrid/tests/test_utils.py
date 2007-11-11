from numpy.testing import *

import numpy
import scipy
from numpy import matrix,array,diag,zeros,sqrt
from scipy.sparse import csr_matrix


set_package_path()
from scipy.sandbox.multigrid.utils import infinity_norm, diag_sparse, \
                                          symmetric_rescaling, \
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


    def check_symmetric_rescaling(self):
        cases = []
        cases.append( diag_sparse(array([1,2,3,4])) )
        cases.append( diag_sparse(array([1,0,3,4])) )
        
        A = array([ [ 5.5,  3.5,  4.8],
                    [ 2. ,  9.9,  0.5],
                    [ 6.5,  2.6,  5.7]])
        A = csr_matrix( A )
        cases.append(A)
        P = diag_sparse([1,0,1])
        cases.append( P*A*P )
        P = diag_sparse([0,1,0])
        cases.append( P*A*P )
        P = diag_sparse([1,-1,1])
        cases.append( P*A*P )

        for A in cases:
            D_sqrt,D_sqrt_inv,DAD = symmetric_rescaling(A)

            assert_almost_equal( diag_sparse(A) > 0, diag_sparse(DAD) )
            assert_almost_equal( diag_sparse(DAD), D_sqrt*D_sqrt_inv )

            D_sqrt,D_sqrt_inv = diag_sparse(D_sqrt),diag_sparse(D_sqrt_inv)
            assert_almost_equal((D_sqrt_inv*A*D_sqrt_inv).todense(), DAD.todense())

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
