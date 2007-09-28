from numpy.testing import *

import numpy
import scipy
from scipy import matrix,array,diag
from scipy.sparse import csr_matrix


set_package_path()
from scipy.sandbox.multigrid.utils import infinity_norm,diag_sparse
restore_path()


class test_utils(NumpyTestCase):
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




if __name__ == '__main__':
    NumpyTest().run()

