"""test sparse matrix construction functions"""

import numpy
from numpy import array
from numpy.testing import *

set_package_path()
from scipy.sparse import csc_matrix, csr_matrix, dok_matrix, coo_matrix, \
     spidentity, speye, spkron, extract_diagonal, lil_matrix, lil_eye, \
     lil_diags, spdiags
from scipy.linsolve import splu
restore_path()

class TestConstructUtils(NumpyTestCase):
    def check_identity(self):
        a = spidentity(3)
        b = array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype='d')
        assert_array_equal(a.toarray(), b)

    def check_eye(self):
        a = speye(2, 3 )
#        print a, a.__repr__
        b = array([[1, 0, 0], [0, 1, 0]], dtype='d')
        assert_array_equal(a.toarray(), b)

        a = speye(3, 2)
#        print a, a.__repr__
        b = array([[1, 0], [0, 1], [0, 0]], dtype='d')
        assert_array_equal( a.toarray(), b)

        a = speye(3, 3)
#        print a, a.__repr__
        b = array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype='d')
        assert_array_equal(a.toarray(), b)

    def check_spkron(self):
        from numpy import kron

        cases = []

        cases.append(array([[ 0]]))
        cases.append(array([[-1]]))
        cases.append(array([[ 4]]))
        cases.append(array([[10]]))
        cases.append(array([[0],[0]]))
        cases.append(array([[0,0]]))
        cases.append(array([[1,2],[3,4]]))
        cases.append(array([[0,2],[5,0]]))
        cases.append(array([[0,2,-6],[8,0,14]]))
        cases.append(array([[5,4],[0,0],[6,0]]))
        cases.append(array([[5,4,4],[1,0,0],[6,0,8]]))
        cases.append(array([[0,1,0,2,0,5,8]]))
        cases.append(array([[0.5,0.125,0,3.25],[0,2.5,0,0]]))
        
        for a in cases:
            for b in cases:
                result = spkron(csr_matrix(a),csr_matrix(b)).todense()
                expected = kron(a,b)

                assert_array_equal(result,expected)

    def check_lil_diags(self):
        assert_array_equal(lil_diags([[1,2,3],[4,5],[6]],
                                     [0,1,2],(3,3)).todense(),
                           [[1,4,6],
                            [0,2,5],
                            [0,0,3]])

        assert_array_equal(lil_diags([[6],[4,5],[1,2,3]],
                                     [2,1,0],(3,3)).todense(),
                           [[1,4,6],
                            [0,2,5],
                            [0,0,3]])

        assert_array_equal(lil_diags([[6,7,8],[4,5],[1,2,3]],
                                     [2,1,0],(3,3)).todense(),
                           [[1,4,6],
                            [0,2,5],
                            [0,0,3]])

        assert_array_equal(lil_diags([[1,2,3],[4,5],[6]],
                                     [0,-1,-2],(3,3)).todense(),
                           [[1,0,0],
                            [4,2,0],
                            [6,5,3]])

        assert_array_equal(lil_diags([[6,7,8],[4,5]],
                                     [-2,-1],(3,3)).todense(),
                           [[0,0,0],
                            [4,0,0],
                            [6,5,0]])

if __name__ == "__main__":
    NumpyTest().run()

