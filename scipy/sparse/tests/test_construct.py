"""test sparse matrix construction functions"""

from numpy import array, kron
from numpy.testing import *

set_package_path()
from scipy.sparse import csr_matrix, \
     spidentity, speye, spkron, spdiags, \
     lil_eye, lil_diags
restore_path()

class TestConstructUtils(NumpyTestCase):
    def check_spdiags(self):
        diags1 = array( [[ 1, 2, 3, 4, 5]] )
        diags2 = array( [[ 1, 2, 3, 4, 5],
                         [ 6, 7, 8, 9,10]] )
        diags3 = array( [[ 1, 2, 3, 4, 5],
                         [ 6, 7, 8, 9,10],
                         [11,12,13,14,15]] )
       
        cases = []
        cases.append( (diags1,  0,  1, 1, [[1]]) )
        cases.append( (diags1, [0], 1, 1, [[1]]) )
        cases.append( (diags1, [0], 2, 1, [[1],[0]]) )
        cases.append( (diags1, [0], 1, 2, [[1,0]]) )
        cases.append( (diags1, [1], 1, 2, [[0,2]]) )
        cases.append( (diags1,[-1], 1, 2, [[0,0]]) )
        cases.append( (diags1, [0], 2, 2, [[1,0],[0,2]]) )
        cases.append( (diags1,[-1], 2, 2, [[0,0],[1,0]]) )
        cases.append( (diags1, [3], 2, 2, [[0,0],[0,0]]) )
        cases.append( (diags1, [0], 3, 4, [[1,0,0,0],[0,2,0,0],[0,0,3,0]]) )
        cases.append( (diags1, [1], 3, 4, [[0,2,0,0],[0,0,3,0],[0,0,0,4]]) )
        cases.append( (diags1, [2], 3, 5, [[0,0,3,0,0],[0,0,0,4,0],[0,0,0,0,5]]) )

        cases.append( (diags2, [0,2], 3, 3, [[1,0,8],[0,2,0],[0,0,3]]) )
        cases.append( (diags2, [-1,0], 3, 4, [[6,0,0,0],[1,7,0,0],[0,2,8,0]]) )
        cases.append( (diags2, [2,-3], 6, 6, [[0,0,3,0,0,0],
                                              [0,0,0,4,0,0],
                                              [0,0,0,0,5,0],
                                              [6,0,0,0,0,0],
                                              [0,7,0,0,0,0],
                                              [0,0,8,0,0,0]]) )

        cases.append( (diags3, [-1,0,1], 6, 6, [[ 6,12, 0, 0, 0, 0],
                                                [ 1, 7,13, 0, 0, 0],
                                                [ 0, 2, 8,14, 0, 0],
                                                [ 0, 0, 3, 9,15, 0],
                                                [ 0, 0, 0, 4,10, 0],
                                                [ 0, 0, 0, 0, 5, 0]]) )
        cases.append( (diags3, [-4,2,-1], 6, 5, [[ 0, 0, 8, 0, 0],
                                                 [11, 0, 0, 9, 0],
                                                 [ 0,12, 0, 0,10],
                                                 [ 0, 0,13, 0, 0],
                                                 [ 1, 0, 0,14, 0],
                                                 [ 0, 2, 0, 0,15]]) )

        for d,o,m,n,result in cases:
            assert_equal( spdiags(d,o,m,n).todense(), result )
        
           
    def check_identity(self):
        a = spidentity(3)
        b = array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype='d')
        assert_array_equal(a.toarray(), b)

    def check_eye(self):
        a = speye(2, 3 )
        b = array([[1, 0, 0], [0, 1, 0]], dtype='d')
        assert_array_equal(a.toarray(), b)

        a = speye(3, 2)
        b = array([[1, 0], [0, 1], [0, 0]], dtype='d')
        assert_array_equal( a.toarray(), b)

        a = speye(3, 3)
        b = array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype='d')
        assert_array_equal(a.toarray(), b)

    def check_spkron(self):
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

