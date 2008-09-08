"""test sparse matrix construction functions"""

import numpy
from numpy import array, matrix
from numpy.testing import *


from scipy.sparse import csr_matrix, coo_matrix

from scipy.sparse.construct import *

sparse_formats = ['csr','csc','coo','bsr','dia','lil','dok']

#TODO check whether format=XXX is respected

class TestConstructUtils(TestCase):
    def test_spdiags(self):
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


    def test_identity(self):
        assert_equal(identity(1).toarray(), [[1]])
        assert_equal(identity(2).toarray(), [[1,0],[0,1]])

        I = identity(3, dtype='int8', format='dia')
        assert_equal( I.dtype, 'int8' )
        assert_equal( I.format, 'dia' )

        for fmt in sparse_formats:
            I = identity( 3, format=fmt )
            assert_equal( I.format, fmt )
            assert_equal( I.toarray(), [[1,0,0],[0,1,0],[0,0,1]])

    def test_eye(self):
        a = eye(2, 3 )
        b = array([[1, 0, 0], [0, 1, 0]], dtype='d')
        assert_array_equal(a.toarray(), b)

        a = eye(3, 2)
        b = array([[1, 0], [0, 1], [0, 0]], dtype='d')
        assert_array_equal( a.toarray(), b)

        a = eye(3, 3)
        b = array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype='d')
        assert_array_equal(a.toarray(), b)

    def test_kron(self):
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
                result   = kron(csr_matrix(a),csr_matrix(b)).todense()
                expected = numpy.kron(a,b)
                assert_array_equal(result,expected)

    def test_kronsum(self):
        cases = []

        cases.append(array([[ 0]]))
        cases.append(array([[-1]]))
        cases.append(array([[ 4]]))
        cases.append(array([[10]]))
        cases.append(array([[1,2],[3,4]]))
        cases.append(array([[0,2],[5,0]]))
        cases.append(array([[0,2,-6],[8,0,14],[0,3,0]]))
        cases.append(array([[1,0,0],[0,5,-1],[4,-2,8]]))

        for a in cases:
            for b in cases:
                result   = kronsum(csr_matrix(a),csr_matrix(b)).todense()
                expected = numpy.kron(numpy.eye(len(b)), a) + \
                        numpy.kron(b, numpy.eye(len(a)))
                assert_array_equal(result,expected)

    def test_vstack(self):

        A = coo_matrix([[1,2],[3,4]])
        B = coo_matrix([[5,6]])

        expected = matrix([[1, 2],
                           [3, 4],
                           [5, 6]])
        assert_equal( vstack( [A,B] ).todense(), expected )

    def test_hstack(self):

        A = coo_matrix([[1,2],[3,4]])
        B = coo_matrix([[5],[6]])

        expected = matrix([[1, 2, 5],
                           [3, 4, 6]])
        assert_equal( hstack( [A,B] ).todense(), expected )

    def test_bmat(self):

        A = coo_matrix([[1,2],[3,4]])
        B = coo_matrix([[5],[6]])
        C = coo_matrix([[7]])

        expected = matrix([[1, 2, 5],
                           [3, 4, 6],
                           [0, 0, 7]])
        assert_equal( bmat( [[A,B],[None,C]] ).todense(), expected )


        expected = matrix([[1, 2, 0],
                           [3, 4, 0],
                           [0, 0, 7]])
        assert_equal( bmat( [[A,None],[None,C]] ).todense(), expected )

        expected = matrix([[0, 5],
                           [0, 6],
                           [7, 0]])
        assert_equal( bmat( [[None,B],[C,None]] ).todense(), expected )

        #TODO test failure cases

    def test_lil_diags(self):
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
    run_module_suite()
