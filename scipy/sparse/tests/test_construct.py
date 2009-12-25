"""test sparse matrix construction functions"""

import numpy as np
from numpy import array, matrix
from numpy.testing import *


from scipy.sparse import csr_matrix, coo_matrix

from scipy.sparse.construct import *
from scipy.sparse.construct import rand as sprand

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
        assert_equal( I.dtype, np.dtype('int8') )
        assert_equal( I.format, 'dia' )

        for fmt in sparse_formats:
            I = identity( 3, format=fmt )
            assert_equal( I.format, fmt )
            assert_equal( I.toarray(), [[1,0,0],[0,1,0],[0,0,1]])

    def test_eye(self):
        assert_equal(eye(1,1).toarray(), [[1]])
        assert_equal(eye(2,3).toarray(), [[1,0,0],[0,1,0]])
        assert_equal(eye(3,2).toarray(), [[1,0],[0,1],[0,0]])
        assert_equal(eye(3,3).toarray(), [[1,0,0],[0,1,0],[0,0,1]])

        assert_equal(eye(3,3,dtype='int16').dtype, np.dtype('int16'))

        for m in [3, 5]:
            for n in [3, 5]:
                for k in range(-5,6):
                    assert_equal(eye(m, n, k=k).toarray(), np.eye(m, n, k=k))

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
                expected = np.kron(a,b)
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
                expected = np.kron(np.eye(len(b)), a) + \
                        np.kron(b, np.eye(len(a)))
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

    def test_rand(self):
        # Simple sanity checks for sparse.rand
        for t in [np.float32, np.float64, np.longdouble]:
            x = sprand(5, 10, density=0.1, dtype=t)
            assert_equal(x.dtype, t)
            assert_equal(x.shape, (5, 10))
            assert_equal(x.nonzero()[0].size, 5)

        x = sprand(5, 10, density=0.1)
        assert_equal(x.dtype, np.double)

        for fmt in ['coo', 'csc', 'csr', 'lil']:
            x = sprand(5, 10, format=fmt)
            assert_equal(x.format, fmt)

        assert_raises(ValueError, lambda: sprand(5, 10, 1.1))
        assert_raises(ValueError, lambda: sprand(5, 10, -0.1))

if __name__ == "__main__":
    run_module_suite()
