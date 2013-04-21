"""test sparse matrix construction functions"""

from __future__ import division, print_function, absolute_import

import numpy as np
from numpy import array, matrix
from numpy.testing import TestCase, run_module_suite, assert_equal, \
        assert_array_equal, assert_raises, assert_array_almost_equal_nulp

from scipy.sparse import csr_matrix, coo_matrix

from scipy.sparse import construct
from scipy.sparse.construct import rand as sprand

sparse_formats = ['csr','csc','coo','bsr','dia','lil','dok']

#TODO check whether format=XXX is respected


class TestConstructUtils(TestCase):
    def test_spdiags(self):
        diags1 = array([[1, 2, 3, 4, 5]])
        diags2 = array([[1, 2, 3, 4, 5],
                         [6, 7, 8, 9,10]])
        diags3 = array([[1, 2, 3, 4, 5],
                         [6, 7, 8, 9,10],
                         [11,12,13,14,15]])

        cases = []
        cases.append((diags1, 0, 1, 1, [[1]]))
        cases.append((diags1, [0], 1, 1, [[1]]))
        cases.append((diags1, [0], 2, 1, [[1],[0]]))
        cases.append((diags1, [0], 1, 2, [[1,0]]))
        cases.append((diags1, [1], 1, 2, [[0,2]]))
        cases.append((diags1,[-1], 1, 2, [[0,0]]))
        cases.append((diags1, [0], 2, 2, [[1,0],[0,2]]))
        cases.append((diags1,[-1], 2, 2, [[0,0],[1,0]]))
        cases.append((diags1, [3], 2, 2, [[0,0],[0,0]]))
        cases.append((diags1, [0], 3, 4, [[1,0,0,0],[0,2,0,0],[0,0,3,0]]))
        cases.append((diags1, [1], 3, 4, [[0,2,0,0],[0,0,3,0],[0,0,0,4]]))
        cases.append((diags1, [2], 3, 5, [[0,0,3,0,0],[0,0,0,4,0],[0,0,0,0,5]]))

        cases.append((diags2, [0,2], 3, 3, [[1,0,8],[0,2,0],[0,0,3]]))
        cases.append((diags2, [-1,0], 3, 4, [[6,0,0,0],[1,7,0,0],[0,2,8,0]]))
        cases.append((diags2, [2,-3], 6, 6, [[0,0,3,0,0,0],
                                              [0,0,0,4,0,0],
                                              [0,0,0,0,5,0],
                                              [6,0,0,0,0,0],
                                              [0,7,0,0,0,0],
                                              [0,0,8,0,0,0]]))

        cases.append((diags3, [-1,0,1], 6, 6, [[6,12, 0, 0, 0, 0],
                                                [1, 7,13, 0, 0, 0],
                                                [0, 2, 8,14, 0, 0],
                                                [0, 0, 3, 9,15, 0],
                                                [0, 0, 0, 4,10, 0],
                                                [0, 0, 0, 0, 5, 0]]))
        cases.append((diags3, [-4,2,-1], 6, 5, [[0, 0, 8, 0, 0],
                                                 [11, 0, 0, 9, 0],
                                                 [0,12, 0, 0,10],
                                                 [0, 0,13, 0, 0],
                                                 [1, 0, 0,14, 0],
                                                 [0, 2, 0, 0,15]]))

        for d,o,m,n,result in cases:
            assert_equal(construct.spdiags(d,o,m,n).todense(), result)

    def test_diags(self):
        a = array([1, 2, 3, 4, 5])
        b = array([6, 7, 8, 9, 10])
        c = array([11, 12, 13, 14, 15])

        cases = []
        cases.append((a[:1], 0, (1, 1), [[1]]))
        cases.append(([a[:1]], [0], (1, 1), [[1]]))
        cases.append(([a[:1]], [0], (2, 1), [[1],[0]]))
        cases.append(([a[:1]], [0], (1, 2), [[1,0]]))
        cases.append(([a[:1]], [1], (1, 2), [[0,1]]))
        cases.append(([a[:2]], [0], (2, 2), [[1,0],[0,2]]))
        cases.append(([a[:1]],[-1], (2, 2), [[0,0],[1,0]]))
        cases.append(([a[:3]], [0], (3, 4), [[1,0,0,0],[0,2,0,0],[0,0,3,0]]))
        cases.append(([a[:3]], [1], (3, 4), [[0,1,0,0],[0,0,2,0],[0,0,0,3]]))
        cases.append(([a[:3]], [2], (3, 5), [[0,0,1,0,0],[0,0,0,2,0],[0,0,0,0,3]]))

        cases.append(([a[:3],b[:1]], [0,2], (3, 3), [[1,0,6],[0,2,0],[0,0,3]]))
        cases.append(([a[:2],b[:3]], [-1,0], (3, 4), [[6,0,0,0],[1,7,0,0],[0,2,8,0]]))
        cases.append(([a[:4],b[:3]], [2,-3], (6, 6), [[0,0,1,0,0,0],
                                                     [0,0,0,2,0,0],
                                                     [0,0,0,0,3,0],
                                                     [6,0,0,0,0,4],
                                                     [0,7,0,0,0,0],
                                                     [0,0,8,0,0,0]]))

        cases.append(([a[:4],b,c[:4]], [-1,0,1], (5, 5), [[6,11, 0, 0, 0],
                                                            [1, 7,12, 0, 0],
                                                            [0, 2, 8,13, 0],
                                                            [0, 0, 3, 9,14],
                                                            [0, 0, 0, 4,10]]))
        cases.append(([a[:2],b[:3],c], [-4,2,-1], (6, 5), [[0, 0, 6, 0, 0],
                                                          [11, 0, 0, 7, 0],
                                                          [0,12, 0, 0, 8],
                                                          [0, 0,13, 0, 0],
                                                          [1, 0, 0,14, 0],
                                                          [0, 2, 0, 0,15]]))

        # scalar case: broadcasting
        cases.append(([1,-2,1], [1,0,-1], (3, 3), [[-2, 1, 0],
                                                    [1, -2, 1],
                                                    [0, 1, -2]]))

        for d, o, shape, result in cases:
            try:
                assert_equal(construct.diags(d, o, shape=shape).todense(),
                             result)

                if shape[0] == shape[1] and hasattr(d[0], '__len__'):
                    # should be able to find the shape automatically
                    assert_equal(construct.diags(d, o).todense(), result)
            except:
                print("%r %r %r" % (d, o, shape))
                raise

    def test_diags_bad(self):
        a = array([1, 2, 3, 4, 5])
        b = array([6, 7, 8, 9, 10])
        c = array([11, 12, 13, 14, 15])

        cases = []
        cases.append(([a[:0]], 0, (1, 1)))
        cases.append(([a], [0], (1, 1)))
        cases.append(([a[:3],b], [0,2], (3, 3)))
        cases.append(([a[:4],b,c[:3]], [-1,0,1], (5, 5)))
        cases.append(([a[:2],c,b[:3]], [-4,2,-1], (6, 5)))
        cases.append(([a[:2],c,b[:3]], [-4,2,-1], None))
        cases.append(([], [-4,2,-1], None))
        cases.append(([1], [-4], (4, 4)))
        cases.append(([a[:0]], [-1], (1, 2)))
        cases.append(([a], 0, None))

        for d, o, shape in cases:
            try:
                assert_raises(ValueError, construct.diags, d, o, shape)
            except:
                print("%r %r %r" % (d, o, shape))
                raise

        assert_raises(TypeError, construct.diags, [[None]], [0])

    def test_diags_vs_diag(self):
        # Check that
        #
        #    diags([a, b, ...], [i, j, ...]) == diag(a, i) + diag(b, j) + ...
        #

        np.random.seed(1234)

        for n_diags in [1, 2, 3, 4, 5, 10]:
            n = 1 + n_diags//2 + np.random.randint(0, 10)

            offsets = np.arange(-n+1, n-1)
            np.random.shuffle(offsets)
            offsets = offsets[:n_diags]

            diagonals = [np.random.rand(n - abs(q)) for q in offsets]

            mat = construct.diags(diagonals, offsets)
            dense_mat = sum([np.diag(x, j) for x, j in zip(diagonals, offsets)])

            assert_array_almost_equal_nulp(mat.todense(), dense_mat)

            if len(offsets) == 1:
                mat = construct.diags(diagonals[0], offsets[0])
                dense_mat = np.diag(diagonals[0], offsets[0])
                assert_array_almost_equal_nulp(mat.todense(), dense_mat)

    def test_diags_dtype(self):
        x = construct.diags([2.2], [0], shape=(2, 2), dtype=int)
        assert_equal(x.dtype, int)
        assert_equal(x.todense(), [[2, 0], [0, 2]])

    def test_diags_one_diagonal(self):
        d = list(range(5))
        for k in range(-5, 6):
            assert_equal(construct.diags(d, k).toarray(),
                         construct.diags([d], [k]).toarray())

    def test_identity(self):
        assert_equal(construct.identity(1).toarray(), [[1]])
        assert_equal(construct.identity(2).toarray(), [[1,0],[0,1]])

        I = construct.identity(3, dtype='int8', format='dia')
        assert_equal(I.dtype, np.dtype('int8'))
        assert_equal(I.format, 'dia')

        for fmt in sparse_formats:
            I = construct.identity(3, format=fmt)
            assert_equal(I.format, fmt)
            assert_equal(I.toarray(), [[1,0,0],[0,1,0],[0,0,1]])

    def test_eye(self):
        assert_equal(construct.eye(1,1).toarray(), [[1]])
        assert_equal(construct.eye(2,3).toarray(), [[1,0,0],[0,1,0]])
        assert_equal(construct.eye(3,2).toarray(), [[1,0],[0,1],[0,0]])
        assert_equal(construct.eye(3,3).toarray(), [[1,0,0],[0,1,0],[0,0,1]])

        assert_equal(construct.eye(3,3,dtype='int16').dtype, np.dtype('int16'))

        for m in [3, 5]:
            for n in [3, 5]:
                for k in range(-5,6):
                    assert_equal(construct.eye(m, n, k=k).toarray(), np.eye(m, n, k=k))
                    if m == n:
                        assert_equal(construct.eye(m, k=k).toarray(), np.eye(m, n, k=k))

    def test_eye_one(self):
        assert_equal(construct.eye(1).toarray(), [[1]])
        assert_equal(construct.eye(2).toarray(), [[1,0],[0,1]])

        I = construct.eye(3, dtype='int8', format='dia')
        assert_equal(I.dtype, np.dtype('int8'))
        assert_equal(I.format, 'dia')

        for fmt in sparse_formats:
            I = construct.eye(3, format=fmt)
            assert_equal(I.format, fmt)
            assert_equal(I.toarray(), [[1,0,0],[0,1,0],[0,0,1]])

    def test_kron(self):
        cases = []

        cases.append(array([[0]]))
        cases.append(array([[-1]]))
        cases.append(array([[4]]))
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
                result = construct.kron(csr_matrix(a),csr_matrix(b)).todense()
                expected = np.kron(a,b)
                assert_array_equal(result,expected)

    def test_kronsum(self):
        cases = []

        cases.append(array([[0]]))
        cases.append(array([[-1]]))
        cases.append(array([[4]]))
        cases.append(array([[10]]))
        cases.append(array([[1,2],[3,4]]))
        cases.append(array([[0,2],[5,0]]))
        cases.append(array([[0,2,-6],[8,0,14],[0,3,0]]))
        cases.append(array([[1,0,0],[0,5,-1],[4,-2,8]]))

        for a in cases:
            for b in cases:
                result = construct.kronsum(csr_matrix(a),csr_matrix(b)).todense()
                expected = np.kron(np.eye(len(b)), a) + \
                        np.kron(b, np.eye(len(a)))
                assert_array_equal(result,expected)

    def test_vstack(self):

        A = coo_matrix([[1,2],[3,4]])
        B = coo_matrix([[5,6]])

        expected = matrix([[1, 2],
                           [3, 4],
                           [5, 6]])
        assert_equal(construct.vstack([A,B]).todense(), expected)

    def test_hstack(self):

        A = coo_matrix([[1,2],[3,4]])
        B = coo_matrix([[5],[6]])

        expected = matrix([[1, 2, 5],
                           [3, 4, 6]])
        assert_equal(construct.hstack([A,B]).todense(), expected)

    def test_bmat(self):

        A = coo_matrix([[1,2],[3,4]])
        B = coo_matrix([[5],[6]])
        C = coo_matrix([[7]])

        expected = matrix([[1, 2, 5],
                           [3, 4, 6],
                           [0, 0, 7]])
        assert_equal(construct.bmat([[A,B],[None,C]]).todense(), expected)

        expected = matrix([[1, 2, 0],
                           [3, 4, 0],
                           [0, 0, 7]])
        assert_equal(construct.bmat([[A,None],[None,C]]).todense(), expected)

        expected = matrix([[0, 5],
                           [0, 6],
                           [7, 0]])
        assert_equal(construct.bmat([[None,B],[C,None]]).todense(), expected)

        #TODO test failure cases

    def test_block_diag_basic(self):
        """ basic test for block_diag """
        A = coo_matrix([[1,2],[3,4]])
        B = coo_matrix([[5],[6]])
        C = coo_matrix([[7]])

        expected = matrix([[1, 2, 0, 0],
                           [3, 4, 0, 0],
                           [0, 0, 5, 0],
                           [0, 0, 6, 0],
                           [0, 0, 0, 7]])

        assert_equal(construct.block_diag((A, B, C)).todense(), expected)

    def test_block_diag_scalar_1d_args(self):
        """ block_diag with scalar and 1d arguments """
        # one 1d matrix and a scalar
        assert_array_equal(construct.block_diag([[2,3], 4]).toarray(),
                           [[2, 3, 0], [0, 0, 4]])

    def test_block_diag_1(self):
        """ block_diag with one matrix """
        assert_equal(construct.block_diag([[1, 0]]).todense(),
                     matrix([[1, 0]]))
        assert_equal(construct.block_diag([[[1, 0]]]).todense(),
                     matrix([[1, 0]]))
        assert_equal(construct.block_diag([[[1], [0]]]).todense(),
                     matrix([[1], [0]]))
        # just on scalar
        assert_equal(construct.block_diag([1]).todense(),
                     matrix([[1]]))

    def test_rand(self):
        # Simple sanity checks for sparse.rand
        for t in [np.float32, np.float64, np.longdouble]:
            x = sprand(5, 10, density=0.1, dtype=t)
            assert_equal(x.dtype, t)
            assert_equal(x.shape, (5, 10))
            assert_equal(x.nonzero()[0].size, 5)

        x1 = sprand(5, 10, density=0.1,
                    random_state=4321)
        assert_equal(x1.dtype, np.double)

        x2 = sprand(5, 10, density=0.1,
                    random_state=np.random.RandomState(4321))

        assert_array_equal(x1.data, x2.data)
        assert_array_equal(x1.row, x2.row)
        assert_array_equal(x1.col, x2.col)

        for fmt in ['coo', 'csc', 'csr', 'lil']:
            x = sprand(5, 10, format=fmt)
            assert_equal(x.format, fmt)

        assert_raises(ValueError, lambda: sprand(5, 10, 1.1))
        assert_raises(ValueError, lambda: sprand(5, 10, -0.1))


if __name__ == "__main__":
    run_module_suite()
