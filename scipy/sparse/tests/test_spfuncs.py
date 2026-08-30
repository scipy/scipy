import numpy as np
from numpy.testing import assert_equal

from scipy.sparse import _spfuncs as spfuncs
from scipy.sparse import csr_array, csc_array, bsr_array
from scipy.sparse._sparsetools import (csr_scale_rows, csr_scale_columns,
                                       bsr_scale_rows, bsr_scale_columns)


class TestSparseFunctions:
    def test_scale_rows_and_cols(self):
        D = np.array([[1, 0, 0, 2, 3],
                   [0, 4, 0, 5, 0],
                   [0, 0, 6, 7, 0]])

        #TODO expose through function
        S = csr_array(D)
        v = np.array([1,2,3])
        csr_scale_rows(3,5,S.indptr,S.indices,S.data,v)
        assert_equal(S.toarray(), np.diag(v) @ D)

        S = csr_array(D)
        v = np.array([1,2,3,4,5])
        csr_scale_columns(3,5,S.indptr,S.indices,S.data,v)
        assert_equal(S.toarray(), D @ np.diag(v))

        # blocks
        E = np.kron(D,[[1,2],[3,4]])
        S = bsr_array(E,blocksize=(2,2))
        v = np.array([1,2,3,4,5,6])
        bsr_scale_rows(3,5,2,2,S.indptr,S.indices,S.data,v)
        assert_equal(S.toarray(), np.diag(v) @ E)

        S = bsr_array(E,blocksize=(2,2))
        v = np.array([1,2,3,4,5,6,7,8,9,10])
        bsr_scale_columns(3,5,2,2,S.indptr,S.indices,S.data,v)
        assert_equal(S.toarray(), E @ np.diag(v))

        E = np.kron(D,[[1,2,3],[4,5,6]])
        S = bsr_array(E,blocksize=(2,3))
        v = np.array([1,2,3,4,5,6])
        bsr_scale_rows(3,5,2,3,S.indptr,S.indices,S.data,v)
        assert_equal(S.toarray(), np.diag(v) @ E)

        S = bsr_array(E,blocksize=(2,3))
        v = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
        bsr_scale_columns(3,5,2,3,S.indptr,S.indices,S.data,v)
        assert_equal(S.toarray(), E @ np.diag(v))

    def test_estimate_blocksize(self):
        mats = []
        mats.append([[0,1],[1,0]])
        mats.append([[1,1,0],[0,0,1],[1,0,1]])
        mats.append([[0],[0],[1]])
        mats = [np.array(x) for x in mats]

        blks = []
        blks.append([[1]])
        blks.append([[1,1],[1,1]])
        blks.append([[1,1],[0,1]])
        blks.append([[1,1,0],[1,0,1],[1,1,1]])
        blks = [np.array(x) for x in blks]

        for A in mats:
            for B in blks:
                X = np.kron(A,B)
                r,c = spfuncs.estimate_blocksize(X)
                assert r >= B.shape[0]
                assert c >= B.shape[1]

    def test_count_blocks(self):
        def gold(A,bs):
            R,C = bs
            I,J = A.nonzero()
            return len(set(zip(I//R,J//C)))

        mats = []
        mats.append([[0]])
        mats.append([[1]])
        mats.append([[1,0]])
        mats.append([[1,1]])
        mats.append([[0,1],[1,0]])
        mats.append([[1,1,0],[0,0,1],[1,0,1]])
        mats.append([[0],[0],[1]])

        for A in mats:
            for B in mats:
                X = np.kron(A,B)
                Y = csr_array(X)
                for R in range(1,6):
                    for C in range(1,6):
                        assert_equal(spfuncs.count_blocks(Y, (R, C)), gold(X, (R, C)))

        X = np.kron([[1,1,0],[0,0,1],[1,0,1]],[[1,1]])
        Y = csc_array(X)
        assert_equal(spfuncs.count_blocks(X, (1, 2)), gold(X, (1, 2)))
        assert_equal(spfuncs.count_blocks(Y, (1, 2)), gold(X, (1, 2)))
