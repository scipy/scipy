from numpy import array, kron, diag, matrix
from scipy.testing import *

from scipy.sparse.spfuncs import *
from scipy.sparse import csr_matrix, csc_matrix, bsr_matrix
from scipy.sparse.sparsetools import csr_scale_rows, csr_scale_columns, \
        bsr_scale_rows, bsr_scale_columns

class TestSparseFunctions(TestCase):
    def test_scale_rows_and_cols(self):
        D = matrix([[1,0,0,2,3],
                    [0,4,0,5,0],
                    [0,0,6,7,0]])
        
        
        #TODO expose through function
        S = csr_matrix(D)
        v = array([1,2,3])
        csr_scale_rows(3,5,S.indptr,S.indices,S.data,v)
        assert_equal(S.todense(), diag(v)*D )

        S = csr_matrix(D)
        v = array([1,2,3,4,5])
        csr_scale_columns(3,5,S.indptr,S.indices,S.data,v)
        assert_equal(S.todense(), D*diag(v) )

        # blocks
        E = kron(D,[[1,2],[3,4]])
        S = bsr_matrix(E,blocksize=(2,2))
        v = array([1,2,3,4,5,6])
        bsr_scale_rows(3,5,2,2,S.indptr,S.indices,S.data,v)
        assert_equal(S.todense(), diag(v)*E )

        S = bsr_matrix(E,blocksize=(2,2))
        v = array([1,2,3,4,5,6,7,8,9,10])
        bsr_scale_columns(3,5,2,2,S.indptr,S.indices,S.data,v)
        assert_equal(S.todense(), E*diag(v) )
        
        E = kron(D,[[1,2,3],[4,5,6]])
        S = bsr_matrix(E,blocksize=(2,3))
        v = array([1,2,3,4,5,6])
        bsr_scale_rows(3,5,2,3,S.indptr,S.indices,S.data,v)
        assert_equal(S.todense(), diag(v)*E )

        S = bsr_matrix(E,blocksize=(2,3))
        v = array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
        bsr_scale_columns(3,5,2,3,S.indptr,S.indices,S.data,v)
        assert_equal(S.todense(), E*diag(v) )



    def test_extract_diagonal(self):
        mats = []
        mats.append( [[1,0,2]] )
        mats.append( [[1],[0],[2]] )
        mats.append( [[0,1],[0,2],[0,3]] )
        mats.append( [[0,0,1],[0,0,2],[0,3,0]] )

        mats.append( kron(mats[0],[[1,2]]) )
        mats.append( kron(mats[0],[[1],[2]]) )
        mats.append( kron(mats[1],[[1,2],[3,4]]) )
        mats.append( kron(mats[2],[[1,2],[3,4]]) )
        mats.append( kron(mats[3],[[1,2],[3,4]]) )

        for m in mats:
            expected = diag(m)
            assert_equal(extract_diagonal(m),expected)
            assert_equal(extract_diagonal(csr_matrix(m)),expected)
            assert_equal(extract_diagonal(csc_matrix(m)),expected)
        
        for m in mats:
            m = array(m)
            M,N = m.shape
            expected = diag(m)
            for R in range(1,M+1):
                for C in range(1,N+1):
                    if M % R == 0 and N % C == 0:
                        result = extract_diagonal( bsr_matrix(m,blocksize=(R,C)) )
                        assert_equal(result,expected)
            

    def test_estimate_blocksize(self):
        mats = []
        mats.append( [[0,1],[1,0]] )
        mats.append( [[1,1,0],[0,0,1],[1,0,1]] )
        mats.append( [[0],[0],[1]] )
        mats = [array(x) for x in mats]
    
        blks = []
        blks.append( [[1]] )
        blks.append( [[1,1],[1,1]] )
        blks.append( [[1,1],[0,1]] )
        blks.append( [[1,1,0],[1,0,1],[1,1,1]] )
        blks = [array(x) for x in blks]

        for A in mats:
            for B in blks:
                X = kron(A,B)
                r,c = estimate_blocksize(X)
                assert(r >= B.shape[0])
                assert(c >= B.shape[1])

    def test_count_blocks(self):
        def gold(A,bs):
            R,C = bs
            I,J = A.nonzero()
            return len( set( zip(I/R,J/C) ) )
        
        mats = []
        mats.append( [[0]] ) 
        mats.append( [[1]] ) 
        mats.append( [[1,0]] ) 
        mats.append( [[1,1]] ) 
        mats.append( [[0,1],[1,0]] )
        mats.append( [[1,1,0],[0,0,1],[1,0,1]] )
        mats.append( [[0],[0],[1]] )

        for A in mats:
            for B in mats:
                X = kron(A,B)
                Y = csr_matrix(X)
                for R in range(1,6):
                    for C in range(1,6):
                        assert_equal(count_blocks(Y,(R,C)),gold(X,(R,C)))
        
        X = kron([[1,1,0],[0,0,1],[1,0,1]],[[1,1]])
        Y = csc_matrix(X)
        assert_equal(count_blocks(X,(1,2)),gold(X,(1,2)))
        assert_equal(count_blocks(Y,(1,2)),gold(X,(1,2)))


if __name__ == "__main__":
    unittests.main()

