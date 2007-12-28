from numpy import array, kron
from numpy.testing import *

set_package_path()
from scipy.sparse.spfuncs import *
from scipy.sparse import csr_matrix, csc_matrix
restore_path()

class TestSparseFunctions(NumpyTestCase):
    def check_estimate_blocksize(self):

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

    def check_count_blocks(self):
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
    NumpyTest().run()

