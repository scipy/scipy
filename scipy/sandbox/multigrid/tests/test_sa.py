from numpy.testing import *

from numpy import sqrt,empty,ones,arange,array_split
from scipy import rand
from scipy.sparse import spdiags,csr_matrix,lil_matrix
import numpy

set_package_path()
import scipy.sandbox.multigrid
from scipy.sandbox.multigrid.sa import sa_strong_connections, sa_constant_interpolation, \
                                        sa_interpolation, sa_fit_candidates
from scipy.sandbox.multigrid.multilevel import poisson_problem1D,poisson_problem2D
restore_path()


def reference_sa_strong_connections(A,epsilon):
    A_coo = A.tocoo()
    S = lil_matrix(A.shape)

    for (i,j,v) in zip(A_coo.row,A_coo.col,A_coo.data):
        if i == j: continue #skip diagonal

        if abs(A[i,j]) >= epsilon*sqrt(abs(A[i,i])*abs(A[j,j])):
            S[i,j] = v

    return S.tocsr()

class TestSAStrongConnections(NumpyTestCase):
    def check_simple(self):
        N = 4
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).tocsr()
        S = spdiags([ -ones(N),-ones(N)],[-1,1],N,N).tocsr()
        assert_array_equal(sa_strong_connections(A,0.50).todense(),S.todense())   #all connections are strong
        assert_array_equal(sa_strong_connections(A,0.51).todense(),0*S.todense()) #no connections are strong
       
        N = 100
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).tocsr()
        S = spdiags([ -ones(N),-ones(N)],[-1,1],N,N).tocsr()
        assert_array_equal(sa_strong_connections(A,0.50).todense(),S.todense())   #all connections are strong
        assert_array_equal(sa_strong_connections(A,0.51).todense(),0*S.todense()) #no connections are strong

    def check_random(self):
        numpy.random.seed(0)

        for N in [2,3,5]:
            A = csr_matrix(rand(N,N))
            for epsilon in [0.0,0.1,0.5,1.0,10.0]:
                S_result = sa_strong_connections(A,epsilon)
                S_expected = reference_sa_strong_connections(A,epsilon)
                assert_array_equal(S_result.todense(),S_expected.todense())

    def check_poisson1D(self):
        for N in [2,3,5,7,10,11,19]:
            A = poisson_problem1D(N)
            for epsilon in [0.0,0.1,0.5,1.0]:
                S_result   = sa_strong_connections(A,epsilon)
                S_expected = reference_sa_strong_connections(A,epsilon)
                assert_array_equal(S_result.todense(),S_expected.todense())

    def check_poisson2D(self):
        for N in [2,3,5,7,10,11]:
            A = poisson_problem2D(N)
            for epsilon in [0.0,0.1,0.5,1.0]:
                S_result   = sa_strong_connections(A,epsilon)
                S_expected = reference_sa_strong_connections(A,epsilon)
                assert_array_equal(S_result.todense(),S_expected.todense())


# note that this method only tests the current implementation, not
# all possible implementations
def reference_sa_constant_interpolation(A,epsilon):
    S = sa_strong_connections(A,epsilon)
    S = array_split(S.indices,S.indptr[1:-1])

    n = A.shape[0]

    R = set(range(n))
    j = 0

    aggregates = empty(n,dtype=A.indices.dtype)
    aggregates[:] = -1

    # Pass #1
    for i,row in enumerate(S):
        Ni = set(row) | set([i])

        if Ni.issubset(R):
            R -= Ni
            for x in Ni:
                aggregates[x] = j
            j += 1

    # Pass #2
    Old_R = R.copy()
    for i,row in enumerate(S):
        if i not in R: continue
        
        for x in row:
            if x not in Old_R:
                aggregates[i] = aggregates[x]
                R.remove(i)
                break

    # Pass #3
    for i,row in enumerate(S):
        if i not in R: continue
        Ni = set(row) | set([i])

        for x in Ni:
            if x in R:
                aggregates[x] = j
            j += 1

    assert(len(R) == 0)

    Pj = aggregates
    Pp = arange(n+1)
    Px = ones(n)
 
    return csr_matrix((Px,Pj,Pp))

class TestSAConstantInterpolation(NumpyTestCase):
    def check_random(self):
        numpy.random.seed(0)
        for N in [2,3,5,10]:
            A = csr_matrix(rand(N,N))
            for epsilon in [0.0,0.1,0.5,1.0]:
                S_result   = sa_constant_interpolation(A,epsilon)
                S_expected = reference_sa_constant_interpolation(A,epsilon)
                assert_array_equal(S_result.todense(),S_expected.todense())

    def check_poisson1D(self):
        for N in [2,3,5,7,10,11,20,21,29,30]:
            A = poisson_problem1D(N)
            for epsilon in [0.0,0.1,0.5,1.0]:
                S_result   = sa_constant_interpolation(A,epsilon)
                S_expected = reference_sa_constant_interpolation(A,epsilon)
                assert_array_equal(S_result.todense(),S_expected.todense())

    def check_poisson2D(self):
        for N in [2,3,5,7,10,11]:
            A = poisson_problem2D(N)
            for epsilon in [0.0,0.1,0.5,1.0]:
                S_result   = sa_constant_interpolation(A,epsilon)
                S_expected = reference_sa_constant_interpolation(A,epsilon)
                assert_array_equal(S_result.todense(),S_expected.todense())

##    def check_sample_data(self):
##        from examples import all_examples,read_matrix
##
##        for filename in all_examples:
##            A = read_matrix(filename)            
##            for epsilon in [0.0,0.08,0.51,1.0]:
##                S_result   = sa_constant_interpolation(A,epsilon)
##                S_expected = reference_sa_constant_interpolation(A,epsilon)
##                assert_array_equal((S_result - S_expected).nnz,0)

if __name__ == '__main__':
    NumpyTest().run()

