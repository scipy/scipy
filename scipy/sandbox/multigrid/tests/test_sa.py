try:
    from sets import Set
    set = Set
except:
    pass

from scipy.testing import *
from numpy import sqrt,empty,ones,arange,array_split,eye,array, \
                  zeros,diag,zeros_like,diff,matrix,hstack,vstack
from numpy.linalg import norm
from scipy import rand
from scipy.sparse import spdiags,csr_matrix,lil_matrix, \
                         isspmatrix_csr,isspmatrix_csc,isspmatrix_coo, \
                         isspmatrix_lil
import numpy


import scipy.sandbox.multigrid
from scipy.sandbox.multigrid.sa import *
from scipy.sandbox.multigrid.utils import diag_sparse
from scipy.sandbox.multigrid.gallery import poisson, linear_elasticity

#def sparsity(A):
#    A = A.copy()
#
#    if isspmatrix_csr(A) or isspmatrix_csc(A) or isspmatrix_coo(A):
#        A.data[:] = 1
#    elif isspmatrix_lil:
#        for row in A.data:
#            row[:] = [1]*len(row)
#    else:
#        raise ValueError,'expected csr,csc,coo, or lil'
#
#    return A

class TestSA(TestCase):
    def setUp(self):
        self.cases = []

        # random matrices
        numpy.random.seed(0)
        for N in [2,3,5]:
            self.cases.append( csr_matrix(rand(N,N)) )

        # poisson problems in 1D and 2D
        for N in [2,3,5,7,10,11,19]:
            self.cases.append( poisson( (N,), format='csr') )
        for N in [2,3,5,7,10,11]:
            self.cases.append( poisson( (N,N), format='csr') )


    def test_sa_strong_connections(self):
        for A in self.cases:
            for epsilon in [0.0,0.1,0.5,1.0,10.0]:
                S_expected = reference_sa_strong_connections(A,epsilon)
                S_result = sa_strong_connections(A,epsilon)
                assert_almost_equal(S_result.todense(),S_expected.todense())
                #assert_array_equal(sparsity(S_result).todense(),sparsity(S_expected).todense())

        ##check simple block examples
        #A = csr_matrix(arange(16).reshape(4,4))
        #A = A + A.T
        #A = A.tobsr(blocksize=(2,2))

        #S_result   = sa_standard_aggregation(A,epsilon=0.0)
        #S_expected = matrix([1,1]).T
        #assert_array_equal(S_result.todense(),S_expected)

        #S_result   = sa_standard_aggregation(A,epsilon=0.5)
        #S_expected = matrix([1,1]).T
        #assert_array_equal(S_result.todense(),S_expected)

        #S_result   = sa_standard_aggregation(A,epsilon=2.0)
        #S_expected = matrix([[1,0],[0,1]])
        #assert_array_equal(S_result.todense(),S_expected)

    def test_sa_standard_aggregation(self):
        for C in self.cases:
            S_expected = reference_sa_standard_aggregation(C)

            S_result   = sa_standard_aggregation(C)
            assert_array_equal(S_result.todense(),S_expected.todense())



#    def test_user_aggregation(self):
#        """check that the sa_interpolation accepts user-defined aggregates"""
#
#        user_cases = []
#
#        #simple 1d example w/ two aggregates
#        A = poisson( (6,), format='csr')
#        AggOp = csr_matrix((ones(6),array([0,0,0,1,1,1]),arange(7)),shape=(6,2))
#        candidates = ones((6,1))
#        user_cases.append((A,AggOp,candidates))
#
#        #simple 1d example w/ two aggregates (not all nodes are aggregated)
#        A = poisson( (6,), format='csr')
#        AggOp = csr_matrix((ones(4),array([0,0,1,1]),array([0,1,1,2,3,3,4])),shape=(6,2))
#        candidates = ones((6,1))
#        user_cases.append((A,AggOp,candidates))
#
#        for A,AggOp,candidates in user_cases:
#            T,coarse_candidates_result = sa_fit_candidates(AggOp,candidates)
#
#            P_result = sa_interpolation(A,candidates,omega=4.0/3.0,AggOp=AggOp)[0]
#            P_expected = sa_smoothed_prolongator(A,T,epsilon=0.0,omega=4.0/3.0)
#
#            assert_almost_equal(P_result.todense(),P_expected.todense())



class TestFitCandidates(TestCase):
    def setUp(self):
        self.cases = []

        ### tests where AggOp includes all DOFs
        # one candidate
        self.cases.append((csr_matrix((ones(5),array([0,0,0,1,1]),arange(6)),shape=(5,2)), ones((5,1)) ))
        self.cases.append((csr_matrix((ones(5),array([1,1,0,0,0]),arange(6)),shape=(5,2)), ones((5,1)) ))
        self.cases.append((csr_matrix((ones(9),array([0,0,0,1,1,1,2,2,2]),arange(10)),shape=(9,3)), ones((9,1)) ))
        self.cases.append((csr_matrix((ones(9),array([2,1,0,0,1,2,1,0,2]),arange(10)),shape=(9,3)), arange(9).reshape(9,1) ))
        # two candidates
        self.cases.append((csr_matrix((ones(4),array([0,0,1,1]),arange(5)),shape=(4,2)), vstack((ones(4),arange(4))).T ))
        self.cases.append((csr_matrix((ones(9),array([0,0,0,1,1,1,2,2,2]),arange(10)),shape=(9,3)), vstack((ones(9),arange(9))).T ))
        self.cases.append((csr_matrix((ones(9),array([0,0,1,1,2,2,3,3,3]),arange(10)),shape=(9,4)), vstack((ones(9),arange(9))).T ))
       
        # block aggregates, one candidate
        self.cases.append((csr_matrix((ones(3),array([0,1,1]),arange(4)),shape=(3,2)), ones((6,1)) ))
        self.cases.append((csr_matrix((ones(3),array([0,1,1]),arange(4)),shape=(3,2)), ones((9,1)) ))
        self.cases.append((csr_matrix((ones(5),array([2,0,2,1,1]),arange(6)),shape=(5,3)), ones((10,1)) ))
        
        # block aggregates, two candidates
        self.cases.append((csr_matrix((ones(3),array([0,1,1]),arange(4)),shape=(3,2)), vstack((ones(6),arange(6))).T ))
        self.cases.append((csr_matrix((ones(3),array([0,1,1]),arange(4)),shape=(3,2)), vstack((ones(9),arange(9))).T ))
        self.cases.append((csr_matrix((ones(5),array([2,0,2,1,1]),arange(6)),shape=(5,3)), vstack((ones(10),arange(10))).T ))

        ### tests where AggOp excludes some DOFs
        # one candidate
        self.cases.append((csr_matrix((ones(4),array([0,0,1,1]),array([0,1,2,2,3,4])),shape=(5,2)), ones((5,1)) ))
        self.cases.append((csr_matrix((ones(4),array([0,0,1,1]),array([0,1,2,2,3,4])),shape=(5,2)), vstack((ones(5),arange(5))).T ))

        # overdetermined blocks
        self.cases.append((csr_matrix((ones(4),array([0,0,1,1]),array([0,1,2,2,3,4])),shape=(5,2)), vstack((ones(5),arange(5),arange(5)**2)).T  ))
        self.cases.append((csr_matrix((ones(6),array([1,3,0,2,1,0]),array([0,0,1,2,2,3,4,5,5,6])),shape=(9,4)), vstack((ones(9),arange(9),arange(9)**2)).T ))
        self.cases.append((csr_matrix((ones(6),array([1,3,0,2,1,0]),array([0,0,1,2,2,3,4,5,5,6])),shape=(9,4)), vstack((ones(9),arange(9))).T ))

    def test_all_cases(self):
        """Test case where aggregation includes all fine nodes"""

        def mask_candidate(AggOp,candidates):
            #mask out all DOFs that are not included in the aggregation
            candidates[diff(AggOp.indptr) == 0,:] = 0

        for AggOp,fine_candidates in self.cases:
            mask_candidate(AggOp,fine_candidates)

            Q,coarse_candidates = sa_fit_candidates(AggOp,fine_candidates)

            #each fine level candidate should be fit exactly
            assert_almost_equal(fine_candidates,Q*coarse_candidates)
            assert_almost_equal(Q*(Q.T*fine_candidates),fine_candidates)


class TestSASolverPerformance(TestCase):
    def setUp(self):
        self.cases = []

        self.cases.append(( poisson( (10000,),  format='csr'), None))
        self.cases.append(( poisson( (100,100), format='csr'), None))
        self.cases.append( linear_elasticity( (100,100), format='bsr') )
        # TODO add unstructured tests


    def test_basic(self):
        """check that method converges at a reasonable rate"""

        for A,B in self.cases:
            ml = smoothed_aggregation_solver(A,B,max_coarse=10,max_levels=10)

            numpy.random.seed(0) #make tests repeatable

            x = rand(A.shape[0])
            b = A*rand(A.shape[0]) #zeros_like(x)

            x_sol,residuals = ml.solve(b,x0=x,maxiter=20,tol=1e-10,return_residuals=True)

            avg_convergence_ratio = (residuals[-1]/residuals[0])**(1.0/len(residuals))
            
            assert(avg_convergence_ratio < 0.3)

    def test_DAD(self):
        A = poisson( (200,200), format='csr' )        

        x = rand(A.shape[0])
        b = rand(A.shape[0])
 
        D     = diag_sparse(1.0/sqrt(10**(12*rand(A.shape[0])-6))).tocsr()
        D_inv = diag_sparse(1.0/D.data)
 
        DAD   = D*A*D
 
        B = ones((A.shape[0],1))
 
        #TODO force 2 level method and check that result is the same
 
        sa = smoothed_aggregation_solver(D*A*D, D_inv * B, max_levels=2)
 
        x_sol,residuals = sa.solve(b,x0=x,maxiter=10,tol=1e-12,return_residuals=True)
 
        avg_convergence_ratio = (residuals[-1]/residuals[0])**(1.0/len(residuals))
        
        assert(avg_convergence_ratio < 0.25)




################################################
##   reference implementations for unittests  ##
################################################
def reference_sa_strong_connections(A,epsilon):
    A_coo = A.tocoo()
    S = lil_matrix(A.shape)

    for (i,j,v) in zip(A_coo.row,A_coo.col,A_coo.data):
        if i == j or abs(v) >= epsilon*sqrt(abs(A[i,i])*abs(A[j,j])):
            S[i,j] += v
        else:
            S[i,i] += v

    return S


# note that this method only tests the current implementation, not
# all possible implementations
def reference_sa_standard_aggregation(C):
    S = array_split(C.indices,C.indptr[1:-1])

    n = C.shape[0]

    R = set(range(n))
    j = 0

    aggregates = empty(n,dtype=C.indices.dtype)
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



if __name__ == '__main__':
    nose.run(argv=['', __file__])
