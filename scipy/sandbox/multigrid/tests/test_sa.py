from numpy.testing import *

from numpy import sqrt,empty,ones,arange,array_split,eye,array, \
                  zeros,diag,zeros_like,diff,matrix
from numpy.linalg import norm                  
from scipy import rand
from scipy.sparse import spdiags,csr_matrix,lil_matrix, \
                         isspmatrix_csr,isspmatrix_csc,isspmatrix_coo, \
                         isspmatrix_lil
import numpy

set_package_path()
import scipy.sandbox.multigrid
from scipy.sandbox.multigrid.sa import sa_strong_connections, sa_constant_interpolation, \
                                        sa_interpolation, sa_fit_candidates, \
                                        sa_smoothed_prolongator
from scipy.sandbox.multigrid.multilevel import poisson_problem1D,poisson_problem2D, \
                                        smoothed_aggregation_solver
from scipy.sandbox.multigrid.utils import diag_sparse
restore_path()


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

class TestSA(NumpyTestCase):
    def setUp(self):
        self.cases = []

        # random matrices
        numpy.random.seed(0)        
        for N in [2,3,5]:
            self.cases.append( csr_matrix(rand(N,N)) ) 
        
        # poisson problems in 1D and 2D
        for N in [2,3,5,7,10,11,19]:
            self.cases.append( poisson_problem1D(N) )
        for N in [2,3,5,7,10,11]:
            self.cases.append( poisson_problem2D(N) )


    def check_sa_strong_connections(self):
        for A in self.cases:
            for epsilon in [0.0,0.1,0.5,1.0,10.0]:
                S_expected = reference_sa_strong_connections(A,epsilon)
                S_result = sa_strong_connections(A,epsilon)
                assert_almost_equal(S_result.todense(),S_expected.todense())
                #assert_array_equal(sparsity(S_result).todense(),sparsity(S_expected).todense())

    def check_sa_constant_interpolation(self):
        for A in self.cases:
            for epsilon in [0.0,0.1,0.5,1.0]:
                S_expected = reference_sa_constant_interpolation(A,epsilon)
                
                S_result   = sa_constant_interpolation(A,epsilon,blocks=None)
                assert_array_equal(S_result.todense(),S_expected.todense())
               
                #blocks=1...N should be the same as blocks=None
                S_result   = sa_constant_interpolation(A,epsilon,blocks=arange(A.shape[0]))
                assert_array_equal(S_result.todense(),S_expected.todense())

        # two aggregates in 1D
        A = poisson_problem1D(6)
        AggOp = csr_matrix((ones(6),array([0,0,0,1,1,1]),arange(7)),dims=(6,2))
        candidates = [ones(6)]

        T_result,coarse_candidates_result = sa_fit_candidates(AggOp,candidates)
        T_expected = csr_matrix((sqrt(1.0/3.0)*ones(6),array([0,0,0,1,1,1]),arange(7)),dims=(6,2))
        assert_almost_equal(T_result.todense(),T_expected.todense())

        #check simple block examples
        A = csr_matrix(arange(16).reshape(4,4))
        A = A + A.T
        blocks = array([0,0,1,1])

        S_result   = sa_constant_interpolation(A,epsilon=0.0,blocks=blocks)
        S_expected = matrix([1,1,1,1]).T
        assert_array_equal(S_result.todense(),S_expected)

        S_result   = sa_constant_interpolation(A,epsilon=0.5,blocks=blocks)
        S_expected = matrix([1,1,1,1]).T
        assert_array_equal(S_result.todense(),S_expected)

        S_result   = sa_constant_interpolation(A,epsilon=2.0,blocks=blocks)
        S_expected = matrix([[1,0],[1,0],[0,1],[0,1]])
        assert_array_equal(S_result.todense(),S_expected)


    def check_user_aggregation(self):
        """check that the sa_interpolation accepts user-defined aggregates"""

        user_cases = []

        #simple 1d example w/ two aggregates
        A = poisson_problem1D(6)
        AggOp = csr_matrix((ones(6),array([0,0,0,1,1,1]),arange(7)),dims=(6,2))
        candidates = [ones(6)]
        user_cases.append((A,AggOp,candidates))

        #simple 1d example w/ two aggregates (not all nodes are aggregated)
        A = poisson_problem1D(6)
        AggOp = csr_matrix((ones(4),array([0,0,1,1]),array([0,1,1,2,3,3,4])),dims=(6,2))
        candidates = [ones(6)]
        user_cases.append((A,AggOp,candidates))

        for A,AggOp,candidates in user_cases:
            T,coarse_candidates_result = sa_fit_candidates(AggOp,candidates)

            P_result = sa_interpolation(A,candidates,omega=4.0/3.0,AggOp=AggOp)[0]
            P_expected = sa_smoothed_prolongator(A,T,epsilon=0.0,omega=4.0/3.0)

            assert_almost_equal(P_result.todense(),P_expected.todense())

                  

class TestFitCandidates(NumpyTestCase):
    def setUp(self):
        self.cases = []

        ## tests where AggOp includes all DOFs
        #one candidate
        self.cases.append((csr_matrix((ones(5),array([0,0,0,1,1]),arange(6)),dims=(5,2)),[ones(5)]))
        self.cases.append((csr_matrix((ones(5),array([1,1,0,0,0]),arange(6)),dims=(5,2)),[ones(5)]))
        self.cases.append((csr_matrix((ones(9),array([0,0,0,1,1,1,2,2,2]),arange(10)),dims=(9,3)),[ones(9)]))
        self.cases.append((csr_matrix((ones(9),array([2,1,0,0,1,2,1,0,2]),arange(10)),dims=(9,3)),[arange(9)]))
        #two candidates
        self.cases.append((csr_matrix((ones(4),array([0,0,1,1]),arange(5)),dims=(4,2)),[ones(4),arange(4)]))
        self.cases.append((csr_matrix((ones(9),array([0,0,0,1,1,1,2,2,2]),arange(10)),dims=(9,3)),[ones(9),arange(9)]))
        self.cases.append((csr_matrix((ones(9),array([0,0,1,1,2,2,3,3,3]),arange(10)),dims=(9,4)),[ones(9),arange(9)]))
        #block candidates
        self.cases.append((csr_matrix((ones(9),array([0,0,0,1,1,1,2,2,2]),arange(10)),dims=(9,3)),[array([1]*9 + [0]*9),arange(2*9)]))
        
        ## tests where AggOp excludes some DOFs
        self.cases.append((csr_matrix((ones(4),array([0,0,1,1]),array([0,1,2,2,3,4])),dims=(5,2)),[ones(5)]))
        self.cases.append((csr_matrix((ones(4),array([0,0,1,1]),array([0,1,2,2,3,4])),dims=(5,2)),[ones(5),arange(5)]))
        self.cases.append((csr_matrix((ones(6),array([1,3,0,2,1,0]),array([0,0,1,2,2,3,4,5,5,6])),dims=(9,4)),[ones(9),arange(9)]))

    def check_all_cases(self):
        """Test case where aggregation includes all fine nodes"""
        
        def mask_candidate(AggOp,candidates):
            #mask out all DOFs that are not included in the aggregation
            for c in candidates:
                c[diff(AggOp.indptr) == 0] = 0

        for AggOp,fine_candidates in self.cases:

            mask_candidate(AggOp,fine_candidates)

            Q,coarse_candidates = sa_fit_candidates(AggOp,fine_candidates)

            assert_equal(len(coarse_candidates),len(fine_candidates))

            #each fine level candidate should be fit exactly
            for fine,coarse in zip(fine_candidates,coarse_candidates):
                assert_almost_equal(fine,Q*coarse)
                assert_almost_equal(Q*(Q.T*fine),fine)

            #aggregate one more level (to a single aggregate)
            K = len(coarse_candidates)
            N = K*AggOp.shape[1]
            AggOp = csr_matrix((ones(N),zeros(N),arange(N + 1)),dims=(N,1)) #aggregate to a single point
            fine_candidates = coarse_candidates
            
            #mask_candidate(AggOp,fine_candidates)  #not needed
            
            #now check the coarser problem
            Q,coarse_candidates = sa_fit_candidates(AggOp,fine_candidates)

            assert_equal(len(coarse_candidates),len(fine_candidates))

            for fine,coarse in zip(fine_candidates,coarse_candidates):
                assert_almost_equal(fine,Q*coarse)
                assert_almost_equal(Q*(Q.T*fine),fine)



class TestSASolverPerformance(NumpyTestCase):
    def setUp(self):
        self.cases = []

        self.cases.append((poisson_problem1D(100),None))
        self.cases.append((poisson_problem2D(50),None))
        # TODO add unstructured tests
       

    def check_basic(self):
        """check that method converges at a reasonable rate"""

        for A,candidates in self.cases:
            ml = smoothed_aggregation_solver(A,candidates,max_coarse=10,max_levels=10)

            numpy.random.seed(0) #make tests repeatable
            
            x = rand(A.shape[0])
            b = A*rand(A.shape[0]) #zeros_like(x)
            
            x_sol,residuals = ml.solve(b,x0=x,maxiter=20,tol=1e-12,return_residuals=True)

            avg_convergence_ratio = (residuals[-1]/residuals[0])**(1.0/len(residuals))
            
            assert(avg_convergence_ratio < 0.5)

    def check_DAD(self):
        for A,candidates in self.cases:

            x = rand(A.shape[0])
            b = A*rand(A.shape[0]) 
            
            D     = diag_sparse(1.0/sqrt(10**(12*rand(A.shape[0])-6))).tocsr()
            D_inv = diag_sparse(1.0/D.data)

            DAD   = D*A*D
           
            if candidates is None:
                candidates = [ ones(A.shape[0]) ]
           
            DAD_candidates = [ (D_inv * c) for c in candidates ]
           
            #TODO force 2 level method and check that result is the same

            #ml = smoothed_aggregation_solver(A,candidates,max_coarse=1,max_levels=2)

            ml = smoothed_aggregation_solver(DAD,DAD_candidates,max_coarse=100,max_levels=2)

            #print (D_inv*ml.Ps[0]).todense()
            
            x_sol,residuals = ml.solve(b,x0=x,maxiter=10,tol=1e-12,return_residuals=True)

            avg_convergence_ratio = (residuals[-1]/residuals[0])**(1.0/len(residuals))
            #print avg_convergence_ratio

            assert(avg_convergence_ratio < 0.5)




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



if __name__ == '__main__':
    NumpyTest().run()

