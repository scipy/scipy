from numpy.testing import *

from numpy import sqrt,empty,ones,arange,array_split,eye,array, \
                  zeros,diag,zeros_like
from numpy.linalg import norm                  
from scipy import rand
from scipy.sparse import spdiags,csr_matrix,lil_matrix, \
                         isspmatrix_csr,isspmatrix_csc,isspmatrix_coo, \
                         isspmatrix_lil
import numpy

set_package_path()
import scipy.sandbox.multigrid
from scipy.sandbox.multigrid.sa import sa_strong_connections, sa_constant_interpolation, \
                                        sa_interpolation, sa_fit_candidates
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

def reference_sa_strong_connections(A,epsilon):
    A_coo = A.tocoo()
    S = lil_matrix(A.shape)

    for (i,j,v) in zip(A_coo.row,A_coo.col,A_coo.data):
        if i == j or abs(v) >= epsilon*sqrt(abs(A[i,i])*abs(A[j,j])):
            S[i,j] += v
        else:
            S[i,i] += v

    return S

class TestSAStrongConnections(NumpyTestCase):
#    def check_simple(self):
#        N = 4
#        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).tocsr()
#        assert_array_equal(sa_strong_connections(A,0.50).todense(),A.todense())   #all connections are strong
#        assert_array_equal(sa_strong_connections(A,0.51).todense(),0*A.todense()) #no connections are strong
#       
#        N = 100
#        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).tocsr()
#        assert_array_equal(sa_strong_connections(A,0.50).todense(),A.todense())   #all connections are strong
#        assert_array_equal(sa_strong_connections(A,0.51).todense(),0*A.todense()) #no connections are strong

    def check_random(self):
        numpy.random.seed(0)

        for N in [2,3,5]:
            A = csr_matrix(rand(N,N))
            for epsilon in [0.0,0.1,0.5,1.0,10.0]:
                S_result = sa_strong_connections(A,epsilon)
                S_expected = reference_sa_strong_connections(A,epsilon)
                assert_almost_equal(S_result.todense(),S_expected.todense())
                #assert_array_equal(sparsity(S_result).todense(),sparsity(S_expected).todense())

    def check_poisson1D(self):
        for N in [2,3,5,7,10,11,19]:
            A = poisson_problem1D(N)
            for epsilon in [0.0,0.1,0.5,1.0]:
                S_result   = sa_strong_connections(A,epsilon)
                S_expected = reference_sa_strong_connections(A,epsilon)
                assert_array_equal(S_result.todense(),S_expected.todense())
                #assert_array_equal(sparsity(S_result).todense(),sparsity(S_expected).todense())

    def check_poisson2D(self):
        for N in [2,3,5,7,10,11]:
            A = poisson_problem2D(N)
            for epsilon in [0.0,0.1,0.5,1.0]:
                S_result   = sa_strong_connections(A,epsilon)
                S_expected = reference_sa_strong_connections(A,epsilon)
                assert_array_equal(S_result.todense(),S_expected.todense())
                #assert_array_equal(sparsity(S_result).todense(),sparsity(S_expected).todense())




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

class TestFitCandidates(NumpyTestCase):
    def setUp(self):
        self.normal_cases = []

        #one candidate
        self.normal_cases.append((csr_matrix((ones(5),array([0,0,0,1,1]),arange(6)),dims=(5,2)),[ones(5)]))
        self.normal_cases.append((csr_matrix((ones(5),array([1,1,0,0,0]),arange(6)),dims=(5,2)),[ones(5)]))
        self.normal_cases.append((csr_matrix((ones(9),array([0,0,0,1,1,1,2,2,2]),arange(10)),dims=(9,3)),[ones(9)]))
        self.normal_cases.append((csr_matrix((ones(9),array([2,1,0,0,1,2,1,0,2]),arange(10)),dims=(9,3)),[arange(9)]))

        #two candidates
        self.normal_cases.append((csr_matrix((ones(4),array([0,0,1,1]),arange(5)),dims=(4,2)),[ones(4),arange(4)]))
        self.normal_cases.append((csr_matrix((ones(9),array([0,0,0,1,1,1,2,2,2]),arange(10)),dims=(9,3)),[ones(9),arange(9)]))
        self.normal_cases.append((csr_matrix((ones(9),array([0,0,1,1,2,2,3,3,3]),arange(10)),dims=(9,4)),[ones(9),arange(9)]))

        #block candidates
        self.normal_cases.append((csr_matrix((ones(9),array([0,0,0,1,1,1,2,2,2]),arange(10)),dims=(9,3)),[array([1]*9 + [0]*9),arange(2*9)]))

        #TODO add test case where aggregation operator has holes

    def check_normal(self):
        """Test case where aggregation includes all fine nodes"""
        
        for AggOp,fine_candidates in self.normal_cases:
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
            
            #now check the coarser problem
            Q,coarse_candidates = sa_fit_candidates(AggOp,fine_candidates)

            assert_equal(len(coarse_candidates),len(fine_candidates))

            for fine,coarse in zip(fine_candidates,coarse_candidates):
                assert_almost_equal(fine,Q*coarse)
                assert_almost_equal(Q*(Q.T*fine),fine)




class TestSASolver(NumpyTestCase):
    def setUp(self):
        self.cases = []

        #self.cases.append((poisson_problem1D(10),None))

        self.cases.append((poisson_problem1D(500),None))
        self.cases.append((poisson_problem2D(50),None))
        
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
            b = A*rand(A.shape[0]) #zeros_like(x)
            
            D = diag_sparse(rand(A.shape[0]))
            D_inv = diag_sparse(1.0/D.data)
            DAD = D*A*D
           
            if candidates is None:
                candidates = [ ones(A.shape[0]) ]
           
            DAD_candidates = [ (D_inv * c) for c in candidates ]
            
            #ml = smoothed_aggregation_solver(A,candidates,max_coarse=1,max_levels=2)

            ml = smoothed_aggregation_solver(DAD,DAD_candidates,max_coarse=100,max_levels=2)

            #print (D_inv*ml.Ps[0]).todense()
            
            x_sol,residuals = ml.solve(b,x0=x,maxiter=10,tol=1e-12,return_residuals=True)

            avg_convergence_ratio = (residuals[-1]/residuals[0])**(1.0/len(residuals))
            print avg_convergence_ratio

            assert(avg_convergence_ratio < 0.5)
        
##    def check_DAD(self):
##        """check that method is invariant to symmetric diagonal scaling (A -> DAD)"""
##
##        for A,A_candidates in self.cases:
##            numpy.random.seed(0) #make tests repeatable
##            
##            x = rand(A.shape[0])
##            b = zeros_like(x)
##
##            D = diag_sparse(rand(A.shape[0]))
##            D_inv = diag_sparse(1.0/D.data)
##            DAD = D*A*D
##            
##
##            if A_candidates is None:
##                A_candidates = [ ones(A.shape[0]) ]
##           
##            DAD_candidates = [ (D_inv * c) for c in A_candidates ]
##            
##            ml_A    = smoothed_aggregation_solver(A,     A_candidates, max_coarse=10, max_levels=10, epsilon=0.0)
##            x_sol_A = ml_A.solve(b, x0=x, maxiter=1, tol=1e-12)
##
##            ml_DAD    = smoothed_aggregation_solver(DAD, DAD_candidates, max_coarse=10, max_levels=10, epsilon=0.0)
##            x_sol_DAD = ml_DAD.solve(b, x0=D*x, maxiter=1, tol=1e-12)
##            
##            assert_almost_equal(x_sol_A, D_inv * x_sol_DAD)


if __name__ == '__main__':
    NumpyTest().run()

