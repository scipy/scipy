from scipy.testing import *

from numpy import arange, ones, zeros, array, eye, vstack, diff
from scipy import rand
from scipy.sparse import csr_matrix


from scipy.sandbox.multigrid.sa import sa_fit_candidates
from scipy.sandbox.multigrid import smoothed_aggregation_solver
#from scipy.sandbox.multigrid.adaptive import augment_candidates

from scipy.sandbox.multigrid.gallery import *
from scipy.sandbox.multigrid.adaptive import *


class TestAdaptiveSA(TestCase):
    def setUp(self):
        from numpy.random import seed
        seed(0)

    def test_poisson(self):
        A = poisson( (100,100), format='csr' )

        asa = adaptive_sa_solver(A, max_candidates = 1)
        sa  = smoothed_aggregation_solver(A, B = ones((A.shape[0],1)) )

        b = rand(A.shape[0])

        sol0,residuals0 = asa.solve(b, maxiter=20, tol=1e-10, return_residuals=True)
        sol1,residuals1 =  sa.solve(b, maxiter=20, tol=1e-10, return_residuals=True)
       
        conv_asa = (residuals0[-1]/residuals0[0])**(1.0/len(residuals0))
        conv_sa  = (residuals1[-1]/residuals1[0])**(1.0/len(residuals1))
        
        assert( conv_asa < 1.1 * conv_sa ) #aSA shouldn't be any worse than SA

#    def test_elasticity(self):
#        A,B = linear_elasticity( (100,100), format='bsr' )
#
#        asa = adaptive_sa_solver(A, max_candidates = 3)
#        sa  = smoothed_aggregation_solver(A, B=B )
#
#        b = rand(A.shape[0])
#
#        sol0,residuals0 = asa.solve(b, maxiter=20, tol=1e-10, return_residuals=True)
#        sol1,residuals1 =  sa.solve(b, maxiter=20, tol=1e-10, return_residuals=True)
#       
#        conv_asa = (residuals0[-1]/residuals0[0])**(1.0/len(residuals0))
#        conv_sa  = (residuals1[-1]/residuals1[0])**(1.0/len(residuals1))
#       
#        print "ASA convergence",conv_asa
#        assert( conv_asa < 1.1 * conv_sa ) #aSA shouldn't be any worse than SA
        
#class TestAugmentCandidates(TestCase):
#    def setUp(self):
#        self.cases = []
#
#        #two candidates
#
#        ##block candidates
#        ##self.cases.append((csr_matrix((ones(9),array([0,0,0,1,1,1,2,2,2]),arange(10)),shape=(9,3)), vstack((array([1]*9 + [0]*9),arange(2*9))).T ))
#
#    def test_first_level(self):
#        cases = []
#
#        ## tests where AggOp includes all DOFs
#        cases.append((csr_matrix((ones(4),array([0,0,1,1]),arange(5)),shape=(4,2)), vstack((ones(4),arange(4))).T ))
#        cases.append((csr_matrix((ones(9),array([0,0,0,1,1,1,2,2,2]),arange(10)),shape=(9,3)), vstack((ones(9),arange(9))).T ))
#        cases.append((csr_matrix((ones(9),array([0,0,1,1,2,2,3,3,3]),arange(10)),shape=(9,4)), vstack((ones(9),arange(9))).T ))
#
#        ## tests where AggOp excludes some DOFs
#        cases.append((csr_matrix((ones(4),array([0,0,1,1]),array([0,1,2,2,3,4])),shape=(5,2)), vstack((ones(5),arange(5))).T ))
#
#        # overdetermined blocks
#        cases.append((csr_matrix((ones(4),array([0,0,1,1]),array([0,1,2,2,3,4])),shape=(5,2)), vstack((ones(5),arange(5),arange(5)**2)).T  ))
#        cases.append((csr_matrix((ones(6),array([1,3,0,2,1,0]),array([0,0,1,2,2,3,4,5,5,6])),shape=(9,4)), vstack((ones(9),arange(9),arange(9)**2)).T ))
#        cases.append((csr_matrix((ones(6),array([1,3,0,2,1,0]),array([0,0,1,2,2,3,4,5,5,6])),shape=(9,4)), vstack((ones(9),arange(9))).T ))
#
#        def mask_candidate(AggOp,candidates):
#            #mask out all DOFs that are not included in the aggregation
#            candidates[diff(AggOp.indptr) == 0,:] = 0
#
#        for AggOp,fine_candidates in cases:
#
#            mask_candidate(AggOp,fine_candidates)
#
#            for i in range(1,fine_candidates.shape[1]):
#                Q_expected,R_expected = sa_fit_candidates(AggOp,fine_candidates[:,:i+1])
#
#                old_Q, old_R = sa_fit_candidates(AggOp,fine_candidates[:,:i])
#
#                Q_result,R_result = augment_candidates(AggOp, old_Q, old_R, fine_candidates[:,[i]])
#
#                # compare against SA method (which is assumed to be correct)
#                assert_almost_equal(Q_expected.todense(),Q_result.todense())
#                assert_almost_equal(R_expected,R_result)
#
#                #each fine level candidate should be fit exactly
#                assert_almost_equal(fine_candidates[:,:i+1],Q_result*R_result)
#                assert_almost_equal(Q_result*(Q_result.T*fine_candidates[:,:i+1]),fine_candidates[:,:i+1])
#

if __name__ == '__main__':
    nose.run(argv=['', __file__])
