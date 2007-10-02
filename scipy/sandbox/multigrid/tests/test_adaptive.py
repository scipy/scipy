from numpy.testing import *

from scipy.sparse import csr_matrix
from scipy import arange,ones,zeros,array,eye

set_package_path()
from scipy.sandbox.multigrid.adaptive import fit_candidates
restore_path()


class TestFitCandidates(NumpyTestCase):
    def setUp(self):
        self.cases = []

        #one candidate
        self.cases.append((csr_matrix((ones(5),array([0,0,0,1,1]),arange(6)),dims=(5,2)),[ones(5)]))
        self.cases.append((csr_matrix((ones(5),array([1,1,0,0,0]),arange(6)),dims=(5,2)),[ones(5)]))
        self.cases.append((csr_matrix((ones(9),array([0,0,0,1,1,1,2,2,2]),arange(10)),dims=(9,3)),[ones(9)]))
        self.cases.append((csr_matrix((ones(9),array([2,1,0,0,1,2,1,0,2]),arange(10)),dims=(9,3)),[arange(9)]))

        #two candidates
        self.cases.append((csr_matrix((ones(4),array([0,0,1,1]),arange(5)),dims=(4,2)),[ones(4),arange(4)]))
        self.cases.append((csr_matrix((ones(9),array([0,0,0,1,1,1,2,2,2]),arange(10)),dims=(9,3)),[ones(9),arange(9)]))
        self.cases.append((csr_matrix((ones(9),array([0,0,1,1,2,2,3,3,3]),arange(10)),dims=(9,4)),[ones(9),arange(9)]))

    def check_all(self):
        for AggOp,fine_candidates in self.cases:
            Q,coarse_candidates = fit_candidates(AggOp,fine_candidates)

            assert_equal(len(coarse_candidates),len(fine_candidates))
            assert_almost_equal((Q.T*Q).todense(),eye(Q.shape[1]))

            for fine,coarse in zip(fine_candidates,coarse_candidates):
                assert_almost_equal(fine,Q*coarse)

            #aggregate one more level (to a single aggregate)
            K = len(coarse_candidates)
            N = K*AggOp.shape[1]
            AggOp = csr_matrix((ones(N),zeros(N),arange(N + 1)),dims=(N,1))
            fine_candidates = coarse_candidates
            
            Q,coarse_candidates = fit_candidates(AggOp,fine_candidates)

            assert_equal(len(coarse_candidates),len(fine_candidates))
            assert_almost_equal((Q.T*Q).todense(),eye(Q.shape[1]))

            for fine,coarse in zip(fine_candidates,coarse_candidates):
                assert_almost_equal(fine,Q*coarse)

if __name__ == '__main__':
    NumpyTest().run()
      

