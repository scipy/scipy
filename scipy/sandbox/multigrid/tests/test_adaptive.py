from numpy.testing import *

from scipy.sparse import csr_matrix
from numpy import arange,ones,zeros,array,eye,vstack,diff

set_package_path()
from scipy.sandbox.multigrid.sa import sa_fit_candidates
from scipy.sandbox.multigrid.adaptive import augment_candidates
restore_path()

#import pdb; pdb.set_trace()

class TestAdaptiveSA(NumpyTestCase):
    def setUp(self):
        pass

class TestAugmentCandidates(NumpyTestCase):
    def setUp(self):
        self.cases = []

        #two candidates

        ##block candidates
        ##self.cases.append((csr_matrix((ones(9),array([0,0,0,1,1,1,2,2,2]),arange(10)),shape=(9,3)), vstack((array([1]*9 + [0]*9),arange(2*9))).T ))



    def check_first_level(self):
        cases = []

        ## tests where AggOp includes all DOFs
        cases.append((csr_matrix((ones(4),array([0,0,1,1]),arange(5)),shape=(4,2)), vstack((ones(4),arange(4))).T ))
        cases.append((csr_matrix((ones(9),array([0,0,0,1,1,1,2,2,2]),arange(10)),shape=(9,3)), vstack((ones(9),arange(9))).T ))
        cases.append((csr_matrix((ones(9),array([0,0,1,1,2,2,3,3,3]),arange(10)),shape=(9,4)), vstack((ones(9),arange(9))).T ))

        ## tests where AggOp excludes some DOFs
        cases.append((csr_matrix((ones(4),array([0,0,1,1]),array([0,1,2,2,3,4])),shape=(5,2)), vstack((ones(5),arange(5))).T ))

        # overdetermined blocks
        cases.append((csr_matrix((ones(4),array([0,0,1,1]),array([0,1,2,2,3,4])),shape=(5,2)), vstack((ones(5),arange(5),arange(5)**2)).T  ))
        cases.append((csr_matrix((ones(6),array([1,3,0,2,1,0]),array([0,0,1,2,2,3,4,5,5,6])),shape=(9,4)), vstack((ones(9),arange(9),arange(9)**2)).T ))
        cases.append((csr_matrix((ones(6),array([1,3,0,2,1,0]),array([0,0,1,2,2,3,4,5,5,6])),shape=(9,4)), vstack((ones(9),arange(9))).T ))

        def mask_candidate(AggOp,candidates):
            #mask out all DOFs that are not included in the aggregation
            candidates[diff(AggOp.indptr) == 0,:] = 0

        for AggOp,fine_candidates in cases:

            mask_candidate(AggOp,fine_candidates)

            for i in range(1,fine_candidates.shape[1]):
                Q_expected,R_expected = sa_fit_candidates(AggOp,fine_candidates[:,:i+1])

                old_Q, old_R = sa_fit_candidates(AggOp,fine_candidates[:,:i])

                Q_result,R_result = augment_candidates(AggOp, old_Q, old_R, fine_candidates[:,[i]])

                # compare against SA method (which is assumed to be correct)
                assert_almost_equal(Q_expected.todense(),Q_result.todense())
                assert_almost_equal(R_expected,R_result)

                #each fine level candidate should be fit exactly
                assert_almost_equal(fine_candidates[:,:i+1],Q_result*R_result)
                assert_almost_equal(Q_result*(Q_result.T*fine_candidates[:,:i+1]),fine_candidates[:,:i+1])


if __name__ == '__main__':
    NumpyTest().run()
