import warnings

from numpy import array, finfo
from numpy.testing import *

from scipy.linalg import norm, inv
from scipy.sparse import spdiags, SparseEfficiencyWarning
from scipy.sparse.linalg.dsolve import spsolve, use_solver

warnings.simplefilter('ignore',SparseEfficiencyWarning)

#TODO add more comprehensive tests
use_solver( useUmfpack = False )

class TestLinsolve(TestCase):
    ## this crashes SuperLU
    #def test_singular(self):
    #    A = csc_matrix( (5,5), dtype='d' )
    #    b = array([1, 2, 3, 4, 5],dtype='d')
    #    x = spsolve(A,b)

    def test_twodiags(self):
        A = spdiags([[1, 2, 3, 4, 5], [6, 5, 8, 9, 10]], [0, 1], 5, 5)
        b = array([1, 2, 3, 4, 5])

        # condition number of A
        cond_A = norm(A.todense(),2) * norm(inv(A.todense()),2)


        for t in ['f','d','F','D']:
            eps = finfo(t).eps #floating point epsilon
            b = b.astype(t)

            for format in ['csc','csr']:
                Asp = A.astype(t).asformat(format)

                x = spsolve(Asp,b)

                assert( norm(b - Asp*x) < 10 * cond_A * eps )


if __name__ == "__main__":
    run_module_suite()
