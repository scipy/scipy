import warnings

from numpy import array, finfo, arange
import numpy.random as random
from numpy.testing import *

from scipy.linalg import norm, inv
from scipy.sparse import spdiags, SparseEfficiencyWarning
from scipy.sparse.linalg.dsolve import spsolve, use_solver, splu, spilu

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


class TestSplu(object):
    def setUp(self):
        n = 40
        d = arange(n) + 1
        self.n = n
        self.A = spdiags((d, 2*d, d[::-1]), (-3, 0, 5), n, n)
        random.seed(1234)

    def test_splu(self):
        x = random.rand(self.n)
        lu = splu(self.A)
        r = self.A*lu.solve(x)
        assert abs(x - r).max() < 1e-13

    def test_spilu(self):
        x = random.rand(self.n)
        lu = spilu(self.A, drop_tol=1e-2, fill_factor=5)
        r = self.A*lu.solve(x)
        assert abs(x - r).max() < 1e-2
        assert abs(x - r).max() > 1e-5

if __name__ == "__main__":
    run_module_suite()
