""" Unit tests for nonlinear solvers
Author: Ondrej Certik
May 2007
"""

from numpy.testing import *

from scipy.optimize import nonlin
from numpy import matrix, diag


def F(x):
    def p3(y):
        return float(y.T*y)*y
    try:
        x=tuple(x.flat)
    except:
        pass
    x=matrix(x).T

    d=matrix(diag([3,2,1.5,1,0.5]))
    c=0.01
    f=-d*x-c*p3(x)

    return tuple(f.flat)

class TestNonlin(TestCase):
    """ Test case for a simple constrained entropy maximization problem
    (the machine translation example of Berger et al in
    Computational Linguistics, vol 22, num 1, pp 39--72, 1996.)
    """
    def setUp(self):
        self.xin=[1,1,1,1,1]


    def test_linearmixing(self):
        x = nonlin.linearmixing(F,self.xin,iter=60,alpha=0.5)
        assert nonlin.norm(x)<1e-7
        assert nonlin.norm(F(x))<1e-7

    def test_broyden1(self):
        x= nonlin.broyden1(F,self.xin,iter=11,alpha=1)
        assert nonlin.norm(x)<1e-9
        assert nonlin.norm(F(x))<1e-9

    def test_broyden2(self):
        x= nonlin.broyden2(F,self.xin,iter=12,alpha=1)
        assert nonlin.norm(x)<1e-9
        assert nonlin.norm(F(x))<1e-9

    def test_broyden3(self):
        x= nonlin.broyden3(F,self.xin,iter=12,alpha=1)
        assert nonlin.norm(x)<1e-9
        assert nonlin.norm(F(x))<1e-9

    def test_exciting(self):
        x= nonlin.excitingmixing(F,self.xin,iter=20,alpha=0.5)
        assert nonlin.norm(x)<1e-5
        assert nonlin.norm(F(x))<1e-5

    def test_anderson(self):
        x= nonlin.anderson(F,self.xin,iter=12,alpha=0.03,M=5)
        assert nonlin.norm(x)<0.33

    def test_anderson2(self):
        x= nonlin.anderson2(F,self.xin,iter=12,alpha=0.6,M=5)
        assert nonlin.norm(x)<0.2

    def test_broydengeneralized(self):
        x= nonlin.broyden_generalized(F,self.xin,iter=60,alpha=0.5,M=0)
        assert nonlin.norm(x)<1e-7
        assert nonlin.norm(F(x))<1e-7
        x= nonlin.broyden_generalized(F,self.xin,iter=61,alpha=0.1,M=1)
        assert nonlin.norm(x)<2e-4
        assert nonlin.norm(F(x))<2e-4

    def xtest_broydenmodified(self):
        x= nonlin.broyden_modified(F,self.xin,iter=12,alpha=1)
        assert nonlin.norm(x)<1e-9
        assert nonlin.norm(F(x))<1e-9

    def test_broyden1modified(self):
        x= nonlin.broyden1_modified(F,self.xin,iter=35,alpha=1)
        assert nonlin.norm(x)<1e-9
        assert nonlin.norm(F(x))<1e-9

    def test_vackar(self):
        x= nonlin.vackar(F,self.xin,iter=11,alpha=1)
        assert nonlin.norm(x)<1e-9
        assert nonlin.norm(F(x))<1e-9


if __name__ == "__main__":
    run_module_suite()
