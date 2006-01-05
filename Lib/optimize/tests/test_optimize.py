""" Unit tests for optimization routines
Author: Ed Schofield
Nov 2005
"""

from numpy.testing import *

set_package_path()
from scipy import optimize
from numpy import array, zeros, float64, dot, log, exp
restore_path()


class test_optimize(ScipyTestCase):
    """ Test case for a simple constrained entropy maximization problem
    (the machine translation example of Berger et al in
    Computational Linguistics, vol 22, num 1, pp 39--72, 1996.)
    """
    def setUp(self):
        self.F = array([[1,1,1],[1,1,0],[1,0,1],[1,0,0],[1,0,0]])
        self.K = array([1., 0.3, 0.5])
        self.startparams = zeros(3, float64)
        self.solution = array([0., -0.524869316, 0.487525860])
        self.maxiter = 1000
        self.funccalls = 0
    

    def func(self, x):
        self.funccalls += 1
        if self.funccalls > 6000:
            raise RuntimeError, "too many iterations in optimization routine"
        log_pdot = dot(self.F, x)
        logZ = log(sum(exp(log_pdot)))
        f = logZ - dot(self.K, x)
        return f
    

    def grad(self, x):
        log_pdot = dot(self.F, x)
        logZ = log(sum(exp(log_pdot)))
        p = exp(log_pdot - logZ)
        return dot(self.F.transpose(), p) - self.K
    

    def check_cg(self):
        """ conjugate gradient optimization routine
        """
        retval = optimize.fmin_cg(self.func, self.startparams, self.grad, (), \
                                  maxiter=self.maxiter, \
                                  full_output=True, disp=False, retall=False)
                
        (params, fopt, func_calls, grad_calls, warnflag) = retval

        err = abs(self.func(params) - self.func(self.solution))
        #print "CG: Difference is: " + str(err)
        assert err < 1e-6


    def check_bfgs(self):
        """ Broyden-Fletcher-Goldfarb-Shanno optimization routine
        """
        retval = optimize.fmin_bfgs(self.func, self.startparams, self.grad, \
                                    args=(), maxiter=self.maxiter, \
                                    full_output=True, disp=False, retall=False)
                
        (params, fopt, gopt, Hopt, func_calls, grad_calls, warnflag) = retval

        err = abs(self.func(params) - self.func(self.solution))
        #print "BFGS: Difference is: " + str(err)
        assert err < 1e-6


    def check_powell(self):
        """ Powell (direction set) optimization routine
        """
        retval = optimize.fmin_powell(self.func, self.startparams, \
                                    args=(), maxiter=self.maxiter, \
                                    full_output=True, disp=False, retall=False)
                
        (params, fopt, direc, numiter, func_calls, warnflag) = retval

        err = abs(self.func(params) - self.func(self.solution))
        #print "Powell: Difference is: " + str(err)
        assert err < 1e-6

    def check_neldermead(self):
        """ Nelder-Mead simplex algorithm
        """
        retval = optimize.fmin(self.func, self.startparams, \
                                    args=(), maxiter=self.maxiter, \
                                    full_output=True, disp=False, retall=False)
                
        (params, fopt, numiter, func_calls, warnflag) = retval

        err = abs(self.func(params) - self.func(self.solution))
        #print "Nelder-Mead: Difference is: " + str(err)
        assert err < 1e-6

    def check_ncg(self):
        """ line-search Newton conjugate gradient optimization routine
        """
        retval = optimize.fmin_ncg(self.func, self.startparams, self.grad,\
                                    args=(), maxiter=self.maxiter, \
                                    full_output=False, disp=False, retall=False)
                
        params = retval

        err = abs(self.func(params) - self.func(self.solution))
        #print "NCG: Difference is: " + str(err)
        assert err < 1e-6


    def check_l_bfgs_b(self):
        """ limited-memory bound-constrained BFGS algorithm
        """
        retval = optimize.fmin_l_bfgs_b(self.func, self.startparams, self.grad,\
                                    args=(), maxfun=self.maxiter)
                
        (params, fopt, d) = retval

        err = abs(self.func(params) - self.func(self.solution))
        #print "LBFGSB: Difference is: " + str(err)
        assert err < 1e-6


if __name__ == "__main__":
    ScipyTest().run()
