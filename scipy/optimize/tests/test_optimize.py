"""
Unit tests for optimization routines from optimize.py and tnc.py

Authors:
   Ed Schofield, Nov 2005
   Andrew Straw, April 2008

To run it in its simplest form::
  nosetests test_optimize.py

"""

from numpy.testing import assert_raises, assert_almost_equal, \
        assert_equal, assert_, TestCase, run_module_suite

from scipy import optimize
from numpy import array, zeros, float64, dot, log, exp, inf, sin, cos
import numpy as np
from scipy.optimize.tnc import RCSTRINGS, MSG_NONE
import numpy.random
from math import pow

class TestOptimize(TestCase):
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
        self.gradcalls = 0
        self.trace = []


    def func(self, x):
        self.funccalls += 1
        if self.funccalls > 6000:
            raise RuntimeError("too many iterations in optimization routine")
        log_pdot = dot(self.F, x)
        logZ = log(sum(exp(log_pdot)))
        f = logZ - dot(self.K, x)
        self.trace.append(x)
        return f


    def grad(self, x):
        self.gradcalls += 1
        log_pdot = dot(self.F, x)
        logZ = log(sum(exp(log_pdot)))
        p = exp(log_pdot - logZ)
        return dot(self.F.transpose(), p) - self.K


    def test_cg(self):
        """ conjugate gradient optimization routine
        """
        retval = optimize.fmin_cg(self.func, self.startparams, self.grad, (), \
                                  maxiter=self.maxiter, \
                                  full_output=True, disp=False, retall=False)

        (params, fopt, func_calls, grad_calls, warnflag) = retval

        err = abs(self.func(params) - self.func(self.solution))
        #print "CG: Difference is: " + str(err)
        assert_(err < 1e-6)

        # Ensure that function call counts are 'known good'; these are from
        # Scipy 0.7.0. Don't allow them to increase.
        assert_(self.funccalls == 9, self.funccalls)
        assert_(self.gradcalls == 7, self.gradcalls)

        # Ensure that the function behaves the same; this is from Scipy 0.7.0
        assert_(np.allclose(self.trace[2:4],
                           [[0, -0.5, 0.5],
                            [0, -5.05700028e-01, 4.95985862e-01]],
                           atol=1e-14, rtol=1e-7), self.trace[2:4])


    def test_bfgs(self):
        """ Broyden-Fletcher-Goldfarb-Shanno optimization routine
        """
        retval = optimize.fmin_bfgs(self.func, self.startparams, self.grad, \
                                    args=(), maxiter=self.maxiter, \
                                    full_output=True, disp=False, retall=False)

        (params, fopt, gopt, Hopt, func_calls, grad_calls, warnflag) = retval

        err = abs(self.func(params) - self.func(self.solution))
        #print "BFGS: Difference is: " + str(err)
        assert_(err < 1e-6)

        # Ensure that function call counts are 'known good'; these are from
        # Scipy 0.7.0. Don't allow them to increase.
        assert_(self.funccalls == 10, self.funccalls)
        assert_(self.gradcalls == 8, self.gradcalls)

        # Ensure that the function behaves the same; this is from Scipy 0.7.0
        assert_(np.allclose(self.trace[6:8],
                           [[0, -5.25060743e-01,   4.87748473e-01],
                            [0, -5.24885582e-01,   4.87530347e-01]],
                           atol=1e-14, rtol=1e-7), self.trace[6:8])


    def test_bfgs_infinite(self):
        """Test corner case where -Inf is the minimum.  See #1494."""
        func = lambda x: -np.e**-x
        fprime = lambda x: -func(x)
        x0 = [0]
        olderr = np.seterr(over='ignore')
        try:
            x = optimize.fmin_bfgs(func, x0, fprime, disp=False)
            assert_(not np.isfinite(func(x)))
        finally:
            np.seterr(**olderr)


    def test_powell(self):
        """ Powell (direction set) optimization routine
        """
        retval = optimize.fmin_powell(self.func, self.startparams, \
                                    args=(), maxiter=self.maxiter, \
                                    full_output=True, disp=False, retall=False)

        (params, fopt, direc, numiter, func_calls, warnflag) = retval

        err = abs(self.func(params) - self.func(self.solution))
        #print "Powell: Difference is: " + str(err)
        assert_(err < 1e-6)

        # Ensure that function call counts are 'known good'; these are from
        # Scipy 0.7.0. Don't allow them to increase.
        #
        # However, some leeway must be added: the exact evaluation
        # count is sensitive to numerical error, and floating-point
        # computations are not bit-for-bit reproducible across
        # machines, and when using e.g. MKL, data alignment
        # etc. affect the rounding error.
        #
        assert_(self.funccalls <= 116 + 20, self.funccalls)
        assert_(self.gradcalls == 0, self.gradcalls)

        # Ensure that the function behaves the same; this is from Scipy 0.7.0
        assert_(np.allclose(self.trace[34:39],
                           [[ 0.72949016, -0.44156936,  0.47100962],
                            [ 0.72949016, -0.44156936,  0.48052496],
                            [ 1.45898031, -0.88313872,  0.95153458],
                            [ 0.72949016, -0.44156936,  0.47576729],
                            [ 1.72949016, -0.44156936,  0.47576729]],
                           atol=1e-14, rtol=1e-7), self.trace[34:39])

    def test_neldermead(self):
        """ Nelder-Mead simplex algorithm
        """
        retval = optimize.fmin(self.func, self.startparams, \
                                    args=(), maxiter=self.maxiter, \
                                    full_output=True, disp=False, retall=False)

        (params, fopt, numiter, func_calls, warnflag) = retval

        err = abs(self.func(params) - self.func(self.solution))
        #print "Nelder-Mead: Difference is: " + str(err)
        assert_(err < 1e-6)

        # Ensure that function call counts are 'known good'; these are from
        # Scipy 0.7.0. Don't allow them to increase.
        assert_(self.funccalls == 167, self.funccalls)
        assert_(self.gradcalls == 0, self.gradcalls)

        # Ensure that the function behaves the same; this is from Scipy 0.7.0
        assert_(np.allclose(self.trace[76:78],
                           [[0.1928968 , -0.62780447,  0.35166118],
                            [0.19572515, -0.63648426,  0.35838135]],
                           atol=1e-14, rtol=1e-7), self.trace[76:78])

    def test_ncg(self):
        """ line-search Newton conjugate gradient optimization routine
        """
        retval = optimize.fmin_ncg(self.func, self.startparams, self.grad,
                                   args=(), maxiter=self.maxiter,
                                   full_output=False, disp=False,
                                   retall=False)

        params = retval

        err = abs(self.func(params) - self.func(self.solution))
        #print "NCG: Difference is: " + str(err)
        assert_(err < 1e-6)

        # Ensure that function call counts are 'known good'; these are from
        # Scipy 0.7.0. Don't allow them to increase.
        assert_(self.funccalls == 7, self.funccalls)
        assert_(self.gradcalls <= 18, self.gradcalls) # 0.9.0
        #assert_(self.gradcalls == 18, self.gradcalls) # 0.8.0
        #assert_(self.gradcalls == 22, self.gradcalls) # 0.7.0

        # Ensure that the function behaves the same; this is from Scipy 0.7.0
        assert_(np.allclose(self.trace[3:5],
                           [[-4.35700753e-07, -5.24869435e-01, 4.87527480e-01],
                            [-4.35700753e-07, -5.24869401e-01, 4.87527774e-01]],
                           atol=1e-6, rtol=1e-7), self.trace[:5])


    def test_l_bfgs_b(self):
        """ limited-memory bound-constrained BFGS algorithm
        """
        retval = optimize.fmin_l_bfgs_b(self.func, self.startparams,
                                        self.grad, args=(),
                                        maxfun=self.maxiter)

        (params, fopt, d) = retval

        err = abs(self.func(params) - self.func(self.solution))
        #print "LBFGSB: Difference is: " + str(err)
        assert_(err < 1e-6)

        # Ensure that function call counts are 'known good'; these are from
        # Scipy 0.7.0. Don't allow them to increase.
        assert_(self.funccalls == 7, self.funccalls)
        assert_(self.gradcalls == 5, self.gradcalls)

        # Ensure that the function behaves the same; this is from Scipy 0.7.0
        assert_(np.allclose(self.trace[3:5],
                           [[0.        , -0.52489628,  0.48753042],
                            [0.        , -0.52489628,  0.48753042]],
                           atol=1e-14, rtol=1e-7), self.trace[3:5])

    def test_brent(self):
        """ brent algorithm
        """
        x = optimize.brent(lambda x: (x-1.5)**2-0.8)
        err1 = abs(x - 1.5)
        x = optimize.brent(lambda x: (x-1.5)**2-0.8, brack = (-3,-2))
        err2 = abs(x - 1.5)
        x = optimize.brent(lambda x: (x-1.5)**2-0.8, full_output=True)
        err3 = abs(x[0] - 1.5)
        x = optimize.brent(lambda x: (x-1.5)**2-0.8, brack = (-15,-1,15))
        err4 = abs(x - 1.5)

        assert_(max((err1,err2,err3,err4)) < 1e-6)


    def test_fminbound(self):
        """Test fminbound
        """
        x = optimize.fminbound(lambda x: (x - 1.5)**2 - 0.8, 0, 1)
        assert_(abs(x - 1) < 1e-5)
        x = optimize.fminbound(lambda x: (x - 1.5)**2 - 0.8, 1, 5)
        assert_(abs(x - 1.5) < 1e-6)
        x = optimize.fminbound(lambda x: (x - 1.5)**2 - 0.8,
                               numpy.array([1]), numpy.array([5]))
        assert_(abs(x - 1.5) < 1e-6)
        assert_raises(ValueError,
                optimize.fminbound, lambda x: (x - 1.5)**2 - 0.8, 5, 1)

    def test_fminbound_scalar(self):
        assert_raises(ValueError,
                      optimize.fminbound, lambda x: (x - 1.5)**2 - 0.8,
                      np.zeros(2), 1)

        assert_almost_equal(
            optimize.fminbound(lambda x: (x - 1.5)**2 - 0.8, 1, np.array(5)),
            1.5)


class TestTnc(TestCase):
    """TNC non-linear optimization.

    These tests are taken from Prof. K. Schittkowski's test examples
    for constrained non-linear programming.

    http://www.uni-bayreuth.de/departments/math/~kschittkowski/home.htm

    """
    tests = []

    def setUp(self):
        def test1fg(x):
            f = 100.0*pow((x[1]-pow(x[0],2)),2)+pow(1.0-x[0],2)
            dif = [0,0]
            dif[1] = 200.0*(x[1]-pow(x[0],2))
            dif[0] = -2.0*(x[0]*(dif[1]-1.0)+1.0)
            return f, dif
        self.tests.append((test1fg, [-2,1], ([-inf,None],[-1.5,None]),
                           [1,1]))
        def test2fg(x):
            f = 100.0*pow((x[1]-pow(x[0],2)),2)+pow(1.0-x[0],2)
            dif = [0,0]
            dif[1] = 200.0*(x[1]-pow(x[0],2))
            dif[0] = -2.0*(x[0]*(dif[1]-1.0)+1.0)
            return f, dif
        self.tests.append((test2fg, [-2,1], [(-inf,None),(1.5,None)],
                           [-1.2210262419616387,1.5]))

        def test3fg(x):
            f = x[1]+pow(x[1]-x[0],2)*1.0e-5
            dif = [0,0]
            dif[0] = -2.0*(x[1]-x[0])*1.0e-5
            dif[1] = 1.0-dif[0]
            return f, dif
        self.tests.append((test3fg, [10,1], [(-inf,None),(0.0, None)],
                           [0,0]))

        def test4fg(x):
            f = pow(x[0]+1.0,3)/3.0+x[1]
            dif = [0,0]
            dif[0] = pow(x[0]+1.0,2)
            dif[1] = 1.0
            return f, dif
        self.tests.append((test4fg, [1.125,0.125], [(1, None),(0, None)],
                           [1,0]))

        def test5fg(x):
            f = sin(x[0]+x[1])+pow(x[0]-x[1],2)-1.5*x[0]+2.5*x[1]+1.0
            dif = [0,0]
            v1 = cos(x[0]+x[1])
            v2 = 2.0*(x[0]-x[1])

            dif[0] = v1+v2-1.5
            dif[1] = v1-v2+2.5
            return f, dif
        self.tests.append((test5fg, [0,0], [(-1.5, 4),(-3,3)],
                           [-0.54719755119659763, -1.5471975511965976]))

        def test38fg(x):
            f = (100.0*pow(x[1]-pow(x[0],2),2) + \
                 pow(1.0-x[0],2)+90.0*pow(x[3]-pow(x[2],2),2) + \
                 pow(1.0-x[2],2)+10.1*(pow(x[1]-1.0,2)+pow(x[3]-1.0,2)) + \
                 19.8*(x[1]-1.0)*(x[3]-1.0))*1.0e-5
            dif = [0,0,0,0]
            dif[0] = (-400.0*x[0]*(x[1]-pow(x[0],2))-2.0*(1.0-x[0]))*1.0e-5
            dif[1] = (200.0*(x[1]-pow(x[0],2))+20.2 \
                      *(x[1]-1.0)+19.8*(x[3]-1.0))*1.0e-5
            dif[2] = (-360.0*x[2]*(x[3]-pow(x[2],2))-2.0\
                      *(1.0-x[2]))*1.0e-5
            dif[3] = (180.0*(x[3]-pow(x[2],2))+20.2\
                      *(x[3]-1.0)+19.8*(x[1]-1.0))*1.0e-5
            return f, dif
        self.tests.append((test38fg, array([-3,-1,-3,-1]), [(-10,10)]*4, [1]*4))

        def test45fg(x):
            f = 2.0-x[0]*x[1]*x[2]*x[3]*x[4]/120.0
            dif = [0]*5
            dif[0] = -x[1]*x[2]*x[3]*x[4]/120.0
            dif[1] = -x[0]*x[2]*x[3]*x[4]/120.0
            dif[2] = -x[0]*x[1]*x[3]*x[4]/120.0
            dif[3] = -x[0]*x[1]*x[2]*x[4]/120.0
            dif[4] = -x[0]*x[1]*x[2]*x[3]/120.0
            return f, dif
        self.tests.append((test45fg, [2]*5, [(0,1),(0,2),(0,3),(0,4),(0,5)],
                           [1,2,3,4,5]))

    def test_tnc(self):
        for fg, x, bounds, xopt in self.tests:
            x, nf, rc = optimize.fmin_tnc(fg, x, bounds=bounds,
                                          messages=MSG_NONE, maxfun=200)
            err = "Failed optimization of %s.\n" \
                  "After %d function evaluations, TNC returned: %s.""" % \
                  (fg.__name__, nf, RCSTRINGS[rc])

        ef = abs(fg(xopt)[0] - fg(x)[0])
        if ef > 1e-8:
            raise err


class TestRosen(TestCase):

    def test_hess(self):
        """Compare rosen_hess(x) times p with rosen_hess_prod(x,p) (ticket #1248)"""
        x = array([3, 4, 5])
        p = array([2, 2, 2])
        hp = optimize.rosen_hess_prod(x, p)
        dothp = np.dot(optimize.rosen_hess(x), p)
        assert_equal(hp, dothp)


if __name__ == "__main__":
    run_module_suite()
