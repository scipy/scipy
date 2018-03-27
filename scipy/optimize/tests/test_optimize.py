"""
Unit tests for optimization routines from optimize.py

Authors:
   Ed Schofield, Nov 2005
   Andrew Straw, April 2008

To run it in its simplest form::
  nosetests test_optimize.py

"""
from __future__ import division, print_function, absolute_import

import itertools

import numpy as np
from numpy.testing import (assert_allclose, assert_equal,
                           assert_,
                           assert_almost_equal, assert_warns,
                           assert_array_less)
import pytest
from pytest import raises as assert_raises

from scipy._lib._numpy_compat import suppress_warnings
from scipy import optimize
from scipy.optimize.optimize import (Function, NelderMead, Optimizer, BFGS)
from scipy.optimize.lbfgsb import LBFGSB


def logit(x):
    return 1 / (1 + np.exp(-x))


def der_logit(x):
    return np.exp(-x) / (1 + np.exp(-x)) ** 2


def test_check_grad():
    # Verify if check_grad is able to estimate the derivative of the
    # logistic function.
    x0 = np.array([1.5])

    r = optimize.check_grad(logit, der_logit, x0)
    assert_almost_equal(r, 0)

    r = optimize.check_grad(logit, der_logit, x0, epsilon=1e-6)
    assert_almost_equal(r, 0)

    # Check if the epsilon parameter is being considered.
    r = abs(optimize.check_grad(logit, der_logit, x0, epsilon=1e-1) - 0)
    assert_(r > 1e-7)


class CheckOptimize(object):
    """ Base test case for a simple constrained entropy maximization problem
    (the machine translation example of Berger et al in
    Computational Linguistics, vol 22, num 1, pp 39--72, 1996.)
    """
    def setup_method(self):
        self.F = np.array([[1,1,1],[1,1,0],[1,0,1],[1,0,0],[1,0,0]])
        self.K = np.array([1., 0.3, 0.5])
        self.startparams = np.zeros(3, np.float64)
        self.solution = np.array([0., -0.524869316, 0.487525860])
        self.maxiter = 1000
        self.funccalls = 0
        self.gradcalls = 0
        self.trace = []

    def func(self, x):
        self.funccalls += 1
        if self.funccalls > 6000:
            raise RuntimeError("too many iterations in optimization routine")
        log_pdot = np.dot(self.F, x)
        logZ = np.log(sum(np.exp(log_pdot)))
        f = logZ - np.dot(self.K, x)
        self.trace.append(x)
        return f

    def grad(self, x):
        self.gradcalls += 1
        log_pdot = np.dot(self.F, x)
        logZ = np.log(sum(np.exp(log_pdot)))
        p = np.exp(log_pdot - logZ)
        return np.dot(self.F.transpose(), p) - self.K

    def hess(self, x):
        log_pdot = np.dot(self.F, x)
        logZ = np.log(sum(np.exp(log_pdot)))
        p = np.exp(log_pdot - logZ)
        return np.dot(self.F.T,
                      np.dot(np.diag(p), self.F - np.dot(self.F.T, p)))

    def hessp(self, x, p):
        return np.dot(self.hess(x), p)


class CheckOptimizeParameterized(CheckOptimize):

    def test_cg(self):
        # conjugate gradient optimization routine
        if self.use_wrapper:
            opts = {'maxiter': self.maxiter, 'disp': self.disp,
                    'return_all': False}
            res = optimize.minimize(self.func, self.startparams, args=(),
                                    method='CG', jac=self.grad,
                                    options=opts)
            params, fopt, func_calls, grad_calls, warnflag = \
                res['x'], res['fun'], res['nfev'], res['njev'], res['status']
        else:
            retval = optimize.fmin_cg(self.func, self.startparams,
                                      self.grad, (), maxiter=self.maxiter,
                                      full_output=True, disp=self.disp,
                                      retall=False)
            (params, fopt, func_calls, grad_calls, warnflag) = retval

        assert_allclose(self.func(params), self.func(self.solution),
                        atol=1e-6)

        # Ensure that function call counts are 'known good'; these are from
        # Scipy 0.7.0. Don't allow them to increase.
        assert_(self.funccalls == 9, self.funccalls)
        assert_(self.gradcalls == 7, self.gradcalls)

        # Ensure that the function behaves the same; this is from Scipy 0.7.0
        assert_allclose(self.trace[2:4],
                        [[0, -0.5, 0.5],
                         [0, -5.05700028e-01, 4.95985862e-01]],
                        atol=1e-14, rtol=1e-7)

    def test_cg_cornercase(self):
        def f(r):
            return 2.5 * (1 - np.exp(-1.5*(r - 0.5)))**2

        # Check several initial guesses. (Too far away from the
        # minimum, the function ends up in the flat region of exp.)
        for x0 in np.linspace(-0.75, 3, 71):
            sol = optimize.minimize(f, [x0], method='CG')
            assert_(sol.success)
            assert_allclose(sol.x, [0.5], rtol=1e-5)

    def test_bfgs(self):
        # Broyden-Fletcher-Goldfarb-Shanno optimization routine
        if self.use_wrapper:
            opts = {'maxiter': self.maxiter, 'disp': self.disp,
                    'return_all': False}
            res = optimize.minimize(self.func, self.startparams,
                                    jac=self.grad, method='BFGS', args=(),
                                    options=opts)

            params, fopt, gopt, Hopt, func_calls, grad_calls, warnflag = (
                    res['x'], res['fun'], res['jac'], res['hess_inv'],
                    res['nfev'], res['njev'], res['status'])
        else:
            retval = optimize.fmin_bfgs(self.func, self.startparams, self.grad,
                                        args=(), maxiter=self.maxiter,
                                        full_output=True, disp=self.disp,
                                        retall=False)
            (params, fopt, gopt, Hopt, func_calls, grad_calls, warnflag) = retval

        assert_allclose(self.func(params), self.func(self.solution),
                        atol=1e-6)

        # Ensure that function call counts are 'known good'; these are from
        # Scipy 0.7.0. Don't allow them to increase.
        assert_(self.funccalls == 10, self.funccalls)
        assert_(self.gradcalls == 8, self.gradcalls)

        # Ensure that the function behaves the same; this is from Scipy 0.7.0
        assert_allclose(self.trace[6:8],
                        [[0, -5.25060743e-01, 4.87748473e-01],
                         [0, -5.24885582e-01, 4.87530347e-01]],
                        atol=1e-14, rtol=1e-7)

    def test_bfgs_infinite(self):
        # Test corner case where -Inf is the minimum.  See gh-2019.
        func = lambda x: -np.e**-x
        fprime = lambda x: -func(x)
        x0 = [0]
        olderr = np.seterr(over='ignore')
        try:
            if self.use_wrapper:
                opts = {'disp': self.disp}
                x = optimize.minimize(func, x0, jac=fprime, method='BFGS',
                                      args=(), options=opts)['x']
            else:
                x = optimize.fmin_bfgs(func, x0, fprime, disp=self.disp)
            assert_(not np.isfinite(func(x)))
        finally:
            np.seterr(**olderr)

    def test_powell(self):
        # Powell (direction set) optimization routine
        if self.use_wrapper:
            opts = {'maxiter': self.maxiter, 'disp': self.disp,
                    'return_all': False}
            res = optimize.minimize(self.func, self.startparams, args=(),
                                    method='Powell', options=opts)
            params, fopt, direc, numiter, func_calls, warnflag = (
                    res['x'], res['fun'], res['direc'], res['nit'],
                    res['nfev'], res['status'])
        else:
            retval = optimize.fmin_powell(self.func, self.startparams,
                                          args=(), maxiter=self.maxiter,
                                          full_output=True, disp=self.disp,
                                          retall=False)
            (params, fopt, direc, numiter, func_calls, warnflag) = retval

        assert_allclose(self.func(params), self.func(self.solution),
                        atol=1e-6)

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
        assert_allclose(self.trace[34:39],
                        [[0.72949016, -0.44156936, 0.47100962],
                         [0.72949016, -0.44156936, 0.48052496],
                         [1.45898031, -0.88313872, 0.95153458],
                         [0.72949016, -0.44156936, 0.47576729],
                         [1.72949016, -0.44156936, 0.47576729]],
                        atol=1e-14, rtol=1e-7)

    def test_neldermead(self):
        # Nelder-Mead simplex algorithm
        if self.use_wrapper:
            opts = {'maxiter': self.maxiter, 'disp': self.disp,
                    'return_all': False}
            res = optimize.minimize(self.func, self.startparams, args=(),
                                    method='Nelder-mead', options=opts)
            params, fopt, numiter, func_calls, warnflag, final_simplex = (
                    res['x'], res['fun'], res['nit'], res['nfev'],
                    res['status'], res['final_simplex'])
        else:
            retval = optimize.fmin(self.func, self.startparams,
                                        args=(), maxiter=self.maxiter,
                                        full_output=True, disp=self.disp,
                                        retall=False)
            (params, fopt, numiter, func_calls, warnflag) = retval

        assert_allclose(self.func(params), self.func(self.solution),
                        atol=1e-6)

        # Ensure that function call counts are 'known good'; these are from
        # Scipy 0.7.0. Don't allow them to increase.
        assert_(self.funccalls == 167, self.funccalls)
        assert_(self.gradcalls == 0, self.gradcalls)

        # Ensure that the function behaves the same; this is from Scipy 0.7.0
        assert_allclose(self.trace[76:78],
                        [[0.1928968, -0.62780447, 0.35166118],
                         [0.19572515, -0.63648426, 0.35838135]],
                        atol=1e-14, rtol=1e-7)

    def test_neldermead_initial_simplex(self):
        # Nelder-Mead simplex algorithm
        simplex = np.zeros((4, 3))
        simplex[...] = self.startparams
        for j in range(3):
            simplex[j+1,j] += 0.1

        if self.use_wrapper:
            opts = {'maxiter': self.maxiter, 'disp': False,
                    'return_all': True, 'initial_simplex': simplex}
            res = optimize.minimize(self.func, self.startparams, args=(),
                                    method='Nelder-mead', options=opts)
            params, fopt, numiter, func_calls, warnflag = \
                    res['x'], res['fun'], res['nit'], res['nfev'], \
                    res['status']
            assert_allclose(res['allvecs'][0], simplex[0])
        else:
            retval = optimize.fmin(self.func, self.startparams,
                                        args=(), maxiter=self.maxiter,
                                        full_output=True, disp=False, retall=False,
                                        initial_simplex=simplex)

            (params, fopt, numiter, func_calls, warnflag) = retval

        assert_allclose(self.func(params), self.func(self.solution),
                        atol=1e-6)

        # Ensure that function call counts are 'known good'; these are from
        # Scipy 0.17.0. Don't allow them to increase.
        assert_(self.funccalls == 100, self.funccalls)
        assert_(self.gradcalls == 0, self.gradcalls)

        # Ensure that the function behaves the same; this is from Scipy 0.15.0
        assert_allclose(self.trace[50:52],
                        [[0.14687474, -0.5103282, 0.48252111],
                         [0.14474003, -0.5282084, 0.48743951]],
                        atol=1e-14, rtol=1e-7)

    def test_neldermead_initial_simplex_bad(self):
        # Check it fails with a bad simplices
        bad_simplices = []

        simplex = np.zeros((3, 2))
        simplex[...] = self.startparams[:2]
        for j in range(2):
            simplex[j+1,j] += 0.1
        bad_simplices.append(simplex)

        simplex = np.zeros((3, 3))
        bad_simplices.append(simplex)

        for simplex in bad_simplices:
            if self.use_wrapper:
                opts = {'maxiter': self.maxiter, 'disp': False,
                        'return_all': False, 'initial_simplex': simplex}
                assert_raises(ValueError,
                              optimize.minimize, self.func, self.startparams, args=(),
                              method='Nelder-mead', options=opts)
            else:
                assert_raises(ValueError, optimize.fmin, self.func, self.startparams,
                              args=(), maxiter=self.maxiter,
                              full_output=True, disp=False, retall=False,
                              initial_simplex=simplex)

    def test_ncg_negative_maxiter(self):
        # Regression test for gh-8241
        opts = {'maxiter': -1}
        result = optimize.minimize(self.func, self.startparams,
                                   method='Newton-CG', jac=self.grad,
                                   args=(), options=opts)
        assert_(result.status == 1)

    def test_ncg(self):
        # line-search Newton conjugate gradient optimization routine
        if self.use_wrapper:
            opts = {'maxiter': self.maxiter, 'disp': self.disp,
                    'return_all': False}
            retval = optimize.minimize(self.func, self.startparams,
                                       method='Newton-CG', jac=self.grad,
                                       args=(), options=opts)['x']
        else:
            retval = optimize.fmin_ncg(self.func, self.startparams, self.grad,
                                       args=(), maxiter=self.maxiter,
                                       full_output=False, disp=self.disp,
                                       retall=False)

        params = retval

        assert_allclose(self.func(params), self.func(self.solution),
                        atol=1e-6)

        # Ensure that function call counts are 'known good'; these are from
        # Scipy 0.7.0. Don't allow them to increase.
        assert_(self.funccalls == 7, self.funccalls)
        assert_(self.gradcalls <= 22, self.gradcalls)  # 0.13.0
        #assert_(self.gradcalls <= 18, self.gradcalls) # 0.9.0
        #assert_(self.gradcalls == 18, self.gradcalls) # 0.8.0
        #assert_(self.gradcalls == 22, self.gradcalls) # 0.7.0

        # Ensure that the function behaves the same; this is from Scipy 0.7.0
        assert_allclose(self.trace[3:5],
                        [[-4.35700753e-07, -5.24869435e-01, 4.87527480e-01],
                         [-4.35700753e-07, -5.24869401e-01, 4.87527774e-01]],
                        atol=1e-6, rtol=1e-7)

    def test_ncg_hess(self):
        # Newton conjugate gradient with Hessian
        if self.use_wrapper:
            opts = {'maxiter': self.maxiter, 'disp': self.disp,
                    'return_all': False}
            retval = optimize.minimize(self.func, self.startparams,
                                       method='Newton-CG', jac=self.grad,
                                       hess=self.hess,
                                       args=(), options=opts)['x']
        else:
            retval = optimize.fmin_ncg(self.func, self.startparams, self.grad,
                                       fhess=self.hess,
                                       args=(), maxiter=self.maxiter,
                                       full_output=False, disp=self.disp,
                                       retall=False)

        params = retval

        assert_allclose(self.func(params), self.func(self.solution),
                        atol=1e-6)

        # Ensure that function call counts are 'known good'; these are from
        # Scipy 0.7.0. Don't allow them to increase.
        assert_(self.funccalls == 7, self.funccalls)
        assert_(self.gradcalls <= 18, self.gradcalls)  # 0.9.0
        # assert_(self.gradcalls == 18, self.gradcalls) # 0.8.0
        # assert_(self.gradcalls == 22, self.gradcalls) # 0.7.0

        # Ensure that the function behaves the same; this is from Scipy 0.7.0
        assert_allclose(self.trace[3:5],
                        [[-4.35700753e-07, -5.24869435e-01, 4.87527480e-01],
                         [-4.35700753e-07, -5.24869401e-01, 4.87527774e-01]],
                        atol=1e-6, rtol=1e-7)

    def test_ncg_hessp(self):
        # Newton conjugate gradient with Hessian times a vector p.
        if self.use_wrapper:
            opts = {'maxiter': self.maxiter, 'disp': self.disp,
                    'return_all': False}
            retval = optimize.minimize(self.func, self.startparams,
                                       method='Newton-CG', jac=self.grad,
                                       hessp=self.hessp,
                                       args=(), options=opts)['x']
        else:
            retval = optimize.fmin_ncg(self.func, self.startparams, self.grad,
                                       fhess_p=self.hessp,
                                       args=(), maxiter=self.maxiter,
                                       full_output=False, disp=self.disp,
                                       retall=False)

        params = retval

        assert_allclose(self.func(params), self.func(self.solution),
                        atol=1e-6)

        # Ensure that function call counts are 'known good'; these are from
        # Scipy 0.7.0. Don't allow them to increase.
        assert_(self.funccalls == 7, self.funccalls)
        assert_(self.gradcalls <= 18, self.gradcalls)  # 0.9.0
        # assert_(self.gradcalls == 18, self.gradcalls) # 0.8.0
        # assert_(self.gradcalls == 22, self.gradcalls) # 0.7.0

        # Ensure that the function behaves the same; this is from Scipy 0.7.0
        assert_allclose(self.trace[3:5],
                        [[-4.35700753e-07, -5.24869435e-01, 4.87527480e-01],
                         [-4.35700753e-07, -5.24869401e-01, 4.87527774e-01]],
                        atol=1e-6, rtol=1e-7)


def test_neldermead_xatol_fatol():
    # gh4484
    # test we can call with fatol, xatol specified
    func = lambda x: x[0]**2 + x[1]**2

    optimize._minimize._minimize_neldermead(func, [1, 1], maxiter=2,
                                            xatol=1e-3, fatol=1e-3)
    assert_warns(DeprecationWarning,
                 optimize._minimize._minimize_neldermead,
                 func, [1, 1], xtol=1e-3, ftol=1e-3, maxiter=2)

def test_neldermead_adaptive():
    func = lambda x: np.sum(x**2)
    p0 = [0.15746215, 0.48087031, 0.44519198, 0.4223638, 0.61505159, 0.32308456,
      0.9692297, 0.4471682, 0.77411992, 0.80441652, 0.35994957, 0.75487856,
      0.99973421, 0.65063887, 0.09626474]
   
    res = optimize.minimize(func, p0, method='Nelder-Mead')
    assert_equal(res.success, False)
 
    res = optimize.minimize(func, p0, method='Nelder-Mead',
                    options={'adaptive':True})
    assert_equal(res.success, True)

class TestOptimizeWrapperDisp(CheckOptimizeParameterized):
    use_wrapper = True
    disp = True


class TestOptimizeWrapperNoDisp(CheckOptimizeParameterized):
    use_wrapper = True
    disp = False


class TestOptimizeNoWrapperDisp(CheckOptimizeParameterized):
    use_wrapper = False
    disp = True


class TestOptimizeNoWrapperNoDisp(CheckOptimizeParameterized):
    use_wrapper = False
    disp = False


class TestOptimizeSimple(CheckOptimize):

    def test_bfgs_nan(self):
        # Test corner case where nan is fed to optimizer.  See gh-2067.
        func = lambda x: x
        fprime = lambda x: np.ones_like(x)
        x0 = [np.nan]
        with np.errstate(over='ignore', invalid='ignore'):
            x = optimize.fmin_bfgs(func, x0, fprime, disp=False)
            assert_(np.isnan(func(x)))

    def test_bfgs_nan_return(self):
        # Test corner cases where fun returns NaN. See gh-4793.

        # First case: NaN from first call.
        func = lambda x: np.nan
        with np.errstate(invalid='ignore'):
            result = optimize.minimize(func, 0)

        assert_(np.isnan(result['fun']))
        assert_(result['success'] is False)

        # Second case: NaN from second call.
        func = lambda x: 0 if x == 0 else np.nan
        fprime = lambda x: np.ones_like(x)  # Steer away from zero.
        with np.errstate(invalid='ignore'):
            result = optimize.minimize(func, 0, jac=fprime)

        assert_(np.isnan(result['fun']))
        assert_(result['success'] is False)

    def test_bfgs_numerical_jacobian(self):
        # BFGS with numerical jacobian and a vector epsilon parameter.
        # define the epsilon parameter using a random vector
        epsilon = np.sqrt(np.finfo(float).eps) * np.random.rand(len(self.solution))

        params = optimize.fmin_bfgs(self.func, self.startparams,
                                    epsilon=epsilon, args=(),
                                    maxiter=self.maxiter, disp=False)

        assert_allclose(self.func(params), self.func(self.solution),
                        atol=1e-6)

    def test_bfgs_gh_2169(self):
        def f(x):
            if x < 0:
                return 1.79769313e+308
            else:
                return x + 1./x
        xs = optimize.fmin_bfgs(f, [10.], disp=False)
        assert_allclose(xs, 1.0, rtol=1e-4, atol=1e-4)

    def test_l_bfgs_b(self):
        # limited-memory bound-constrained BFGS algorithm
        retval = optimize.fmin_l_bfgs_b(self.func, self.startparams,
                                        self.grad, args=(),
                                        maxiter=self.maxiter)

        (params, fopt, d) = retval

        assert_allclose(self.func(params), self.func(self.solution),
                        atol=1e-6)

        # Ensure that function call counts are 'known good'; these are from
        # Scipy 0.7.0. Don't allow them to increase.
        assert_(self.funccalls == 7, self.funccalls)
        assert_(self.gradcalls == 5, self.gradcalls)

        # Ensure that the function behaves the same; this is from Scipy 0.7.0
        assert_allclose(self.trace[3:5],
                        [[0., -0.52489628, 0.48753042],
                         [0., -0.52489628, 0.48753042]],
                        atol=1e-14, rtol=1e-7)

    def test_l_bfgs_b_numjac(self):
        # L-BFGS-B with numerical jacobian
        retval = optimize.fmin_l_bfgs_b(self.func, self.startparams,
                                        approx_grad=True,
                                        maxiter=self.maxiter)

        (params, fopt, d) = retval

        assert_allclose(self.func(params), self.func(self.solution),
                        atol=1e-6)

    def test_l_bfgs_b_funjac(self):
        # L-BFGS-B with combined objective function and jacobian
        def fun(x):
            return self.func(x), self.grad(x)

        retval = optimize.fmin_l_bfgs_b(fun, self.startparams,
                                        maxiter=self.maxiter)

        (params, fopt, d) = retval

        assert_allclose(self.func(params), self.func(self.solution),
                        atol=1e-6)

    def test_l_bfgs_b_maxiter(self):
        # gh7854
        # Ensure that not more than maxiters are ever run.
        class Callback(object):
            def __init__(self):
                self.nit = 0
                self.fun = None
                self.x = None

            def __call__(self, x):
                self.x = x
                self.fun = optimize.rosen(x)
                self.nit += 1

        c = Callback()
        res = optimize.minimize(optimize.rosen, [0., 0.], method='l-bfgs-b',
                                callback=c, options={'maxiter': 5})

        assert_equal(res.nit, 5)
        assert_equal(res.nit, c.nit)
        assert_almost_equal(res.x, c.x)
        assert_almost_equal(res.fun, c.fun)
        assert_equal(res.status, 1)
        assert_(res.success is False)
        assert_equal(res.message, 'STOP: TOTAL NO. of ITERATIONS REACHED LIMIT')

    def test_minimize_l_bfgs_b(self):
        # Minimize with L-BFGS-B method
        opts = {'disp': False, 'maxiter': self.maxiter}
        r = optimize.minimize(self.func, self.startparams,
                              method='L-BFGS-B', jac=self.grad,
                              options=opts)
        assert_allclose(self.func(r.x), self.func(self.solution),
                        atol=1e-6)
        # approximate jacobian
        ra = optimize.minimize(self.func, self.startparams,
                               method='L-BFGS-B', options=opts)
        assert_allclose(self.func(ra.x), self.func(self.solution),
                        atol=1e-6)
        # check that function evaluations in approximate jacobian are counted
        assert_(ra.nfev > r.nfev)

    def test_minimize_l_bfgs_b_ftol(self):
        # Check that the `ftol` parameter in l_bfgs_b works as expected
        v0 = None
        for tol in [1e-1, 1e-4, 1e-7, 1e-10]:
            opts = {'disp': False, 'maxiter': self.maxiter, 'ftol': tol}
            sol = optimize.minimize(self.func, self.startparams,
                                    method='L-BFGS-B', jac=self.grad,
                                    options=opts)
            v = self.func(sol.x)

            if v0 is None:
                v0 = v
            else:
                assert_(v < v0)

            assert_allclose(v, self.func(self.solution), rtol=tol)

    def test_minimize_l_bfgs_maxls(self):
        # check that the maxls is passed down to the Fortran routine
        sol = optimize.minimize(optimize.rosen, np.array([-1.2,1.0]),
                                method='L-BFGS-B', jac=optimize.rosen_der,
                                options={'disp': False, 'maxls': 1})
        assert_(not sol.success)

    def test_minimize_l_bfgs_b_maxfun_interruption(self):
        # gh-6162
        f = optimize.rosen
        g = optimize.rosen_der
        values = []
        x0 = np.ones(7) * 1000

        def objfun(x):
            value = f(x)
            values.append(value)
            return value

        # Look for an interesting test case.
        # Request a maxfun that stops at a particularly bad function
        # evaluation somewhere between 100 and 300 evaluations.
        low, medium, high = 30, 100, 300
        optimize.fmin_l_bfgs_b(objfun, x0, fprime=g, maxfun=high)
        v, k = max((y, i) for i, y in enumerate(values[medium:]))
        maxfun = medium + k
        # If the minimization strategy is reasonable,
        # the minimize() result should not be worse than the best
        # of the first 30 function evaluations.
        target = min(values[:low])
        xmin, fmin, d = optimize.fmin_l_bfgs_b(f, x0, fprime=g, maxfun=maxfun)
        assert_array_less(fmin, target)

    def test_custom(self):
        # This function comes from the documentation example.
        def custmin(fun, x0, args=(), maxfev=None, stepsize=0.1,
                maxiter=100, callback=None, **options):
            bestx = x0
            besty = fun(x0)
            funcalls = 1
            niter = 0
            improved = True
            stop = False

            while improved and not stop and niter < maxiter:
                improved = False
                niter += 1
                for dim in range(np.size(x0)):
                    for s in [bestx[dim] - stepsize, bestx[dim] + stepsize]:
                        testx = np.copy(bestx)
                        testx[dim] = s
                        testy = fun(testx, *args)
                        funcalls += 1
                        if testy < besty:
                            besty = testy
                            bestx = testx
                            improved = True
                    if callback is not None:
                        callback(bestx)
                    if maxfev is not None and funcalls >= maxfev:
                        stop = True
                        break

            return optimize.OptimizeResult(fun=besty, x=bestx, nit=niter,
                                           nfev=funcalls, success=(niter > 1))

        x0 = [1.35, 0.9, 0.8, 1.1, 1.2]
        res = optimize.minimize(optimize.rosen, x0, method=custmin,
                                options=dict(stepsize=0.05))
        assert_allclose(res.x, 1.0, rtol=1e-4, atol=1e-4)

    def test_minimize_tol_parameter(self):
        # Check that the minimize() tol= argument does something
        def func(z):
            x, y = z
            return x**2*y**2 + x**4 + 1

        def dfunc(z):
            x, y = z
            return np.array([2*x*y**2 + 4*x**3, 2*x**2*y])

        for method in ['nelder-mead', 'powell', 'cg', 'bfgs',
                       'newton-cg', 'l-bfgs-b', 'tnc',
                       'cobyla', 'slsqp']:
            if method in ('nelder-mead', 'powell', 'cobyla'):
                jac = None
            else:
                jac = dfunc

            sol1 = optimize.minimize(func, [1, 1], jac=jac, tol=1e-10,
                                     method=method)
            sol2 = optimize.minimize(func, [1, 1], jac=jac, tol=1.0,
                                     method=method)
            assert_(func(sol1.x) < func(sol2.x),
                    "%s: %s vs. %s" % (method, func(sol1.x), func(sol2.x)))

    @pytest.mark.parametrize('method', ['nelder-mead', 'powell', 'cg', 'bfgs', 'newton-cg',
                              'l-bfgs-b', 'tnc', 'cobyla', 'slsqp'])
    def test_no_increase(self, method):
        # Check that the solver doesn't return a value worse than the
        # initial point.

        def func(x):
            return (x - 1)**2

        def bad_grad(x):
            # purposefully invalid gradient function, simulates a case
            # where line searches start failing
            return 2*(x - 1) * (-1) - 2

        x0 = np.array([2.0])
        f0 = func(x0)
        jac = bad_grad
        if method in ['nelder-mead', 'powell', 'cobyla']:
            jac = None
        sol = optimize.minimize(func, x0, jac=jac, method=method,
                                options=dict(maxiter=20))
        assert_equal(func(sol.x), sol.fun)

        if method == 'slsqp':
            pytest.xfail("SLSQP returns slightly worse")
        assert_(func(sol.x) <= f0)

    def test_slsqp_respect_bounds(self):
        # Regression test for gh-3108
        def f(x):
            return sum((x - np.array([1., 2., 3., 4.]))**2)

        def cons(x):
            a = np.array([[-1, -1, -1, -1], [-3, -3, -2, -1]])
            return np.concatenate([np.dot(a, x) + np.array([5, 10]), x])

        x0 = np.array([0.5, 1., 1.5, 2.])
        res = optimize.minimize(f, x0, method='slsqp',
                                constraints={'type': 'ineq', 'fun': cons})
        assert_allclose(res.x, np.array([0., 2, 5, 8])/3, atol=1e-12)

    def test_minimize_automethod(self):
        def f(x):
            return x**2

        def cons(x):
            return x - 2

        x0 = np.array([10.])
        sol_0 = optimize.minimize(f, x0)
        sol_1 = optimize.minimize(f, x0, constraints=[{'type': 'ineq', 'fun': cons}])
        sol_2 = optimize.minimize(f, x0, bounds=[(5, 10)])
        sol_3 = optimize.minimize(f, x0, constraints=[{'type': 'ineq', 'fun': cons}], bounds=[(5, 10)])
        sol_4 = optimize.minimize(f, x0, constraints=[{'type': 'ineq', 'fun': cons}], bounds=[(1, 10)])
        for sol in [sol_0, sol_1, sol_2, sol_3, sol_4]:
            assert_(sol.success)
        assert_allclose(sol_0.x, 0, atol=1e-7)
        assert_allclose(sol_1.x, 2, atol=1e-7)
        assert_allclose(sol_2.x, 5, atol=1e-7)
        assert_allclose(sol_3.x, 5, atol=1e-7)
        assert_allclose(sol_4.x, 2, atol=1e-7)

    def test_minimize_coerce_args_param(self):
        # Regression test for gh-3503
        def Y(x, c):
            return np.sum((x-c)**2)

        def dY_dx(x, c=None):
            return 2*(x-c)

        c = np.array([3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5])
        xinit = np.random.randn(len(c))
        optimize.minimize(Y, xinit, jac=dY_dx, args=(c), method="BFGS")

    def test_initial_step_scaling(self):
        # Check that optimizer initial step is not huge even if the
        # function and gradients are

        scales = [1e-50, 1, 1e50]
        methods = ['CG', 'BFGS', 'L-BFGS-B', 'Newton-CG']

        def f(x):
            if first_step_size[0] is None and x[0] != x0[0]:
                first_step_size[0] = abs(x[0] - x0[0])
            if abs(x).max() > 1e4:
                raise AssertionError("Optimization stepped far away!")
            return scale*(x[0] - 1)**2

        def g(x):
            return np.array([scale*(x[0] - 1)])

        for scale, method in itertools.product(scales, methods):
            if method in ('CG', 'BFGS'):
                options = dict(gtol=scale*1e-8)
            else:
                options = dict()

            if scale < 1e-10 and method in ('L-BFGS-B', 'Newton-CG'):
                # XXX: return initial point if they see small gradient
                continue

            x0 = [-1.0]
            first_step_size = [None]
            res = optimize.minimize(f, x0, jac=g, method=method,
                                    options=options)

            err_msg = "{0} {1}: {2}: {3}".format(method, scale, first_step_size,
                                                 res)

            assert_(res.success, err_msg)
            assert_allclose(res.x, [1.0], err_msg=err_msg)
            assert_(res.nit <= 3, err_msg)

            if scale > 1e-10:
                if method in ('CG', 'BFGS'):
                    assert_allclose(first_step_size[0], 1.01, err_msg=err_msg)
                else:
                    # Newton-CG and L-BFGS-B use different logic for the first step,
                    # but are both scaling invariant with step sizes ~ 1
                    assert_(first_step_size[0] > 0.5 and first_step_size[0] < 3,
                            err_msg)
            else:
                # step size has upper bound of ||grad||, so line
                # search makes many small steps
                pass


class TestLBFGSBBounds(object):
    def setup_method(self):
        self.bounds = ((1, None), (None, None))
        self.solution = (1, 0)

    def fun(self, x, p=2.0):
        return 1.0 / p * (x[0]**p + x[1]**p)

    def jac(self, x, p=2.0):
        return x**(p - 1)

    def fj(self, x, p=2.0):
        return self.fun(x, p), self.jac(x, p)

    def test_l_bfgs_b_bounds(self):
        x, f, d = optimize.fmin_l_bfgs_b(self.fun, [0, -1],
                                         fprime=self.jac,
                                         bounds=self.bounds)
        assert_(d['warnflag'] == 0, d['task'])
        assert_allclose(x, self.solution, atol=1e-6)

    def test_l_bfgs_b_funjac(self):
        # L-BFGS-B with fun and jac combined and extra arguments
        x, f, d = optimize.fmin_l_bfgs_b(self.fj, [0, -1], args=(2.0, ),
                                         bounds=self.bounds)
        assert_(d['warnflag'] == 0, d['task'])
        assert_allclose(x, self.solution, atol=1e-6)

    def test_minimize_l_bfgs_b_bounds(self):
        # Minimize with method='L-BFGS-B' with bounds
        res = optimize.minimize(self.fun, [0, -1], method='L-BFGS-B',
                                jac=self.jac, bounds=self.bounds)
        assert_(res['success'], res['message'])
        assert_allclose(res.x, self.solution, atol=1e-6)


class TestOptimizeScalar(object):
    def setup_method(self):
        self.solution = 1.5

    def fun(self, x, a=1.5):
        """Objective function"""
        return (x - a)**2 - 0.8

    def test_brent(self):
        x = optimize.brent(self.fun)
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.brent(self.fun, brack=(-3, -2))
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.brent(self.fun, full_output=True)
        assert_allclose(x[0], self.solution, atol=1e-6)

        x = optimize.brent(self.fun, brack=(-15, -1, 15))
        assert_allclose(x, self.solution, atol=1e-6)

    def test_golden(self):
        x = optimize.golden(self.fun)
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.golden(self.fun, brack=(-3, -2))
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.golden(self.fun, full_output=True)
        assert_allclose(x[0], self.solution, atol=1e-6)

        x = optimize.golden(self.fun, brack=(-15, -1, 15))
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.golden(self.fun, tol=0)
        assert_allclose(x, self.solution)

        maxiter_test_cases = [0, 1, 5]
        for maxiter in maxiter_test_cases:
            x0 = optimize.golden(self.fun, maxiter=0, full_output=True)
            x = optimize.golden(self.fun, maxiter=maxiter, full_output=True)
            nfev0, nfev = x0[2], x[2]
            assert_equal(nfev - nfev0, maxiter)

    def test_fminbound(self):
        x = optimize.fminbound(self.fun, 0, 1)
        assert_allclose(x, 1, atol=1e-4)

        x = optimize.fminbound(self.fun, 1, 5)
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.fminbound(self.fun, np.array([1]), np.array([5]))
        assert_allclose(x, self.solution, atol=1e-6)
        assert_raises(ValueError, optimize.fminbound, self.fun, 5, 1)

    def test_fminbound_scalar(self):
        try:
            optimize.fminbound(self.fun, np.zeros((1, 2)), 1)
            self.fail("exception not raised")
        except ValueError as e:
            assert_('must be scalar' in str(e))

        x = optimize.fminbound(self.fun, 1, np.array(5))
        assert_allclose(x, self.solution, atol=1e-6)

    def test_minimize_scalar(self):
        # combine all tests above for the minimize_scalar wrapper
        x = optimize.minimize_scalar(self.fun).x
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.minimize_scalar(self.fun, method='Brent')
        assert_(x.success)

        x = optimize.minimize_scalar(self.fun, method='Brent',
                                     options=dict(maxiter=3))
        assert_(not x.success)

        x = optimize.minimize_scalar(self.fun, bracket=(-3, -2),
                                    args=(1.5, ), method='Brent').x
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.minimize_scalar(self.fun, method='Brent',
                                    args=(1.5,)).x
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.minimize_scalar(self.fun, bracket=(-15, -1, 15),
                                    args=(1.5, ), method='Brent').x
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.minimize_scalar(self.fun, bracket=(-3, -2),
                                     args=(1.5, ), method='golden').x
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.minimize_scalar(self.fun, method='golden',
                                     args=(1.5,)).x
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.minimize_scalar(self.fun, bracket=(-15, -1, 15),
                                     args=(1.5, ), method='golden').x
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.minimize_scalar(self.fun, bounds=(0, 1), args=(1.5,),
                                     method='Bounded').x
        assert_allclose(x, 1, atol=1e-4)

        x = optimize.minimize_scalar(self.fun, bounds=(1, 5), args=(1.5, ),
                                    method='bounded').x
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.minimize_scalar(self.fun, bounds=(np.array([1]),
                                                      np.array([5])),
                                    args=(np.array([1.5]), ),
                                    method='bounded').x
        assert_allclose(x, self.solution, atol=1e-6)

        assert_raises(ValueError, optimize.minimize_scalar, self.fun,
                      bounds=(5, 1), method='bounded', args=(1.5, ))

        assert_raises(ValueError, optimize.minimize_scalar, self.fun,
                      bounds=(np.zeros(2), 1), method='bounded', args=(1.5, ))

        x = optimize.minimize_scalar(self.fun, bounds=(1, np.array(5)),
                                     method='bounded').x
        assert_allclose(x, self.solution, atol=1e-6)

    def test_minimize_scalar_custom(self):
        # This function comes from the documentation example.
        def custmin(fun, bracket, args=(), maxfev=None, stepsize=0.1,
                maxiter=100, callback=None, **options):
            bestx = (bracket[1] + bracket[0]) / 2.0
            besty = fun(bestx)
            funcalls = 1
            niter = 0
            improved = True
            stop = False

            while improved and not stop and niter < maxiter:
                improved = False
                niter += 1
                for testx in [bestx - stepsize, bestx + stepsize]:
                    testy = fun(testx, *args)
                    funcalls += 1
                    if testy < besty:
                        besty = testy
                        bestx = testx
                        improved = True
                if callback is not None:
                    callback(bestx)
                if maxfev is not None and funcalls >= maxfev:
                    stop = True
                    break

            return optimize.OptimizeResult(fun=besty, x=bestx, nit=niter,
                                           nfev=funcalls, success=(niter > 1))

        res = optimize.minimize_scalar(self.fun, bracket=(0, 4), method=custmin,
                                       options=dict(stepsize=0.05))
        assert_allclose(res.x, self.solution, atol=1e-6)

    def test_minimize_scalar_coerce_args_param(self):
        # Regression test for gh-3503
        optimize.minimize_scalar(self.fun, args=1.5)


def test_brent_negative_tolerance():
    assert_raises(ValueError, optimize.brent, np.cos, tol=-.01)


class TestNewtonCg(object):
    def test_rosenbrock(self):
        x0 = np.array([-1.2, 1.0])
        sol = optimize.minimize(optimize.rosen, x0,
                                jac=optimize.rosen_der,
                                hess=optimize.rosen_hess,
                                tol=1e-5,
                                method='Newton-CG')
        assert_(sol.success, sol.message)
        assert_allclose(sol.x, np.array([1, 1]), rtol=1e-4)

    def test_himmelblau(self):
        x0 = np.array(himmelblau_x0)
        sol = optimize.minimize(himmelblau,
                                x0,
                                jac=himmelblau_grad,
                                hess=himmelblau_hess,
                                method='Newton-CG',
                                tol=1e-6)
        assert_(sol.success, sol.message)
        assert_allclose(sol.x, himmelblau_xopt, rtol=1e-4)
        assert_allclose(sol.fun, himmelblau_min, atol=1e-4)


class TestRosen(object):

    def test_hess(self):
        # Compare rosen_hess(x) times p with rosen_hess_prod(x,p).  See gh-1775
        x = np.array([3, 4, 5])
        p = np.array([2, 2, 2])
        hp = optimize.rosen_hess_prod(x, p)
        dothp = np.dot(optimize.rosen_hess(x), p)
        assert_equal(hp, dothp)


def himmelblau(p):
    """
    R^2 -> R^1 test function for optimization.  The function has four local
    minima where himmelblau(xopt) == 0.
    """
    x, y = p
    a = x*x + y - 11
    b = x + y*y - 7
    return a*a + b*b


def himmelblau_grad(p):
    x, y = p
    return np.array([4*x**3 + 4*x*y - 42*x + 2*y**2 - 14,
                     2*x**2 + 4*x*y + 4*y**3 - 26*y - 22])


def himmelblau_hess(p):
    x, y = p
    return np.array([[12*x**2 + 4*y - 42, 4*x + 4*y],
                     [4*x + 4*y, 4*x + 12*y**2 - 26]])


himmelblau_x0 = [-0.27, -0.9]
himmelblau_xopt = [3, 2]
himmelblau_min = 0.0


def test_minimize_multiple_constraints():
    # Regression test for gh-4240.
    def func(x):
        return np.array([25 - 0.2 * x[0] - 0.4 * x[1] - 0.33 * x[2]])

    def func1(x):
        return np.array([x[1]])

    def func2(x):
        return np.array([x[2]])

    cons = ({'type': 'ineq', 'fun': func},
            {'type': 'ineq', 'fun': func1},
            {'type': 'ineq', 'fun': func2})

    f = lambda x: -1 * (x[0] + x[1] + x[2])

    res = optimize.minimize(f, [0, 0, 0], method='SLSQP', constraints=cons)
    assert_allclose(res.x, [125, 0, 0], atol=1e-10)


class TestOptimizeResultAttributes(object):
    # Test that all minimizers return an OptimizeResult containing
    # all the OptimizeResult attributes
    def setup_method(self):
        self.x0 = [5, 5]
        self.func = optimize.rosen
        self.jac = optimize.rosen_der
        self.hess = optimize.rosen_hess
        self.hessp = optimize.rosen_hess_prod
        self.bounds = [(0., 10.), (0., 10.)]

    def test_attributes_present(self):
        methods = ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'Newton-CG',
                   'L-BFGS-B', 'TNC', 'COBYLA', 'SLSQP', 'dogleg',
                   'trust-ncg']
        attributes = ['nit', 'nfev', 'x', 'success', 'status', 'fun',
                      'message']
        skip = {'COBYLA': ['nit']}
        for method in methods:
            with suppress_warnings() as sup:
                sup.filter(RuntimeWarning,
                           "Method .+ does not use (gradient|Hessian.*) information")
                res = optimize.minimize(self.func, self.x0, method=method,
                                        jac=self.jac, hess=self.hess,
                                        hessp=self.hessp)
            for attribute in attributes:
                if method in skip and attribute in skip[method]:
                    continue

                assert_(hasattr(res, attribute))
                assert_(attribute in dir(res))


class TestBrute:
    # Test the "brute force" method
    def setup_method(self):
        self.params = (2, 3, 7, 8, 9, 10, 44, -1, 2, 26, 1, -2, 0.5)
        self.rranges = (slice(-4, 4, 0.25), slice(-4, 4, 0.25))
        self.solution = np.array([-1.05665192, 1.80834843])

    def f1(self, z, *params):
        x, y = z
        a, b, c, d, e, f, g, h, i, j, k, l, scale = params
        return (a * x**2 + b * x * y + c * y**2 + d*x + e*y + f)

    def f2(self, z, *params):
        x, y = z
        a, b, c, d, e, f, g, h, i, j, k, l, scale = params
        return (-g*np.exp(-((x-h)**2 + (y-i)**2) / scale))

    def f3(self, z, *params):
        x, y = z
        a, b, c, d, e, f, g, h, i, j, k, l, scale = params
        return (-j*np.exp(-((x-k)**2 + (y-l)**2) / scale))

    def func(self, z, *params):
        return self.f1(z, *params) + self.f2(z, *params) + self.f3(z, *params)

    def test_brute(self):
        # test fmin
        resbrute = optimize.brute(self.func, self.rranges, args=self.params,
                                  full_output=True, finish=optimize.fmin)
        assert_allclose(resbrute[0], self.solution, atol=1e-3)
        assert_allclose(resbrute[1], self.func(self.solution, *self.params),
                        atol=1e-3)

        # test minimize
        resbrute = optimize.brute(self.func, self.rranges, args=self.params,
                                  full_output=True,
                                  finish=optimize.minimize)
        assert_allclose(resbrute[0], self.solution, atol=1e-3)
        assert_allclose(resbrute[1], self.func(self.solution, *self.params),
                        atol=1e-3)

    def test_1D(self):
        # test that for a 1D problem the test function is passed an array,
        # not a scalar.
        def f(x):
            assert_(len(x.shape) == 1)
            assert_(x.shape[0] == 1)
            return x ** 2

        optimize.brute(f, [(-1, 1)], Ns=3, finish=None)


class TestIterationLimits(object):
    # Tests that optimisation does not give up before trying requested
    # number of iterations or evaluations. And that it does not succeed
    # by exceeding the limits.
    def setup_method(self):
        self.funcalls = 0

    def slow_func(self, v):
        self.funcalls += 1
        r,t = np.sqrt(v[0]**2+v[1]**2), np.arctan2(v[0],v[1])
        return np.sin(r*20 + t)+r*0.5

    def test_neldermead_limit(self):
        self.check_limits("Nelder-Mead", 200)

    def test_powell_limit(self):
        self.check_limits("powell", 1000)

    def check_limits(self, method, default_iters):
        for start_v in [[0.1,0.1], [1,1], [2,2]]:
            for mfev in [50, 500, 5000]:
                self.funcalls = 0
                res = optimize.minimize(self.slow_func, start_v,
                      method=method, options={"maxfev":mfev})
                assert_(self.funcalls == res["nfev"])
                if res["success"]:
                    assert_(res["nfev"] <= mfev)
                else:
                    assert_(res["nfev"] >= mfev)
            for mit in [50, 500,5000]:
                res = optimize.minimize(self.slow_func, start_v,
                      method=method, options={"maxiter":mit})
                if res["success"]:
                    assert_(res["nit"] <= mit)
                else:
                    assert_(res["nit"] >= mit)
            for mfev,mit in [[50,50], [5000,5000],[5000,np.inf]]:
                self.funcalls = 0
                res = optimize.minimize(self.slow_func, start_v,
                      method=method, options={"maxiter":mit, "maxfev":mfev})
                assert_(self.funcalls == res["nfev"])
                if res["success"]:
                    assert_(res["nfev"] < mfev and res["nit"] <= mit)
                else:
                    assert_(res["nfev"] >= mfev or res["nit"] >= mit)
            for mfev,mit in [[np.inf,None], [None,np.inf]]:
                self.funcalls = 0
                res = optimize.minimize(self.slow_func, start_v,
                      method=method, options={"maxiter":mit, "maxfev":mfev})
                assert_(self.funcalls == res["nfev"])
                if res["success"]:
                    if mfev is None:
                        assert_(res["nfev"] < default_iters*2)
                    else:
                        assert_(res["nit"] <= default_iters*2)
                else:
                    assert_(res["nfev"] >= default_iters*2 or
                        res["nit"] >= default_iters*2)


class ExpSin(Function):
    def __init__(self):
        super(ExpSin, self).__init__()
        self.func_calls = 0
        self.grad_calls = 0
        self.hess_calls = 0

    def func(self, x):
        self.func_calls += 1
        return np.exp(-x[0]) + np.sin(x[1])

    def grad(self, x):
        self.grad_calls += 1
        return np.array([-np.exp(-x[0]), np.cos(x[1])])

    def hess(self, x):
        self.hess_calls += 1
        return np.exp(-x[0]) - np.sin(x[1])


def expsin_f(x):
    return np.exp(-x[0]) + np.sin(x[1])


def expsin_g(x):
    return np.array([-np.exp(-x[0]), np.cos(x[1])])


def expsin_h(x):
    return np.exp(-x[0]) - np.sin(x[1])


class TestFunction(object):
    def setup_method(self):
        pass

    def test_f(self):
        x = [0.5, 0.5]
        es = ExpSin()
        f = es(x)

        # evaluated Function is correct compared to
        # function
        assert_equal(f, expsin_f(x))
        assert_almost_equal(f, 1.0859561983168364)

        func = Function(func=expsin_f, grad=expsin_g)
        # both func.__call__ and func.func work.
        assert_equal(func(x), f)
        assert_equal(func.f_calls, 1)
        assert_equal(func.func(x), f)
        assert_equal(func.f_calls, 1)
        assert_equal(func.g_calls, 0)
        assert_equal(func.h_calls, 0)

        # func_and_grad works, with correct number of function
        # calls
        f_fg, g_fg = func.func_and_grad(x)
        assert_equal(f_fg, f)
        assert_equal(func.f_calls, 1)
        assert_equal(func.g_calls, 1)

        # test that function calls are correct for func_and_grad
        # (numerical diff)
        es = ExpSin()
        func = Function(func=es.func)
        f, g = func.func_and_grad(x)
        assert_equal(func.f_calls, es.func_calls)

    def test_g(self):
        x = [0.5, 0.5]
        es = ExpSin()

        g_ana = es.grad(x)
        # evaluated g is the same analytically and numerically
        # function
        func = Function(es.func, fd_method='2-point')
        assert_allclose(func.grad(x), g_ana)

        func = Function(es.func, fd_method='3-point')
        assert_allclose(func.grad(x), g_ana)

        func = Function(es.func, fd_method='abs')
        assert_allclose(func.grad(x), g_ana)

        # number of function/grad calls is correct
        es = ExpSin()
        g = es.grad(x)
        assert_equal(es.func_calls, 0)

        func = Function(es.func, fd_method='abs')
        es.func_calls = 0
        func.grad(x)
        assert_equal(func.f_calls, es.func_calls)
        assert_equal(func.g_calls, 1)

        func = Function(es.func, fd_method='2-point')
        es.func_calls = 0
        func.grad(x)
        assert_equal(func.f_calls, es.func_calls)
        assert_equal(func.g_calls, 1)

        func = Function(es.func, fd_method='3-point')
        es.func_calls = 0
        func.grad(x)
        assert_equal(func.f_calls, es.func_calls)
        assert_equal(func.g_calls, 1)

    def test_grad_logit(self):
        # Verify grad of logistic function
        x0 = np.array([1.5])

        func = Function(func=logit, grad=der_logit)
        assert_allclose(func.grad(x0), der_logit(x0))

        func = Function(func=logit, fd_method='3-point')
        assert_allclose(func.grad(x0), der_logit(x0))

        func = Function(func=logit, fd_method='2-point')
        assert_allclose(func.grad(x0), der_logit(x0))

        func = Function(func=logit, fd_method='abs')
        assert_allclose(func.grad(x0), der_logit(x0))

        # hand check absolute step finite difference
        func = Function(func=logit, fd_method='abs', step=1e-1)
        f0 = func(x0)
        f1 = func(x0 + 0.1)
        der = (f1 - f0) / 0.1
        assert_allclose(func.grad(x0), der)

    def test_fgh_berger(self):
        berger = CheckOptimize()
        berger.setup_method()

        function = Function(func=berger.func, grad=berger.grad,
                            hess=berger.hess)

        solution = np.array([0., -0.524869316, 0.487525860])

        assert_allclose(function.func(solution), berger.func(solution))
        assert_allclose(function.grad(solution), berger.grad(solution))
        assert_allclose(function.hess(solution), berger.hess(solution))

        # check numdiff grad on Berger
        function = Function(func=berger.func, fd_method='3-point')
        assert_almost_equal(function.grad(solution), berger.grad(solution))

        function = Function(func=berger.func, fd_method='2-point')
        assert_almost_equal(function.grad(solution), berger.grad(solution))

        function = Function(func=berger.func, fd_method='abs')
        assert_almost_equal(function.grad(solution), berger.grad(solution))

        # check Hessian on Berger
        # numdiff from grad
        function = Function(func=berger.func, grad=berger.grad,
                            fd_method='3-point')
        assert_allclose(function.hess(solution), berger.hess(solution),
                        atol=1e-8)

        # TODO THIS DOESN'T work with 'abs', and is not close enough with
        # 3-point
        # numdiff from func
        function = Function(func=berger.func, fd_method='3-point')
        assert_allclose(function.hess(solution+0.01),
                            berger.hess(solution+0.01),
                            atol=1e-5)

    def test_args_kwargs(self):
        # we want args and kwargs to be passed to the function we're trying
        # to calculate.

        # define new logit with arg and kwarg. Set the default to 2, so
        # we know if the kwarg is set by the caller.
        def logit(x, num, den=2):
            return num / (den + np.exp(-x))

        x0 = np.array([1.5])

        # make a Function with the arg and kwarg
        func = Function(func=logit, args=(1, ),kwargs={'den': 1},
                        fd_method='3-point')

        # check that the func and grad eval are correct.
        assert_allclose(func.func(x0), logit(x0, 1, den=1))
        assert_allclose(func.grad(x0), der_logit(x0))
        f, g = func.func_and_grad(x0)
        assert_allclose(f, logit(x0, 1, den=1))
        assert_allclose(g, der_logit(x0))


class Callback():
    def __init__(self, limit_it=None, error_it=None):
        self.counts = 0
        self.x = []
        # use limit_it to see if we should stop iteration early.
        self.limit_it = limit_it
        self.error_it = error_it

    def __call__(self, res):
        # res should be an intermediate OptimizeResult
        assert_(isinstance(res, optimize.OptimizeResult))
        self.counts += 1
        self.x.append(res.x)
        if self.counts == self.limit_it:
            raise StopIteration
        if self.error_it == self.counts:
            raise RuntimeError


@pytest.mark.parametrize("opt,args,kwargs", [
    (Optimizer, (), {}),
    (NelderMead, (), {'x0': [1.5]}),
    (BFGS, (), {'x0': [1.5], 'gtol':1e-4}),
    (LBFGSB, (), {'x0': [1.5]}),
])
class Test_Optimizer(object):
    # Optimizer base class functionality

    def setup_method(self):
        pass

    def test_initalisation(self, opt, args, kwargs):
        # optimizer object must fail __init__ if not passed Function.
        with pytest.raises(ValueError):
            opt(None, *args, **kwargs)

    def test_func(self, opt, args, kwargs):
        # test that the optimizer itself can calculate a func, gradient
        # and with the correct number of function calls.
        func = Function(func=logit, fd_method='3-point')
        x0 = np.array([1.5])
        optimizer = opt(func, *args, **kwargs)

        assert_equal(optimizer.func(x0), logit((x0)))
        assert_equal(optimizer.nfev, 1)
        assert_allclose(optimizer.grad(x0), der_logit(x0))
        # a grad call with 3 point central difference should require 3 func
        # calls
        assert_equal(optimizer.nfev, 4)
        assert_equal(optimizer.njev, 1)
        _ = optimizer.func_and_grad(x0)
        assert_equal(optimizer.nfev, 7)
        assert_equal(optimizer.njev, 2)

        func = Function(func=logit, grad=der_logit)
        optimizer = opt(func, *args, **kwargs)
        _ = optimizer.grad(x0)
        assert_equal(optimizer.nfev, 0)
        assert_equal(optimizer.njev, 1)

    def test_attributes(self, opt, args, kwargs):
        # this is also a loose test to see if the super class init has been
        # called, as all these attributes should be present.
        optimizer = opt(Function(logit), *args, **kwargs)

        required_attr = ['x', 'fun', 'status', 'message', 'warn_flag', 'nit',
                         'options', '_hyper']

        assert_(np.all([hasattr(optimizer, attr) for attr in required_attr]))

    def test_smoke(self, opt, args, kwargs):
        # smoke test for Optimizers with rosen to check if thing runs
        if opt is Optimizer:
            return

        func = Function(optimize.rosen)
        kwargs['x0'] = np.array([1.1, 1.1])

        # TODO: tighten gtol for BFGS, investigate whether warn_flag can be
        # amended to allow success if line search fails to find better solution
        callee = Callback()
        optimizer = opt(func, *args, **kwargs)
        res = optimizer.solve(maxiter=1000, callback=callee)
        assert_equal(res.nit, callee.counts)
        assert_equal(optimizer.nit, callee.counts)
        assert_(res.nit <= 1000)
        # all the optimizers should pass on Rosen, especially when we're so
        # close to the solution.
        assert_(optimizer.converged())
        assert_(res.success)
        assert_(optimizer.warn_flag == 0)
        assert_equal(optimizer.fun, res.fun)
        assert_equal(res.x, callee.x[-1])
        assert_equal(optimizer.x, callee.x[-1])

        # examine if the problem dimensionality is correct
        assert_equal(optimizer.N, 2)

        # solve shouldn't run if maxiter is 0.
        callee = Callback()
        optimizer = opt(func, *args, **kwargs)
        res = optimizer.solve(maxiter=0, callback=callee)
        assert_equal(callee.counts, 0)
        assert_equal(optimizer.nit, 0)
        assert_equal(res.nit, optimizer.nit)

        # solve shouldn't run if maxfun is 0.
        # NOTE optimizer.nfev will be 1, because as a minimum the optimizer
        # has to have function evaluated for the result.
        callee = Callback()
        optimizer = opt(func, *args, **kwargs)
        res = optimizer.solve(maxfun=0, callback=callee)
        assert_equal(callee.counts, 0)
        assert_equal(optimizer.nfev, 1)
        assert_equal(optimizer.fun, func(res.x))

        # test that stepwise iteration is possible. This will check that
        # __next__ is overridden otherwise there will be a NotImplemented Error
        optimizer = opt(func, *args, **kwargs)
        x, f = next(optimizer)
        assert_equal(optimizer.nit, 1)

        # tests that the result from the iterator is a copy. This is because
        # users may not expect arrays to be passed by reference.
        assert_(id(x) is not id(optimizer.x))

        # test __call__ method
        callee = Callback()
        optimizer = opt(func, *args, **kwargs)
        res = optimizer(5, callback=callee)
        assert_equal(callee.counts, 5)
        assert_equal(optimizer.nit, 5)
        # we only asked for a small amount of iterations, so we should hit
        # limit
        assert_(optimizer.warn_flag != 0)

        # test that iteration halts if callback raises StopIteration
        callee = Callback(limit_it=2)
        optimizer = opt(func, *args, **kwargs)
        res = optimizer.solve(callback=callee)
        assert_equal(callee.counts, 2)
        assert_equal(optimizer.nit, 2)

        # test that iteration halts if function raises StopIteration
        # with only 3 function evaluations it shouldn't have converged
        # (unless it's a super solver).
        def fwrapper(function):
            calls = [0]

            def wrapped_fun(x):
                val = function(x)
                calls[0] += 1
                if calls[0] == 3:
                    raise StopIteration
                return val

            return wrapped_fun

        func = Function(fwrapper(optimize.rosen))
        optimizer = opt(func, *args, **kwargs)
        res = optimizer.solve()
        assert_(not res.success)

        # test that if callback raises other errors they are still raised.
        optimizer = opt(func, *args, **kwargs)
        with pytest.raises(RuntimeError):
            callee = Callback(error_it=2)
            res = optimizer.solve(callback=callee)
        assert_equal(optimizer.nit, 2)

        # test that the context manager works
        with opt(func, *args, **kwargs) as optimizer:
            # more detailed _finish_up behaviour should be tested for each
            # Optimizer
            res = optimizer.solve()
