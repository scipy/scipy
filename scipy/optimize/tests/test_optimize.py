"""
Unit tests for optimization routines from optimize.py

Authors:
   Ed Schofield, Nov 2005
   Andrew Straw, April 2008

To run it in its simplest form::
  nosetests test_optimize.py

"""
from __future__ import division, print_function, absolute_import

from numpy.testing import assert_raises, assert_allclose, \
        assert_equal, assert_, TestCase, run_module_suite, dec

from scipy import optimize
import numpy as np


class TestOptimize(object):
    """ Test case for a simple constrained entropy maximization problem
    (the machine translation example of Berger et al in
    Computational Linguistics, vol 22, num 1, pp 39--72, 1996.)
    """
    def setUp(self):
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

    def test_cg(self, use_wrapper=False):
        """ conjugate gradient optimization routine """
        if use_wrapper:
            opts = {'maxiter': self.maxiter, 'disp': False,
                    'return_all': False}
            res = optimize.minimize(self.func, self.startparams, args=(),
                                    method='CG', jac=self.grad,
                                    options=opts)

            params, fopt, func_calls, grad_calls, warnflag = \
                res['x'], res['fun'], res['nfev'], res['njev'], res['status']
        else:
            retval = optimize.fmin_cg(self.func, self.startparams, self.grad, (),
                                      maxiter=self.maxiter,
                                      full_output=True, disp=False, retall=False)

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

    def test_bfgs(self, use_wrapper=False):
        """ Broyden-Fletcher-Goldfarb-Shanno optimization routine """
        if use_wrapper:
            opts = {'maxiter': self.maxiter, 'disp': False,
                    'return_all': False}
            res = optimize.minimize(self.func, self.startparams,
                                    jac=self.grad, method='BFGS', args=(),
                                    options=opts)

            params, fopt, gopt, Hopt, func_calls, grad_calls, warnflag = \
                    res['x'], res['fun'], res['jac'], res['hess_inv'], \
                    res['nfev'], res['njev'], res['status']
        else:
            retval = optimize.fmin_bfgs(self.func, self.startparams, self.grad,
                                        args=(), maxiter=self.maxiter,
                                        full_output=True, disp=False, retall=False)

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

    def test_bfgs_nan(self):
        """Test corner case where nan is fed to optimizer.  See #1542."""
        func = lambda x: x
        fprime = lambda x: np.ones_like(x)
        x0 = [np.nan]
        olderr = np.seterr(over='ignore')
        try:
            x = optimize.fmin_bfgs(func, x0, fprime, disp=False)
            assert_(np.isnan(func(x)))
        finally:
            np.seterr(**olderr)

    def test_bfgs_numerical_jacobian(self):
        """ BFGS with numerical jacobian and a vector epsilon parameter """
        # define the epsilon parameter using a random vector
        epsilon = np.sqrt(np.finfo(float).eps) * np.random.rand(len(self.solution))

        params = optimize.fmin_bfgs(self.func, self.startparams,
                                    epsilon=epsilon, args=(),
                                    maxiter=self.maxiter, disp=False)

        assert_allclose(self.func(params), self.func(self.solution),
                        atol=1e-6)

    def test_bfgs_infinite(self, use_wrapper=False):
        """Test corner case where -Inf is the minimum.  See #1494."""
        func = lambda x: -np.e**-x
        fprime = lambda x: -func(x)
        x0 = [0]
        olderr = np.seterr(over='ignore')
        try:
            if use_wrapper:
                opts = {'disp': False}
                x = optimize.minimize(func, x0, jac=fprime, method='BFGS',
                                      args=(), options=opts)['x']
            else:
                x = optimize.fmin_bfgs(func, x0, fprime, disp=False)
            assert_(not np.isfinite(func(x)))
        finally:
            np.seterr(**olderr)

    def test_powell(self, use_wrapper=False):
        """ Powell (direction set) optimization routine
        """
        if use_wrapper:
            opts = {'maxiter': self.maxiter, 'disp': False,
                    'return_all': False}
            res = optimize.minimize(self.func, self.startparams, args=(),
                                    method='Powell', options=opts)
            params, fopt, direc, numiter, func_calls, warnflag = \
                    res['x'], res['fun'], res['direc'], res['nit'], \
                    res['nfev'], res['status']
        else:
            retval = optimize.fmin_powell(self.func, self.startparams,
                                        args=(), maxiter=self.maxiter,
                                        full_output=True, disp=False, retall=False)

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

    def test_neldermead(self, use_wrapper=False):
        """ Nelder-Mead simplex algorithm
        """
        if use_wrapper:
            opts = {'maxiter': self.maxiter, 'disp': False,
                    'return_all': False}
            res = optimize.minimize(self.func, self.startparams, args=(),
                                    method='Nelder-mead', options=opts)
            params, fopt, numiter, func_calls, warnflag = \
                    res['x'], res['fun'], res['nit'], res['nfev'], \
                    res['status']
        else:
            retval = optimize.fmin(self.func, self.startparams,
                                        args=(), maxiter=self.maxiter,
                                        full_output=True, disp=False, retall=False)

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

    def test_ncg(self, use_wrapper=False):
        """ line-search Newton conjugate gradient optimization routine
        """
        if use_wrapper:
            opts = {'maxiter': self.maxiter, 'disp': False,
                    'return_all': False}
            retval = optimize.minimize(self.func, self.startparams,
                                       method='Newton-CG', jac=self.grad,
                                       args=(), options=opts)['x']
        else:
            retval = optimize.fmin_ncg(self.func, self.startparams, self.grad,
                                       args=(), maxiter=self.maxiter,
                                       full_output=False, disp=False,
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

    def test_ncg_hess(self, use_wrapper=False):
        """ Newton conjugate gradient with Hessian """
        if use_wrapper:
            opts = {'maxiter': self.maxiter, 'disp': False,
                    'return_all': False}
            retval = optimize.minimize(self.func, self.startparams,
                                       method='Newton-CG', jac=self.grad,
                                       hess=self.hess,
                                       args=(), options=opts)['x']
        else:
            retval = optimize.fmin_ncg(self.func, self.startparams, self.grad,
                                       fhess=self.hess,
                                       args=(), maxiter=self.maxiter,
                                       full_output=False, disp=False,
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

    def test_ncg_hessp(self, use_wrapper=False):
        """ Newton conjugate gradient with Hessian times a vector p """
        if use_wrapper:
            opts = {'maxiter': self.maxiter, 'disp': False,
                    'return_all': False}
            retval = optimize.minimize(self.func, self.startparams,
                                       method='Newton-CG', jac=self.grad,
                                       hessp=self.hessp,
                                       args=(), options=opts)['x']
        else:
            retval = optimize.fmin_ncg(self.func, self.startparams, self.grad,
                                       fhess_p=self.hessp,
                                       args=(), maxiter=self.maxiter,
                                       full_output=False, disp=False,
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

    def test_l_bfgs_b(self):
        """ limited-memory bound-constrained BFGS algorithm
        """
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
        """ L-BFGS-B with numerical jacobian """
        retval = optimize.fmin_l_bfgs_b(self.func, self.startparams,
                                        approx_grad=True,
                                        maxiter=self.maxiter)

        (params, fopt, d) = retval

        assert_allclose(self.func(params), self.func(self.solution),
                        atol=1e-6)

    def test_l_bfgs_b_funjac(self):
        """ L-BFGS-B with combined objective function and jacobian """
        def fun(x):
            return self.func(x), self.grad(x)

        retval = optimize.fmin_l_bfgs_b(fun, self.startparams,
                                        maxiter=self.maxiter)

        (params, fopt, d) = retval

        assert_allclose(self.func(params), self.func(self.solution),
                        atol=1e-6)

    def test_minimize_l_bfgs_b(self):
        """ Minimize with L-BFGS-B method """
        opts = {'disp': False, 'maxiter': self.maxiter}
        x = optimize.minimize(self.func, self.startparams,
                              method='L-BFGS-B', jac=self.grad,
                              options=opts)['x']
        assert_allclose(self.func(x), self.func(self.solution),
                        atol=1e-6)

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

    def test_minimize(self):
        """Tests for the minimize wrapper."""
        self.setUp()
        self.test_bfgs(True)
        self.setUp()
        self.test_bfgs_infinite(True)
        self.setUp()
        self.test_cg(True)
        self.setUp()
        self.test_ncg(True)
        self.setUp()
        self.test_ncg_hess(True)
        self.setUp()
        self.test_ncg_hessp(True)
        self.setUp()
        self.test_neldermead(True)
        self.setUp()
        self.test_powell(True)

    def test_minimize_tol_parameter(self):
        # Check that the minimize() tol= argument does something
        def func(z):
            x, y = z
            return x**2*y**2 + x**4 + 1

        def dfunc(z):
            x, y = z
            return np.array([2*x*y**2 + 4*x**3, 2*x**2*y])

        for method in ['nelder-mead', 'powell', 'cg', 'bfgs',
                       'newton-cg', 'anneal', 'l-bfgs-b', 'tnc',
                       'cobyla', 'slsqp']:
            if method in ('nelder-mead', 'powell', 'anneal', 'cobyla'):
                jac = None
            else:
                jac = dfunc
            sol1 = optimize.minimize(func, [1,1], jac=jac, tol=1e-10,
                                     method=method)
            sol2 = optimize.minimize(func, [1,1], jac=jac, tol=1.0,
                                     method=method)
            assert_(func(sol1.x) < func(sol2.x),
                    "%s: %s vs. %s" % (method, func(sol1.x), func(sol2.x)))

    def test_no_increase(self):
        # Check that the solver doesn't return a value worse than the
        # initial point.

        def func(x):
            return (x - 1)**2

        def bad_grad(x):
            # purposefully invalid gradient function, simulates a case
            # where line searches start failing
            return 2*(x - 1) * (-1) - 2

        def check(method):
            x0 = np.array([2.0])
            f0 = func(x0)
            jac = bad_grad
            if method in ['nelder-mead', 'powell', 'anneal', 'cobyla']:
                jac = None
            sol = optimize.minimize(func, x0, jac=jac, method=method,
                                    options=dict(maxiter=20))
            assert_equal(func(sol.x), sol.fun)

            dec.knownfailureif(method == 'slsqp', "SLSQP returns slightly worse")(lambda: None)()
            assert_(func(sol.x) <= f0)

        for method in ['nelder-mead', 'powell', 'cg', 'bfgs',
                       'newton-cg', 'anneal', 'l-bfgs-b', 'tnc',
                       'cobyla', 'slsqp']:
            yield check, method


class TestLBFGSBBounds(TestCase):
    """ Tests for L-BFGS-B with bounds """
    def setUp(self):
        self.bounds = ((1, None), (None, None))
        self.solution = (1, 0)

    def fun(self, x, p=2.0):
        return 1.0 / p * (x[0]**p + x[1]**p)

    def jac(self, x, p=2.0):
        return x**(p - 1)

    def fj(self, x, p=2.0):
        return self.fun(x, p), self.jac(x, p)

    def test_l_bfgs_b_bounds(self):
        """ L-BFGS-B with bounds """
        x, f, d = optimize.fmin_l_bfgs_b(self.fun, [0, -1],
                                         fprime=self.jac,
                                         bounds=self.bounds)
        assert_(d['warnflag'] == 0, d['task'])
        assert_allclose(x, self.solution, atol=1e-6)

    def test_l_bfgs_b_funjac(self):
        """ L-BFGS-B with fun and jac combined and extra arguments """
        x, f, d = optimize.fmin_l_bfgs_b(self.fj, [0, -1], args=(2.0, ),
                                         bounds=self.bounds)
        assert_(d['warnflag'] == 0, d['task'])
        assert_allclose(x, self.solution, atol=1e-6)

    def test_minimize_l_bfgs_b_bounds(self):
        """ Minimize with method='L-BFGS-B' with bounds """
        res = optimize.minimize(self.fun, [0, -1], method='L-BFGS-B',
                                jac=self.jac, bounds=self.bounds)
        assert_(res['success'], res['message'])
        assert_allclose(res.x, self.solution, atol=1e-6)


class TestOptimizeScalar(TestCase):
    """Tests for scalar optimizers"""
    def setUp(self):
        self.solution = 1.5

    def fun(self, x, a=1.5):
        """Objective function"""
        return (x - a)**2 - 0.8

    def test_brent(self):
        """ brent algorithm """
        x = optimize.brent(self.fun)
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.brent(self.fun, brack=(-3, -2))
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.brent(self.fun, full_output=True)
        assert_allclose(x[0], self.solution, atol=1e-6)

        x = optimize.brent(self.fun, brack=(-15, -1, 15))
        assert_allclose(x, self.solution, atol=1e-6)

    def test_golden(self):
        """ golden algorithm """
        x = optimize.golden(self.fun)
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.golden(self.fun, brack=(-3, -2))
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.golden(self.fun, full_output=True)
        assert_allclose(x[0], self.solution, atol=1e-6)

        x = optimize.golden(self.fun, brack=(-15, -1, 15))
        assert_allclose(x, self.solution, atol=1e-6)

    def test_fminbound(self):
        """Test fminbound """
        x = optimize.fminbound(self.fun, 0, 1)
        assert_allclose(x, 1, atol=1e-4)

        x = optimize.fminbound(self.fun, 1, 5)
        assert_allclose(x, self.solution, atol=1e-6)

        x = optimize.fminbound(self.fun, np.array([1]), np.array([5]))
        assert_allclose(x, self.solution, atol=1e-6)
        assert_raises(ValueError, optimize.fminbound, self.fun, 5, 1)

    def test_fminbound_scalar(self):
        assert_raises(ValueError, optimize.fminbound, self.fun,
                      np.zeros(2), 1)

        x = optimize.fminbound(self.fun, 1, np.array(5))
        assert_allclose(x, self.solution, atol=1e-6)

    def test_minimize_scalar(self):
        # combine all tests above for the minimize_scalar wrapper
        x = optimize.minimize_scalar(self.fun).x
        assert_allclose(x, self.solution, atol=1e-6)

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


class TestRosen(TestCase):

    def test_hess(self):
        """Compare rosen_hess(x) times p with rosen_hess_prod(x,p) (ticket #1248)"""
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

if __name__ == "__main__":
    run_module_suite()
