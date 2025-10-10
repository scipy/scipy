"""
Unit test for SLSQP optimization.
"""
from numpy.testing import (assert_, assert_array_almost_equal,
                           assert_allclose, assert_equal)
from pytest import raises as assert_raises
import pytest
import numpy as np
import scipy

from scipy.optimize import (fmin_slsqp, minimize, Bounds, NonlinearConstraint,
                            OptimizeResult, AnalysisError)


class MyCallBack:
    """pass a custom callback function

    This makes sure it's being used.
    """
    def __init__(self):
        self.been_called = False
        self.ncalls = 0

    def __call__(self, x):
        assert not isinstance(x, OptimizeResult)
        self.been_called = True
        self.ncalls += 1

    def callback2(self, intermediate_result):
        assert isinstance(intermediate_result, OptimizeResult)
        self.been_called = True
        self.ncalls += 1

    def callback3(self, intermediate_result):
        assert isinstance(intermediate_result, OptimizeResult)
        raise StopIteration


class TestSLSQP:
    """
    Test SLSQP algorithm using Example 14.4 from Numerical Methods for
    Engineers by Steven Chapra and Raymond Canale.
    This example maximizes the function f(x) = 2*x*y + 2*x - x**2 - 2*y**2,
    which has a maximum at x=2, y=1.
    """
    def setup_method(self):
        self.opts = {'disp': False}

    def fun(self, d, sign=1.0):
        """
        Arguments:
        d     - A list of two elements, where d[0] represents x and d[1] represents y
                 in the following equation.
        sign - A multiplier for f. Since we want to optimize it, and the SciPy
               optimizers can only minimize functions, we need to multiply it by
               -1 to achieve the desired solution
        Returns:
        2*x*y + 2*x - x**2 - 2*y**2

        """
        x = d[0]
        y = d[1]
        return sign*(2*x*y + 2*x - x**2 - 2*y**2)

    def jac(self, d, sign=1.0):
        """
        This is the derivative of fun, returning a NumPy array
        representing df/dx and df/dy.

        """
        x = d[0]
        y = d[1]
        dfdx = sign*(-2*x + 2*y + 2)
        dfdy = sign*(2*x - 4*y)
        return np.array([dfdx, dfdy], float)

    def fun_and_jac(self, d, sign=1.0):
        return self.fun(d, sign), self.jac(d, sign)

    def f_eqcon(self, x, sign=1.0):
        """ Equality constraint """
        return np.array([x[0] - x[1]])

    def fprime_eqcon(self, x, sign=1.0):
        """ Equality constraint, derivative """
        return np.array([[1, -1]])

    def f_eqcon_scalar(self, x, sign=1.0):
        """ Scalar equality constraint """
        return self.f_eqcon(x, sign)[0]

    def fprime_eqcon_scalar(self, x, sign=1.0):
        """ Scalar equality constraint, derivative """
        return self.fprime_eqcon(x, sign)[0].tolist()

    def f_ieqcon(self, x, sign=1.0):
        """ Inequality constraint """
        return np.array([x[0] - x[1] - 1.0])

    def fprime_ieqcon(self, x, sign=1.0):
        """ Inequality constraint, derivative """
        return np.array([[1, -1]])

    def f_ieqcon2(self, x):
        """ Vector inequality constraint """
        return np.asarray(x)

    def fprime_ieqcon2(self, x):
        """ Vector inequality constraint, derivative """
        return np.identity(x.shape[0])

    # minimize
    def test_minimize_unbounded_approximated(self):
        # Minimize, method='SLSQP': unbounded, approximated jacobian.
        jacs = [None, False, '2-point', '3-point']
        for jac in jacs:
            res = minimize(self.fun, [-1.0, 1.0], args=(-1.0, ),
                           jac=jac, method='SLSQP',
                           options=self.opts)
            assert_(res['success'], res['message'])
            assert_allclose(res.x, [2, 1])

    def test_minimize_unbounded_given(self):
        # Minimize, method='SLSQP': unbounded, given Jacobian.
        res = minimize(self.fun, [-1.0, 1.0], args=(-1.0, ),
                       jac=self.jac, method='SLSQP', options=self.opts)
        assert_(res['success'], res['message'])
        assert_allclose(res.x, [2, 1])

    def test_minimize_bounded_approximated(self):
        # Minimize, method='SLSQP': bounded, approximated jacobian.
        jacs = [None, False, '2-point', '3-point']
        for jac in jacs:
            with np.errstate(invalid='ignore'):
                res = minimize(self.fun, [-1.0, 1.0], args=(-1.0, ),
                               jac=jac,
                               bounds=((2.5, None), (None, 0.5)),
                               method='SLSQP', options=self.opts)
            assert_(res['success'], res['message'])
            assert_allclose(res.x, [2.5, 0.5])
            assert_(2.5 <= res.x[0])
            assert_(res.x[1] <= 0.5)

    def test_minimize_unbounded_combined(self):
        # Minimize, method='SLSQP': unbounded, combined function and Jacobian.
        res = minimize(self.fun_and_jac, [-1.0, 1.0], args=(-1.0, ),
                       jac=True, method='SLSQP', options=self.opts)
        assert_(res['success'], res['message'])
        assert_allclose(res.x, [2, 1])

    def test_minimize_equality_approximated(self):
        # Minimize with method='SLSQP': equality constraint, approx. jacobian.
        jacs = [None, False, '2-point', '3-point']
        for jac in jacs:
            res = minimize(self.fun, [-1.0, 1.0], args=(-1.0, ),
                           jac=jac,
                           constraints={'type': 'eq',
                                        'fun': self.f_eqcon,
                                        'args': (-1.0, )},
                           method='SLSQP', options=self.opts)
            assert_(res['success'], res['message'])
            assert_allclose(res.x, [1, 1])

    def test_minimize_equality_given(self):
        # Minimize with method='SLSQP': equality constraint, given Jacobian.
        res = minimize(self.fun, [-1.0, 1.0], jac=self.jac,
                       method='SLSQP', args=(-1.0,),
                       constraints={'type': 'eq', 'fun':self.f_eqcon,
                                    'args': (-1.0, )},
                       options=self.opts)
        assert_(res['success'], res['message'])
        assert_allclose(res.x, [1, 1])

    def test_minimize_equality_given2(self):
        # Minimize with method='SLSQP': equality constraint, given Jacobian
        # for fun and const.
        res = minimize(self.fun, [-1.0, 1.0], method='SLSQP',
                       jac=self.jac, args=(-1.0,),
                       constraints={'type': 'eq',
                                    'fun': self.f_eqcon,
                                    'args': (-1.0, ),
                                    'jac': self.fprime_eqcon},
                       options=self.opts)
        assert_(res['success'], res['message'])
        assert_allclose(res.x, [1, 1])

    def test_minimize_equality_given_cons_scalar(self):
        # Minimize with method='SLSQP': scalar equality constraint, given
        # Jacobian for fun and const.
        res = minimize(self.fun, [-1.0, 1.0], method='SLSQP',
                       jac=self.jac, args=(-1.0,),
                       constraints={'type': 'eq',
                                    'fun': self.f_eqcon_scalar,
                                    'args': (-1.0, ),
                                    'jac': self.fprime_eqcon_scalar},
                       options=self.opts)
        assert_(res['success'], res['message'])
        assert_allclose(res.x, [1, 1])

    def test_minimize_inequality_given(self):
        # Minimize with method='SLSQP': inequality constraint, given Jacobian.
        res = minimize(self.fun, [-1.0, 1.0], method='SLSQP',
                       jac=self.jac, args=(-1.0, ),
                       constraints={'type': 'ineq',
                                    'fun': self.f_ieqcon,
                                    'args': (-1.0, )},
                       options=self.opts)
        assert_(res['success'], res['message'])
        assert_allclose(res.x, [2, 1], atol=1e-3)

    def test_minimize_inequality_given_vector_constraints(self):
        # Minimize with method='SLSQP': vector inequality constraint, given
        # Jacobian.
        res = minimize(self.fun, [-1.0, 1.0], jac=self.jac,
                       method='SLSQP', args=(-1.0,),
                       constraints={'type': 'ineq',
                                    'fun': self.f_ieqcon2,
                                    'jac': self.fprime_ieqcon2},
                       options=self.opts)
        assert_(res['success'], res['message'])
        assert_allclose(res.x, [2, 1])

    def test_minimize_bounded_constraint(self):
        # when the constraint makes the solver go up against a parameter
        # bound make sure that the numerical differentiation of the
        # jacobian doesn't try to exceed that bound using a finite difference.
        # gh11403
        def c(x):
            assert 0 <= x[0] <= 1 and 0 <= x[1] <= 1, x
            return x[0] ** 0.5 + x[1]

        def f(x):
            assert 0 <= x[0] <= 1 and 0 <= x[1] <= 1, x
            return -x[0] ** 2 + x[1] ** 2

        cns = [NonlinearConstraint(c, 0, 1.5)]
        x0 = np.asarray([0.9, 0.5])
        bnd = Bounds([0., 0.], [1.0, 1.0])
        minimize(f, x0, method='SLSQP', bounds=bnd, constraints=cns)

    def test_minimize_bound_equality_given2(self):
        # Minimize with method='SLSQP': bounds, eq. const., given jac. for
        # fun. and const.
        res = minimize(self.fun, [-1.0, 1.0], method='SLSQP',
                       jac=self.jac, args=(-1.0, ),
                       bounds=[(-0.8, 1.), (-1, 0.8)],
                       constraints={'type': 'eq',
                                    'fun': self.f_eqcon,
                                    'args': (-1.0, ),
                                    'jac': self.fprime_eqcon},
                       options=self.opts)
        assert_(res['success'], res['message'])
        assert_allclose(res.x, [0.8, 0.8], atol=1e-3)
        assert_(-0.8 <= res.x[0] <= 1)
        assert_(-1 <= res.x[1] <= 0.8)

    # fmin_slsqp
    def test_unbounded_approximated(self):
        # SLSQP: unbounded, approximated Jacobian.
        res = fmin_slsqp(self.fun, [-1.0, 1.0], args=(-1.0, ),
                         iprint = 0, full_output = 1)
        x, fx, its, imode, smode = res
        assert_(imode == 0, imode)
        assert_array_almost_equal(x, [2, 1])

    def test_unbounded_given(self):
        # SLSQP: unbounded, given Jacobian.
        res = fmin_slsqp(self.fun, [-1.0, 1.0], args=(-1.0, ),
                         fprime = self.jac, iprint = 0,
                         full_output = 1)
        x, fx, its, imode, smode = res
        assert_(imode == 0, imode)
        assert_array_almost_equal(x, [2, 1])

    def test_equality_approximated(self):
        # SLSQP: equality constraint, approximated Jacobian.
        res = fmin_slsqp(self.fun,[-1.0,1.0], args=(-1.0,),
                         eqcons = [self.f_eqcon],
                         iprint = 0, full_output = 1)
        x, fx, its, imode, smode = res
        assert_(imode == 0, imode)
        assert_array_almost_equal(x, [1, 1])

    def test_equality_given(self):
        # SLSQP: equality constraint, given Jacobian.
        res = fmin_slsqp(self.fun, [-1.0, 1.0],
                         fprime=self.jac, args=(-1.0,),
                         eqcons = [self.f_eqcon], iprint = 0,
                         full_output = 1)
        x, fx, its, imode, smode = res
        assert_(imode == 0, imode)
        assert_array_almost_equal(x, [1, 1])

    def test_equality_given2(self):
        # SLSQP: equality constraint, given Jacobian for fun and const.
        res = fmin_slsqp(self.fun, [-1.0, 1.0],
                         fprime=self.jac, args=(-1.0,),
                         f_eqcons = self.f_eqcon,
                         fprime_eqcons = self.fprime_eqcon,
                         iprint = 0,
                         full_output = 1)
        x, fx, its, imode, smode = res
        assert_(imode == 0, imode)
        assert_array_almost_equal(x, [1, 1])

    def test_inequality_given(self):
        # SLSQP: inequality constraint, given Jacobian.
        res = fmin_slsqp(self.fun, [-1.0, 1.0],
                         fprime=self.jac, args=(-1.0, ),
                         ieqcons = [self.f_ieqcon],
                         iprint = 0, full_output = 1)
        x, fx, its, imode, smode = res
        assert_(imode == 0, imode)
        assert_array_almost_equal(x, [2, 1], decimal=3)

    def test_bound_equality_given2(self):
        # SLSQP: bounds, eq. const., given jac. for fun. and const.
        res = fmin_slsqp(self.fun, [-1.0, 1.0],
                         fprime=self.jac, args=(-1.0, ),
                         bounds = [(-0.8, 1.), (-1, 0.8)],
                         f_eqcons = self.f_eqcon,
                         fprime_eqcons = self.fprime_eqcon,
                         iprint = 0, full_output = 1)
        x, fx, its, imode, smode = res
        assert_(imode == 0, imode)
        assert_array_almost_equal(x, [0.8, 0.8], decimal=3)
        assert_(-0.8 <= x[0] <= 1)
        assert_(-1 <= x[1] <= 0.8)

    def test_scalar_constraints(self):
        # Regression test for gh-2182
        x = fmin_slsqp(lambda z: z**2, [3.],
                       ieqcons=[lambda z: z[0] - 1],
                       iprint=0)
        assert_array_almost_equal(x, [1.])

        x = fmin_slsqp(lambda z: z**2, [3.],
                       f_ieqcons=lambda z: [z[0] - 1],
                       iprint=0)
        assert_array_almost_equal(x, [1.])

    def test_integer_bounds(self):
        # This should not raise an exception
        fmin_slsqp(lambda z: z**2 - 1, [0], bounds=[[0, 1]], iprint=0)

    def test_array_bounds(self):
        # NumPy used to treat n-dimensional 1-element arrays as scalars
        # in some cases.  The handling of `bounds` by `fmin_slsqp` still
        # supports this behavior.
        bounds = [(-np.inf, np.inf), (np.array([2]), np.array([3]))]
        x = fmin_slsqp(lambda z: np.sum(z**2 - 1), [2.5, 2.5], bounds=bounds,
                       iprint=0)
        assert_array_almost_equal(x, [0, 2])

    def test_obj_must_return_scalar(self):
        # Regression test for Github Issue #5433
        # If objective function does not return a scalar, raises ValueError
        with assert_raises(ValueError):
            fmin_slsqp(lambda x: [0, 1], [1, 2, 3])

    def test_obj_returns_scalar_in_list(self):
        # Test for Github Issue #5433 and PR #6691
        # Objective function should be able to return length-1 Python list
        #  containing the scalar
        fmin_slsqp(lambda x: [0], [1, 2, 3], iprint=0)

    def test_callback(self):
        # Minimize, method='SLSQP': unbounded, approximated jacobian. Check for callback
        callback = MyCallBack()
        res = minimize(self.fun, [-1.0, 1.0], args=(-1.0, ),
                       method='SLSQP', callback=callback, options=self.opts)
        assert res.success
        assert res.message
        assert callback.been_called
        assert_equal(callback.ncalls, res['nit'])

        res = minimize(
            self.fun,
            [-1.0, 1.0],
            args=(-1.0, ),
            method='SLSQP',
            callback=callback.callback2,
            options=self.opts
        )
        assert res.success
        assert callback.been_called

        res = minimize(
            self.fun,
            [-1.0, 1.0],
            args=(-1.0, ),
            method='SLSQP',
            callback=callback.callback3,
            options=self.opts
        )
        assert not res.success
        assert res.message.startswith("`callback` raised `StopIteration`")

    def test_inconsistent_linearization(self):
        # SLSQP must be able to solve this problem, even if the
        # linearized problem at the starting point is infeasible.

        # Linearized constraints are
        #
        #    2*x0[0]*x[0] >= 1
        #
        # At x0 = [0, 1], the second constraint is clearly infeasible.
        # This triggers a call with n2==1 in the LSQ subroutine.
        x = [0, 1]
        def f1(x):
            return x[0] + x[1] - 2
        def f2(x):
            return x[0] ** 2 - 1
        sol = minimize(
            lambda x: x[0]**2 + x[1]**2,
            x,
            constraints=({'type':'eq','fun': f1},
                         {'type':'ineq','fun': f2}),
            bounds=((0,None), (0,None)),
            method='SLSQP')
        x = sol.x

        assert_allclose(f1(x), 0, atol=1e-8)
        assert_(f2(x) >= -1e-8)
        assert_(sol.success, sol)

    def test_regression_5743(self):
        # SLSQP must not indicate success for this problem,
        # which is infeasible.
        x = [1, 2]
        sol = minimize(
            lambda x: x[0]**2 + x[1]**2,
            x,
            constraints=({'type':'eq','fun': lambda x: x[0]+x[1]-1},
                         {'type':'ineq','fun': lambda x: x[0]-2}),
            bounds=((0,None), (0,None)),
            method='SLSQP')
        assert_(not sol.success, sol)

    def test_gh_6676(self):
        def func(x):
            return (x[0] - 1)**2 + 2*(x[1] - 1)**2 + 0.5*(x[2] - 1)**2

        sol = minimize(func, [0, 0, 0], method='SLSQP')
        assert_(sol.jac.shape == (3,))

    def test_invalid_bounds(self):
        # Raise correct error when lower bound is greater than upper bound.
        # See Github issue 6875.
        bounds_list = [
            ((1, 2), (2, 1)),
            ((2, 1), (1, 2)),
            ((2, 1), (2, 1)),
            ((np.inf, 0), (np.inf, 0)),
            ((1, -np.inf), (0, 1)),
        ]
        for bounds in bounds_list:
            with assert_raises(ValueError):
                minimize(self.fun, [-1.0, 1.0], bounds=bounds, method='SLSQP')

    def test_bounds_clipping(self):
        #
        # SLSQP returns bogus results for initial guess out of bounds, gh-6859
        #
        def f(x):
            return (x[0] - 1)**2

        sol = minimize(f, [10], method='slsqp', bounds=[(None, 0)])
        assert_(sol.success)
        assert_allclose(sol.x, 0, atol=1e-10)

        sol = minimize(f, [-10], method='slsqp', bounds=[(2, None)])
        assert_(sol.success)
        assert_allclose(sol.x, 2, atol=1e-10)

        sol = minimize(f, [-10], method='slsqp', bounds=[(None, 0)])
        assert_(sol.success)
        assert_allclose(sol.x, 0, atol=1e-10)

        sol = minimize(f, [10], method='slsqp', bounds=[(2, None)])
        assert_(sol.success)
        assert_allclose(sol.x, 2, atol=1e-10)

        sol = minimize(f, [-0.5], method='slsqp', bounds=[(-1, 0)])
        assert_(sol.success)
        assert_allclose(sol.x, 0, atol=1e-10)

        sol = minimize(f, [10], method='slsqp', bounds=[(-1, 0)])
        assert_(sol.success)
        assert_allclose(sol.x, 0, atol=1e-10)

    def test_infeasible_initial(self):
        # Check SLSQP behavior with infeasible initial point
        def f(x):
            x, = x
            return x*x - 2*x + 1

        cons_u = [{'type': 'ineq', 'fun': lambda x: 0 - x}]
        cons_l = [{'type': 'ineq', 'fun': lambda x: x - 2}]
        cons_ul = [{'type': 'ineq', 'fun': lambda x: 0 - x},
                   {'type': 'ineq', 'fun': lambda x: x + 1}]

        sol = minimize(f, [10], method='slsqp', constraints=cons_u)
        assert_(sol.success)
        assert_allclose(sol.x, 0, atol=1e-10)

        sol = minimize(f, [-10], method='slsqp', constraints=cons_l)
        assert_(sol.success)
        assert_allclose(sol.x, 2, atol=1e-10)

        sol = minimize(f, [-10], method='slsqp', constraints=cons_u)
        assert_(sol.success)
        assert_allclose(sol.x, 0, atol=1e-10)

        sol = minimize(f, [10], method='slsqp', constraints=cons_l)
        assert_(sol.success)
        assert_allclose(sol.x, 2, atol=1e-10)

        sol = minimize(f, [-0.5], method='slsqp', constraints=cons_ul)
        assert_(sol.success)
        assert_allclose(sol.x, 0, atol=1e-10)

        sol = minimize(f, [10], method='slsqp', constraints=cons_ul)
        assert_(sol.success)
        assert_allclose(sol.x, 0, atol=1e-10)

    @pytest.mark.xfail(scipy.show_config(mode='dicts')['Compilers']['fortran']['name']
                       == "intel-llvm",
                       reason="Runtime warning due to floating point issues, not logic")
    def test_inconsistent_inequalities(self):
        # gh-7618

        def cost(x):
            return -1 * x[0] + 4 * x[1]

        def ineqcons1(x):
            return x[1] - x[0] - 1

        def ineqcons2(x):
            return x[0] - x[1]

        # The inequalities are inconsistent, so no solution can exist:
        #
        # x1 >= x0 + 1
        # x0 >= x1

        x0 = (1,5)
        bounds = ((-5, 5), (-5, 5))
        cons = (dict(type='ineq', fun=ineqcons1), dict(type='ineq', fun=ineqcons2))
        res = minimize(cost, x0, method='SLSQP', bounds=bounds, constraints=cons)

        assert_(not res.success)

    def test_new_bounds_type(self):
        def f(x):
            return x[0] ** 2 + x[1] ** 2
        bounds = Bounds([1, 0], [np.inf, np.inf])
        sol = minimize(f, [0, 0], method='slsqp', bounds=bounds)
        assert_(sol.success)
        assert_allclose(sol.x, [1, 0])

    def test_nested_minimization(self):

        class NestedProblem:

            def __init__(self):
                self.F_outer_count = 0

            def F_outer(self, x):
                self.F_outer_count += 1
                if self.F_outer_count > 1000:
                    raise Exception("Nested minimization failed to terminate.")
                inner_res = minimize(self.F_inner, (3, 4), method="SLSQP")
                assert_(inner_res.success)
                assert_allclose(inner_res.x, [1, 1])
                return x[0]**2 + x[1]**2 + x[2]**2

            def F_inner(self, x):
                return (x[0] - 1)**2 + (x[1] - 1)**2

            def solve(self):
                outer_res = minimize(self.F_outer, (5, 5, 5), method="SLSQP")
                assert_(outer_res.success)
                assert_allclose(outer_res.x, [0, 0, 0])

        problem = NestedProblem()
        problem.solve()

    def test_gh1758(self):
        # the test suggested in gh1758
        # https://nlopt.readthedocs.io/en/latest/NLopt_Tutorial/
        # implement two equality constraints, in R^2.
        def fun(x):
            return np.sqrt(x[1])

        def f_eqcon(x):
            """ Equality constraint """
            return x[1] - (2 * x[0]) ** 3

        def f_eqcon2(x):
            """ Equality constraint """
            return x[1] - (-x[0] + 1) ** 3

        c1 = {'type': 'eq', 'fun': f_eqcon}
        c2 = {'type': 'eq', 'fun': f_eqcon2}

        res = minimize(fun, [8, 0.25], method='SLSQP',
                       constraints=[c1, c2], bounds=[(-0.5, 1), (0, 8)])

        np.testing.assert_allclose(res.fun, 0.5443310539518)
        np.testing.assert_allclose(res.x, [0.33333333, 0.2962963])
        assert res.success

    def test_gh9640(self):
        cons = ({'type': 'ineq', 'fun': lambda x: -x[0] - x[1] - 3},
                {'type': 'ineq', 'fun': lambda x: x[1] + x[2] - 2})
        bnds = ((-2, 2), (-2, 2), (-2, 2))

        def target(x):
            return 1
        x0 = [-1.8869783504471584, -0.640096352696244, -0.8174212253407696]
        res = minimize(target, x0, method='SLSQP', bounds=bnds, constraints=cons,
                       options={'disp':False, 'maxiter':10000})

        # The problem is infeasible, so it cannot succeed
        assert not res.success

    def test_parameters_stay_within_bounds(self):
        # gh11403. For some problems the SLSQP Fortran code suggests a step
        # outside one of the lower/upper bounds. When this happens
        # approx_derivative complains because it's being asked to evaluate
        # a gradient outside its domain.

        # gh21872, removal of random initial position, replacing with specific
        # starting point, because success/fail depends on the seed used.
        bounds = Bounds(np.array([0.1]), np.array([1.0]))
        x0 = np.array(bounds.lb + (bounds.ub - bounds.lb) *
                      0.417022004702574)

        def f(x):
            assert (x >= bounds.lb).all()
            return np.linalg.norm(x)
        # The following should not raise any warnings which was the case, with the
        # old Fortran code.
        res = minimize(f, x0, method='SLSQP', bounds=bounds)
        assert res.success
    
    def test_fnc_with_simple_invalid_region(self):
        # Test function where slsqp is encouraged to evalutate the square
        # root of a negative design variable to demonstrate recoverability
        # by raising an AnalysisError.

        # Tracks the number of invalid steps to ensure the test is valid
        invalid_step_count = [0]

        def fun(x):
            return (x - 0.5)**4 + 1

        def jac_fun(x):
            return 4*(x - 0.5)**3

        def con(x):
            if x < 0:
                invalid_step_count[0] += 1
                raise AnalysisError()
            return np.sqrt(x)

        def jac_con(x):
            return 0.5*x**(-0.5)

        constraints = [
            NonlinearConstraint(
                con,
                lb=0.0,
                ub=1.0,
                jac=jac_con,
            )
        ]

        res = minimize(
            fun,
            x0=10.0,
            jac=jac_fun,
            method='SLSQP',
            constraints=constraints,
            tol=1e-9
            )
        
        # Ensures the optimization converges
        assert res.success
        np.testing.assert_allclose(res.fun, 1.0, rtol=1e-4)
        np.testing.assert_allclose(res.x, [0.5], rtol=2e-2)

        # Ensure the test is valid (the function entered the invalid region)
        assert invalid_step_count[0] > 0

    def test_func_with_invalid_obj_region(self):
        # Test where the objective function has regions where evaluation
        # is not possible. The objective function is based on the Bean
        # function defined in Engineering Design Optimization by Dr.
        # Martins and Dr. Ning. The invalid regions is a periodic set of
        # holes in the design space. The minima is in the valid design
        # space.

        # Tracks the number of invalid steps to ensure the test is valid
        invalid_step_count = [0]

        def fun(x):
            if np.cos(2.0*x[0]) + np.cos(2.0*x[1]) > 0.5:
                invalid_step_count[0] += 1
                raise AnalysisError("Invalid Region")
            return (1-x[0])**2 + (1-x[1])**2 + 0.5*(2.0*x[1]-x[0]**2)**2
        
        res = minimize(fun, [2.0, 2.0], method="SLSQP")

        # Ensures the optimization converges
        assert res.success
        np.testing.assert_allclose(res.fun, 0.09194, rtol=1e-3)
        np.testing.assert_allclose(res.x, [1.2134, 0.82414], rtol=1e-3)


        # Ensure the test is valid (the function entered the invalid region)
        assert invalid_step_count[0] > 0

    def test_func_with_value_error_fails(self):
        # Test to ensure ValueError causes `minimize` to fail

        def fun(x):
            if x < 10:
                raise ValueError()
            return x**2
        
        with pytest.raises(ValueError):
            minimize(fun, [10.0], method="SLSQP")

    def test_func_with_invalid_con_region(self):
        # Test where the constraint function has regions where evaluation
        # is not possible. The objective function is based on the Bean
        # function defined in Engineering Design Optimization by Dr.
        # Martins and Dr. Ning. The invalid regions is a periodic set of
        # holes in the design space. The minima is in the valid design
        # space. The constraint is intentionally inactive.

        # Tracks the number of invalid steps to ensure the test is valid
        invalid_step_count = [0]

        def fun(x):
            return (1-x[0])**2 + (1-x[1])**2 + 0.5*(2.0*x[1]-x[0]**2)**2
        
        def con(x):
            if np.cos(2.0*x[0]) + np.cos(2.0*x[1]) > 0.5:
                invalid_step_count[0] += 1
                raise AnalysisError("Invalid Region")
            return x[0] + 5
        
        res = minimize(fun, [2.0, 2.0], method="SLSQP", constraints=[
            {"type": "ineq", "fun": con}
        ])

        # Ensures the optimization converges
        assert res.success
        np.testing.assert_allclose(res.fun, 0.09194, rtol=1e-3)
        np.testing.assert_allclose(res.x, [1.2134, 0.82414], rtol=1e-3)

        # Ensure the test is valid (the function entered the invalid region)
        assert invalid_step_count[0] > 0

    def test_func_with_minima_in_invalid_region(self):
        # Test where the objective function has regions where evaluation
        # is not possible and the minima lies within this region. In this
        # case, the optimization will fail and exit to the last point where
        # the function evaluation was possible. The objective function is 
        # based on the Bean function defined in Engineering Design
        # Optimization by Dr. Martins and Dr. Ning.

        def fun(x):
            if x[0] < 2.0:
                raise AnalysisError("Invalid Region")
            return (1-x[0])**2 + (1-x[1])**2 + 0.5*(2.0*x[1]-x[0]**2)**2

        res = minimize(fun, [5.0, 10.0], method="SLSQP")

        # Ensure optimization fails
        assert not res.success
        assert res.status == 9  # Iteration limit failure

        # Check final point
        assert res.x[0] >= 2  # Does not end in invalid region
        fun(res.x)  # Ensure func can be evaluated without error at final point

        # Check that the optimization exits near the border of the invalid region
        assert res.x[0] < 2.1

    def test_divergent_simulation_problem(self):
        # Test problem involving a numerical simulation with scaling
        # such that it is intentionally divergent. Tests that SLSQP
        # is able to recover from the divence via raising an AnalysisError

        prob = DivergentSimulationProblem()
        constraints = [{'type': 'ineq', 'fun': prob.cons, 'jac': prob.cons_jac}]

        x0 = 1.0
        res = minimize(
            prob.obj,
            x0/prob.SCALE_X,
            method='SLSQP',
            jac=prob.obj_jac,
            constraints=constraints,
            tol=1e-6
        )

        assert res.success
        np.testing.assert_allclose(res.fun + 1, prob.fglob + 1.0, rtol=1e-3)
        np.testing.assert_allclose(res.x, [prob.xglob], rtol=1e-3)

        # Ensure the test is valid (the function entered the invalid region)
        assert prob.divergent_simulations >= 1

    def test_callbacks_in_invalid_region(self):
        # Tests that AnalysisError support does not interfere with callbacks

        all_x = []
        def callback(x):
            all_x.append(x)

        def fun(x):
            if x[0] < 2.0:
                raise AnalysisError("Invalid Region")
            return x**2

        # Optimization will fail here, but that's not relevant
        minimize(fun, 5.0, method="SLSQP", callback=callback)

        # Check that the callback still triggered even in the invalid region
        assert any([x < 2.0 for x in all_x])

def test_slsqp_segfault_wrong_workspace_computation():
    # See gh-14915
    # This problem is not well-defined, however should not cause a segfault.
    # The previous F77 workspace computation did not handle only equality-
    # constrained problems correctly.
    rng = np.random.default_rng(1742651087222879)
    x = rng.uniform(size=[22,365])
    target = np.linspace(0.9, 4.0, 50)

    def metric(v, weights):
        return [[0, 0],[1, 1]]

    def efficient_metric(v, target):
        def metric_a(weights):
            return metric(v, weights)[1][0]

        def metric_b(weights, v):
            return metric(v, weights)[0][0]

        constraints = ({'type': 'eq', 'fun': lambda x: metric_a(x) - target},
                       {'type': 'eq', 'fun': lambda x: np.sum(x) - 1})
        weights = np.array([len(v)*[1./len(v)]])[0]
        result = minimize(metric_b,
                          weights,
                          args=(v,),
                          method='SLSQP',
                          constraints=constraints)
        return result

    efficient_metric(x, target)


class DivergentSimulationProblem:
    """
    Problem that optimizes a logistic-map parameter so that the
    final population density is 0.25, while enforcing monotonic increase in
    population over time via inequality constraints.

    The logistic map is iterated for 51 steps starting from ``y0 = 0.1``.
    The objective minimizes ``(y_T - 0.25)**2`` and the constraints require
    ``np.diff(y) >= 0``.

    The logistic map is defined as:

    .. math::

        y_{t+1} = x * y_t * (1 - y_t)

    where ``x`` is the parameter to be optimized.

    Notes
    -----
    * The simulation raises ``AnalysisError`` when the state or its
      sensitivities diverge to infinity or overflow is encountered
    * An intentionally poor scaling factor ``SCALE_X`` is used to drive
    the optimizer towards the divergent region

    References
    ----------
    .. [1] R. M. May, "Simple mathematical models with very complicated dynamics,"
           Nature, 261, 459-467 (1976). https://doi.org/10.1038/261459a0
    """

    def __init__(self):

        # Scale factor to encourage large steps
        self.SCALE_X = 10.0

        # Analytic Solution
        self.fglob = 0.0
        self.xglob = (4.0/3.0) / self.SCALE_X  # Won't be exact as we only do 51 steps

        # Track function calls where simulation diverges
        self.divergent_simulations = 0

        # Sim Cache
        self.cache_x = None
        self.cache_res = None

    # Avoid excessive calls to simulation
    def get_sim_results(self, x):
        x = float(x[0])
        if self.cache_x is None or x != self.cache_x:
            self.cache_res = self.simulate(x)
            self.cache_x = x
        return self.cache_res

    def scale_x(self, x):
        return self.SCALE_X * x

    def simulate(self, x):
        num_iter = 51
        y, dy_dx = [0.0]*num_iter, [0.0]*num_iter
        y[0] = 0.1

        # Run simulation and propagate gradients
        for i in range(1, num_iter):
            y[i] = x*y[i-1]*(1.0 - y[i-1])
            dy_dx[i] = y[i-1]*(1.0 - y[i-1]) + x*(1.0 - 2.0*y[i-1])*dy_dx[i-1]

        # Check for divergent simulation
        if any([np.any(np.isinf(y)), np.any(np.isinf(dy_dx))]):
            self.divergent_simulations += 1
            raise AnalysisError()

        return y, dy_dx

    def obj(self, x):
        x_scaled = self.scale_x(x)
        y, _ = self.get_sim_results(x_scaled)
        try:
            obj = (y[-1] - 0.25)**2
        except OverflowError:
            self.divergent_simulations += 1
            raise AnalysisError()
        return obj

    def obj_jac(self, x):
        x_scaled = self.scale_x(x)
        y, dy_dx = self.get_sim_results(x_scaled)
        dobj_dx = 2.0 * (y[-1] - 0.25) * dy_dx[-1]
        dobj_dxscaled =  dobj_dx * self.SCALE_X
        return dobj_dxscaled

    def cons(self, x):
        x_scaled = self.scale_x(x)
        y, _ = self.get_sim_results(x_scaled)
        return np.diff(y)

    def cons_jac(self, x):
        x_scaled = self.scale_x(x)
        _, dy_dx = self.get_sim_results(x_scaled)
        dcon_dx = np.diff(dy_dx).reshape(-1, 1)
        dcon_dxscaled = dcon_dx * self.SCALE_X
        return dcon_dxscaled
