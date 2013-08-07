# Authors: Nils Wagner, Ed Schofield, Pauli Virtanen, John Travers
"""
Tests for numerical integration.
"""
from __future__ import division, print_function, absolute_import

import numpy
from numpy import arange, zeros, array, dot, sqrt, cos, sin, eye, pi, exp, \
                  allclose

from scipy.lib.six.moves import xrange

from numpy.testing import assert_, TestCase, run_module_suite, \
        assert_array_almost_equal, assert_raises, assert_allclose, \
        assert_array_equal, assert_equal
from scipy.integrate import odeint, ode, complex_ode

#------------------------------------------------------------------------------
# Test ODE integrators
#------------------------------------------------------------------------------


class TestOdeint(TestCase):
    # Check integrate.odeint
    def _do_problem(self, problem):
        t = arange(0.0, problem.stop_t, 0.05)
        z, infodict = odeint(problem.f, problem.z0, t, full_output=True)
        assert_(problem.verify(z, t))

    def test_odeint(self):
        for problem_cls in PROBLEMS:
            problem = problem_cls()
            if problem.cmplx:
                continue
            self._do_problem(problem)


class TestOde(TestCase):
    # Check integrate.ode
    def _do_problem(self, problem, integrator, method='adams'):

        # ode has callback arguments in different order than odeint
        f = lambda t, z: problem.f(z, t)
        jac = None
        if hasattr(problem, 'jac'):
            jac = lambda t, z: problem.jac(z, t)

        ig = ode(f, jac)
        ig.set_integrator(integrator,
                          atol=problem.atol/10,
                          rtol=problem.rtol/10,
                          method=method)
        ig.set_initial_value(problem.z0, t=0.0)
        z = ig.integrate(problem.stop_t)

        assert_(ig.successful(), (problem, method))
        assert_(problem.verify(array([z]), problem.stop_t), (problem, method))

    def test_vode(self):
        # Check the vode solver
        for problem_cls in PROBLEMS:
            problem = problem_cls()
            if problem.cmplx:
                continue
            if not problem.stiff:
                self._do_problem(problem, 'vode', 'adams')
            self._do_problem(problem, 'vode', 'bdf')

    def test_zvode(self):
        # Check the zvode solver
        for problem_cls in PROBLEMS:
            problem = problem_cls()
            if not problem.stiff:
                self._do_problem(problem, 'zvode', 'adams')
            self._do_problem(problem, 'zvode', 'bdf')

    def test_lsoda(self):
        # Check the lsoda solver
        for problem_cls in PROBLEMS:
            problem = problem_cls()
            if problem.cmplx:
                continue
            self._do_problem(problem, 'lsoda')

    def test_dopri5(self):
        # Check the dopri5 solver
        for problem_cls in PROBLEMS:
            problem = problem_cls()
            if problem.cmplx:
                continue
            if problem.stiff:
                continue
            if hasattr(problem, 'jac'):
                continue
            self._do_problem(problem, 'dopri5')

    def test_dop853(self):
        # Check the dop853 solver
        for problem_cls in PROBLEMS:
            problem = problem_cls()
            if problem.cmplx:
                continue
            if problem.stiff:
                continue
            if hasattr(problem, 'jac'):
                continue
            self._do_problem(problem, 'dop853')

    def test_concurrent_fail(self):
        for sol in ('vode', 'zvode', 'lsoda'):
            f = lambda t, y: 1.0

            r = ode(f).set_integrator(sol)
            r.set_initial_value(0, 0)

            r2 = ode(f).set_integrator(sol)
            r2.set_initial_value(0, 0)

            r.integrate(r.t + 0.1)
            r2.integrate(r2.t + 0.1)

            assert_raises(RuntimeError, r.integrate, r.t + 0.1)

    def test_concurrent_ok(self):
        f = lambda t, y: 1.0

        for k in xrange(3):
            for sol in ('vode', 'zvode', 'lsoda', 'dopri5', 'dop853'):
                r = ode(f).set_integrator(sol)
                r.set_initial_value(0, 0)

                r2 = ode(f).set_integrator(sol)
                r2.set_initial_value(0, 0)

                r.integrate(r.t + 0.1)
                r2.integrate(r2.t + 0.1)
                r2.integrate(r2.t + 0.1)

                assert_allclose(r.y, 0.1)
                assert_allclose(r2.y, 0.2)

            for sol in ('dopri5', 'dop853'):
                r = ode(f).set_integrator(sol)
                r.set_initial_value(0, 0)

                r2 = ode(f).set_integrator(sol)
                r2.set_initial_value(0, 0)

                r.integrate(r.t + 0.1)
                r.integrate(r.t + 0.1)
                r2.integrate(r2.t + 0.1)
                r.integrate(r.t + 0.1)
                r2.integrate(r2.t + 0.1)

                assert_allclose(r.y, 0.3)
                assert_allclose(r2.y, 0.2)


class TestComplexOde(TestCase):
    # Check integrate.complex_ode
    def _do_problem(self, problem, integrator, method='adams'):

        # ode has callback arguments in different order than odeint
        f = lambda t, z: problem.f(z, t)
        jac = None
        if hasattr(problem, 'jac'):
            jac = lambda t, z: problem.jac(z, t)
        ig = complex_ode(f, jac)
        ig.set_integrator(integrator,
                          atol=problem.atol/10,
                          rtol=problem.rtol/10,
                          method=method)
        ig.set_initial_value(problem.z0, t=0.0)
        z = ig.integrate(problem.stop_t)
        z2 = ig.y

        assert_array_equal(z, z2)
        assert_(ig.successful(), (problem, method))
        assert_(problem.verify(array([z]), problem.stop_t), (problem, method))

    def test_vode(self):
        # Check the vode solver
        for problem_cls in PROBLEMS:
            problem = problem_cls()
            if not problem.stiff:
                self._do_problem(problem, 'vode', 'adams')
            else:
                self._do_problem(problem, 'vode', 'bdf')

    def test_lsoda(self):
        # Check the lsoda solver
        for problem_cls in PROBLEMS:
            problem = problem_cls()
            self._do_problem(problem, 'lsoda')

    def test_dopri5(self):
        # Check the dopri5 solver
        for problem_cls in PROBLEMS:
            problem = problem_cls()
            if problem.stiff:
                continue
            if hasattr(problem, 'jac'):
                continue
            self._do_problem(problem, 'dopri5')

    def test_dop853(self):
        # Check the dop853 solver
        for problem_cls in PROBLEMS:
            problem = problem_cls()
            if problem.stiff:
                continue
            if hasattr(problem, 'jac'):
                continue
            self._do_problem(problem, 'dop853')

class TestSolout(TestCase):
    # Check integrate.ode correctly handles solout for dopri5 and dop853
    def _run_solout_test(self, integrator):
        # Check correct usage of solout
        ts = []
        ys = []
        t0 = 0.0
        tend = 10.0
        y0 = [1.0, 2.0]
        def solout(t, y):
            ts.append(t)
            ys.append(y.copy())
        def rhs(t, y):
            return [y[0] + y[1], -y[1]**2]
        ig = ode(rhs).set_integrator(integrator)
        ig.set_solout(solout)
        ig.set_initial_value(y0, t0)
        ret = ig.integrate(tend)
        assert_array_equal(ys[0], y0)
        assert_array_equal(ys[-1], ret)
        assert_equal(ts[0], t0)
        assert_equal(ts[-1], tend)

    def test_solout(self):
        for integrator in ('dopri5', 'dop853'):
            self._run_solout_test(integrator)

    def _run_solout_break_test(self, integrator):
        # Check correct usage of stopping via solout
        ts = []
        ys = []
        t0 = 0.0
        tend = 10.0
        y0 = [1.0, 2.0]
        def solout(t, y):
            ts.append(t)
            ys.append(y.copy())
            if t > tend/2.0:
                return -1
        def rhs(t, y):
            return [y[0] + y[1], -y[1]**2]
        ig = ode(rhs).set_integrator(integrator)
        ig.set_solout(solout)
        ig.set_initial_value(y0, t0)
        ret = ig.integrate(tend)
        assert_array_equal(ys[0], y0)
        assert_array_equal(ys[-1], ret)
        assert_equal(ts[0], t0)
        assert_(ts[-1] > tend/2.0)
        assert_(ts[-1] < tend)

    def test_solout_break(self):
        for integrator in ('dopri5', 'dop853'):
            self._run_solout_break_test(integrator)


class TestComplexSolout(TestCase):
    # Check integrate.ode correctly handles solout for dopri5 and dop853
    def _run_solout_test(self, integrator):
        # Check correct usage of solout
        ts = []
        ys = []
        t0 = 0.0
        tend = 20.0
        y0 = [0.0]
        def solout(t, y):
            ts.append(t)
            ys.append(y.copy())
        def rhs(t, y):
            return [1.0/(t - 10.0 - 1j)]
        ig = complex_ode(rhs).set_integrator(integrator)
        ig.set_solout(solout)
        ig.set_initial_value(y0, t0)
        ret = ig.integrate(tend)
        assert_array_equal(ys[0], y0)
        assert_array_equal(ys[-1], ret)
        assert_equal(ts[0], t0)
        assert_equal(ts[-1], tend)

    def test_solout(self):
        for integrator in ('dopri5', 'dop853'):
            self._run_solout_test(integrator)

    def _run_solout_break_test(self, integrator):
        # Check correct usage of stopping via solout
        ts = []
        ys = []
        t0 = 0.0
        tend = 20.0
        y0 = [0.0]
        def solout(t, y):
            ts.append(t)
            ys.append(y.copy())
            if t > tend/2.0:
                return -1
        def rhs(t, y):
            return [1.0/(t - 10.0 - 1j)]
        ig = complex_ode(rhs).set_integrator(integrator)
        ig.set_solout(solout)
        ig.set_initial_value(y0, t0)
        ret = ig.integrate(tend)
        assert_array_equal(ys[0], y0)
        assert_array_equal(ys[-1], ret)
        assert_equal(ts[0], t0)
        assert_(ts[-1] > tend/2.0)
        assert_(ts[-1] < tend)

    def test_solout_break(self):
        for integrator in ('dopri5', 'dop853'):
            self._run_solout_break_test(integrator)


#------------------------------------------------------------------------------
# Test problems
#------------------------------------------------------------------------------


class ODE:
    """
    ODE problem
    """
    stiff = False
    cmplx = False
    stop_t = 1
    z0 = []

    atol = 1e-6
    rtol = 1e-5


class SimpleOscillator(ODE):
    r"""
    Free vibration of a simple oscillator::
        m \ddot{u} + k u = 0, u(0) = u_0 \dot{u}(0) \dot{u}_0
    Solution::
        u(t) = u_0*cos(sqrt(k/m)*t)+\dot{u}_0*sin(sqrt(k/m)*t)/sqrt(k/m)
    """
    stop_t = 1 + 0.09
    z0 = array([1.0, 0.1], float)

    k = 4.0
    m = 1.0

    def f(self, z, t):
        tmp = zeros((2,2), float)
        tmp[0,1] = 1.0
        tmp[1,0] = -self.k / self.m
        return dot(tmp, z)

    def verify(self, zs, t):
        omega = sqrt(self.k / self.m)
        u = self.z0[0]*cos(omega*t)+self.z0[1]*sin(omega*t)/omega
        return allclose(u, zs[:,0], atol=self.atol, rtol=self.rtol)


class ComplexExp(ODE):
    r"""The equation :lm:`\dot u = i u`"""
    stop_t = 1.23*pi
    z0 = exp([1j,2j,3j,4j,5j])
    cmplx = True

    def f(self, z, t):
        return 1j*z

    def jac(self, z, t):
        return 1j*eye(5)

    def verify(self, zs, t):
        u = self.z0 * exp(1j*t)
        return allclose(u, zs, atol=self.atol, rtol=self.rtol)


class Pi(ODE):
    r"""Integrate 1/(t + 1j) from t=-10 to t=10"""
    stop_t = 20
    z0 = [0]
    cmplx = True

    def f(self, z, t):
        return array([1./(t - 10 + 1j)])

    def verify(self, zs, t):
        u = -2j*numpy.arctan(10)
        return allclose(u, zs[-1,:], atol=self.atol, rtol=self.rtol)

PROBLEMS = [SimpleOscillator, ComplexExp, Pi]

#------------------------------------------------------------------------------


def f(t, x):
    dxdt = [x[1], -x[0]]
    return dxdt


def jac(t, x):
    j = array([[0.0, 1.0],
               [-1.0, 0.0]])
    return j


def f1(t, x, omega):
    dxdt = [omega*x[1], -omega*x[0]]
    return dxdt


def jac1(t, x, omega):
    j = array([[0.0, omega],
               [-omega, 0.0]])
    return j


def f2(t, x, omega1, omega2):
    dxdt = [omega1*x[1], -omega2*x[0]]
    return dxdt


def jac2(t, x, omega1, omega2):
    j = array([[0.0, omega1],
               [-omega2, 0.0]])
    return j


def fv(t, x, omega):
    dxdt = [omega[0]*x[1], -omega[1]*x[0]]
    return dxdt


def jacv(t, x, omega):
    j = array([[0.0, omega[0]],
               [-omega[1], 0.0]])
    return j


class ODECheckParameterUse(object):
    """Call an ode-class solver with several cases of parameter use."""

    # This class is intentionally not a TestCase subclass.
    # solver_name must be set before tests can be run with this class.

    # Set these in subclasses.
    solver_name = ''
    solver_uses_jac = False

    def _get_solver(self, f, jac):
        solver = ode(f, jac)
        if self.solver_uses_jac:
            solver.set_integrator(self.solver_name, atol=1e-9, rtol=1e-7,
                                  with_jacobian=self.solver_uses_jac)
        else:
            # XXX Shouldn't set_integrator *always* accept the keyword arg
            # 'with_jacobian', and perhaps raise an exception if it is set
            # to True if the solver can't actually use it?
            solver.set_integrator(self.solver_name, atol=1e-9, rtol=1e-7)
        return solver

    def _check_solver(self, solver):
        ic = [1.0, 0.0]
        solver.set_initial_value(ic, 0.0)
        solver.integrate(pi)
        assert_array_almost_equal(solver.y, [-1.0, 0.0])

    def test_no_params(self):
        solver = self._get_solver(f, jac)
        self._check_solver(solver)

    def test_one_scalar_param(self):
        solver = self._get_solver(f1, jac1)
        omega = 1.0
        solver.set_f_params(omega)
        if self.solver_uses_jac:
            solver.set_jac_params(omega)
        self._check_solver(solver)

    def test_two_scalar_params(self):
        solver = self._get_solver(f2, jac2)
        omega1 = 1.0
        omega2 = 1.0
        solver.set_f_params(omega1, omega2)
        if self.solver_uses_jac:
            solver.set_jac_params(omega1, omega2)
        self._check_solver(solver)

    def test_vector_param(self):
        solver = self._get_solver(fv, jacv)
        omega = [1.0, 1.0]
        solver.set_f_params(omega)
        if self.solver_uses_jac:
            solver.set_jac_params(omega)
        self._check_solver(solver)


class DOPRI5CheckParameterUse(ODECheckParameterUse, TestCase):
    solver_name = 'dopri5'
    solver_uses_jac = False


class DOP853CheckParameterUse(ODECheckParameterUse, TestCase):
    solver_name = 'dop853'
    solver_uses_jac = False


class VODECheckParameterUse(ODECheckParameterUse, TestCase):
    solver_name = 'vode'
    solver_uses_jac = True


class ZVODECheckParameterUse(ODECheckParameterUse, TestCase):
    solver_name = 'zvode'
    solver_uses_jac = True


class LSODACheckParameterUse(ODECheckParameterUse, TestCase):
    solver_name = 'lsoda'
    solver_uses_jac = True


if __name__ == "__main__":
    run_module_suite()
