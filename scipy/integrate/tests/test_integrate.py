# Authors: Nils Wagner, Ed Schofield, Pauli Virtanen, John Travers
"""
Tests for numerical integration.
"""

import numpy
from numpy import arange, zeros, array, dot, sqrt, cos, sin, eye, pi, exp, \
                  allclose

from numpy.testing import *
from scipy.integrate import odeint, ode, complex_ode

#------------------------------------------------------------------------------
# Test ODE integrators
#------------------------------------------------------------------------------

class TestOdeint(TestCase):
    """
    Check integrate.odeint
    """
    def _do_problem(self, problem):
        t = arange(0.0, problem.stop_t, 0.05)
        z, infodict = odeint(problem.f, problem.z0, t, full_output=True)
        assert problem.verify(z, t)

    def test_odeint(self):
        for problem_cls in PROBLEMS:
            problem = problem_cls()
            if problem.cmplx: continue
            self._do_problem(problem)

class TestOde(TestCase):
    """
    Check integrate.ode
    """
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

        assert ig.successful(), (problem, method)
        assert problem.verify(array([z]), problem.stop_t), (problem, method)

    def test_vode(self):
        """Check the vode solver"""
        for problem_cls in PROBLEMS:
            problem = problem_cls()
            if problem.cmplx: continue
            if not problem.stiff:
                self._do_problem(problem, 'vode', 'adams')
            self._do_problem(problem, 'vode', 'bdf')

    def test_zvode(self):
        """Check the zvode solver"""
        for problem_cls in PROBLEMS:
            problem = problem_cls()
            if not problem.stiff:
                self._do_problem(problem, 'zvode', 'adams')
            self._do_problem(problem, 'zvode', 'bdf')

    def test_dopri5(self):
        """Check the dopri5 solver"""
        for problem_cls in PROBLEMS:
            problem = problem_cls()
            if problem.cmplx: continue
            if problem.stiff: continue
            if hasattr(problem, 'jac'): continue
            self._do_problem(problem, 'dopri5')
            
    def test_dop853(self):
        """Check the dop853 solver"""
        for problem_cls in PROBLEMS:
            problem = problem_cls()
            if problem.cmplx: continue
            if problem.stiff: continue
            if hasattr(problem, 'jac'): continue
            self._do_problem(problem, 'dop853')

class TestComplexOde(TestCase):
    """
    Check integrate.complex_ode
    """
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

        assert ig.successful(), (problem, method)
        assert problem.verify(array([z]), problem.stop_t), (problem, method)

    def test_vode(self):
        """Check the vode solver"""
        for problem_cls in PROBLEMS:
            problem = problem_cls()
            if not problem.stiff:
                self._do_problem(problem, 'vode', 'adams')
            else:
                self._do_problem(problem, 'vode', 'bdf')

    def test_dopri5(self):
        """Check the dopri5 solver"""
        for problem_cls in PROBLEMS:
            problem = problem_cls()
            if problem.stiff: continue
            if hasattr(problem, 'jac'): continue
            self._do_problem(problem, 'dopri5')
            
    def test_dop853(self):
        """Check the dop853 solver"""
        for problem_cls in PROBLEMS:
            problem = problem_cls()
            if problem.stiff: continue
            if hasattr(problem, 'jac'): continue
            self._do_problem(problem, 'dop853')

#------------------------------------------------------------------------------
# Test problems
#------------------------------------------------------------------------------

class ODE:
    """
    ODE problem
    """
    stiff   = False
    cmplx   = False
    stop_t  = 1
    z0      = []

    atol    = 1e-6
    rtol    = 1e-5

class SimpleOscillator(ODE):
    r"""
    Free vibration of a simple oscillator::
        m \ddot{u} + k u = 0, u(0) = u_0 \dot{u}(0) \dot{u}_0
    Solution::
        u(t) = u_0*cos(sqrt(k/m)*t)+\dot{u}_0*sin(sqrt(k/m)*t)/sqrt(k/m)
    """
    stop_t  = 1 + 0.09
    z0      = array([1.0, 0.1], float)

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
    stop_t  = 1.23*pi
    z0      = exp([1j,2j,3j,4j,5j])
    cmplx   = True

    def f(self, z, t):
        return 1j*z

    def jac(self, z, t):
        return 1j*eye(5)

    def verify(self, zs, t):
        u = self.z0 * exp(1j*t)
        return allclose(u, zs, atol=self.atol, rtol=self.rtol)

class Pi(ODE):
    r"""Integrate 1/(t + 1j) from t=-10 to t=10"""
    stop_t  = 20
    z0      = [0]
    cmplx   = True

    def f(self, z, t):
        return array([1./(t - 10 + 1j)])
    def verify(self, zs, t):
        u = -2j*numpy.arctan(10)
        return allclose(u, zs[-1,:], atol=self.atol, rtol=self.rtol)

PROBLEMS = [SimpleOscillator, ComplexExp, Pi]

#------------------------------------------------------------------------------

if __name__ == "__main__":
    run_module_suite()
