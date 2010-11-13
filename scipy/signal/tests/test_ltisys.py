
import warnings

import numpy as np
from numpy.testing import assert_almost_equal, assert_equal, run_module_suite

from scipy.signal.ltisys import ss2tf, lsim2, impulse2, step2, lti
from scipy.signal.filter_design import BadCoefficients

class TestSS2TF:
    def tst_matrix_shapes(self, p, q, r):
        ss2tf(np.zeros((p, p)),
              np.zeros((p, q)),
              np.zeros((r, p)),
              np.zeros((r, q)), 0)

    def test_basic(self):
        for p, q, r in [
            (3, 3, 3),
            (1, 3, 3),
            (1, 1, 1)]:
            yield self.tst_matrix_shapes, p, q, r


class Test_lsim2(object):

    def test_01(self):
        t = np.linspace(0,10,1001)
        u = np.zeros_like(t)
        # First order system: x'(t) + x(t) = u(t), x(0) = 1.
        # Exact solution is x(t) = exp(-t).
        system = ([1.0],[1.0,1.0])
        tout, y, x = lsim2(system, u, t, X0=[1.0])
        expected_x = np.exp(-tout)
        assert_almost_equal(x[:,0], expected_x)

    def test_02(self):
        t = np.array([0.0, 1.0, 1.0, 3.0])
        u = np.array([0.0, 0.0, 1.0, 1.0])
        # Simple integrator: x'(t) = u(t)
        system = ([1.0],[1.0,0.0])
        tout, y, x = lsim2(system, u, t, X0=[1.0])
        expected_x = np.maximum(1.0, tout)
        assert_almost_equal(x[:,0], expected_x)

    def test_03(self):
        t = np.array([0.0, 1.0, 1.0, 1.1, 1.1, 2.0])
        u = np.array([0.0, 0.0, 1.0, 1.0, 0.0, 0.0])
        # Simple integrator:  x'(t) = u(t)
        system = ([1.0],[1.0, 0.0])
        tout, y, x = lsim2(system, u, t, hmax=0.01)
        expected_x = np.array([0.0, 0.0, 0.0, 0.1, 0.1, 0.1])
        assert_almost_equal(x[:,0], expected_x)

    def test_04(self):
        t = np.linspace(0, 10, 1001)
        u = np.zeros_like(t)
        # Second order system with a repeated root: x''(t) + 2*x(t) + x(t) = 0.
        # With initial conditions x(0)=1.0 and x'(t)=0.0, the exact solution
        # is (1-t)*exp(-t).
        system = ([1.0], [1.0, 2.0, 1.0])
        tout, y, x = lsim2(system, u, t, X0=[1.0, 0.0])
        expected_x = (1.0 - tout) * np.exp(-tout)
        assert_almost_equal(x[:,0], expected_x)

    def test_05(self):
        # The call to lsim2 triggers a "BadCoefficients" warning from
        # scipy.signal.filter_design, but the test passes.  I think the warning
        # is related to the incomplete handling of multi-input systems in
        # scipy.signal.

        # A system with two state variables, two inputs, and one output.
        A = np.array([[-1.0, 0.0], [0.0, -2.0]])
        B = np.array([[1.0, 0.0], [0.0, 1.0]])
        C = np.array([1.0, 0.0])
        D = np.zeros((1,2))

        t = np.linspace(0, 10.0, 101)
        warnings.simplefilter("ignore", BadCoefficients)
        try:
            tout, y, x = lsim2((A,B,C,D), T=t, X0=[1.0, 1.0])
        finally:
            del warnings.filters[0]
        expected_y = np.exp(-tout)
        expected_x0 = np.exp(-tout)
        expected_x1 = np.exp(-2.0*tout)
        assert_almost_equal(y, expected_y)
        assert_almost_equal(x[:,0], expected_x0)
        assert_almost_equal(x[:,1], expected_x1)

    def test_06(self):
        """Test use of the default values of the arguments `T` and `U`."""
        # Second order system with a repeated root: x''(t) + 2*x(t) + x(t) = 0.
        # With initial conditions x(0)=1.0 and x'(t)=0.0, the exact solution
        # is (1-t)*exp(-t).
        system = ([1.0], [1.0, 2.0, 1.0])
        tout, y, x = lsim2(system, X0=[1.0, 0.0])
        expected_x = (1.0 - tout) * np.exp(-tout)
        assert_almost_equal(x[:,0], expected_x)

class Test_impulse2(object):

    def test_01(self):
        # First order system: x'(t) + x(t) = u(t)
        # Exact impulse response is x(t) = exp(-t).
        system = ([1.0],[1.0,1.0])
        tout, y = impulse2(system)
        expected_y = np.exp(-tout)
        assert_almost_equal(y, expected_y)

    def test_02(self):
        """Specify the desired time values for the output."""

        # First order system: x'(t) + x(t) = u(t)
        # Exact impulse response is x(t) = exp(-t).
        system = ([1.0],[1.0,1.0])
        n = 21
        t = np.linspace(0, 2.0, n)
        tout, y = impulse2(system, T=t)
        assert_equal(tout.shape, (n,))
        assert_almost_equal(tout, t)
        expected_y = np.exp(-t)
        assert_almost_equal(y, expected_y)

    def test_03(self):
        """Specify an initial condition as a scalar."""

        # First order system: x'(t) + x(t) = u(t), x(0)=3.0
        # Exact impulse response is x(t) = 4*exp(-t).
        system = ([1.0],[1.0,1.0])
        tout, y = impulse2(system, X0=3.0)
        expected_y = 4.0*np.exp(-tout)
        assert_almost_equal(y, expected_y)

    def test_04(self):
        """Specify an initial condition as a list."""

        # First order system: x'(t) + x(t) = u(t), x(0)=3.0
        # Exact impulse response is x(t) = 4*exp(-t).
        system = ([1.0],[1.0,1.0])
        tout, y = impulse2(system, X0=[3.0])
        expected_y = 4.0*np.exp(-tout)
        assert_almost_equal(y, expected_y)

    def test_05(self):
        # Simple integrator: x'(t) = u(t)
        system = ([1.0],[1.0,0.0])
        tout, y = impulse2(system)
        expected_y = np.ones_like(tout)
        assert_almost_equal(y, expected_y)

    def test_06(self):
        # Second order system with a repeated root: x''(t) + 2*x(t) + x(t) = u(t)
        # The exact impulse response is t*exp(-t).
        system = ([1.0], [1.0, 2.0, 1.0])
        tout, y = impulse2(system)
        expected_y = tout * np.exp(-tout)
        assert_almost_equal(y, expected_y)

class Test_step2(object):

    def test_01(self):
        # First order system: x'(t) + x(t) = u(t)
        # Exact step response is x(t) = 1 - exp(-t).
        system = ([1.0],[1.0,1.0])
        tout, y = step2(system)
        expected_y = 1.0 - np.exp(-tout)
        assert_almost_equal(y, expected_y)

    def test_02(self):
        """Specify the desired time values for the output."""

        # First order system: x'(t) + x(t) = u(t)
        # Exact step response is x(t) = 1 - exp(-t).
        system = ([1.0],[1.0,1.0])
        n = 21
        t = np.linspace(0, 2.0, n)
        tout, y = step2(system, T=t)
        assert_equal(tout.shape, (n,))
        assert_almost_equal(tout, t)
        expected_y = 1 - np.exp(-t)
        assert_almost_equal(y, expected_y)

    def test_03(self):
        """Specify an initial condition as a scalar."""

        # First order system: x'(t) + x(t) = u(t), x(0)=3.0
        # Exact step response is x(t) = 1 + 2*exp(-t).
        system = ([1.0],[1.0,1.0])
        tout, y = step2(system, X0=3.0)
        expected_y = 1 + 2.0*np.exp(-tout)
        assert_almost_equal(y, expected_y)

    def test_04(self):
        """Specify an initial condition as a list."""

        # First order system: x'(t) + x(t) = u(t), x(0)=3.0
        # Exact step response is x(t) = 1 + 2*exp(-t).
        system = ([1.0],[1.0,1.0])
        tout, y = step2(system, X0=[3.0])
        expected_y = 1 + 2.0*np.exp(-tout)
        assert_almost_equal(y, expected_y)

    def test_05(self):
        # Simple integrator: x'(t) = u(t)
        # Exact step response is x(t) = t.
        system = ([1.0],[1.0,0.0])
        tout, y = step2(system, atol=1e-10, rtol=1e-8)
        expected_y = tout
        assert_almost_equal(y, expected_y)

    def test_06(self):
        # Second order system with a repeated root: x''(t) + 2*x(t) + x(t) = u(t)
        # The exact step response is 1 - (1 + t)*exp(-t).
        system = ([1.0], [1.0, 2.0, 1.0])
        tout, y = step2(system, atol=1e-10, rtol=1e-8)
        expected_y = 1 - (1 + tout) * np.exp(-tout)
        assert_almost_equal(y, expected_y)

if __name__ == "__main__":
    run_module_suite()
