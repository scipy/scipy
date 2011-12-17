"""
Unit tests for optimization routines from minpack.py.
"""

from numpy.testing import assert_, assert_almost_equal, assert_array_equal, \
        assert_array_almost_equal, TestCase, run_module_suite, assert_raises
import numpy as np
from numpy import array, float64

from scipy import optimize
from scipy.optimize.minpack import leastsq, curve_fit, fixed_point


class ReturnShape(object):
    """This class exists to create a callable that does not have a 'func_name' attribute.

    __init__ takes the argument 'shape', which should be a tuple of ints.  When an instance
    it called with a single argument 'x', it returns numpy.ones(shape).
    """
    def __init__(self, shape):
        self.shape = shape

    def __call__(self, x):
        return np.ones(self.shape)

def dummy_func(x, shape):
    """A function that returns an array of ones of the given shape.
    `x` is ignored.
    """
    return np.ones(shape)


class TestFSolve(object):
    def pressure_network(self, flow_rates, Qtot, k):
        """Evaluate non-linear equation system representing
        the pressures and flows in a system of n parallel pipes::

            f_i = P_i - P_0, for i = 1..n
            f_0 = sum(Q_i) - Qtot

        Where Q_i is the flow rate in pipe i and P_i the pressure in that pipe.
        Pressure is modeled as a P=kQ**2 where k is a valve coefficient and
        Q is the flow rate.

        Parameters
        ----------
        flow_rates : float
            A 1D array of n flow rates [kg/s].
        k : float
            A 1D array of n valve coefficients [1/kg m].
        Qtot : float
            A scalar, the total input flow rate [kg/s].

        Returns
        -------
        F : float
            A 1D array, F[i] == f_i.

        """
        P = k * flow_rates**2
        F = np.hstack((P[1:] - P[0], flow_rates.sum() - Qtot))
        return F

    def pressure_network_jacobian(self, flow_rates, Qtot, k):
        """Return the jacobian of the equation system F(flow_rates)
        computed by `pressure_network` with respect to
        *flow_rates*. See `pressure_network` for the detailed
        description of parrameters.

        Returns
        -------
        jac : float
            *n* by *n* matrix ``df_i/dQ_i`` where ``n = len(flow_rates)``
            and *f_i* and *Q_i* are described in the doc for `pressure_network`
        """
        n = len(flow_rates)
        pdiff = np.diag(flow_rates[1:] * 2 * k[1:] - 2 * flow_rates[0] * k[0])

        jac = np.empty((n, n))
        jac[:n-1, :n-1] = pdiff
        jac[:n-1, n-1] = 0
        jac[n-1, :] = np.ones(n)

        return jac

    def test_pressure_network_no_gradient(self):
        """fsolve without gradient, equal pipes -> equal flows"""
        k = np.ones(4) * 0.5
        Qtot = 4
        initial_guess = array([2., 0., 2., 0.])
        final_flows = optimize.fsolve(
            self.pressure_network, initial_guess, args=(Qtot, k))
        assert_array_almost_equal(final_flows, np.ones(4))

    def test_pressure_network_with_gradient(self):
        """fsolve with gradient, equal pipes -> equal flows"""
        k = np.ones(4) * 0.5
        Qtot = 4
        initial_guess = array([2., 0., 2., 0.])
        final_flows = optimize.fsolve(
            self.pressure_network, initial_guess, args=(Qtot, k),
            fprime=self.pressure_network_jacobian)
        assert_array_almost_equal(final_flows, np.ones(4))

    def test_wrong_shape_func_callable(self):
        """The callable 'func' has no 'func_name' attribute."""
        func = ReturnShape(1)
        # x0 is a list of two elements, but func will return an array with
        # length 1, so this should result in a TypeError.
        x0 = [1.5, 2.0]
        assert_raises(TypeError, optimize.fsolve, func, x0)

    def test_wrong_shape_func_function(self):
        # x0 is a list of two elements, but func will return an array with
        # length 1, so this should result in a TypeError.
        x0 = [1.5, 2.0]
        assert_raises(TypeError, optimize.fsolve, dummy_func, x0, args=((1,),))

    def test_wrong_shape_fprime_callable(self):
        """The callables 'func' and 'deriv_func' have no 'func_name' attribute."""
        func = ReturnShape(1)
        deriv_func = ReturnShape((2,2))
        assert_raises(TypeError, optimize.fsolve, func, x0=[0,1], fprime=deriv_func)

    def test_wrong_shape_fprime_function(self):
        func = lambda x: dummy_func(x, (2,))
        deriv_func = lambda x: dummy_func(x, (3,3))
        assert_raises(TypeError, optimize.fsolve, func, x0=[0,1], fprime=deriv_func)


class TestLeastSq(TestCase):
    def setUp(self):
        x = np.linspace(0, 10, 40)
        a,b,c = 3.1, 42, -304.2
        self.x = x
        self.abc = a,b,c
        y_true = a*x**2 + b*x + c
        np.random.seed(0)
        self.y_meas = y_true + 0.01*np.random.standard_normal(y_true.shape)

    def residuals(self, p, y, x):
        a,b,c = p
        err = y-(a*x**2 + b*x + c)
        return err

    def test_basic(self):
        p0 = array([0,0,0])
        params_fit, ier = leastsq(self.residuals, p0,
                                  args=(self.y_meas, self.x))
        assert_(ier in (1,2,3,4), 'solution not found (ier=%d)'%ier)
        # low precision due to random
        assert_array_almost_equal(params_fit, self.abc, decimal=2)

    def test_full_output(self):
        p0 = array([0,0,0])
        full_output = leastsq(self.residuals, p0,
                              args=(self.y_meas, self.x),
                              full_output=True)
        params_fit, cov_x, infodict, mesg, ier = full_output
        assert_(ier in (1,2,3,4), 'solution not found: %s'%mesg)

    def test_input_untouched(self):
        p0 = array([0,0,0],dtype=float64)
        p0_copy = array(p0, copy=True)
        full_output = leastsq(self.residuals, p0,
                              args=(self.y_meas, self.x),
                              full_output=True)
        params_fit, cov_x, infodict, mesg, ier = full_output
        assert_(ier in (1,2,3,4), 'solution not found: %s'%mesg)
        assert_array_equal(p0, p0_copy)

    def test_wrong_shape_func_callable(self):
        """The callable 'func' has no 'func_name' attribute."""
        func = ReturnShape(1)
        # x0 is a list of two elements, but func will return an array with
        # length 1, so this should result in a TypeError.
        x0 = [1.5, 2.0]
        assert_raises(TypeError, optimize.leastsq, func, x0)

    def test_wrong_shape_func_function(self):
        # x0 is a list of two elements, but func will return an array with
        # length 1, so this should result in a TypeError.
        x0 = [1.5, 2.0]
        assert_raises(TypeError, optimize.leastsq, dummy_func, x0, args=((1,),))

    def test_wrong_shape_Dfun_callable(self):
        """The callables 'func' and 'deriv_func' have no 'func_name' attribute."""
        func = ReturnShape(1)
        deriv_func = ReturnShape((2,2))
        assert_raises(TypeError, optimize.leastsq, func, x0=[0,1], Dfun=deriv_func)

    def test_wrong_shape_Dfun_function(self):
        func = lambda x: dummy_func(x, (2,))
        deriv_func = lambda x: dummy_func(x, (3,3))
        assert_raises(TypeError, optimize.leastsq, func, x0=[0,1], Dfun=deriv_func)


class TestCurveFit(TestCase):
    def setUp(self):
        self.y = array([1.0, 3.2, 9.5, 13.7])
        self.x = array([1.0, 2.0, 3.0, 4.0])

    def test_one_argument(self):
        def func(x,a):
            return x**a
        popt, pcov = curve_fit(func, self.x, self.y)
        assert_(len(popt) == 1)
        assert_(pcov.shape == (1,1))
        assert_almost_equal(popt[0], 1.9149, decimal=4)
        assert_almost_equal(pcov[0,0], 0.0016, decimal=4)

        # Test if we get the same with full_output. Regression test for #1415.
        res = curve_fit(func, self.x, self.y, full_output=1)
        (popt2, pcov2, infodict, errmsg, ier) = res
        assert_array_almost_equal(popt, popt2)

    def test_two_argument(self):
        def func(x, a, b):
            return b*x**a
        popt, pcov = curve_fit(func, self.x, self.y)
        assert_(len(popt) == 2)
        assert_(pcov.shape == (2,2))
        assert_array_almost_equal(popt, [1.7989, 1.1642], decimal=4)
        assert_array_almost_equal(pcov, [[0.0852, -0.1260],[-0.1260, 0.1912]],
                                  decimal=4)

    def test_func_is_classmethod(self):
        class test_self(object):
            """This class tests if curve_fit passes the correct number of
               arguments when the model function is a class instance method.
            """
            def func(self, x, a, b):
                return b * x**a

        test_self_inst = test_self()
        popt, pcov = curve_fit(test_self_inst.func, self.x, self.y)
        assert_(pcov.shape == (2,2))
        assert_array_almost_equal(popt, [1.7989, 1.1642], decimal=4)
        assert_array_almost_equal(pcov, [[0.0852, -0.1260], [-0.1260, 0.1912]],
                                  decimal=4)


class TestFixedPoint(TestCase):

    def test_scalar_trivial(self):
        """f(x) = 2x; fixed point should be x=0"""
        def func(x):
            return 2.0*x
        x0 = 1.0
        x = fixed_point(func, x0)
        assert_almost_equal(x, 0.0)

    def test_scalar_basic1(self):
        """f(x) = x**2; x0=1.05; fixed point should be x=1"""
        def func(x):
            return x**2
        x0 = 1.05
        x = fixed_point(func, x0)
        assert_almost_equal(x, 1.0)

    def test_scalar_basic2(self):
        """f(x) = x**0.5; x0=1.05; fixed point should be x=1"""
        def func(x):
            return x**0.5
        x0 = 1.05
        x = fixed_point(func, x0)
        assert_almost_equal(x, 1.0)

    def test_array_trivial(self):
        def func(x):
            return 2.0*x
        x0 = [0.3, 0.15]
        olderr = np.seterr(all='ignore')
        try:
            x = fixed_point(func, x0)
        finally:
            np.seterr(**olderr)
        assert_almost_equal(x, [0.0, 0.0])

    def test_array_basic1(self):
        """f(x) = c * x**2; fixed point should be x=1/c"""
        def func(x, c):
            return c * x**2
        c = array([0.75, 1.0, 1.25])
        x0 = [1.1, 1.15, 0.9]
        olderr = np.seterr(all='ignore')
        try:
            x = fixed_point(func, x0, args=(c,))
        finally:
            np.seterr(**olderr)
        assert_almost_equal(x, 1.0/c)

    def test_array_basic2(self):
        """f(x) = c * x**0.5; fixed point should be x=c**2"""
        def func(x, c):
            return c * x**0.5
        c = array([0.75, 1.0, 1.25])
        x0 = [0.8, 1.1, 1.1]
        x = fixed_point(func, x0, args=(c,))
        assert_almost_equal(x, c**2)


if __name__ == "__main__":
    run_module_suite()
