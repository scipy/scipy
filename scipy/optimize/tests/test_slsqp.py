from numpy.testing import assert_array_almost_equal, TestCase, run_module_suite
import numpy as np

from scipy.optimize import  fmin_slsqp


class TestSLSQP(TestCase):
    """Test fmin_slsqp using Example 14.4 from Numerical Methods for
    Engineers by Steven Chapra and Raymond Canale.  This example
    maximizes the function f(x) = 2*x*y + 2*x - x**2 - 2*y**2, which
    has a maximum at x=2,y=1.

    """

    def _testfunc(self,d,*args):
        """
        Arguments:
        d     - A list of two elements, where d[0] represents x and d[1] represents y
                 in the following equation.
        sign - A multiplier for f.  Since we want to optimize it, and the scipy
               optimizers can only minimize functions, we need to multiply it by
               -1 to achieve the desired solution
        Returns:
        2*x*y + 2*x - x**2 - 2*y**2

        """
        try:
            sign = args[0]
        except:
            sign = 1.0
        x = d[0]
        y = d[1]
        return sign*(2*x*y + 2*x - x**2 - 2*y**2)

    def _testfunc_deriv(self,d,*args):
        """
        This is the derivative of testfunc, returning a numpy array
        representing df/dx and df/dy.

        """
        try:
            sign = args[0]
        except:
            sign = 1.0
        x = d[0]
        y = d[1]
        dfdx = sign*(-2*x + 2*y + 2)
        dfdy = sign*(2*x - 4*y)
        return np.array([ dfdx, dfdy ],float)

    def test_unbounded_approximated(self):
        res =  fmin_slsqp(self._testfunc, [-1.0,1.0], args = (-1.0,),
                          iprint = 0, full_output = 1)
        x,fx,its,imode,smode = res
        assert_array_almost_equal(x,[2,1])

    def test_unbounded_given(self):
        res = fmin_slsqp(self._testfunc,[-1.0,1.0], args = (-1.0,),
                       iprint = 0, full_output = 1)
        x,fx,its,imode,smode = res
        assert_array_almost_equal(x,[2,1])

    def test_bound_approximated(self):
        res = fmin_slsqp(self._testfunc,[-1.0,1.0], args = (-1.0,),
                         eqcons = [lambda x, y: x[0]-x[1] ],
                         iprint = 0, full_output = 1)
        x,fx,its,imode,smode = res
        assert_array_almost_equal(x,[1,1])

    def test_bound_equality_given(self):
        res = fmin_slsqp(self._testfunc,[-1.0,1.0],
                         fprime = self._testfunc_deriv,
                         args = (-1.0,), eqcons = [lambda x, y: x[0]-x[1] ],
                         iprint = 0, full_output = 1)
        x,fx,its,imode,smode = res
        assert_array_almost_equal(x,[1,1])

    def test_bound_equality_inequality_given(self):
        res = fmin_slsqp(self._testfunc,[-1.0,1.0],
                         fprime = self._testfunc_deriv,
                         args = (-1.0,),
                         ieqcons = [lambda x, y: x[0]-x[1]-1.0],
                         iprint=0, full_output=1)
        x,fx,its,imode,smode = res
        assert_array_almost_equal(x,[2,1],decimal=3)

if __name__ == "__main__":
    run_module_suite()
