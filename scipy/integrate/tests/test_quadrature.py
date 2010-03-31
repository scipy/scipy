
import numpy
from numpy import cos, sin, pi
from numpy.testing import *

from scipy.integrate import quadrature, romberg, romb, newton_cotes

class TestQuadrature(TestCase):
    def quad(self, x, a, b, args):
        raise NotImplementedError

    def test_quadrature(self):
        # Typical function with two extra arguments:
        def myfunc(x,n,z):       # Bessel function integrand
            return cos(n*x-z*sin(x))/pi
        val, err = quadrature(myfunc,0,pi,(2,1.8))
        table_val = 0.30614353532540296487
        assert_almost_equal(val, table_val, decimal=7)

    def test_romberg(self):
        # Typical function with two extra arguments:
        def myfunc(x, n, z):       # Bessel function integrand
            return cos(n*x-z*sin(x))/pi
        val = romberg(myfunc,0,pi, args=(2, 1.8))
        table_val = 0.30614353532540296487
        assert_almost_equal(val, table_val, decimal=7)

    def test_romb(self):
        assert_equal(romb(numpy.arange(17)),128)

    def test_non_dtype(self):
        # Check that we work fine with functions returning float
        import math
        valmath = romberg(math.sin, 0, 1)
        expected_val = 0.45969769413185085
        assert_almost_equal(valmath, expected_val, decimal=7)

    def test_newton_cotes(self):
        """Test the first few degrees, for evenly spaced points."""
        n = 1
        wts, errcoff = newton_cotes(n, 1)
        assert_equal(wts, n*numpy.array([0.5, 0.5]))
        assert_almost_equal(errcoff, -n**3/12.0)

        n = 2
        wts, errcoff = newton_cotes(n, 1)
        assert_almost_equal(wts, n*numpy.array([1.0, 4.0, 1.0])/6.0)
        assert_almost_equal(errcoff, -n**5/2880.0)

        n = 3
        wts, errcoff = newton_cotes(n, 1)
        assert_almost_equal(wts, n*numpy.array([1.0, 3.0, 3.0, 1.0])/8.0)
        assert_almost_equal(errcoff, -n**5/6480.0)

        n = 4
        wts, errcoff = newton_cotes(n, 1)
        assert_almost_equal(wts, n*numpy.array([7.0, 32.0, 12.0, 32.0, 7.0])/90.0)
        assert_almost_equal(errcoff, -n**7/1935360.0)

    def test_newton_cotes2(self):
        """Test newton_cotes with points that are not evenly spaced."""
        
        x = numpy.array([0.0, 1.5, 2.0])
        y = x**2
        wts, errcoff = newton_cotes(x)
        exact_integral = 8.0/3
        numeric_integral = numpy.dot(wts, y)
        assert_almost_equal(numeric_integral, exact_integral)
        
        x = numpy.array([0.0, 1.4, 2.1, 3.0])
        y = x**2
        wts, errcoff = newton_cotes(x)
        exact_integral = 9.0
        numeric_integral = numpy.dot(wts, y)
        assert_almost_equal(numeric_integral, exact_integral)        


if __name__ == "__main__":
    run_module_suite()
