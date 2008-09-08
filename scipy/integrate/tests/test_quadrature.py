
import numpy
from numpy import cos, sin, pi
from numpy.testing import *

from scipy.integrate import quadrature, romberg, romb

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

if __name__ == "__main__":
    run_module_suite()
