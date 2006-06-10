
import numpy
from numpy import cos, sin, pi
from numpy.testing import *

set_package_path()
from scipy.integrate import quadrature, romberg
restore_path()

def assert_quad((value, err), tabledValue, errTol=1.5e-8):
    assert abs(value-tabledValue) < err, (value, tabledValue, err)
    if errTol is not None:
        assert err < errTol, (err, errTol)

class test_quadrature(ScipyTestCase):
    def quad(self, x, a, b, args):
        raise NotImplementedError

    def check_quadrature(self):
        # Typical function with two extra arguments:
        def myfunc(x,n,z):       # Bessel function integrand
            return cos(n*x-z*sin(x))/pi
        assert_quad(quadrature(myfunc,0,pi,(2,1.8)), 0.30614353532540296487)

    def check_romberg(self):
        # Typical function with two extra arguments:
        def myfunc(x, n, z):       # Bessel function integrand
            return cos(n*x-z*sin(x))/pi
        assert_quad(romberg(myfunc,0,pi, args=(2, 1.8)), 0.30614353532540296487)

if __name__ == "__main__":
    ScipyTest().run()
