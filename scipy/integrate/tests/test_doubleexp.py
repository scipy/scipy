from numpy import exp, sin, cos
from numpy.testing import *
from scipy.integrate import de_integrate

def assert_de(value, expected, errTol=1e-10):
    assert abs(value - expected) < errTol, (value, expected, errTol)

class TestDE(TestCase):
    def test_typical(self):
        def f(x):
            return exp(-x/5.0) * (2 + sin(2 * x))
        def symbolic(x):
            return -5.0 / 101.0 * exp(-x / 5.0) * (10 * cos(2 * x) + sin(2 * x) + 202)
        assert_de(de_integrate(f, 0, 10, target_error=1e-10), symbolic(10) - symbolic(0), errTol=1e-10)

    def test_singularity(self):
        def f(x):
            def _pow(a,b):
                if a==0.0:
                    return 0.0
                return a**b
            return _pow(1-x,5) * _pow(x,-1.0/3)
        assert_de(de_integrate(f, 0, 1, target_error=1e-10), 2187.0/5236, errTol=1e-10)

if __name__ == "__main__":
    run_module_suite()
