import math

from numpy.testing import *

from scipy.optimize import cobyla as co

class TestCobyla(TestCase):
    def test_simple(self):

        function = lambda x: x[0]**2 + abs(x[1])**3
        con1 = lambda x: x[0]**2 + x[1]**2 - 25
        con2 = lambda x: -con1(x)

        x = co.fmin_cobyla(function, [4.95,0.66], [con1, con2], rhobeg=1,
                           rhoend=1e-5, iprint=0, maxfun=100)
        x1 = 2.0/3
        x0 = math.sqrt(25-x1*x1)
        assert_almost_equal(x, [x0, x1], decimal=5)

if __name__ == "__main__":
    run_module_suite()
