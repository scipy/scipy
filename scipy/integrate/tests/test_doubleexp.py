import numpy as np
from numpy import exp, sin, cos, log, sqrt, pi
from numpy.testing import *
from scipy.integrate import quad_de, quad
import scipy.integrate.doubleexp

def assert_de(value, expected, tol=1e-10):
    assert abs(value - expected) < tol, (value, expected, tol)

class TestDE(TestCase):
    def test_typical(self):
        def f(x):
            return exp(-x/5.0) * (2 + sin(2*x))
        def symbolic(x):
            return -5.0/101.0 * exp(-x/5.0) * (10*cos(2*x) + sin(2*x) + 202)
        assert_de(quad_de(f, 0, 10, tol=1e-10)[0], symbolic(10) - symbolic(0),
                  tol=1e-10)

    def test_singularity(self):
        def f(x):
            return (1-x)**5 * x**(-1.0/3)
        assert_de(quad_de(f, 0, 1, tol=1e-10)[0], 2187.0/5236, tol=1e-10)

    def test_constant(self):
        def f(x):
            return 2
        assert_de(quad_de(f, 0, 1, tol=1e-10)[0], 2, tol=1e-10)

    def _check_convergence(self, f, a, b, error_formula, info="", plot=False,
                           exact=None):
        num_levels = len(scipy.integrate.doubleexp._abscissas)

        res = []
        N = []
        for level in xrange(num_levels):
            res.append(quad_de(f, a, b, max_level=level, tol=1e-99)[0])
            N.append(len(scipy.integrate.doubleexp._abscissas[level]))
        res = np.array(res)
        N = np.array(N)
        if exact is None:
            exact = res[-1]
        err = abs(res - exact)

        formula, err_a, err_b = error_formula
        expected = eval(formula)

        # XXX: remove this plotting from the final version
        if plot:
            import matplotlib.pyplot as plt
            plt.semilogy(N, expected, '-', N, err, '.')
            plt.legend(('expected error from [TM]', 'actual error'))
            plt.xlim(0, err_b)
            plt.ylim(1e-15, 1)
            plt.title(info)
            plt.show()

        # check convergence
        j = (N >= err_a)
        assert np.all(err[j] <= expected[j])

    def test_convergence_fig6(self):
        # Compare convergence against Fig. 6 in [TM];
        # the formula below is a line fitted to the figure, in given N interval
        dI_fitted = ("10*exp(-2.95*N/log(N))", 15, 68)
        def f1(x):
            return 1.0/((x-2) * (1-x)**(.25) * (1+x)**(.75))
        self._check_convergence(f1, -1, 1, dI_fitted, info="Fig 6", plot=True,
                                exact=-sqrt(2)*pi*3**(-3./4))

    def test_convergence_fig7(self):
        # Compare convergence against Fig. 7 in [TM];
        dI_fitted = ("5e4*exp(-3.4*N/log(N))", 15, 70)
        def f1(x):
            return cos(pi*x)/sqrt(1-x)
        self._check_convergence(f1, -1, 1, dI_fitted, info="Fig 7", plot=True,
                                exact=-0.69049458874660501715279861110)


if __name__ == "__main__":
    run_module_suite()
