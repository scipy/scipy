import numpy as np
from numpy import exp, sin, cos, log, sqrt, pi
from numpy.testing import *
from scipy.integrate import quad_de, quad
from scipy.integrate import _doubleexp

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

    def test_vs_quad(self):
        cases = [
            (lambda x: exp(-x), 0, np.inf),
            (lambda x: exp(-x), 1.23, np.inf),
            (lambda x: exp(-x), np.inf, 1.23),
            (lambda x: exp(-x**2), -np.inf, np.inf),
            (lambda x: exp(-x**2), np.inf, -np.inf),
        ]

        for j, (f, a, b) in enumerate(cases):
            c1 = quad_de(f, a, b, tol=1e-10)[0]
            if a < b:
                c2 = quad(f, a, b, epsabs=1e-12)[0]
            else:
                c2 = -quad(f, b, a, epsabs=1e-12)[0]

            assert_allclose(c1, c2, atol=1e-10, rtol=0, err_msg=repr(j))


    # Check that the convergence really is double exponential.
    # Ideally, we'd like to compare to the figures in [TM], but this
    # doesn't work so well --- in [TM], the sampling step h is tuned
    # for the best possible result for each N. We cannot do this in
    # quad_de, as the abscissas are precomputed, and an adaptive
    # scheme is used.

    def _check_convergence(self, f, a, b, error_formula, info="", plot=False,
                           tol=0, exact=None):

        if a == np.inf and b == np.inf:
            abscissas = _doubleexp._abscissas_minf_inf
        elif b == np.inf:
            abscissas = _doubleexp._abscissas_0_inf
        else:
            abscissas = _doubleexp._abscissas_m1_1

        num_levels = len(abscissas)

        res = []
        N = [0]
        for level in xrange(num_levels):
            res.append(quad_de(f, a, b, max_level=level, tol=1e-99)[0])
            N.append(N[-1] + len(abscissas[level]))

        res = np.array(res)
        N = np.array(N[1:])
        if exact is None:
            exact = res[-1]
        err = abs(res - exact)

        formula, err_a, err_b = error_formula
        expected = eval(formula)

        if plot:
            import matplotlib.pyplot as plt
            plt.semilogy(N, expected, '-', N, err, '.')
            plt.legend(('expected error', 'actual error'))
            plt.xlim(0, max(N))
            plt.ylim(1e-15, 1)
            plt.title(info)
            plt.show()

        # check convergence
        j = (N >= err_a)
        assert np.all(err[j] <= tol + expected[j]), (N[j], err[j],
                                                     tol+expected[j])

    def test_double_exp_convergence_m1_1(self):
        dI_fitted = ("1e4*exp(-1.3*N/log(N))", 15, 70)
        def f1(x):
            return 1/(1 + x**2) * cos(pi*x)
        self._check_convergence(f1, -1, 1, dI_fitted, info="-1,1",
                                tol=1e-15,
                                exact=0.226109298772811084590)

    def test_double_exp_convergence_0_inf(self):
        dI_fitted = ("4*exp(-0.3*N/log(N))", 15, 70)
        def f1(x):
            return exp(-x) * cos(x)
        self._check_convergence(f1, 0, np.inf, dI_fitted, info="0,inf",
                                tol=1e-15,
                                exact=0.5)

    def test_double_exp_convergence_minf_inf(self):
        dI_fitted = ("1e3*exp(-0.9*N/log(N))", 15, 70)
        def f1(x):
            return 1/(1 + x**4)
        self._check_convergence(f1, -np.inf, np.inf, dI_fitted, info="-inf,inf",
                                tol=1e-15,
                                exact=pi/sqrt(2))


if __name__ == "__main__":
    run_module_suite()
