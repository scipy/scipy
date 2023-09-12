import numpy as np
from numpy.testing import (assert_array_equal, assert_equal,
        assert_array_almost_equal_nulp, assert_almost_equal, assert_allclose)
from pytest import raises as assert_raises
import pytest
from scipy.special import gammaln, multigammaln, betaln, multivariate_betaln


class TestMultiGammaLn:

    def test1(self):
        # A test of the identity
        #     Gamma_1(a) = Gamma(a)
        np.random.seed(1234)
        a = np.abs(np.random.randn())
        assert_array_equal(multigammaln(a, 1), gammaln(a))

    def test2(self):
        # A test of the identity
        #     Gamma_2(a) = sqrt(pi) * Gamma(a) * Gamma(a - 0.5)
        a = np.array([2.5, 10.0])
        result = multigammaln(a, 2)
        expected = np.log(np.sqrt(np.pi)) + gammaln(a) + gammaln(a - 0.5)
        assert_almost_equal(result, expected)

    def test_bararg(self):
        assert_raises(ValueError, multigammaln, 0.5, 1.2)


def _check_multigammaln_array_result(a, d):
    # Test that the shape of the array returned by multigammaln
    # matches the input shape, and that all the values match
    # the value computed when multigammaln is called with a scalar.
    result = multigammaln(a, d)
    assert_array_equal(a.shape, result.shape)
    a1 = a.ravel()
    result1 = result.ravel()
    for i in range(a.size):
        assert_array_almost_equal_nulp(result1[i], multigammaln(a1[i], d))


def test_multigammaln_array_arg():
    # Check that the array returned by multigammaln has the correct
    # shape and contains the correct values.  The cases have arrays
    # with several differnent shapes.
    # The cases include a regression test for ticket #1849
    # (a = np.array([2.0]), an array with a single element).
    np.random.seed(1234)

    cases = [
        # a, d
        (np.abs(np.random.randn(3, 2)) + 5, 5),
        (np.abs(np.random.randn(1, 2)) + 5, 5),
        (np.arange(10.0, 18.0).reshape(2, 2, 2), 3),
        (np.array([2.0]), 3),
        (np.float64(2.0), 3),
    ]

    for a, d in cases:
        _check_multigammaln_array_result(a, d)

class TestMultivariateBetaln:
    def test_1d(self):
        # test that multivariate beta and one dimensional beta functions yield
        # same result in 1D
        a = 1
        b = 2
        assert_allclose(betaln(a, b), multivariate_betaln([a, b]), rtol=1e-16)

    # reference values were computed via mpmath
    # from mpmath import mp
    # mp.dps = 500
    # def multibetaln(x):
    #     x = [mp.mpf(coeff) for coeff in x]
    #     first_term = mp.fsum([mp.loggamma(coeff) for coeff in x])
    #     second_term = mp.loggamma(mp.fsum(x))
    #     return float(first_term - second_term)

    @pytest.mark.parametrize('alpha, ref',
                             [(np.logspace(-100, 100, 1000),
                               -4.828190914524598e+100),
                              (np.logspace(-300, -100, 100),
                               45821.45294191622),
                              (np.logspace(150, 300, 200),
                               -6.864052935635053e+299)])
    def test_accuracy(self, alpha, ref):
        assert_allclose(multivariate_betaln(alpha), ref, rtol=1e-13)


