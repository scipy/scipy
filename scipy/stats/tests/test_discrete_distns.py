from __future__ import division, print_function, absolute_import

from scipy.stats import hypergeom, bernoulli, boltzmann, yulesimon
import numpy as np
from numpy.testing import assert_almost_equal, assert_equal, assert_allclose, assert_array_equal


def test_hypergeom_logpmf():
    # symmetries test
    # f(k,N,K,n) = f(n-k,N,N-K,n) = f(K-k,N,K,N-n) = f(k,N,n,K)
    k = 5
    N = 50
    K = 10
    n = 5
    logpmf1 = hypergeom.logpmf(k, N, K, n)
    logpmf2 = hypergeom.logpmf(n - k, N, N - K, n)
    logpmf3 = hypergeom.logpmf(K - k, N, K, N - n)
    logpmf4 = hypergeom.logpmf(k, N, n, K)
    assert_almost_equal(logpmf1, logpmf2, decimal=12)
    assert_almost_equal(logpmf1, logpmf3, decimal=12)
    assert_almost_equal(logpmf1, logpmf4, decimal=12)

    # test related distribution
    # Bernoulli distribution if n = 1
    k = 1
    N = 10
    K = 7
    n = 1
    hypergeom_logpmf = hypergeom.logpmf(k, N, K, n)
    bernoulli_logpmf = bernoulli.logpmf(k, K/N)
    assert_almost_equal(hypergeom_logpmf, bernoulli_logpmf, decimal=12)


def test_boltzmann_upper_bound():
    k = np.arange(-3, 5)

    N = 1
    p = boltzmann.pmf(k, 0.123, N)
    expected = k == 0
    assert_equal(p, expected)

    lam = np.log(2)
    N = 3
    p = boltzmann.pmf(k, lam, N)
    expected = [0, 0, 0, 4/7, 2/7, 1/7, 0, 0]
    assert_allclose(p, expected, rtol=1e-13)

    c = boltzmann.cdf(k, lam, N)
    expected = [0, 0, 0, 4/7, 6/7, 1, 1, 1]
    assert_allclose(c, expected, rtol=1e-13)

def test_yule_simon_pmf():
    # test pmf implementation
    alpha = 5
    rv = yulesimon(alpha)
    ys_pmf = rv.pmf([0, 1, 2, 3, 4])
    expected = [0.0, 0.83333333, 0.11904762, 0.0297619, 0.00992063]
    assert_allclose(ys_pmf, expected, rtol=1e8)
    # test if negative parameter results is nan pmf
    assert_array_equal(yulesimon(-1).pmf(1), np.nan)
    # test log_pmf function
    assert_almost_equal(yulesimon(1).logpmf(1), -0.6931471805599453, decimal=12)
    # test cdf implementation
    assert_almost_equal(yulesimon(1).cdf(1), 0.5, decimal=12)
    # test sf implementation
    assert_almost_equal(yulesimon(1).sf(1), 0.5, decimal=12)
    # test logsf implementation
    assert_almost_equal(yulesimon(1).logsf(1), -0.6931471805599453, decimal=12)
    #test stats implementation
    assert_allclose(yulesimon(11).stats('m'), (1.0999999999),rtol=1e-8)
