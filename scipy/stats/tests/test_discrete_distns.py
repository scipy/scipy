from scipy.stats import (betabinom, hypergeom, nhypergeom, bernoulli,
                         boltzmann, skellam, fnch, randint)

import numpy as np
from numpy.testing import assert_almost_equal, assert_equal, assert_allclose
from scipy.special import binom


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


def test_nhypergeom_pmf():
    # test with hypergeom
    M, n, r = 45, 13, 8
    k = 6
    NHG = nhypergeom.pmf(k, M, n, r)
    HG = hypergeom.pmf(k, M, n, k+r-1) * (M - n - (r-1)) / (M - (k+r-1))
    assert_allclose(HG, NHG, rtol=1e-10)


def test_nhypergeom_pmfcdf():
    # test pmf and cdf with arbitrary values.
    M = 8
    n = 3
    r = 4
    support = np.arange(n+1)
    pmf = nhypergeom.pmf(support, M, n, r)
    cdf = nhypergeom.cdf(support, M, n, r)
    assert_allclose(pmf, [1/14, 3/14, 5/14, 5/14], rtol=1e-13)
    assert_allclose(cdf, [1/14, 4/14, 9/14, 1.0], rtol=1e-13)


def test_nhypergeom_r0():
    # test with `r = 0`.
    M = 10
    n = 3
    r = 0
    pmf = nhypergeom.pmf([[0, 1, 2, 0], [1, 2, 0, 3]], M, n, r)
    assert_allclose(pmf, [[1, 0, 0, 1], [0, 0, 1, 0]], rtol=1e-13)


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


def test_betabinom_a_and_b_unity():
    # test limiting case that betabinom(n, 1, 1) is a discrete uniform
    # distribution from 0 to n
    n = 20
    k = np.arange(n + 1)
    p = betabinom(n, 1, 1).pmf(k)
    expected = np.repeat(1 / (n + 1), n + 1)
    assert_almost_equal(p, expected)


def test_betabinom_bernoulli():
    # test limiting case that betabinom(1, a, b) = bernoulli(a / (a + b))
    a = 2.3
    b = 0.63
    k = np.arange(2)
    p = betabinom(1, a, b).pmf(k)
    expected = bernoulli(a / (a + b)).pmf(k)
    assert_almost_equal(p, expected)


def test_skellam_gh11474():
    # test issue reported in gh-11474 caused by `cdfchn`
    mu = [1, 10, 100, 1000, 5000, 5050, 5100, 5250, 6000]
    cdf = skellam.cdf(0, mu, mu)
    # generated in R
    # library(skellam)
    # options(digits = 16)
    # mu = c(1, 10, 100, 1000, 5000, 5050, 5100, 5250, 6000)
    # pskellam(0, mu, mu, TRUE)
    cdf_expected = [0.6542541612768356, 0.5448901559424127, 0.5141135799745580,
                    0.5044605891382528, 0.5019947363350450, 0.5019848365953181,
                    0.5019750827993392, 0.5019466621805060, 0.5018209330219539]
    assert_allclose(cdf, cdf_expected)


class TestFNCH():
    np.random.seed(2)  # seeds 0 and 1 had some xl = xu; randint failed
    shape = (2, 4, 3)
    max_m = 100
    m1 = np.random.randint(1, max_m, size=shape)    # red balls
    m2 = np.random.randint(1, max_m, size=shape)    # white balls
    N = m1 + m2                                     # total balls
    n = randint.rvs(0, N, size=N.shape)             # number of draws
    xl = np.maximum(0, n-m2)                        # lower bound of support
    xu = np.minimum(n, m1)                          # upper bound of support
    x = randint.rvs(xl, xu, size=xl.shape)

    def test_fnch_hypergeom(self):
        # Fisher's noncentral hypergeometric distribution reduces to the
        # hypergeometric distribution when odds = 1
        x, N, m1, n = self.x, self.N, self.m1, self.n
        assert_allclose(fnch.pmf(x, N, m1, n, odds=1),
                        hypergeom.pmf(x, N, m1, n))

    def test_fnch_naive(self):
        # test against a very simple implementation
        x, N, m1, n = self.x, self.N, self.m1, self.n
        odds = np.random.rand(*x.shape)*2

        @np.vectorize
        def pmf(x, N, m1, n, w):
            # simple implementation of FNCH pmf
            m2 = N - m1
            xl = np.maximum(0, n-m2)
            xu = np.minimum(n, m1)

            def f(x):
                t1 = binom(m1, x)
                t2 = binom(m2, n - x)
                return t1 * t2 * w**x

            P0 = sum((f(y) for y in range(xl, xu + 1)))
            return f(x) / P0

        assert_allclose(fnch.pmf(x, N, m1, n, odds), pmf(x, N, m1, n, odds))
