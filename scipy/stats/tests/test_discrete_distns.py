import pytest
from scipy.stats import (betabinom, hypergeom, nhypergeom, bernoulli,
                         boltzmann, skellam, binom, nbinom)

import numpy as np
from numpy.testing import assert_almost_equal, assert_equal, assert_allclose


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


def test_issue_10317():
    alpha, n, p = 0.9, 10, 1
    assert_equal(nbinom.interval(alpha=alpha, n=n, p=p), (0, 0))


def test_issue_11134():
    alpha, n, p = 0.95, 10, 0
    assert_equal(binom.interval(alpha=alpha, n=n, p=p), (0, 0))


def test_issue_7406():
    np.random.seed(0)
    assert_equal(binom.ppf(np.random.rand(10), 0, 0.5), 0)

    # Also check that endpoints (q=0, q=1) are correct
    assert_equal(binom.ppf(0, 0, 0.5), -1)
    assert_equal(binom.ppf(1, 0, 0.5), 0)


def test_issue_5122():
    p = 0
    n = np.random.randint(100, size=10)

    x = 0
    ppf = binom.ppf(x, n, p)
    assert_equal(ppf, -1)

    x = np.linspace(0.01, 0.99, 10)
    ppf = binom.ppf(x, n, p)
    assert_equal(ppf, 0)

    x = 1
    ppf = binom.ppf(x, n, p)
    assert_equal(ppf, n)


def test_issue_1603():
    assert_equal(binom(1000, np.logspace(-3, -100)).ppf(0.01), 0)


def test_issue_5503():
    p = 0.5
    x = np.logspace(3, 14, 12)
    assert_allclose(binom.cdf(x, 2*x, p), 0.5, atol=1e-2)


@pytest.mark.parametrize('x, n, p, cdf_desired', [
    (300, 1000, 3/10, 0.51559351981411995636),
    (3000, 10000, 3/10, 0.50493298381929698016),
    (30000, 100000, 3/10, 0.50156000591726422864),
    (300000, 1000000, 3/10, 0.50049331906666960038),
    (3000000, 10000000, 3/10, 0.50015600124585261196),
    (30000000, 100000000, 3/10, 0.50004933192735230102),
    (30010000, 100000000, 3/10, 0.98545384016570790717),
    (29990000, 100000000, 3/10, 0.01455017177985268670),
    (29950000, 100000000, 3/10, 5.02250963487432024943e-28),
])
def test_issue_5503pt2(x, n, p, cdf_desired):
    assert_allclose(binom.cdf(x, n, p), cdf_desired)


def test_issue_5503pt3():
    # From Wolfram Alpha: CDF[BinomialDistribution[1e12, 1e-12], 2]
    assert_allclose(binom.cdf(2, 10**12, 10**-12), 0.91969860292869777384)


def test_issue_6682():
    # Reference value from R:
    # options(digits=16)
    # print(pnbinom(250, 50, 32/63, lower.tail=FALSE))
    assert_allclose(nbinom.sf(250, 50, 32./63.), 1.460458510976452e-35)


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
