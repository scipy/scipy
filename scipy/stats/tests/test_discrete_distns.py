from scipy.stats import (betabinom, hypergeom, nhypergeom, bernoulli,
                         boltzmann, skellam, zipf, zipfian, nbinom)

import numpy as np
from numpy.testing import assert_almost_equal, assert_equal, assert_allclose
import pytest


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


class TestZipfian:
    def test_zipfian_asymptotic(self):
        # test limiting case that zipfian(a, n) -> zipf(a) as n-> oo
        a = 6.5
        N = 10000000
        k = np.arange(1, 21)
        assert_allclose(zipfian.pmf(k, a, N), zipf.pmf(k, a))
        assert_allclose(zipfian.cdf(k, a, N), zipf.cdf(k, a))
        assert_allclose(zipfian.sf(k, a, N), zipf.sf(k, a))
        assert_allclose(zipfian.stats(a, N, moments='msvk'),
                        zipf.stats(a, moments='msvk'))

    def test_zipfian_continuity(self):
        # test that zipfian(0.999999, n) ~ zipfian(1.000001, n)
        # (a = 1 switches between methods of calculating harmonic sum)
        alt1, agt1 = 0.99999999, 1.00000001
        N = 30
        k = np.arange(1, N  + 1)
        assert_allclose(zipfian.pmf(k, alt1, N), zipfian.pmf(k, agt1, N))
        assert_allclose(zipfian.cdf(k, alt1, N), zipfian.cdf(k, agt1, N))
        assert_allclose(zipfian.sf(k, alt1, N), zipfian.sf(k, agt1, N))
        assert_allclose(zipfian.stats(alt1, N, moments='msvk'),
                        zipfian.stats(agt1, N, moments='msvk'), rtol=2e-7)

    def test_zipfian_R(self):
        # test against R VGAM package
        # library(VGAM)
        # k <- c(13, 16,  1,  4,  4,  8, 10, 19,  5,  7)
        # a <- c(1.56712977, 3.72656295, 5.77665117, 9.12168729, 5.79977172,
        #        4.92784796, 9.36078764, 4.3739616 , 7.48171872, 4.6824154)
        # n <- c(70, 80, 48, 65, 83, 89, 50, 30, 20, 20)
        # pmf <- dzipf(k, N = n, shape = a)
        # cdf <- pzipf(k, N = n, shape = a)
        # print(pmf)
        # print(cdf)
        np.random.seed(0)
        k = np.random.randint(1, 20, size=10)
        a = np.random.rand(10)*10 + 1
        n = np.random.randint(1, 100, size=10)
        pmf = [8.076972e-03, 2.950214e-05, 9.799333e-01, 3.216601e-06,
               3.158895e-04, 3.412497e-05, 4.350472e-10, 2.405773e-06,
               5.860662e-06, 1.053948e-04]
        cdf = [0.8964133, 0.9998666, 0.9799333, 0.9999995, 0.9998584,
               0.9999458, 1.0000000, 0.9999920, 0.9999977, 0.9998498]
        # skip the first point; zipUC is not accurate for low a, n
        assert_allclose(zipfian.pmf(k, a, n)[1:], pmf[1:], rtol=1e-6)
        assert_allclose(zipfian.cdf(k, a, n)[1:], cdf[1:], rtol=5e-5)

    np.random.seed(0)
    naive_tests = np.vstack((np.logspace(-2, 1, 10),
                             np.random.randint(2, 40, 10))).T
    @pytest.mark.parametrize("a, n", naive_tests)
    def test_zipfian_naive(self, a, n):
        # test against bare-bones implementation

        @np.vectorize
        def Hns(n, s):
            """Naive implementation of harmonic sum"""
            return (1/np.arange(1, n+1)**s).sum()

        @np.vectorize
        def pzip(k, a, n):
            """Naive implementation of zipfian pmf"""
            if k < 1 or k > n:
                return 0.
            else:
                return 1 / k**a / Hns(n, a)

        k = np.arange(n+1)
        pmf = pzip(k, a, n)
        cdf = np.cumsum(pmf)
        mean = np.average(k, weights=pmf)
        var = np.average((k - mean)**2, weights=pmf)
        std = var**0.5
        skew = np.average(((k-mean)/std)**3, weights=pmf)
        kurtosis = np.average(((k-mean)/std)**4, weights=pmf) - 3
        assert_allclose(zipfian.pmf(k, a, n), pmf)
        assert_allclose(zipfian.cdf(k, a, n), cdf)
        assert_allclose(zipfian.stats(a, n, moments="mvsk"),
                        [mean, var, skew, kurtosis])


@pytest.mark.parametrize("mu, q, expected",
                         [[10, 120, -1.240089881791596e-38],
                          [1500, 0, -86.61466680572661]])
def test_nbinom_11465(mu, q, expected):
    # test nbinom.logcdf at extreme tails
    size = 20
    n, p = size, size/(size+mu)
    # In R:
    # options(digits=16)
    # pnbinom(mu=10, size=20, q=120, log.p=TRUE)
    assert_allclose(nbinom.logcdf(q, n, p), expected)
