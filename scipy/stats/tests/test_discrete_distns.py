from scipy.stats import (betabinom, hypergeom, nhypergeom, bernoulli,
                         boltzmann, skellam, fnch, wnch, randint, nbinom)

import numpy as np
from numpy.testing import assert_almost_equal, assert_equal, assert_allclose
from scipy.special import binom
import pytest
from scipy.optimize import root_scalar
from scipy.integrate import quad


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


class TestNCH():
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
    odds = np.random.rand(*x.shape)*2

    # test output is more readable when function names (strings) are passed
    @pytest.mark.parametrize('dist_name', ['fnch', 'wnch'])
    def test_nch_hypergeom(self, dist_name):
        # Both noncentral hypergeometric distributions reduce to the
        # hypergeometric distribution when odds = 1
        dists = {'fnch': fnch, 'wnch': wnch}
        dist = dists[dist_name]
        x, N, m1, n = self.x, self.N, self.m1, self.n
        assert_allclose(dist.pmf(x, N, m1, n, odds=1),
                        hypergeom.pmf(x, N, m1, n))

    def test_fnch_naive(self):
        # test against a very simple implementation
        x, N, m1, n, odds = self.x, self.N, self.m1, self.n, self.odds

        @np.vectorize
        def pmf_mean_var(x, N, m1, n, w):
            # simple implementation of FNCH pmf
            m2 = N - m1
            xl = np.maximum(0, n-m2)
            xu = np.minimum(n, m1)

            def f(x):
                t1 = binom(m1, x)
                t2 = binom(m2, n - x)
                return t1 * t2 * w**x

            def P(k):
                return sum((f(y)*y**k for y in range(xl, xu + 1)))

            P0 = P(0)
            P1 = P(1)
            P2 = P(2)
            pmf = f(x) / P0
            mean = P1 / P0
            var = P2 / P0 - (P1 / P0)**2
            return pmf, mean, var

        pmf, mean, var = pmf_mean_var(x, N, m1, n, odds)
        assert_allclose(fnch.pmf(x, N, m1, n, odds), pmf)
        assert_allclose(fnch.stats(N, m1, n, odds, moments='m'), mean)
        assert_allclose(fnch.stats(N, m1, n, odds, moments='v'), var)

    def test_wnch_naive(self):
        # test against a very simple implementation

        np.random.seed(2)
        shape = (2, 4, 3)
        max_m = 100
        m1 = np.random.randint(1, max_m, size=shape)
        m2 = np.random.randint(1, max_m, size=shape)
        N = m1 + m2
        n = randint.rvs(0, N, size=N.shape)
        xl = np.maximum(0, n-m2)
        xu = np.minimum(n, m1)
        x = randint.rvs(xl, xu, size=xl.shape)
        w = np.random.rand(*x.shape)*2

        def support(N, m1, n, w):
            m2 = N - m1
            xl = np.maximum(0, n-m2)
            xu = np.minimum(n, m1)
            return xl, xu

        @np.vectorize
        def mean(N, m1, n, w):
            m2 = N - m1
            xl, xu = support(N, m1, n, w)

            def fun(u):
                return u/m1 + (1 - (n-u)/m2)**w - 1

            return root_scalar(fun, bracket=(xl, xu)).root

        assert_allclose(wnch.mean(N, m1, n, w), mean(N, m1, n, w), rtol=2e-2)

        @np.vectorize
        def variance(N, m1, n, w):
            m2 = N - m1
            u = mean(N, m1, n, w)
            a = u * (m1 - u)
            b = (n-u)*(u + m2 - n)
            return N*a*b / ((N-1) * (m1*b + m2*a))

        assert_allclose(wnch.stats(N, m1, n, w, moments='v'),
                        variance(N, m1, n, w), rtol=5e-2)

        @np.vectorize
        def pmf(x, N, m1, n, w):
            m2 = N - m1
            xl, xu = support(N, m1, n, w)

            def integrand(t):
                D = w*(m1 - x) + (m2 - (n-x))
                res = (1-t**(w/D))**x * (1-t**(1/D))**(n-x)
                return res

            def f(x):
                t1 = binom(m1, x)
                t2 = binom(m2, n - x)
                the_integral = quad(integrand, 0, 1,
                                    epsrel=1e-16, epsabs=1e-16)
                return t1 * t2 * the_integral[0]

            return f(x)

        pmf0 = pmf(x, N, m1, n, w)
        pmf1 = wnch.pmf(x, N, m1, n, w)

        atol, rtol = 1e-6, 1e-6
        i = np.abs(pmf1 - pmf0) < atol + rtol*np.abs(pmf0)
        assert(i.sum() > np.prod(shape) / 2)  # works at least half the time

        # for those that fail, discredit the naive implementation
        for N, m1, n, w in zip(N[~i], m1[~i], n[~i], w[~i]):
            # get the support
            m2 = N - m1
            xl, xu = support(N, m1, n, w)
            x = np.arange(xl, xu + 1)

            # calculate sum of pmf over the support
            # the naive implementation is very wrong in these cases
            assert pmf(x, N, m1, n, w).sum() < .5
            assert_allclose(wnch.pmf(x, N, m1, n, w).sum(), 1)

    def test_wallenius_against_mpmath(self):
        # precompute data with mpmath since niave implementation above
        # is not reliable. See source code in gh-13330.
        M = 50
        n = 30
        N = 20
        odds = 2.25
        # Expected results, computed with mpmath.
        sup = np.arange(21)
        pmf = np.array([3.699003068656875e-20,
                        5.89398584245431e-17,
                        2.1594437742911123e-14,
                        3.221458044649955e-12,
                        2.4658279241205077e-10,
                        1.0965862603981212e-08,
                        3.057890479665704e-07,
                        5.622818831643761e-06,
                        7.056482841531681e-05,
                        0.000618899425358671,
                        0.003854172932571669,
                        0.01720592676256026,
                        0.05528844897093792,
                        0.12772363313574242,
                        0.21065898367825722,
                        0.24465958845359234,
                        0.1955114898110033,
                        0.10355390084949237,
                        0.03414490375225675,
                        0.006231989845775931,
                        0.0004715577304677075])
        mean = 14.808018384813426
        var = 2.6085975877923717

        # wnch.pmf returns 0 for pmf(0) and pmf(1), and pmf(2)
        # has only three digits of accuracy (~ 2.1511e-14).
        assert_allclose(wnch.pmf(sup, M, n, N, odds), pmf,
                        rtol=1e-13, atol=1e-13)
        assert_allclose(wnch.mean(M, n, N, odds), mean, rtol=1e-13)
        assert_allclose(wnch.var(M, n, N, odds), var, rtol=1e-11)

    @pytest.mark.parametrize('dist_name', ['fnch', 'wnch'])
    def test_rvs_shape(self, dist_name):
        # Check that when given a size with more dimensions than the
        # dimensions of the broadcast parameters, rvs returns an array
        # with the correct shape.
        dists = {'fnch': fnch, 'wnch': wnch}
        dist = dists[dist_name]
        x = dist.rvs(50, 30, [[10], [20]], [0.5, 1.0, 2.0], size=(5, 1, 2, 3))
        assert x.shape == (5, 1, 2, 3)


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
