""" Test functions for stats module

"""
from __future__ import division, print_function, absolute_import

import warnings
import re
import sys

from numpy.testing import (TestCase, run_module_suite, assert_equal,
    assert_array_equal, assert_almost_equal, assert_array_almost_equal,
    assert_allclose, assert_, assert_raises, rand, dec)
from nose import SkipTest

import numpy
import numpy as np
from numpy import typecodes, array
from scipy.lib._version import NumpyVersion
from scipy import special
import scipy.stats as stats
from scipy.stats._distn_infrastructure import argsreduce
import scipy.stats.distributions

from scipy.special import xlogy


# python -OO strips docstrings
DOCSTRINGS_STRIPPED = sys.flags.optimize > 1


# Generate test cases to test cdf and distribution consistency.
# Note that this list does not include all distributions.
dists = ['uniform','norm','lognorm','expon','beta',
         'powerlaw','bradford','burr','fisk','cauchy','halfcauchy',
         'foldcauchy','gamma','gengamma','loggamma',
         'alpha','anglit','arcsine','betaprime',
         'dgamma','exponweib','exponpow','frechet_l','frechet_r',
         'gilbrat','f','ncf','chi2','chi','nakagami','genpareto',
         'genextreme','genhalflogistic','pareto','lomax','halfnorm',
         'halflogistic','fatiguelife','foldnorm','ncx2','t','nct',
         'weibull_min','weibull_max','dweibull','maxwell','rayleigh',
         'genlogistic', 'logistic','gumbel_l','gumbel_r','gompertz',
         'hypsecant', 'laplace', 'reciprocal','triang','tukeylambda',
         'vonmises', 'vonmises_line', 'pearson3']


def _assert_hasattr(a, b, msg=None):
    if msg is None:
        msg = '%s does not have attribute %s' % (a, b)
    assert_(hasattr(a, b), msg=msg)


def test_api_regression():
    # https://github.com/scipy/scipy/issues/3802
    _assert_hasattr(scipy.stats.distributions, 'f_gen')


# check function for test generator


def check_distribution(dist, args, alpha):
    D,pval = stats.kstest(dist,'', args=args, N=1000)
    if (pval < alpha):
        D,pval = stats.kstest(dist,'',args=args, N=1000)
        # if (pval < alpha):
        #    D,pval = stats.kstest(dist,'',args=args, N=1000)
        assert_(pval > alpha, msg="D = " + str(D) + "; pval = " + str(pval) +
               "; alpha = " + str(alpha) + "\nargs = " + str(args))

# nose test generator


def test_all_distributions():
    for dist in dists:
        distfunc = getattr(stats, dist)
        nargs = distfunc.numargs
        alpha = 0.01
        if dist == 'fatiguelife':
            alpha = 0.001

        if dist == 'frechet':
            args = tuple(2*rand(1))+(0,)+tuple(2*rand(2))
        elif dist == 'triang':
            args = tuple(rand(nargs))
        elif dist == 'reciprocal':
            vals = rand(nargs)
            vals[1] = vals[0] + 1.0
            args = tuple(vals)
        elif dist == 'vonmises':
            yield check_distribution, dist, (10,), alpha
            yield check_distribution, dist, (101,), alpha
            args = tuple(1.0+rand(nargs))
        else:
            args = tuple(1.0+rand(nargs))

        yield check_distribution, dist, args, alpha


def check_vonmises_pdf_periodic(k,l,s,x):
    vm = stats.vonmises(k,loc=l,scale=s)
    assert_almost_equal(vm.pdf(x),vm.pdf(x % (2*numpy.pi*s)))


def check_vonmises_cdf_periodic(k,l,s,x):
    vm = stats.vonmises(k,loc=l,scale=s)
    assert_almost_equal(vm.cdf(x) % 1,vm.cdf(x % (2*numpy.pi*s)) % 1)


def test_vonmises_pdf_periodic():
    for k in [0.1, 1, 101]:
        for x in [0,1,numpy.pi,10,100]:
            yield check_vonmises_pdf_periodic, k, 0, 1, x
            yield check_vonmises_pdf_periodic, k, 1, 1, x
            yield check_vonmises_pdf_periodic, k, 0, 10, x

            yield check_vonmises_cdf_periodic, k, 0, 1, x
            yield check_vonmises_cdf_periodic, k, 1, 1, x
            yield check_vonmises_cdf_periodic, k, 0, 10, x


def test_vonmises_line_support():
    assert_equal(stats.vonmises_line.a, -np.pi)
    assert_equal(stats.vonmises_line.b, np.pi)


class TestRandInt(TestCase):
    def test_rvs(self):
        vals = stats.randint.rvs(5,30,size=100)
        assert_(numpy.all(vals < 30) & numpy.all(vals >= 5))
        assert_(len(vals) == 100)
        vals = stats.randint.rvs(5,30,size=(2,50))
        assert_(numpy.shape(vals) == (2,50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.randint.rvs(15,46)
        assert_((val >= 15) & (val < 46))
        assert_(isinstance(val, numpy.ScalarType), msg=repr(type(val)))
        val = stats.randint(15,46).rvs(3)
        assert_(val.dtype.char in typecodes['AllInteger'])

    def test_pdf(self):
        k = numpy.r_[0:36]
        out = numpy.where((k >= 5) & (k < 30), 1.0/(30-5), 0)
        vals = stats.randint.pmf(k,5,30)
        assert_array_almost_equal(vals,out)

    def test_cdf(self):
        x = numpy.r_[0:36:100j]
        k = numpy.floor(x)
        out = numpy.select([k >= 30,k >= 5],[1.0,(k-5.0+1)/(30-5.0)],0)
        vals = stats.randint.cdf(x,5,30)
        assert_array_almost_equal(vals, out, decimal=12)


class TestBinom(TestCase):
    def test_rvs(self):
        vals = stats.binom.rvs(10, 0.75, size=(2, 50))
        assert_(numpy.all(vals >= 0) & numpy.all(vals <= 10))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.binom.rvs(10, 0.75)
        assert_(isinstance(val, int))
        val = stats.binom(10, 0.75).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllInteger'])

    def test_pmf(self):
        # regression test for Ticket #1842
        vals1 = stats.binom.pmf(100, 100,1)
        vals2 = stats.binom.pmf(0, 100,0)
        assert_allclose(vals1, 1.0, rtol=1e-15, atol=0)
        assert_allclose(vals2, 1.0, rtol=1e-15, atol=0)

    def test_entropy(self):
        # Basic entropy tests.
        b = stats.binom(2, 0.5)
        expected_p = np.array([0.25, 0.5, 0.25])
        expected_h = -sum(xlogy(expected_p, expected_p))
        h = b.entropy()
        assert_allclose(h, expected_h)

        b = stats.binom(2, 0.0)
        h = b.entropy()
        assert_equal(h, 0.0)

        b = stats.binom(2, 1.0)
        h = b.entropy()
        assert_equal(h, 0.0)

    def test_warns_p0(self):
        # no spurious warnigns are generated for p=0; gh-3817
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            assert_equal(stats.binom(n=2, p=0).mean(), 0)
            assert_equal(stats.binom(n=2, p=0).std(), 0)


class TestBernoulli(TestCase):
    def test_rvs(self):
        vals = stats.bernoulli.rvs(0.75, size=(2, 50))
        assert_(numpy.all(vals >= 0) & numpy.all(vals <= 1))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.bernoulli.rvs(0.75)
        assert_(isinstance(val, int))
        val = stats.bernoulli(0.75).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllInteger'])

    def test_entropy(self):
        # Simple tests of entropy.
        b = stats.bernoulli(0.25)
        expected_h = -0.25*np.log(0.25) - 0.75*np.log(0.75)
        h = b.entropy()
        assert_allclose(h, expected_h)

        b = stats.bernoulli(0.0)
        h = b.entropy()
        assert_equal(h, 0.0)

        b = stats.bernoulli(1.0)
        h = b.entropy()
        assert_equal(h, 0.0)


class TestNBinom(TestCase):
    def test_rvs(self):
        vals = stats.nbinom.rvs(10, 0.75, size=(2, 50))
        assert_(numpy.all(vals >= 0))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.nbinom.rvs(10, 0.75)
        assert_(isinstance(val, int))
        val = stats.nbinom(10, 0.75).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllInteger'])

    def test_pmf(self):
        # regression test for ticket 1779
        assert_allclose(np.exp(stats.nbinom.logpmf(700, 721, 0.52)),
                        stats.nbinom.pmf(700, 721, 0.52))
        # logpmf(0,1,1) shouldn't return nan (regression test for gh-4029)
        val = scipy.stats.nbinom.logpmf(0,1,1)
        assert_equal(val,0)


class TestGeom(TestCase):
    def test_rvs(self):
        vals = stats.geom.rvs(0.75, size=(2, 50))
        assert_(numpy.all(vals >= 0))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.geom.rvs(0.75)
        assert_(isinstance(val, int))
        val = stats.geom(0.75).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllInteger'])

    def test_pmf(self):
        vals = stats.geom.pmf([1,2,3],0.5)
        assert_array_almost_equal(vals,[0.5,0.25,0.125])

    def test_logpmf(self):
        # regression test for ticket 1793
        vals1 = np.log(stats.geom.pmf([1,2,3], 0.5))
        vals2 = stats.geom.logpmf([1,2,3], 0.5)
        assert_allclose(vals1, vals2, rtol=1e-15, atol=0)

        # regression test for gh-4028
        val = stats.geom.logpmf(1, 1)
        assert_equal(val, 0.0)

    def test_cdf_sf(self):
        vals = stats.geom.cdf([1, 2, 3], 0.5)
        vals_sf = stats.geom.sf([1, 2, 3], 0.5)
        expected = array([0.5, 0.75, 0.875])
        assert_array_almost_equal(vals, expected)
        assert_array_almost_equal(vals_sf, 1-expected)

    def test_logcdf_logsf(self):
        vals = stats.geom.logcdf([1, 2, 3], 0.5)
        vals_sf = stats.geom.logsf([1, 2, 3], 0.5)
        expected = array([0.5, 0.75, 0.875])
        assert_array_almost_equal(vals, np.log(expected))
        assert_array_almost_equal(vals_sf, np.log1p(-expected))

    def test_ppf(self):
        vals = stats.geom.ppf([0.5, 0.75, 0.875], 0.5)
        expected = array([1.0, 2.0, 3.0])
        assert_array_almost_equal(vals, expected)


class TestTruncnorm(TestCase):
    def test_ppf_ticket1131(self):
        vals = stats.truncnorm.ppf([-0.5,0,1e-4,0.5, 1-1e-4,1,2], -1., 1.,
                               loc=[3]*7, scale=2)
        expected = np.array([np.nan, 1, 1.00056419, 3, 4.99943581, 5, np.nan])
        assert_array_almost_equal(vals, expected)

    def test_isf_ticket1131(self):
        vals = stats.truncnorm.isf([-0.5,0,1e-4,0.5, 1-1e-4,1,2], -1., 1.,
                                   loc=[3]*7, scale=2)
        expected = np.array([np.nan, 5, 4.99943581, 3, 1.00056419, 1, np.nan])
        assert_array_almost_equal(vals, expected)

    def test_gh_2477_small_values(self):
        # Check a case that worked in the original issue.
        low, high = -11, -10
        x = stats.truncnorm.rvs(low, high, 0, 1, size=10)
        assert_(low < x.min() < x.max() < high)
        # Check a case that failed in the original issue.
        low, high = 10, 11
        x = stats.truncnorm.rvs(low, high, 0, 1, size=10)
        assert_(low < x.min() < x.max() < high)

    def test_gh_2477_large_values(self):
        # Check a case that fails because of extreme tailness.
        raise SkipTest('truncnorm rvs is know to fail at extreme tails')
        low, high = 100, 101
        x = stats.truncnorm.rvs(low, high, 0, 1, size=10)
        assert_(low < x.min() < x.max() < high)

    def test_gh_1489_trac_962_rvs(self):
        # Check the original example.
        low, high = 10, 15
        x = stats.truncnorm.rvs(low, high, 0, 1, size=10)
        assert_(low < x.min() < x.max() < high)


class TestHypergeom(TestCase):
    def test_rvs(self):
        vals = stats.hypergeom.rvs(20, 10, 3, size=(2, 50))
        assert_(numpy.all(vals >= 0) &
               numpy.all(vals <= 3))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.hypergeom.rvs(20, 3, 10)
        assert_(isinstance(val, int))
        val = stats.hypergeom(20, 3, 10).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllInteger'])

    def test_precision(self):
        # comparison number from mpmath
        M = 2500
        n = 50
        N = 500
        tot = M
        good = n
        hgpmf = stats.hypergeom.pmf(2, tot, good, N)
        assert_almost_equal(hgpmf, 0.0010114963068932233, 11)

    def test_cdf_above_one(self):
        # for some values of parameters, hypergeom cdf was >1, see gh-2238
        assert_(0 <= stats.hypergeom.cdf(30, 13397950, 4363, 12390) <= 1.0)

    def test_precision2(self):
        # Test hypergeom precision for large numbers.  See #1218.
        # Results compared with those from R.
        oranges = 9.9e4
        pears = 1.1e5
        fruits_eaten = np.array([3, 3.8, 3.9, 4, 4.1, 4.2, 5]) * 1e4
        quantile = 2e4
        res = []
        for eaten in fruits_eaten:
            res.append(stats.hypergeom.sf(quantile, oranges + pears, oranges, eaten))
        expected = np.array([0, 1.904153e-114, 2.752693e-66, 4.931217e-32,
                             8.265601e-11, 0.1237904, 1])
        assert_allclose(res, expected, atol=0, rtol=5e-7)

        # Test with array_like first argument
        quantiles = [1.9e4, 2e4, 2.1e4, 2.15e4]
        res2 = stats.hypergeom.sf(quantiles, oranges + pears, oranges, 4.2e4)
        expected2 = [1, 0.1237904, 6.511452e-34, 3.277667e-69]
        assert_allclose(res2, expected2, atol=0, rtol=5e-7)

    def test_entropy(self):
        # Simple tests of entropy.
        hg = stats.hypergeom(4, 1, 1)
        h = hg.entropy()
        expected_p = np.array([0.75, 0.25])
        expected_h = -np.sum(xlogy(expected_p, expected_p))
        assert_allclose(h, expected_h)

        hg = stats.hypergeom(1, 1, 1)
        h = hg.entropy()
        assert_equal(h, 0.0)


class TestLoggamma(TestCase):

    def test_stats(self):
        # The following precomputed values are from the table in section 2.2
        # of "A Statistical Study of Log-Gamma Distribution", by Ping Shing
        # Chan (thesis, McMaster University, 1993).
        table = np.array([
                # c,    mean,   var,    skew,    exc. kurt.
                0.5, -1.9635, 4.9348, -1.5351, 4.0000,
                1.0, -0.5772, 1.6449, -1.1395, 2.4000,
                12.0, 2.4427, 0.0869, -0.2946, 0.1735,
            ]).reshape(-1, 5)
        for c, mean, var, skew, kurt in table:
            computed = stats.loggamma.stats(c, moments='msvk')
            assert_array_almost_equal(computed, [mean, var, skew, kurt],
                                      decimal=4)


class TestLogser(TestCase):
    def test_rvs(self):
        vals = stats.logser.rvs(0.75, size=(2, 50))
        assert_(numpy.all(vals >= 1))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.logser.rvs(0.75)
        assert_(isinstance(val, int))
        val = stats.logser(0.75).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllInteger'])


class TestPareto(TestCase):
    def test_stats(self):
        # Check the stats() method with some simple values. Also check
        # that the calculations do not trigger RuntimeWarnings.
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)

            m, v, s, k = stats.pareto.stats(0.5, moments='mvsk')
            assert_equal(m, np.inf)
            assert_equal(v, np.inf)
            assert_equal(s, np.nan)
            assert_equal(k, np.nan)

            m, v, s, k = stats.pareto.stats(1.0, moments='mvsk')
            assert_equal(m, np.inf)
            assert_equal(v, np.inf)
            assert_equal(s, np.nan)
            assert_equal(k, np.nan)

            m, v, s, k = stats.pareto.stats(1.5, moments='mvsk')
            assert_equal(m, 3.0)
            assert_equal(v, np.inf)
            assert_equal(s, np.nan)
            assert_equal(k, np.nan)

            m, v, s, k = stats.pareto.stats(2.0, moments='mvsk')
            assert_equal(m, 2.0)
            assert_equal(v, np.inf)
            assert_equal(s, np.nan)
            assert_equal(k, np.nan)

            m, v, s, k = stats.pareto.stats(2.5, moments='mvsk')
            assert_allclose(m, 2.5 / 1.5)
            assert_allclose(v, 2.5 / (1.5*1.5*0.5))
            assert_equal(s, np.nan)
            assert_equal(k, np.nan)

            m, v, s, k = stats.pareto.stats(3.0, moments='mvsk')
            assert_allclose(m, 1.5)
            assert_allclose(v, 0.75)
            assert_equal(s, np.nan)
            assert_equal(k, np.nan)

            m, v, s, k = stats.pareto.stats(3.5, moments='mvsk')
            assert_allclose(m, 3.5 / 2.5)
            assert_allclose(v, 3.5 / (2.5*2.5*1.5))
            assert_allclose(s, (2*4.5/0.5)*np.sqrt(1.5/3.5))
            assert_equal(k, np.nan)

            m, v, s, k = stats.pareto.stats(4.0, moments='mvsk')
            assert_allclose(m, 4.0 / 3.0)
            assert_allclose(v, 4.0 / 18.0)
            assert_allclose(s, 2*(1+4.0)/(4.0-3) * np.sqrt((4.0-2)/4.0))
            assert_equal(k, np.nan)

            m, v, s, k = stats.pareto.stats(4.5, moments='mvsk')
            assert_allclose(m, 4.5 / 3.5)
            assert_allclose(v, 4.5 / (3.5*3.5*2.5))
            assert_allclose(s, (2*5.5/1.5) * np.sqrt(2.5/4.5))
            assert_allclose(k, 6*(4.5**3 + 4.5**2 - 6*4.5 - 2)/(4.5*1.5*0.5))


class TestGenpareto(TestCase):
    def test_ab(self):
        # c >= 0: a, b = [0, inf]
        for c in [1., 0.]:
            c = np.asarray(c)
            stats.genpareto._argcheck(c)  # ugh
            assert_equal(stats.genpareto.a, 0.)
            assert_(np.isposinf(stats.genpareto.b))

        # c < 0: a=0, b=1/|c|
        c = np.asarray(-2.)
        stats.genpareto._argcheck(c)
        assert_allclose([stats.genpareto.a, stats.genpareto.b], [0., 0.5])

    def test_c0(self):
        # with c=0, genpareto reduces to the exponential distribution
        rv = stats.genpareto(c=0.)
        x = np.linspace(0, 10., 30)
        assert_allclose(rv.pdf(x), stats.expon.pdf(x))
        assert_allclose(rv.cdf(x), stats.expon.cdf(x))
        assert_allclose(rv.sf(x), stats.expon.sf(x))

        q = np.linspace(0., 1., 10)
        assert_allclose(rv.ppf(q), stats.expon.ppf(q))

    def test_cm1(self):
        # with c=-1, genpareto reduces to the uniform distr on [0, 1]
        rv = stats.genpareto(c=-1.)
        x = np.linspace(0, 10., 30)
        assert_allclose(rv.pdf(x), stats.uniform.pdf(x))
        assert_allclose(rv.cdf(x), stats.uniform.cdf(x))
        assert_allclose(rv.sf(x), stats.uniform.sf(x))

        q = np.linspace(0., 1., 10)
        assert_allclose(rv.ppf(q), stats.uniform.ppf(q))

        # logpdf(1., c=-1) should be zero
        assert_allclose(rv.logpdf(1), 0)

    def test_x_inf(self):
        # make sure x=inf is handled gracefully
        rv = stats.genpareto(c=0.1)
        assert_allclose([rv.pdf(np.inf), rv.cdf(np.inf)], [0., 1.])
        assert_(np.isneginf(rv.logpdf(np.inf)))

        rv = stats.genpareto(c=0.)
        assert_allclose([rv.pdf(np.inf), rv.cdf(np.inf)], [0., 1.])
        assert_(np.isneginf(rv.logpdf(np.inf)))

        rv = stats.genpareto(c=-1.)
        assert_allclose([rv.pdf(np.inf), rv.cdf(np.inf)], [0., 1.])
        assert_(np.isneginf(rv.logpdf(np.inf)))

    def test_c_continuity(self):
        # pdf is continuous at c=0, -1
        x = np.linspace(0, 10, 30)
        for c in [0, -1]:
            pdf0 = stats.genpareto.pdf(x, c)
            for dc in [1e-14, -1e-14]:
                pdfc = stats.genpareto.pdf(x, c + dc)
                assert_allclose(pdf0, pdfc, atol=1e-12)

            cdf0 = stats.genpareto.cdf(x, c)
            for dc in [1e-14, 1e-14]:
                cdfc = stats.genpareto.cdf(x, c + dc)
                assert_allclose(cdf0, cdfc, atol=1e-12)

    def test_c_continuity_ppf(self):
        q = np.r_[np.logspace(1e-12, 0.01, base=0.1),
                  np.linspace(0.01, 1, 30, endpoint=False),
                  1. - np.logspace(1e-12, 0.01, base=0.1)]
        for c in [0., -1.]:
            ppf0 = stats.genpareto.ppf(q, c)
            for dc in [1e-14, -1e-14]:
                ppfc = stats.genpareto.ppf(q, c + dc)
                assert_allclose(ppf0, ppfc, atol=1e-12)

    def test_c_continuity_isf(self):
        q = np.r_[np.logspace(1e-12, 0.01, base=0.1),
                  np.linspace(0.01, 1, 30, endpoint=False),
                  1. - np.logspace(1e-12, 0.01, base=0.1)]
        for c in [0., -1.]:
            isf0 = stats.genpareto.isf(q, c)
            for dc in [1e-14, -1e-14]:
                isfc = stats.genpareto.isf(q, c + dc)
                assert_allclose(isf0, isfc, atol=1e-12)

    def test_cdf_ppf_roundtrip(self):
        # this should pass with machine precision. hat tip @pbrod
        q = np.r_[np.logspace(1e-12, 0.01, base=0.1),
                  np.linspace(0.01, 1, 30, endpoint=False),
                  1. - np.logspace(1e-12, 0.01, base=0.1)]
        for c in [1e-8, -1e-18, 1e-15, -1e-15]:
            assert_allclose(stats.genpareto.cdf(stats.genpareto.ppf(q, c), c),
                    q, atol=1e-15)


class TestPearson3(TestCase):
    def test_rvs(self):
        vals = stats.pearson3.rvs(0.1, size=(2, 50))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllFloat'])
        val = stats.pearson3.rvs(0.5)
        assert_(isinstance(val, float))
        val = stats.pearson3(0.5).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllFloat'])
        assert_(len(val) == 3)

    def test_pdf(self):
        vals = stats.pearson3.pdf(2, [0.0, 0.1, 0.2])
        assert_allclose(vals, np.array([0.05399097, 0.05555481, 0.05670246]),
                        atol=1e-6)
        vals = stats.pearson3.pdf(-3, 0.1)
        assert_allclose(vals, np.array([0.00313791]), atol=1e-6)
        vals = stats.pearson3.pdf([-3,-2,-1,0,1], 0.1)
        assert_allclose(vals, np.array([0.00313791, 0.05192304, 0.25028092,
                                        0.39885918, 0.23413173]), atol=1e-6)

    def test_cdf(self):
        vals = stats.pearson3.cdf(2, [0.0, 0.1, 0.2])
        assert_allclose(vals, np.array([0.97724987, 0.97462004, 0.97213626]),
                        atol=1e-6)
        vals = stats.pearson3.cdf(-3, 0.1)
        assert_allclose(vals, [0.00082256], atol=1e-6)
        vals = stats.pearson3.cdf([-3,-2,-1,0,1], 0.1)
        assert_allclose(vals, [8.22563821e-04, 1.99860448e-02, 1.58550710e-01,
                               5.06649130e-01, 8.41442111e-01], atol=1e-6)


class TestPoisson(TestCase):
    def test_rvs(self):
        vals = stats.poisson.rvs(0.5, size=(2, 50))
        assert_(numpy.all(vals >= 0))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.poisson.rvs(0.5)
        assert_(isinstance(val, int))
        val = stats.poisson(0.5).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllInteger'])

    def test_stats(self):
        mu = 16.0
        result = stats.poisson.stats(mu, moments='mvsk')
        assert_allclose(result, [mu, mu, np.sqrt(1.0/mu), 1.0/mu])


class TestZipf(TestCase):
    def test_rvs(self):
        vals = stats.zipf.rvs(1.5, size=(2, 50))
        assert_(numpy.all(vals >= 1))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.zipf.rvs(1.5)
        assert_(isinstance(val, int))
        val = stats.zipf(1.5).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllInteger'])

    def test_moments(self):
        # n-th moment is finite iff a > n + 1
        m, v = stats.zipf.stats(a=2.8)
        assert_(np.isfinite(m))
        assert_equal(v, np.inf)

        s, k = stats.zipf.stats(a=4.8, moments='sk')
        assert_(not np.isfinite([s, k]).all())


class TestDLaplace(TestCase):
    def test_rvs(self):
        vals = stats.dlaplace.rvs(1.5, size=(2, 50))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.dlaplace.rvs(1.5)
        assert_(isinstance(val, int))
        val = stats.dlaplace(1.5).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllInteger'])
        assert_(stats.dlaplace.rvs(0.8) is not None)

    def test_stats(self):
        # compare the explicit formulas w/ direct summation using pmf
        a = 1.
        dl = stats.dlaplace(a)
        m, v, s, k = dl.stats('mvsk')

        N = 37
        xx = np.arange(-N, N+1)
        pp = dl.pmf(xx)
        m2, m4 = np.sum(pp*xx**2), np.sum(pp*xx**4)
        assert_equal((m, s), (0,0))
        assert_allclose((v, k), (m2, m4/m2**2 - 3.), atol=1e-14, rtol=1e-8)

    def test_stats2(self):
        a = np.log(2.)
        dl = stats.dlaplace(a)
        m, v, s, k = dl.stats('mvsk')
        assert_equal((m, s), (0.,0.))
        assert_allclose((v, k), (4., 3.25))


class TestInvGamma(TestCase):
    @dec.skipif(NumpyVersion(np.__version__) < '1.7.0',
                "assert_* funcs broken with inf/nan")
    def test_invgamma_inf_gh_1866(self):
        # invgamma's moments are only finite for a>n
        # specific numbers checked w/ boost 1.54
        with warnings.catch_warnings():
            warnings.simplefilter('error', RuntimeWarning)
            mvsk = stats.invgamma.stats(a=19.31, moments='mvsk')
            assert_allclose(mvsk,
                [0.05461496450, 0.0001723162534, 1.020362676, 2.055616582])

            a = [1.1, 3.1, 5.6]
            mvsk = stats.invgamma.stats(a=a, moments='mvsk')
            expected = ([10., 0.476190476, 0.2173913043],       # mmm
                        [np.inf, 0.2061430632, 0.01312749422],  # vvv
                        [np.nan, 41.95235392, 2.919025532],     # sss
                        [np.nan, np.nan, 24.51923076])          # kkk
            for x, y in zip(mvsk, expected):
                assert_almost_equal(x, y)


class TestF(TestCase):
    def test_f_moments(self):
        # n-th moment of F distributions is only finite for n < dfd / 2
        m, v, s, k = stats.f.stats(11, 6.5, moments='mvsk')
        assert_(np.isfinite(m))
        assert_(np.isfinite(v))
        assert_(np.isfinite(s))
        assert_(not np.isfinite(k))

    def test_moments_warnings(self):
        # no warnings should be generated for dfd = 2, 4, 6, 8 (div by zero)
        with warnings.catch_warnings():
            warnings.simplefilter('error', RuntimeWarning)
            stats.f.stats(dfn=[11]*4, dfd=[2, 4, 6, 8], moments='mvsk')

    @dec.knownfailureif(True, 'f stats does not properly broadcast')
    def test_stats_broadcast(self):
        # stats do not fully broadcast just yet
        mv = stats.f.stats(dfn=11, dfd=[11, 12])


def test_rvgeneric_std():
    # Regression test for #1191
    assert_array_almost_equal(stats.t.std([5, 6]), [1.29099445, 1.22474487])


class TestRvDiscrete(TestCase):
    def test_rvs(self):
        states = [-1,0,1,2,3,4]
        probability = [0.0,0.3,0.4,0.0,0.3,0.0]
        samples = 1000
        r = stats.rv_discrete(name='sample',values=(states,probability))
        x = r.rvs(size=samples)
        assert_(isinstance(x, numpy.ndarray))

        for s,p in zip(states,probability):
            assert_(abs(sum(x == s)/float(samples) - p) < 0.05)

        x = r.rvs()
        assert_(isinstance(x, int))

    def test_entropy(self):
        # Basic tests of entropy.
        pvals = np.array([0.25, 0.45, 0.3])
        p = stats.rv_discrete(values=([0, 1, 2], pvals))
        expected_h = -sum(xlogy(pvals, pvals))
        h = p.entropy()
        assert_allclose(h, expected_h)

        p = stats.rv_discrete(values=([0, 1, 2], [1.0, 0, 0]))
        h = p.entropy()
        assert_equal(h, 0.0)


class TestExpon(TestCase):
    def test_zero(self):
        assert_equal(stats.expon.pdf(0),1)

    def test_tail(self):  # Regression test for ticket 807
        assert_equal(stats.expon.cdf(1e-18), 1e-18)
        assert_equal(stats.expon.isf(stats.expon.sf(40)), 40)


class TestGenExpon(TestCase):
    def test_pdf_unity_area(self):
        from scipy.integrate import simps
        # PDF should integrate to one
        assert_almost_equal(simps(stats.genexpon.pdf(numpy.arange(0,10,0.01),
                                                     0.5, 0.5, 2.0),
                                  dx=0.01), 1, 1)

    def test_cdf_bounds(self):
        # CDF should always be positive
        cdf = stats.genexpon.cdf(numpy.arange(0, 10, 0.01), 0.5, 0.5, 2.0)
        assert_(numpy.all((0 <= cdf) & (cdf <= 1)))


class TestExponpow(TestCase):
    def test_tail(self):
        assert_almost_equal(stats.exponpow.cdf(1e-10, 2.), 1e-20)
        assert_almost_equal(stats.exponpow.isf(stats.exponpow.sf(5, .8), .8), 5)


class TestSkellam(TestCase):
    def test_pmf(self):
        # comparison to R
        k = numpy.arange(-10, 15)
        mu1, mu2 = 10, 5
        skpmfR = numpy.array(
                   [4.2254582961926893e-005, 1.1404838449648488e-004,
                    2.8979625801752660e-004, 6.9177078182101231e-004,
                    1.5480716105844708e-003, 3.2412274963433889e-003,
                    6.3373707175123292e-003, 1.1552351566696643e-002,
                    1.9606152375042644e-002, 3.0947164083410337e-002,
                    4.5401737566767360e-002, 6.1894328166820688e-002,
                    7.8424609500170578e-002, 9.2418812533573133e-002,
                    1.0139793148019728e-001, 1.0371927988298846e-001,
                    9.9076583077406091e-002, 8.8546660073089561e-002,
                    7.4187842052486810e-002, 5.8392772862200251e-002,
                    4.3268692953013159e-002, 3.0248159818374226e-002,
                    1.9991434305603021e-002, 1.2516877303301180e-002,
                    7.4389876226229707e-003])

        assert_almost_equal(stats.skellam.pmf(k, mu1, mu2), skpmfR, decimal=15)

    def test_cdf(self):
        # comparison to R, only 5 decimals
        k = numpy.arange(-10, 15)
        mu1, mu2 = 10, 5
        skcdfR = numpy.array(
                   [6.4061475386192104e-005, 1.7810985988267694e-004,
                    4.6790611790020336e-004, 1.1596768997212152e-003,
                    2.7077485103056847e-003, 5.9489760066490718e-003,
                    1.2286346724161398e-002, 2.3838698290858034e-002,
                    4.3444850665900668e-002, 7.4392014749310995e-002,
                    1.1979375231607835e-001, 1.8168808048289900e-001,
                    2.6011268998306952e-001, 3.5253150251664261e-001,
                    4.5392943399683988e-001, 5.5764871387982828e-001,
                    6.5672529695723436e-001, 7.4527195703032389e-001,
                    8.1945979908281064e-001, 8.7785257194501087e-001,
                    9.2112126489802404e-001, 9.5136942471639818e-001,
                    9.7136085902200120e-001, 9.8387773632530240e-001,
                    9.9131672394792536e-001])

        assert_almost_equal(stats.skellam.cdf(k, mu1, mu2), skcdfR, decimal=5)


class TestLognorm(TestCase):
    def test_pdf(self):
        # Regression test for Ticket #1471: avoid nan with 0/0 situation
        with np.errstate(divide='ignore'):
            pdf = stats.lognorm.pdf([0, 0.5, 1], 1)
            assert_array_almost_equal(pdf, [0.0, 0.62749608, 0.39894228])


class TestBeta(TestCase):
    def test_logpdf(self):
        # Regression test for Ticket #1326: avoid nan with 0*log(0) situation
        logpdf = stats.beta.logpdf(0,1,0.5)
        assert_almost_equal(logpdf, -0.69314718056)
        logpdf = stats.beta.logpdf(0,0.5,1)
        assert_almost_equal(logpdf, np.inf)

    def test_logpdf_ticket_1866(self):
        alpha, beta = 267, 1472
        x = np.array([0.2, 0.5, 0.6])
        b = stats.beta(alpha, beta)
        assert_allclose(b.logpdf(x).sum(), -1201.699061824062)
        assert_allclose(b.pdf(x), np.exp(b.logpdf(x)))


class TestBetaPrime(TestCase):
    def test_logpdf(self):
        alpha, beta = 267, 1472
        x = np.array([0.2, 0.5, 0.6])
        b = stats.betaprime(alpha, beta)
        assert_(np.isfinite(b.logpdf(x)).all())
        assert_allclose(b.pdf(x), np.exp(b.logpdf(x)))

    def test_cdf(self):
        # regression test for issue 4030: Implentation of 
        # scipy.stats.betaprime.cdf()
        x = stats.betaprime.cdf(0,0.2,0.3)
        assert_equal(x,0.0)


class TestGamma(TestCase):
    def test_pdf(self):
        # a few test cases to compare with R
        pdf = stats.gamma.pdf(90, 394, scale=1./5)
        assert_almost_equal(pdf, 0.002312341)

        pdf = stats.gamma.pdf(3, 10, scale=1./5)
        assert_almost_equal(pdf, 0.1620358)

    def test_logpdf(self):
        # Regression test for Ticket #1326: cornercase avoid nan with 0*log(0)
        # situation
        logpdf = stats.gamma.logpdf(0,1)
        assert_almost_equal(logpdf, 0)


class TestChi2(TestCase):
    # regression tests after precision improvements, ticket:1041, not verified
    def test_precision(self):
        assert_almost_equal(stats.chi2.pdf(1000, 1000), 8.919133934753128e-003, 14)
        assert_almost_equal(stats.chi2.pdf(100, 100), 0.028162503162596778, 14)


class TestArrayArgument(TestCase):  # test for ticket:992
    def test_noexception(self):
        rvs = stats.norm.rvs(loc=(np.arange(5)), scale=np.ones(5), size=(10,5))
        assert_equal(rvs.shape, (10,5))


class TestDocstring(TestCase):
    def test_docstrings(self):
        # See ticket #761
        if stats.rayleigh.__doc__ is not None:
            self.assertTrue("rayleigh" in stats.rayleigh.__doc__.lower())
        if stats.bernoulli.__doc__ is not None:
            self.assertTrue("bernoulli" in stats.bernoulli.__doc__.lower())

    def test_no_name_arg(self):
        # If name is not given, construction shouldn't fail.  See #1508.
        stats.rv_continuous()
        stats.rv_discrete()


class TestEntropy(TestCase):
    def test_entropy_positive(self):
        # See ticket #497
        pk = [0.5,0.2,0.3]
        qk = [0.1,0.25,0.65]
        eself = stats.entropy(pk,pk)
        edouble = stats.entropy(pk,qk)
        assert_(0.0 == eself)
        assert_(edouble >= 0.0)

    def test_entropy_base(self):
        pk = np.ones(16, float)
        S = stats.entropy(pk, base=2.)
        assert_(abs(S - 4.) < 1.e-5)

        qk = np.ones(16, float)
        qk[:8] = 2.
        S = stats.entropy(pk, qk)
        S2 = stats.entropy(pk, qk, base=2.)
        assert_(abs(S/S2 - np.log(2.)) < 1.e-5)

    def test_entropy_zero(self):
        # Test for PR-479
        assert_almost_equal(stats.entropy([0, 1, 2]), 0.63651416829481278,
                            decimal=12)

    def test_entropy_2d(self):
        pk = [[0.1, 0.2], [0.6, 0.3], [0.3, 0.5]]
        qk = [[0.2, 0.1], [0.3, 0.6], [0.5, 0.3]]
        assert_array_almost_equal(stats.entropy(pk, qk),
                [0.1933259, 0.18609809])

    @dec.skipif(NumpyVersion(np.__version__) < '1.7.0',
                "assert_* funcs broken with inf/nan")
    def test_entropy_2d_zero(self):
        pk = [[0.1, 0.2], [0.6, 0.3], [0.3, 0.5]]
        qk = [[0.0, 0.1], [0.3, 0.6], [0.5, 0.3]]
        assert_array_almost_equal(stats.entropy(pk, qk),
                [np.inf, 0.18609809])

        pk[0][0] = 0.0
        assert_array_almost_equal(stats.entropy(pk, qk),
                [0.17403988, 0.18609809])


def TestArgsreduce():
    a = array([1,3,2,1,2,3,3])
    b,c = argsreduce(a > 1, a, 2)

    assert_array_equal(b, [3,2,2,3,3])
    assert_array_equal(c, [2,2,2,2,2])

    b,c = argsreduce(2 > 1, a, 2)
    assert_array_equal(b, a[0])
    assert_array_equal(c, [2])

    b,c = argsreduce(a > 0, a, 2)
    assert_array_equal(b, a)
    assert_array_equal(c, [2] * numpy.size(a))


class TestFitMethod(object):
    skip = ['ncf']

    @dec.slow
    def test_fit(self):
        def check(func, dist, args, alpha):
            if dist in self.skip:
                raise SkipTest("%s fit known to fail" % dist)
            distfunc = getattr(stats, dist)
            with np.errstate(all='ignore'):
                res = distfunc.rvs(*args, **{'size':200})
                vals = distfunc.fit(res)
                vals2 = distfunc.fit(res, optimizer='powell')
            # Only check the length of the return
            # FIXME: should check the actual results to see if we are 'close'
            #   to what was created --- but what is 'close' enough
            if dist == 'frechet':
                assert_(len(vals) == len(args))
                assert_(len(vals2) == len(args))
            else:
                assert_(len(vals) == 2+len(args))
                assert_(len(vals2) == 2+len(args))

        for func, dist, args, alpha in test_all_distributions():
            yield check, func, dist, args, alpha

    @dec.slow
    def test_fix_fit(self):
        def check(func, dist, args, alpha):
            # Not sure why 'ncf', and 'beta' are failing
            # frechet has different len(args) than distfunc.numargs
            if dist in self.skip + ['frechet']:
                raise SkipTest("%s fit known to fail" % dist)
            distfunc = getattr(stats, dist)
            with np.errstate(all='ignore'):
                res = distfunc.rvs(*args, **{'size':200})
                vals = distfunc.fit(res,floc=0)
                vals2 = distfunc.fit(res,fscale=1)
                assert_(len(vals) == 2+len(args))
                assert_(vals[-2] == 0)
                assert_(vals2[-1] == 1)
                assert_(len(vals2) == 2+len(args))
                if len(args) > 0:
                    vals3 = distfunc.fit(res, f0=args[0])
                    assert_(len(vals3) == 2+len(args))
                    assert_(vals3[0] == args[0])
                if len(args) > 1:
                    vals4 = distfunc.fit(res, f1=args[1])
                    assert_(len(vals4) == 2+len(args))
                    assert_(vals4[1] == args[1])
                if len(args) > 2:
                    vals5 = distfunc.fit(res, f2=args[2])
                    assert_(len(vals5) == 2+len(args))
                    assert_(vals5[2] == args[2])

        for func, dist, args, alpha in test_all_distributions():
            yield check, func, dist, args, alpha

    def test_fix_fit_2args_lognorm(self):
        # Regression test for #1551.
        np.random.seed(12345)
        with np.errstate(all='ignore'):
            x = stats.lognorm.rvs(0.25, 0., 20.0, size=20)
            assert_allclose(np.array(stats.lognorm.fit(x, floc=0, fscale=20)),
                            [0.25888672, 0, 20], atol=1e-5)

    def test_fix_fit_norm(self):
        x = np.arange(1, 6)

        loc, scale = stats.norm.fit(x)
        assert_almost_equal(loc, 3)
        assert_almost_equal(scale, np.sqrt(2))

        loc, scale = stats.norm.fit(x, floc=2)
        assert_equal(loc, 2)
        assert_equal(scale, np.sqrt(3))

        loc, scale = stats.norm.fit(x, fscale=2)
        assert_almost_equal(loc, 3)
        assert_equal(scale, 2)

    def test_fix_fit_gamma(self):
        x = np.arange(1, 6)
        meanlog = np.log(x).mean()

        # A basic test of gamma.fit with floc=0.
        floc = 0
        a, loc, scale = stats.gamma.fit(x, floc=floc)
        s = np.log(x.mean()) - meanlog
        assert_almost_equal(np.log(a) - special.digamma(a), s, decimal=5)
        assert_equal(loc, floc)
        assert_almost_equal(scale, x.mean()/a, decimal=8)

        # Regression tests for gh-2514.
        # The problem was that if `floc=0` was given, any other fixed
        # parameters were ignored.
        f0 = 1
        floc = 0
        a, loc, scale = stats.gamma.fit(x, f0=f0, floc=floc)
        assert_equal(a, f0)
        assert_equal(loc, floc)
        assert_almost_equal(scale, x.mean()/a, decimal=8)

        f0 = 2
        floc = 0
        a, loc, scale = stats.gamma.fit(x, f0=f0, floc=floc)
        assert_equal(a, f0)
        assert_equal(loc, floc)
        assert_almost_equal(scale, x.mean()/a, decimal=8)

        # loc and scale fixed.
        floc = 0
        fscale = 2
        a, loc, scale = stats.gamma.fit(x, floc=floc, fscale=fscale)
        assert_equal(loc, floc)
        assert_equal(scale, fscale)
        c = meanlog - np.log(fscale)
        assert_almost_equal(special.digamma(a), c)

    def test_fix_fit_beta(self):
        # Test beta.fit when both floc and fscale are given.

        def mlefunc(a, b, x):
            # Zeros of this function are critical points of
            # the maximum likelihood function.
            n = len(x)
            s1 = np.log(x).sum()
            s2 = np.log(1-x).sum()
            psiab = special.psi(a + b)
            func = [s1 - n * (-psiab + special.psi(a)),
                    s2 - n * (-psiab + special.psi(b))]
            return func

        # Basic test with floc and fscale given.
        x = np.array([0.125, 0.25, 0.5])
        a, b, loc, scale = stats.beta.fit(x, floc=0, fscale=1)
        assert_equal(loc, 0)
        assert_equal(scale, 1)
        assert_allclose(mlefunc(a, b, x), [0,0], atol=1e-6)

        # Basic test with f0, floc and fscale given.
        # This is also a regression test for gh-2514.
        x = np.array([0.125, 0.25, 0.5])
        a, b, loc, scale = stats.beta.fit(x, f0=2, floc=0, fscale=1)
        assert_equal(a, 2)
        assert_equal(loc, 0)
        assert_equal(scale, 1)
        da, db = mlefunc(a, b, x)
        assert_allclose(db, 0, atol=1e-5)

        # Same floc and fscale values as above, but reverse the data
        # and fix b (f1).
        x2 = 1 - x
        a2, b2, loc2, scale2 = stats.beta.fit(x2, f1=2, floc=0, fscale=1)
        assert_equal(b2, 2)
        assert_equal(loc2, 0)
        assert_equal(scale2, 1)
        da, db = mlefunc(a2, b2, x2)
        assert_allclose(da, 0, atol=1e-5)
        # a2 of this test should equal b from above.
        assert_almost_equal(a2, b)

        # Check for detection of data out of bounds when floc and fscale
        # are given.
        assert_raises(ValueError, stats.beta.fit, x, floc=0.5, fscale=1)
        y = np.array([0, .5, 1])
        assert_raises(ValueError, stats.beta.fit, y, floc=0, fscale=1)
        assert_raises(ValueError, stats.beta.fit, y, floc=0, fscale=1, f0=2)
        assert_raises(ValueError, stats.beta.fit, y, floc=0, fscale=1, f1=2)

        # Check that attempting to fix all the parameters raises a ValueError.
        assert_raises(ValueError, stats.beta.fit, y, f0=0, f1=1,
                                                     floc=2, fscale=3)


class TestFrozen(TestCase):
    # Test that a frozen distribution gives the same results as the original object.
    #
    # Only tested for the normal distribution (with loc and scale specified)
    # and for the gamma distribution (with a shape parameter specified).
    def test_norm(self):
        dist = stats.norm
        frozen = stats.norm(loc=10.0, scale=3.0)

        result_f = frozen.pdf(20.0)
        result = dist.pdf(20.0, loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.cdf(20.0)
        result = dist.cdf(20.0, loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.ppf(0.25)
        result = dist.ppf(0.25, loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.isf(0.25)
        result = dist.isf(0.25, loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.sf(10.0)
        result = dist.sf(10.0, loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.median()
        result = dist.median(loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.mean()
        result = dist.mean(loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.var()
        result = dist.var(loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.std()
        result = dist.std(loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.entropy()
        result = dist.entropy(loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.moment(2)
        result = dist.moment(2,loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        assert_equal(frozen.a, dist.a)
        assert_equal(frozen.b, dist.b)

    def test_gamma(self):
        a = 2.0
        dist = stats.gamma
        frozen = stats.gamma(a)

        result_f = frozen.pdf(20.0)
        result = dist.pdf(20.0, a)
        assert_equal(result_f, result)

        result_f = frozen.cdf(20.0)
        result = dist.cdf(20.0, a)
        assert_equal(result_f, result)

        result_f = frozen.ppf(0.25)
        result = dist.ppf(0.25, a)
        assert_equal(result_f, result)

        result_f = frozen.isf(0.25)
        result = dist.isf(0.25, a)
        assert_equal(result_f, result)

        result_f = frozen.sf(10.0)
        result = dist.sf(10.0, a)
        assert_equal(result_f, result)

        result_f = frozen.median()
        result = dist.median(a)
        assert_equal(result_f, result)

        result_f = frozen.mean()
        result = dist.mean(a)
        assert_equal(result_f, result)

        result_f = frozen.var()
        result = dist.var(a)
        assert_equal(result_f, result)

        result_f = frozen.std()
        result = dist.std(a)
        assert_equal(result_f, result)

        result_f = frozen.entropy()
        result = dist.entropy(a)
        assert_equal(result_f, result)

        result_f = frozen.moment(2)
        result = dist.moment(2, a)
        assert_equal(result_f, result)

        assert_equal(frozen.a, frozen.dist.a)
        assert_equal(frozen.b, frozen.dist.b)

    def test_regression_ticket_1293(self):
        # Create a frozen distribution.
        frozen = stats.lognorm(1)
        # Call one of its methods that does not take any keyword arguments.
        m1 = frozen.moment(2)
        # Now call a method that takes a keyword argument.
        frozen.stats(moments='mvsk')
        # Call moment(2) again.
        # After calling stats(), the following was raising an exception.
        # So this test passes if the following does not raise an exception.
        m2 = frozen.moment(2)
        # The following should also be true, of course.  But it is not
        # the focus of this test.
        assert_equal(m1, m2)

    def test_ab(self):
        # test that the support of a frozen distribution
        # (i) remains frozen even if it changes for the original one
        # (ii) is actually correct if the shape parameters are such that
        #      the values of [a, b] are not the default [0, inf]
        # take a genpareto as an example where the support
        # depends on the value of the shape parameter:
        # for c > 0: a, b = 0, inf
        # for c < 0: a, b = 0, -1/c
        rv = stats.genpareto(c=-0.1)
        a, b = rv.dist.a, rv.dist.b
        assert_equal([a, b], [0., 10.])
        assert_equal([rv.a, rv.b], [0., 10.])

        stats.genpareto.pdf(0, c=0.1)  # this changes genpareto.b
        assert_equal([rv.dist.a, rv.dist.b], [a, b])
        assert_equal([rv.a, rv.b], [a, b])

        rv1 = stats.genpareto(c=0.1)
        assert_(rv1.dist is not rv.dist)

    def test_rv_frozen_in_namespace(self):
        # Regression test for gh-3522
        assert_(hasattr(stats.distributions, 'rv_frozen'))

    def test_random_state(self):
        # only check that the random_state attribute exists,
        frozen = stats.norm()
        assert_(hasattr(frozen, 'random_state'))

        # ... that it can be set,
        frozen.random_state = 42
        assert_equal(frozen.random_state.get_state(),
                     np.random.RandomState(42).get_state())

        # ... and that .rvs method accepts it as an argument
        rndm = np.random.RandomState(1234)
        frozen.rvs(size=8, random_state=rndm)

    def test_expect(self):
        # smoke test the expect method of the frozen distribution
        # only take a gamma w/loc and scale and poisson with loc specified
        def func(x):
            return x

        gm = stats.gamma(a=2, loc=3, scale=4)
        gm_val = gm.expect(func, lb=1, ub=2, conditional=True)
        gamma_val = stats.gamma.expect(func, args=(2,), loc=3, scale=4,
                                       lb=1, ub=2, conditional=True)
        assert_allclose(gm_val, gamma_val)

        p = stats.poisson(3, loc=4)
        p_val = p.expect(func)
        poisson_val = stats.poisson.expect(func, args=(3,), loc=4)
        assert_allclose(p_val, poisson_val)


class TestExpect(TestCase):
    # Test for expect method.
    #
    # Uses normal distribution and beta distribution for finite bounds, and
    # hypergeom for discrete distribution with finite support
    def test_norm(self):
        v = stats.norm.expect(lambda x: (x-5)*(x-5), loc=5, scale=2)
        assert_almost_equal(v, 4, decimal=14)

        m = stats.norm.expect(lambda x: (x), loc=5, scale=2)
        assert_almost_equal(m, 5, decimal=14)

        lb = stats.norm.ppf(0.05, loc=5, scale=2)
        ub = stats.norm.ppf(0.95, loc=5, scale=2)
        prob90 = stats.norm.expect(lambda x: 1, loc=5, scale=2, lb=lb, ub=ub)
        assert_almost_equal(prob90, 0.9, decimal=14)

        prob90c = stats.norm.expect(lambda x: 1, loc=5, scale=2, lb=lb, ub=ub,
                                    conditional=True)
        assert_almost_equal(prob90c, 1., decimal=14)

    def test_beta(self):
        # case with finite support interval
        v = stats.beta.expect(lambda x: (x-19/3.)*(x-19/3.), args=(10,5),
                              loc=5, scale=2)
        assert_almost_equal(v, 1./18., decimal=13)

        m = stats.beta.expect(lambda x: x, args=(10,5), loc=5., scale=2.)
        assert_almost_equal(m, 19/3., decimal=13)

        ub = stats.beta.ppf(0.95, 10, 10, loc=5, scale=2)
        lb = stats.beta.ppf(0.05, 10, 10, loc=5, scale=2)
        prob90 = stats.beta.expect(lambda x: 1., args=(10,10), loc=5.,
                                   scale=2.,lb=lb, ub=ub, conditional=False)
        assert_almost_equal(prob90, 0.9, decimal=13)

        prob90c = stats.beta.expect(lambda x: 1, args=(10,10), loc=5,
                                    scale=2, lb=lb, ub=ub, conditional=True)
        assert_almost_equal(prob90c, 1., decimal=13)

    def test_hypergeom(self):
        # test case with finite bounds

        # without specifying bounds
        m_true, v_true = stats.hypergeom.stats(20, 10, 8, loc=5.)
        m = stats.hypergeom.expect(lambda x: x, args=(20, 10, 8), loc=5.)
        assert_almost_equal(m, m_true, decimal=13)

        v = stats.hypergeom.expect(lambda x: (x-9.)**2, args=(20, 10, 8),
                                   loc=5.)
        assert_almost_equal(v, v_true, decimal=14)

        # with bounds, bounds equal to shifted support
        v_bounds = stats.hypergeom.expect(lambda x: (x-9.)**2, args=(20, 10, 8),
                                          loc=5., lb=5, ub=13)
        assert_almost_equal(v_bounds, v_true, decimal=14)

        # drop boundary points
        prob_true = 1-stats.hypergeom.pmf([5, 13], 20, 10, 8, loc=5).sum()
        prob_bounds = stats.hypergeom.expect(lambda x: 1, args=(20, 10, 8),
                                          loc=5., lb=6, ub=12)
        assert_almost_equal(prob_bounds, prob_true, decimal=13)

        # conditional
        prob_bc = stats.hypergeom.expect(lambda x: 1, args=(20, 10, 8), loc=5.,
                                           lb=6, ub=12, conditional=True)
        assert_almost_equal(prob_bc, 1, decimal=14)

        # check simple integral
        prob_b = stats.hypergeom.expect(lambda x: 1, args=(20, 10, 8),
                                        lb=0, ub=8)
        assert_almost_equal(prob_b, 1, decimal=13)

    def test_poisson(self):
        # poisson, use lower bound only
        prob_bounds = stats.poisson.expect(lambda x: 1, args=(2,), lb=3,
                                      conditional=False)
        prob_b_true = 1-stats.poisson.cdf(2,2)
        assert_almost_equal(prob_bounds, prob_b_true, decimal=14)

        prob_lb = stats.poisson.expect(lambda x: 1, args=(2,), lb=2,
                                       conditional=True)
        assert_almost_equal(prob_lb, 1, decimal=14)

    def test_genhalflogistic(self):
        # genhalflogistic, changes upper bound of support in _argcheck
        # regression test for gh-2622
        halflog = stats.genhalflogistic
        # check consistency when calling expect twice with the same input
        res1 = halflog.expect(args=(1.5,))
        halflog.expect(args=(0.5,))
        res2 = halflog.expect(args=(1.5,))
        assert_almost_equal(res1, res2, decimal=14)

    def test_rice_overflow(self):
        # rice.pdf(999, 0.74) was inf since special.i0 silentyly overflows
        # check that using i0e fixes it
        assert_(np.isfinite(stats.rice.pdf(999, 0.74)))

        assert_(np.isfinite(stats.rice.expect(lambda x: 1, args=(0.74,))))
        assert_(np.isfinite(stats.rice.expect(lambda x: 2, args=(0.74,))))
        assert_(np.isfinite(stats.rice.expect(lambda x: 3, args=(0.74,))))


class TestNct(TestCase):
    def test_nc_parameter(self):
        # Parameter values c<=0 were not enabled (gh-2402).
        # For negative values c and for c=0 results of rv.cdf(0) below were nan
        rv = stats.nct(5, 0)
        assert_equal(rv.cdf(0), 0.5)
        rv = stats.nct(5, -1)
        assert_almost_equal(rv.cdf(0), 0.841344746069, decimal=10)

    def test_broadcasting(self):
        res = stats.nct.pdf(5, np.arange(4,7)[:,None], np.linspace(0.1, 1, 4))
        expected = array([[0.00321886, 0.00557466, 0.00918418, 0.01442997],
                          [0.00217142, 0.00395366, 0.00683888, 0.01126276],
                          [0.00153078, 0.00291093, 0.00525206, 0.00900815]])
        assert_allclose(res, expected, rtol=1e-5)

    def text_variance_gh_issue_2401(self):
        # Computation of the variance of a non-central t-distribution resulted
        # in a TypeError: ufunc 'isinf' not supported for the input types,
        # and the inputs could not be safely coerced to any supported types
        # according to the casting rule 'safe'
        rv = stats.nct(4, 0)
        assert_equal(rv.var(), 2.0)

    def test_nct_inf_moments(self):
        # n-th moment of nct only exists for df > n
        m, v, s, k = stats.nct.stats(df=1.9, nc=0.3, moments='mvsk')
        assert_(np.isfinite(m))
        assert_equal([v, s, k], [np.inf, np.nan, np.nan])

        m, v, s, k = stats.nct.stats(df=3.1, nc=0.3, moments='mvsk')
        assert_(np.isfinite([m, v, s]).all())
        assert_equal(k, np.nan)


class TestRice(TestCase):
    def test_rice_zero_b(self):
        # rice distribution should work with b=0, cf gh-2164
        x = [0.2, 1., 5.]
        assert_(np.isfinite(stats.rice.pdf(x, b=0.)).all())
        assert_(np.isfinite(stats.rice.logpdf(x, b=0.)).all())
        assert_(np.isfinite(stats.rice.cdf(x, b=0.)).all())
        assert_(np.isfinite(stats.rice.logcdf(x, b=0.)).all())

        q = [0.1, 0.1, 0.5, 0.9]
        assert_(np.isfinite(stats.rice.ppf(q, b=0.)).all())

        mvsk = stats.rice.stats(0, moments='mvsk')
        assert_(np.isfinite(mvsk).all())

        # furthermore, pdf is continuous as b\to 0
        # rice.pdf(x, b\to 0) = x exp(-x^2/2) + O(b^2)
        # see e.g. Abramovich & Stegun 9.6.7 & 9.6.10
        b = 1e-8
        assert_allclose(stats.rice.pdf(x, 0), stats.rice.pdf(x, b),
                atol=b, rtol=0)

    def test_rice_rvs(self):
        rvs = stats.rice.rvs
        assert_equal(rvs(b=3.).size, 1)
        assert_equal(rvs(b=3., size=(3, 5)).shape, (3, 5))


class TestErlang(TestCase):
    def test_erlang_runtimewarning(self):
        # erlang should generate a RuntimeWarning if a non-integer
        # shape parameter is used.
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)

            # The non-integer shape parameter 1.3 should trigger a RuntimeWarning
            assert_raises(RuntimeWarning,
                              stats.erlang.rvs, 1.3, loc=0, scale=1, size=4)

            # Calling the fit method with `f0` set to an integer should
            # *not* trigger a RuntimeWarning.  It should return the same
            # values as gamma.fit(...).
            data = [0.5, 1.0, 2.0, 4.0]
            result_erlang = stats.erlang.fit(data, f0=1)
            result_gamma = stats.gamma.fit(data, f0=1)
            assert_allclose(result_erlang, result_gamma, rtol=1e-3)


class TestExponWeib(TestCase):

    def test_pdf_logpdf(self):
        # Regression test for gh-3508.
        x = 0.1
        a = 1.0
        c = 100.0
        p = stats.exponweib.pdf(x, a, c)
        logp = stats.exponweib.logpdf(x, a, c)
        # Expected values were computed with mpmath.
        assert_allclose([p, logp],
                        [1.0000000000000054e-97, -223.35075402042244])

    def test_a_is_1(self):
        # For issue gh-3508.
        # Check that when a=1, the pdf and logpdf methods of exponweib are the
        # same as those of weibull_min.
        x = np.logspace(-4, -1, 4)
        a = 1
        c = 100

        p = stats.exponweib.pdf(x, a, c)
        expected = stats.weibull_min.pdf(x, c)
        assert_allclose(p, expected)

        logp = stats.exponweib.logpdf(x, a, c)
        expected = stats.weibull_min.logpdf(x, c)
        assert_allclose(logp, expected)

    def test_a_is_1_c_is_1(self):
        # When a = 1 and c = 1, the distribution is exponential.
        x = np.logspace(-8, 1, 10)
        a = 1
        c = 1

        p = stats.exponweib.pdf(x, a, c)
        expected = stats.expon.pdf(x)
        assert_allclose(p, expected)

        logp = stats.exponweib.logpdf(x, a, c)
        expected = stats.expon.logpdf(x)
        assert_allclose(logp, expected)


class TestRdist(TestCase):
    @dec.slow
    def test_rdist_cdf_gh1285(self):
        # check workaround in rdist._cdf for issue gh-1285.
        distfn = stats.rdist
        values = [0.001, 0.5, 0.999]
        assert_almost_equal(distfn.cdf(distfn.ppf(values, 541.0), 541.0),
                                values, decimal=5)


def test_540_567():
    # test for nan returned in tickets 540, 567
    assert_almost_equal(stats.norm.cdf(-1.7624320982),0.03899815971089126,
                            decimal=10, err_msg='test_540_567')
    assert_almost_equal(stats.norm.cdf(-1.7624320983),0.038998159702449846,
                            decimal=10, err_msg='test_540_567')
    assert_almost_equal(stats.norm.cdf(1.38629436112, loc=0.950273420309,
                            scale=0.204423758009),0.98353464004309321,
                            decimal=10, err_msg='test_540_567')


def test_regression_ticket_1316():
    # The following was raising an exception, because _construct_default_doc()
    # did not handle the default keyword extradoc=None.  See ticket #1316.
    g = stats._continuous_distns.gamma_gen(name='gamma')


def test_regression_ticket_1326():
    # adjust to avoid nan with 0*log(0)
    assert_almost_equal(stats.chi2.pdf(0.0, 2), 0.5, 14)


def test_regression_tukey_lambda():
    # Make sure that Tukey-Lambda distribution correctly handles non-positive lambdas.
    x = np.linspace(-5.0, 5.0, 101)

    olderr = np.seterr(divide='ignore')
    try:
        for lam in [0.0, -1.0, -2.0, np.array([[-1.0], [0.0], [-2.0]])]:
            p = stats.tukeylambda.pdf(x, lam)
            assert_((p != 0.0).all())
            assert_(~np.isnan(p).all())

        lam = np.array([[-1.0], [0.0], [2.0]])
        p = stats.tukeylambda.pdf(x, lam)
    finally:
        np.seterr(**olderr)

    assert_(~np.isnan(p).all())
    assert_((p[0] != 0.0).all())
    assert_((p[1] != 0.0).all())
    assert_((p[2] != 0.0).any())
    assert_((p[2] == 0.0).any())


@dec.skipif(DOCSTRINGS_STRIPPED)
def test_regression_ticket_1421():
    assert_('pdf(x, mu, loc=0, scale=1)' not in stats.poisson.__doc__)
    assert_('pmf(x,' in stats.poisson.__doc__)


def test_nan_arguments_gh_issue_1362():
    with np.errstate(invalid='ignore'):
        assert_(np.isnan(stats.t.logcdf(1, np.nan)))
        assert_(np.isnan(stats.t.cdf(1, np.nan)))
        assert_(np.isnan(stats.t.logsf(1, np.nan)))
        assert_(np.isnan(stats.t.sf(1, np.nan)))
        assert_(np.isnan(stats.t.pdf(1, np.nan)))
        assert_(np.isnan(stats.t.logpdf(1, np.nan)))
        assert_(np.isnan(stats.t.ppf(1, np.nan)))
        assert_(np.isnan(stats.t.isf(1, np.nan)))

        assert_(np.isnan(stats.bernoulli.logcdf(np.nan, 0.5)))
        assert_(np.isnan(stats.bernoulli.cdf(np.nan, 0.5)))
        assert_(np.isnan(stats.bernoulli.logsf(np.nan, 0.5)))
        assert_(np.isnan(stats.bernoulli.sf(np.nan, 0.5)))
        assert_(np.isnan(stats.bernoulli.pmf(np.nan, 0.5)))
        assert_(np.isnan(stats.bernoulli.logpmf(np.nan, 0.5)))
        assert_(np.isnan(stats.bernoulli.ppf(np.nan, 0.5)))
        assert_(np.isnan(stats.bernoulli.isf(np.nan, 0.5)))


def test_frozen_fit_ticket_1536():
    np.random.seed(5678)
    true = np.array([0.25, 0., 0.5])
    x = stats.lognorm.rvs(true[0], true[1], true[2], size=100)

    olderr = np.seterr(divide='ignore')
    try:
        params = np.array(stats.lognorm.fit(x, floc=0.))
    finally:
        np.seterr(**olderr)

    assert_almost_equal(params, true, decimal=2)

    params = np.array(stats.lognorm.fit(x, fscale=0.5, loc=0))
    assert_almost_equal(params, true, decimal=2)

    params = np.array(stats.lognorm.fit(x, f0=0.25, loc=0))
    assert_almost_equal(params, true, decimal=2)

    params = np.array(stats.lognorm.fit(x, f0=0.25, floc=0))
    assert_almost_equal(params, true, decimal=2)

    np.random.seed(5678)
    loc = 1
    floc = 0.9
    x = stats.norm.rvs(loc, 2., size=100)
    params = np.array(stats.norm.fit(x, floc=floc))
    expected = np.array([floc, np.sqrt(((x-floc)**2).mean())])
    assert_almost_equal(params, expected, decimal=4)


def test_regression_ticket_1530():
    # Check the starting value works for Cauchy distribution fit.
    np.random.seed(654321)
    rvs = stats.cauchy.rvs(size=100)
    params = stats.cauchy.fit(rvs)
    expected = (0.045, 1.142)
    assert_almost_equal(params, expected, decimal=1)


def test_tukeylambda_stats_ticket_1545():
    # Some test for the variance and kurtosis of the Tukey Lambda distr.
    # See test_tukeylamdba_stats.py for more tests.

    mv = stats.tukeylambda.stats(0, moments='mvsk')
    # Known exact values:
    expected = [0, np.pi**2/3, 0, 1.2]
    assert_almost_equal(mv, expected, decimal=10)

    mv = stats.tukeylambda.stats(3.13, moments='mvsk')
    # 'expected' computed with mpmath.
    expected = [0, 0.0269220858861465102, 0, -0.898062386219224104]
    assert_almost_equal(mv, expected, decimal=10)

    mv = stats.tukeylambda.stats(0.14, moments='mvsk')
    # 'expected' computed with mpmath.
    expected = [0, 2.11029702221450250, 0, -0.02708377353223019456]
    assert_almost_equal(mv, expected, decimal=10)


def test_poisson_logpmf_ticket_1436():
    assert_(np.isfinite(stats.poisson.logpmf(1500, 200)))


def test_powerlaw_stats():
    """Test the powerlaw stats function.

    This unit test is also a regression test for ticket 1548.

    The exact values are:
    mean:
        mu = a / (a + 1)
    variance:
        sigma**2 = a / ((a + 2) * (a + 1) ** 2)
    skewness:
        One formula (see http://en.wikipedia.org/wiki/Skewness) is
            gamma_1 = (E[X**3] - 3*mu*E[X**2] + 2*mu**3) / sigma**3
        A short calculation shows that E[X**k] is a / (a + k), so gamma_1
        can be implemented as
            n = a/(a+3) - 3*(a/(a+1))*a/(a+2) + 2*(a/(a+1))**3
            d = sqrt(a/((a+2)*(a+1)**2)) ** 3
            gamma_1 = n/d
        Either by simplifying, or by a direct calculation of mu_3 / sigma**3,
        one gets the more concise formula:
            gamma_1 = -2.0 * ((a - 1) / (a + 3)) * sqrt((a + 2) / a)
    kurtosis: (See http://en.wikipedia.org/wiki/Kurtosis)
        The excess kurtosis is
            gamma_2 = mu_4 / sigma**4 - 3
        A bit of calculus and algebra (sympy helps) shows that
            mu_4 = 3*a*(3*a**2 - a + 2) / ((a+1)**4 * (a+2) * (a+3) * (a+4))
        so
            gamma_2 = 3*(3*a**2 - a + 2) * (a+2) / (a*(a+3)*(a+4)) - 3
        which can be rearranged to
            gamma_2 = 6 * (a**3 - a**2 - 6*a + 2) / (a*(a+3)*(a+4))
    """
    cases = [(1.0, (0.5, 1./12, 0.0, -1.2)),
             (2.0, (2./3, 2./36, -0.56568542494924734, -0.6))]
    for a, exact_mvsk in cases:
        mvsk = stats.powerlaw.stats(a, moments="mvsk")
        assert_array_almost_equal(mvsk, exact_mvsk)


def test_powerlaw_edge():
    # Regression test for gh-3986.
    p = stats.powerlaw.logpdf(0, 1)
    assert_equal(p, 0.0)


def test_exponpow_edge():
    # Regression test for gh-3982.
    p = stats.exponpow.logpdf(0, 1)
    assert_equal(p, 0.0)

    # Check pdf and logpdf at x = 0 for other values of b.
    p = stats.exponpow.pdf(0, [0.25, 1.0, 1.5])
    assert_equal(p, [np.inf, 1.0, 0.0])
    p = stats.exponpow.logpdf(0, [0.25, 1.0, 1.5])
    assert_equal(p, [np.inf, 0.0, -np.inf])


def test_gengamma_edge():
    # Regression test for gh-3985.
    p = stats.gengamma.pdf(0, 1, 1)
    assert_equal(p, 1.0)


def test_ksone_fit_freeze():
    # Regression test for ticket #1638.
    d = np.array(
        [-0.18879233, 0.15734249, 0.18695107, 0.27908787, -0.248649,
         -0.2171497, 0.12233512, 0.15126419, 0.03119282, 0.4365294,
          0.08930393, -0.23509903, 0.28231224, -0.09974875, -0.25196048,
          0.11102028, 0.1427649, 0.10176452, 0.18754054, 0.25826724,
          0.05988819, 0.0531668, 0.21906056, 0.32106729, 0.2117662,
          0.10886442, 0.09375789, 0.24583286, -0.22968366, -0.07842391,
         -0.31195432, -0.21271196, 0.1114243, -0.13293002, 0.01331725,
         -0.04330977, -0.09485776, -0.28434547, 0.22245721, -0.18518199,
         -0.10943985, -0.35243174, 0.06897665, -0.03553363, -0.0701746,
         -0.06037974, 0.37670779, -0.21684405])

    try:
        olderr = np.seterr(invalid='ignore')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            warnings.simplefilter('ignore', RuntimeWarning)
            stats.ksone.fit(d)
    finally:
        np.seterr(**olderr)


def test_norm_logcdf():
    # Test precision of the logcdf of the normal distribution.
    # This precision was enhanced in ticket 1614.
    x = -np.asarray(list(range(0, 120, 4)))
    # Values from R
    expected = [-0.69314718, -10.36010149, -35.01343716, -75.41067300,
                -131.69539607, -203.91715537, -292.09872100, -396.25241451,
                -516.38564863, -652.50322759, -804.60844201, -972.70364403,
                -1156.79057310, -1356.87055173, -1572.94460885, -1805.01356068,
                -2053.07806561, -2317.13866238, -2597.19579746, -2893.24984493,
                -3205.30112136, -3533.34989701, -3877.39640444, -4237.44084522,
                -4613.48339520, -5005.52420869, -5413.56342187, -5837.60115548,
                -6277.63751711, -6733.67260303]

    olderr = np.seterr(divide='ignore')
    try:
        assert_allclose(stats.norm().logcdf(x), expected, atol=1e-8)
    finally:
        np.seterr(**olderr)


def test_levy_cdf_ppf():
    # Test levy.cdf, including small arguments.
    x = np.array([1000, 1.0, 0.5, 0.1, 0.01, 0.001])

    # Expected values were calculated separately with mpmath.
    # E.g.
    # >>> mpmath.mp.dps = 100
    # >>> x = mpmath.mp.mpf('0.01')
    # >>> cdf = mpmath.erfc(mpmath.sqrt(1/(2*x)))
    expected = np.array([0.9747728793699604,
                         0.3173105078629141,
                         0.1572992070502851,
                         0.0015654022580025495,
                         1.523970604832105e-23,
                         1.795832784800726e-219])

    y = stats.levy.cdf(x)
    assert_allclose(y, expected, rtol=1e-10)

    # ppf(expected) should get us back to x.
    xx = stats.levy.ppf(expected)
    assert_allclose(xx, x, rtol=1e-13)


def test_hypergeom_interval_1802():
    # these two had endless loops
    assert_equal(stats.hypergeom.interval(.95, 187601, 43192, 757),
                 (152.0, 197.0))
    assert_equal(stats.hypergeom.interval(.945, 187601, 43192, 757),
                 (152.0, 197.0))
    # this was working also before
    assert_equal(stats.hypergeom.interval(.94, 187601, 43192, 757),
                 (153.0, 196.0))

    # degenerate case .a == .b
    assert_equal(stats.hypergeom.ppf(0.02, 100, 100, 8), 8)
    assert_equal(stats.hypergeom.ppf(1, 100, 100, 8), 8)


def test_distribution_too_many_args():
    # Check that a TypeError is raised when too many args are given to a method
    # Regression test for ticket 1815.
    x = np.linspace(0.1, 0.7, num=5)
    assert_raises(TypeError, stats.gamma.pdf, x, 2, 3, loc=1.0)
    assert_raises(TypeError, stats.gamma.pdf, x, 2, 3, 4, loc=1.0)
    assert_raises(TypeError, stats.gamma.pdf, x, 2, 3, 4, 5)
    assert_raises(TypeError, stats.gamma.pdf, x, 2, 3, loc=1.0, scale=0.5)
    assert_raises(TypeError, stats.gamma.rvs, 2., 3, loc=1.0, scale=0.5)
    assert_raises(TypeError, stats.gamma.cdf, x, 2., 3, loc=1.0, scale=0.5)
    assert_raises(TypeError, stats.gamma.ppf, x, 2., 3, loc=1.0, scale=0.5)
    assert_raises(TypeError, stats.gamma.stats, 2., 3, loc=1.0, scale=0.5)
    assert_raises(TypeError, stats.gamma.entropy, 2., 3, loc=1.0, scale=0.5)
    assert_raises(TypeError, stats.gamma.fit, x, 2., 3, loc=1.0, scale=0.5)

    # These should not give errors
    stats.gamma.pdf(x, 2, 3)  # loc=3
    stats.gamma.pdf(x, 2, 3, 4)  # loc=3, scale=4
    stats.gamma.stats(2., 3)
    stats.gamma.stats(2., 3, 4)
    stats.gamma.stats(2., 3, 4, 'mv')
    stats.gamma.rvs(2., 3, 4, 5)
    stats.gamma.fit(stats.gamma.rvs(2., size=7), 2.)

    # Also for a discrete distribution
    stats.geom.pmf(x, 2, loc=3)  # no error, loc=3
    assert_raises(TypeError, stats.geom.pmf, x, 2, 3, 4)
    assert_raises(TypeError, stats.geom.pmf, x, 2, 3, loc=4)

    # And for distributions with 0, 2 and 3 args respectively
    assert_raises(TypeError, stats.expon.pdf, x, 3, loc=1.0)
    assert_raises(TypeError, stats.exponweib.pdf, x, 3, 4, 5, loc=1.0)
    assert_raises(TypeError, stats.exponweib.pdf, x, 3, 4, 5, 0.1, 0.1)
    assert_raises(TypeError, stats.ncf.pdf, x, 3, 4, 5, 6, loc=1.0)
    assert_raises(TypeError, stats.ncf.pdf, x, 3, 4, 5, 6, 1.0, scale=0.5)
    stats.ncf.pdf(x, 3, 4, 5, 6, 1.0)  # 3 args, plus loc/scale


def test_ncx2_tails_ticket_955():
    # Trac #955 -- check that the cdf computed by special functions
    # matches the integrated pdf
    a = stats.ncx2.cdf(np.arange(20, 25, 0.2), 2, 1.07458615e+02)
    b = stats.ncx2._cdfvec(np.arange(20, 25, 0.2), 2, 1.07458615e+02)
    assert_allclose(a, b, rtol=1e-3, atol=0)


def test_foldnorm_zero():
    # Parameter value c=0 was not enabled, see gh-2399.
    rv = stats.foldnorm(0, scale=1)
    assert_equal(rv.cdf(0), 0)  # rv.cdf(0) previously resulted in: nan


def test_stats_shapes_argcheck():
    # stats method was failing for vector shapes if some of the values
    # were outside of the allowed range, see gh-2678
    mv3 = stats.invgamma.stats([0.0, 0.5, 1.0], 1, 0.5)  # 0 is not a legal `a`
    mv2 = stats.invgamma.stats([0.5, 1.0], 1, 0.5)
    mv2_augmented = tuple(np.r_[np.nan, _] for _ in mv2)
    assert_equal(mv2_augmented, mv3)

    mv3 = stats.lognorm.stats([2, 2.4, -1])  # -1 is not a legal shape parameter
    mv2 = stats.lognorm.stats([2, 2.4])
    mv2_augmented = tuple(np.r_[_, np.nan] for _ in mv2)
    assert_equal(mv2_augmented, mv3)

    # FIXME: this is only a quick-and-dirty test of a quick-and-dirty bugfix.
    # stats method with multiple shape parameters is not properly vectorized
    # anyway, so some distributions may or may not fail.


## Test subclassing distributions w/ explicit shapes

class _distr_gen(stats.rv_continuous):
    def _pdf(self, x, a):
        return 42


class _distr2_gen(stats.rv_continuous):
    def _cdf(self, x, a):
        return 42 * a + x


class _distr3_gen(stats.rv_continuous):
    def _pdf(self, x, a, b):
        return a + b

    def _cdf(self, x, a):
        # Different # of shape params from _pdf, to be able to check that
        # inspection catches the inconsistency."""
        return 42 * a + x


class _distr6_gen(stats.rv_continuous):
    # Two shape parameters (both _pdf and _cdf defined, consistent shapes.)
    def _pdf(self, x, a, b):
        return a*x + b

    def _cdf(self, x, a, b):
        return 42 * a + x


class TestSubclassingExplicitShapes(TestCase):
    # Construct a distribution w/ explicit shapes parameter and test it.

    def test_correct_shapes(self):
        dummy_distr = _distr_gen(name='dummy', shapes='a')
        assert_equal(dummy_distr.pdf(1, a=1), 42)

    def test_wrong_shapes_1(self):
        dummy_distr = _distr_gen(name='dummy', shapes='A')
        assert_raises(TypeError, dummy_distr.pdf, 1, **dict(a=1))

    def test_wrong_shapes_2(self):
        dummy_distr = _distr_gen(name='dummy', shapes='a, b, c')
        dct = dict(a=1, b=2, c=3)
        assert_raises(TypeError, dummy_distr.pdf, 1, **dct)

    def test_shapes_string(self):
        # shapes must be a string
        dct = dict(name='dummy', shapes=42)
        assert_raises(TypeError, _distr_gen, **dct)

    def test_shapes_identifiers_1(self):
        # shapes must be a comma-separated list of valid python identifiers
        dct = dict(name='dummy', shapes='(!)')
        assert_raises(SyntaxError, _distr_gen, **dct)

    def test_shapes_identifiers_2(self):
        dct = dict(name='dummy', shapes='4chan')
        assert_raises(SyntaxError, _distr_gen, **dct)

    def test_shapes_identifiers_3(self):
        dct = dict(name='dummy', shapes='m(fti)')
        assert_raises(SyntaxError, _distr_gen, **dct)

    def test_shapes_identifiers_nodefaults(self):
        dct = dict(name='dummy', shapes='a=2')
        assert_raises(SyntaxError, _distr_gen, **dct)

    def test_shapes_args(self):
        dct = dict(name='dummy', shapes='*args')
        assert_raises(SyntaxError, _distr_gen, **dct)

    def test_shapes_kwargs(self):
        dct = dict(name='dummy', shapes='**kwargs')
        assert_raises(SyntaxError, _distr_gen, **dct)

    def test_shapes_keywords(self):
        # python keywords cannot be used for shape parameters
        dct = dict(name='dummy', shapes='a, b, c, lambda')
        assert_raises(SyntaxError, _distr_gen, **dct)

    def test_shapes_signature(self):
        # test explicit shapes which agree w/ the signature of _pdf
        class _dist_gen(stats.rv_continuous):
            def _pdf(self, x, a):
                return stats.norm._pdf(x) * a

        dist = _dist_gen(shapes='a')
        assert_equal(dist.pdf(0.5, a=2), stats.norm.pdf(0.5)*2)

    def test_shapes_signature_inconsistent(self):
        # test explicit shapes which do not agree w/ the signature of _pdf
        class _dist_gen(stats.rv_continuous):
            def _pdf(self, x, a):
                return stats.norm._pdf(x) * a

        dist = _dist_gen(shapes='a, b')
        assert_raises(TypeError, dist.pdf, 0.5, **dict(a=1, b=2))

    def test_star_args(self):
        # test _pdf with only starargs
        # NB: **kwargs of pdf will never reach _pdf
        class _dist_gen(stats.rv_continuous):
            def _pdf(self, x, *args):
                extra_kwarg = args[0]
                return stats.norm._pdf(x) * extra_kwarg

        dist = _dist_gen(shapes='extra_kwarg')
        assert_equal(dist.pdf(0.5, extra_kwarg=33), stats.norm.pdf(0.5)*33)
        assert_equal(dist.pdf(0.5, 33), stats.norm.pdf(0.5)*33)
        assert_raises(TypeError, dist.pdf, 0.5, **dict(xxx=33))

    def test_star_args_2(self):
        # test _pdf with named & starargs
        # NB: **kwargs of pdf will never reach _pdf
        class _dist_gen(stats.rv_continuous):
            def _pdf(self, x, offset, *args):
                extra_kwarg = args[0]
                return stats.norm._pdf(x) * extra_kwarg + offset

        dist = _dist_gen(shapes='offset, extra_kwarg')
        assert_equal(dist.pdf(0.5, offset=111, extra_kwarg=33),
                     stats.norm.pdf(0.5)*33 + 111)
        assert_equal(dist.pdf(0.5, 111, 33),
                     stats.norm.pdf(0.5)*33 + 111)

    def test_extra_kwarg(self):
        # **kwargs to _pdf are ignored.
        # this is a limitation of the framework (_pdf(x, *goodargs))
        class _distr_gen(stats.rv_continuous):
            def _pdf(self, x, *args, **kwargs):
                # _pdf should handle *args, **kwargs itself.  Here "handling" is
                # ignoring *args and looking for ``extra_kwarg`` and using that.
                extra_kwarg = kwargs.pop('extra_kwarg', 1)
                return stats.norm._pdf(x) * extra_kwarg

        dist = _distr_gen(shapes='extra_kwarg')
        assert_equal(dist.pdf(1, extra_kwarg=3), stats.norm.pdf(1))

    def shapes_empty_string(self):
        # shapes='' is equivalent to shapes=None
        class _dist_gen(stats.rv_continuous):
            def _pdf(self, x):
                return stats.norm.pdf(x)

        dist = _dist_gen(shapes='')
        assert_equal(dist.pdf(0.5), stats.norm.pdf(0.5))


class TestSubclassingNoShapes(TestCase):
    # Construct a distribution w/o explicit shapes parameter and test it.

    def test_only__pdf(self):
        dummy_distr = _distr_gen(name='dummy')
        assert_equal(dummy_distr.pdf(1, a=1), 42)

    def test_only__cdf(self):
        # _pdf is determined from _cdf by taking numerical derivative
        dummy_distr = _distr2_gen(name='dummy')
        assert_almost_equal(dummy_distr.pdf(1, a=1), 1)

    @dec.skipif(DOCSTRINGS_STRIPPED)
    def test_signature_inspection(self):
        # check that _pdf signature inspection works correctly, and is used in
        # the class docstring
        dummy_distr = _distr_gen(name='dummy')
        assert_equal(dummy_distr.numargs, 1)
        assert_equal(dummy_distr.shapes, 'a')
        res = re.findall('logpdf\(x, a, loc=0, scale=1\)',
                         dummy_distr.__doc__)
        assert_(len(res) == 1)

    @dec.skipif(DOCSTRINGS_STRIPPED)
    def test_signature_inspection_2args(self):
        # same for 2 shape params and both _pdf and _cdf defined
        dummy_distr = _distr6_gen(name='dummy')
        assert_equal(dummy_distr.numargs, 2)
        assert_equal(dummy_distr.shapes, 'a, b')
        res = re.findall('logpdf\(x, a, b, loc=0, scale=1\)',
                         dummy_distr.__doc__)
        assert_(len(res) == 1)

    def test_signature_inspection_2args_incorrect_shapes(self):
        # both _pdf and _cdf defined, but shapes are inconsistent: raises
        try:
            _distr3_gen(name='dummy')
        except TypeError:
            pass
        else:
            raise AssertionError('TypeError not raised.')

    def test_defaults_raise(self):
        # default arguments should raise
        class _dist_gen(stats.rv_continuous):
            def _pdf(self, x, a=42):
                return 42
        assert_raises(TypeError, _dist_gen, **dict(name='dummy'))

    def test_starargs_raise(self):
        # without explicit shapes, *args are not allowed
        class _dist_gen(stats.rv_continuous):
            def _pdf(self, x, a, *args):
                return 42
        assert_raises(TypeError, _dist_gen, **dict(name='dummy'))

    def test_kwargs_raise(self):
        # without explicit shapes, **kwargs are not allowed
        class _dist_gen(stats.rv_continuous):
            def _pdf(self, x, a, **kwargs):
                return 42
        assert_raises(TypeError, _dist_gen, **dict(name='dummy'))


@dec.skipif(DOCSTRINGS_STRIPPED)
def test_docstrings():
    badones = [',\s*,', '\(\s*,', '^\s*:']
    for distname in stats.__all__:
        dist = getattr(stats, distname)
        if isinstance(dist, (stats.rv_discrete, stats.rv_continuous)):
            for regex in badones:
                assert_(re.search(regex, dist.__doc__) is None)


def test_infinite_input():
    assert_almost_equal(stats.skellam.sf(np.inf, 10, 11), 0)
    assert_almost_equal(stats.ncx2._cdf(np.inf, 8, 0.1), 1)


if __name__ == "__main__":
    run_module_suite()
