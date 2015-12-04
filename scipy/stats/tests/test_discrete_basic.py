from __future__ import division, print_function, absolute_import

import numpy.testing as npt
import numpy as np
from scipy._lib.six import xrange

from scipy import stats
from common_tests import (check_normalization, check_moment, check_mean_expect,
                          check_var_expect, check_skew_expect,
                          check_kurt_expect, check_entropy,
                          check_private_entropy, check_edge_support,
                          check_named_args, check_random_state_property,
                          check_pickling)
from scipy.stats._distr_params import distdiscrete
knf = npt.dec.knownfailureif

vals = ([1, 2, 3, 4], [0.1, 0.2, 0.3, 0.4])
distdiscrete += [[stats.rv_discrete(values=vals), ()]]


def test_discrete_basic():
    for distname, arg in distdiscrete:
        try:
            distfn = getattr(stats, distname)
        except TypeError:
            distfn = distname
            distname = 'sample distribution'
        np.random.seed(9765456)
        rvs = distfn.rvs(size=2000, *arg)
        supp = np.unique(rvs)
        m, v = distfn.stats(*arg)
        yield check_cdf_ppf, distfn, arg, supp, distname + ' cdf_ppf'

        yield check_pmf_cdf, distfn, arg, distname
        yield check_oth, distfn, arg, supp, distname + ' oth'
        yield check_edge_support, distfn, arg

        alpha = 0.01
        yield (check_discrete_chisquare, distfn, arg, rvs, alpha,
               distname + ' chisquare')

    seen = set()
    for distname, arg in distdiscrete:
        if distname in seen:
            continue
        seen.add(distname)
        try:
            distfn = getattr(stats, distname)
        except TypeError:
            distfn = distname
            distname = 'sample distribution'
        locscale_defaults = (0,)
        meths = [distfn.pmf, distfn.logpmf, distfn.cdf, distfn.logcdf,
                 distfn.logsf]
        # make sure arguments are within support
        spec_k = {'randint': 11, 'hypergeom': 4, 'bernoulli': 0, }
        k = spec_k.get(distname, 1)
        yield check_named_args, distfn, k, arg, locscale_defaults, meths
        if distname != 'sample distribution':
            yield check_scale_docstring, distfn
        yield check_random_state_property, distfn, arg
        yield check_pickling, distfn, arg

        # Entropy
        yield check_entropy, distfn, arg, distname
        if distfn.__class__._entropy != stats.rv_discrete._entropy:
            yield check_private_entropy, distfn, arg, stats.rv_discrete


def test_moments():
    for distname, arg in distdiscrete:
        try:
            distfn = getattr(stats, distname)
        except TypeError:
            distfn = distname
            distname = 'sample distribution'
        m, v, s, k = distfn.stats(*arg, moments='mvsk')
        yield check_normalization, distfn, arg, distname

        # compare `stats` and `moment` methods
        yield check_moment, distfn, arg, m, v, distname
        yield check_mean_expect, distfn, arg, m, distname
        yield check_var_expect, distfn, arg, m, v, distname
        yield check_skew_expect, distfn, arg, m, v, s, distname

        cond = distname in ['zipf']
        msg = distname + ' fails kurtosis'
        yield knf(cond, msg)(check_kurt_expect), distfn, arg, m, v, k, distname

        # frozen distr moments
        yield check_moment_frozen, distfn, arg, m, 1
        yield check_moment_frozen, distfn, arg, v+m*m, 2


def check_cdf_ppf(distfn, arg, supp, msg):
    # cdf is a step function, and ppf(q) = min{k : cdf(k) >= q, k integer}
    npt.assert_array_equal(distfn.ppf(distfn.cdf(supp, *arg), *arg),
                           supp, msg + '-roundtrip')
    npt.assert_array_equal(distfn.ppf(distfn.cdf(supp, *arg) - 1e-8, *arg),
                           supp, msg + '-roundtrip')
    supp1 = supp[supp < distfn.b]
    npt.assert_array_equal(distfn.ppf(distfn.cdf(supp1, *arg) + 1e-8, *arg),
                           supp1 + distfn.inc, msg + 'ppf-cdf-next')
    # -1e-8 could cause an error if pmf < 1e-8


def check_pmf_cdf(distfn, arg, distname):
    startind = int(distfn.ppf(0.01, *arg) - 1)
    index = list(range(startind, startind + 10))
    cdfs = distfn.cdf(index, *arg)
    pmfs_cum = distfn.pmf(index, *arg).cumsum()

    atol, rtol = 1e-10, 1e-10
    if distname == 'skellam':    # ncx2 accuracy
        atol, rtol = 1e-5, 1e-5
    npt.assert_allclose(cdfs - cdfs[0], pmfs_cum - pmfs_cum[0],
                        atol=atol, rtol=rtol)


def check_moment_frozen(distfn, arg, m, k):
    npt.assert_allclose(distfn(*arg).moment(k), m,
                        atol=1e-10, rtol=1e-10)


def check_oth(distfn, arg, supp, msg):
    # checking other methods of distfn
    npt.assert_allclose(distfn.sf(supp, *arg), 1. - distfn.cdf(supp, *arg),
                        atol=1e-10, rtol=1e-10)

    q = np.linspace(0.01, 0.99, 20)
    npt.assert_allclose(distfn.isf(q, *arg), distfn.ppf(1. - q, *arg),
                        atol=1e-10, rtol=1e-10)

    median_sf = distfn.isf(0.5, *arg)
    npt.assert_(distfn.sf(median_sf - 1, *arg) > 0.5)
    npt.assert_(distfn.cdf(median_sf + 1, *arg) > 0.5)


def check_discrete_chisquare(distfn, arg, rvs, alpha, msg):
    """Perform chisquare test for random sample of a discrete distribution

    Parameters
    ----------
    distname : string
        name of distribution function
    arg : sequence
        parameters of distribution
    alpha : float
        significance level, threshold for p-value

    Returns
    -------
    result : bool
        0 if test passes, 1 if test fails

    uses global variable debug for printing results

    """
    n = len(rvs)
    nsupp = 20
    wsupp = 1.0/nsupp

    # construct intervals with minimum mass 1/nsupp
    # intervals are left-half-open as in a cdf difference
    distsupport = xrange(max(distfn.a, -1000), min(distfn.b, 1000) + 1)
    last = 0
    distsupp = [max(distfn.a, -1000)]
    distmass = []
    for ii in distsupport:
        current = distfn.cdf(ii, *arg)
        if current - last >= wsupp - 1e-14:
            distsupp.append(ii)
            distmass.append(current - last)
            last = current
            if current > (1 - wsupp):
                break
    if distsupp[-1] < distfn.b:
        distsupp.append(distfn.b)
        distmass.append(1 - last)
    distsupp = np.array(distsupp)
    distmass = np.array(distmass)

    # convert intervals to right-half-open as required by histogram
    histsupp = distsupp + 1e-8
    histsupp[0] = distfn.a

    # find sample frequencies and perform chisquare test
    freq, hsupp = np.histogram(rvs, histsupp)
    cdfs = distfn.cdf(distsupp, *arg)
    chis, pval = stats.chisquare(np.array(freq), n*distmass)

    npt.assert_(pval > alpha,
                'chisquare - test for %s at arg = %s with pval = %s' %
                (msg, str(arg), str(pval)))


def check_scale_docstring(distfn):
    if distfn.__doc__ is not None:
        # Docstrings can be stripped if interpreter is run with -OO
        npt.assert_('scale' not in distfn.__doc__)


if __name__ == "__main__":
    npt.run_module_suite()
