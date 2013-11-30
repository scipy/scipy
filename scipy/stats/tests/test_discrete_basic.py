from __future__ import division, print_function, absolute_import

import inspect

import numpy.testing as npt
import numpy as np
from scipy.lib.six import xrange

from scipy import stats
from common_tests import (check_normalization, check_moment, check_mean_expect,
        check_var_expect, check_skew_expect, check_kurt_expect,
        check_entropy, check_private_entropy, check_edge_support)

DECIMAL_meanvar = 0  # 1  # was 0

distdiscrete = [
    ['bernoulli',(0.3,)],
    ['binom', (5, 0.4)],
    ['boltzmann',(1.4, 19)],
    ['dlaplace', (0.8,)],  # 0.5
    ['geom', (0.5,)],
    ['hypergeom',(30, 12, 6)],
    ['hypergeom',(21,3,12)],  # numpy.random (3,18,12) numpy ticket:921
    ['hypergeom',(21,18,11)],  # numpy.random (18,3,11) numpy ticket:921
    ['logser', (0.6,)],  # reenabled, numpy ticket:921
    ['nbinom', (5, 0.5)],
    ['nbinom', (0.4, 0.4)],  # from tickets: 583
    ['planck', (0.51,)],   # 4.1
    ['poisson', (0.6,)],
    ['randint', (7, 31)],
    ['skellam', (15, 8)],
    ['zipf', (6.5,)]
]


def test_discrete_basic():
    for distname, arg in distdiscrete:
        distfn = getattr(stats,distname)
        #npt.assert_(stats.dlaplace.rvs(0.8) is not None)
        np.random.seed(9765456)
        rvs = distfn.rvs(size=2000,*arg)
        supp = np.unique(rvs)
        m, v = distfn.stats(*arg)
        # yield npt.assert_almost_equal(rvs.mean(), m, decimal=4,err_msg='mean')
        # yield npt.assert_almost_equal, rvs.mean(), m, 2, 'mean' # does not work
        yield check_sample_meanvar, rvs.mean(), m, distname + ' sample mean test'
        yield check_sample_meanvar, rvs.var(), v, distname + ' sample var test'
        yield check_cdf_ppf, distfn, arg, distname + ' cdf_ppf'
        yield check_cdf_ppf2, distfn, arg, supp, distname + ' cdf_ppf'
        yield check_pmf_cdf, distfn, arg, distname + ' pmf_cdf'

        yield check_oth, distfn, arg, distname + ' oth'
        skurt = stats.kurtosis(rvs)
        sskew = stats.skew(rvs)
        yield check_sample_skew_kurt, distfn, arg, skurt, sskew, \
                      distname + ' skew_kurt'

        alpha = 0.01
        yield check_discrete_chisquare, distfn, arg, rvs, alpha, \
                      distname + ' chisquare'

        yield check_edge_support, distfn, arg

    seen = set()
    for distname, arg in distdiscrete:
        if distname in seen:
            continue
        seen.add(distname)
        distfn = getattr(stats,distname)
        locscale_defaults = (0,)
        meths = [distfn.pmf, distfn.logpmf, distfn.cdf, distfn.logcdf,
                 distfn.logsf]
        # make sure arguments are within support
        spec_k = {'randint': 11, 'hypergeom': 4, 'bernoulli': 0, }
        k = spec_k.get(distname, 1)
        yield check_named_args, distfn, k, arg, locscale_defaults, meths
        yield check_scale_docstring, distfn

        # Entropy
        yield check_entropy, distfn, arg, distname
        if distfn.__class__._entropy != stats.rv_discrete._entropy:
            yield check_private_entropy, distfn, arg, stats.rv_discrete


def test_moments():
    knf = npt.dec.knownfailureif
    for distname, arg in distdiscrete:
        distfn = getattr(stats,distname)
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


@npt.dec.skipif(True)
def test_discrete_private():
    # testing private methods mostly for debugging
    #   some tests might fail by design,
    #   e.g. incorrect definition of distfn.a and distfn.b
    for distname, arg in distdiscrete:
        distfn = getattr(stats,distname)
        rvs = distfn.rvs(size=10000,*arg)
        m,v = distfn.stats(*arg)

        yield check_ppf_ppf, distfn, arg
        yield check_cdf_ppf_private, distfn, arg, distname
        yield check_generic_moment, distfn, arg, m, 1, 3   # last is decimal
        yield check_generic_moment, distfn, arg, v+m*m, 2, 3  # last is decimal
        yield check_moment_frozen, distfn, arg, m, 1, 3   # last is decimal
        yield check_moment_frozen, distfn, arg, v+m*m, 2, 3  # last is decimal


def check_sample_meanvar(sm,m,msg):
    if not np.isinf(m):
        npt.assert_almost_equal(sm, m, decimal=DECIMAL_meanvar, err_msg=msg +
                                ' - finite moment')
    else:
        npt.assert_(sm > 10000, msg='infinite moment, sm = ' + str(sm))


def check_sample_var(sm,m,msg):
    npt.assert_almost_equal(sm, m, decimal=DECIMAL_meanvar, err_msg=msg + 'var')


def check_cdf_ppf(distfn,arg,msg):
    ppf05 = distfn.ppf(0.5,*arg)
    cdf05 = distfn.cdf(ppf05,*arg)
    npt.assert_almost_equal(distfn.ppf(cdf05-1e-6,*arg),ppf05,
                            err_msg=msg + 'ppf-cdf-median')
    npt.assert_((distfn.ppf(cdf05+1e-4,*arg) > ppf05), msg + 'ppf-cdf-next')


def check_cdf_ppf2(distfn,arg,supp,msg):
    npt.assert_array_equal(distfn.ppf(distfn.cdf(supp,*arg),*arg),
                           supp, msg + '-roundtrip')
    npt.assert_array_equal(distfn.ppf(distfn.cdf(supp,*arg)-1e-8,*arg),
                           supp, msg + '-roundtrip')
    # -1e-8 could cause an error if pmf < 1e-8


def check_cdf_ppf_private(distfn,arg,msg):
    ppf05 = distfn._ppf(0.5,*arg)
    cdf05 = distfn.cdf(ppf05,*arg)
    npt.assert_almost_equal(distfn._ppf(cdf05-1e-6,*arg),ppf05,
                            err_msg=msg + '_ppf-cdf-median ')
    npt.assert_((distfn._ppf(cdf05+1e-4,*arg) > ppf05), msg + '_ppf-cdf-next')


def check_ppf_ppf(distfn, arg):
    npt.assert_(distfn.ppf(0.5,*arg) < np.inf)
    ppfs = distfn.ppf([0.5,0.9],*arg)
    ppf_s = [distfn._ppf(0.5,*arg), distfn._ppf(0.9,*arg)]
    npt.assert_(np.all(ppfs < np.inf))
    npt.assert_(ppf_s[0] == distfn.ppf(0.5,*arg))
    npt.assert_(ppf_s[1] == distfn.ppf(0.9,*arg))
    npt.assert_(ppf_s[0] == ppfs[0])
    npt.assert_(ppf_s[1] == ppfs[1])


def check_pmf_cdf(distfn, arg, msg):
    startind = np.int(distfn._ppf(0.01,*arg)-1)
    index = list(range(startind,startind+10))
    cdfs = distfn.cdf(index,*arg)
    npt.assert_almost_equal(cdfs, distfn.pmf(index, *arg).cumsum() +
                            cdfs[0] - distfn.pmf(index[0],*arg),
                            decimal=4, err_msg=msg + 'pmf-cdf')


def check_generic_moment(distfn, arg, m, k, decim):
    npt.assert_almost_equal(distfn.generic_moment(k,*arg), m, decimal=decim,
                            err_msg=str(distfn) + ' generic moment test')


def check_moment_frozen(distfn, arg, m, k, decim):
    npt.assert_almost_equal(distfn(*arg).moment(k), m, decimal=decim,
                            err_msg=str(distfn) + ' frozen moment test')


def check_oth(distfn, arg, msg):
    # checking other methods of distfn
    meanint = round(float(distfn.stats(*arg)[0]))  # closest integer to mean
    npt.assert_almost_equal(distfn.sf(meanint, *arg), 1 -
                            distfn.cdf(meanint, *arg), decimal=8)
    median_sf = distfn.isf(0.5, *arg)

    npt.assert_(distfn.sf(median_sf - 1, *arg) > 0.5)
    npt.assert_(distfn.cdf(median_sf + 1, *arg) > 0.5)
    npt.assert_equal(distfn.isf(0.5, *arg), distfn.ppf(0.5, *arg))


def check_sample_skew_kurt(distfn, arg, sk, ss, msg):
    k,s = distfn.stats(moments='ks', *arg)
    check_sample_meanvar, sk, k, msg + 'sample skew test'
    check_sample_meanvar, ss, s, msg + 'sample kurtosis test'



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
        current = distfn.cdf(ii,*arg)
        if current - last >= wsupp-1e-14:
            distsupp.append(ii)
            distmass.append(current - last)
            last = current
            if current > (1-wsupp):
                break
    if distsupp[-1] < distfn.b:
        distsupp.append(distfn.b)
        distmass.append(1-last)
    distsupp = np.array(distsupp)
    distmass = np.array(distmass)

    # convert intervals to right-half-open as required by histogram
    histsupp = distsupp+1e-8
    histsupp[0] = distfn.a

    # find sample frequencies and perform chisquare test
    freq,hsupp = np.histogram(rvs,histsupp)
    cdfs = distfn.cdf(distsupp,*arg)
    (chis,pval) = stats.chisquare(np.array(freq),n*distmass)

    npt.assert_(pval > alpha, 'chisquare - test for %s'
           ' at arg = %s with pval = %s' % (msg,str(arg),str(pval)))


def check_named_args(distfn, x, shape_args, defaults, meths):
    """Check calling w/ named arguments."""
    ### This is a copy-paste from test_continuous_basic.py; dedupe?

    # check consistency of shapes, numargs and _parse signature
    signature = inspect.getargspec(distfn._parse_args)
    npt.assert_(signature.varargs is None)
    npt.assert_(signature.keywords is None)
    npt.assert_(signature.defaults == defaults)

    shape_argnames = signature.args[1:-len(defaults)]  # self, a, b, loc=0, scale=1
    if distfn.shapes:
        shapes_ = distfn.shapes.replace(',', ' ').split()
    else:
        shapes_ = ''
    npt.assert_(len(shapes_) == distfn.numargs)
    npt.assert_(len(shapes_) == len(shape_argnames))

    # check calling w/ named arguments
    shape_args = list(shape_args)

    vals = [meth(x, *shape_args) for meth in meths]
    npt.assert_(np.all(np.isfinite(vals)))

    names, a, k = shape_argnames[:], shape_args[:], {}
    while names:
        k.update({names.pop(): a.pop()})
        v = [meth(x, *a, **k) for meth in meths]
        npt.assert_array_equal(vals, v)
        if not 'n' in k.keys():
            # `n` is first parameter of moment(), so can't be used as named arg
            npt.assert_equal(distfn.moment(3, *a, **k),
                             distfn.moment(3, *shape_args))

    # unknown arguments should not go through:
    k.update({'kaboom': 42})
    npt.assert_raises(TypeError, distfn.cdf, x, **k)


def check_scale_docstring(distfn):
    if distfn.__doc__ is not None:
        # Docstrings can be stripped if interpreter is run with -OO
        npt.assert_('scale' not in distfn.__doc__)


if __name__ == "__main__":
    npt.run_module_suite()
