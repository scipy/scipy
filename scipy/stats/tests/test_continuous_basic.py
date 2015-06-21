from __future__ import division, print_function, absolute_import

import warnings

import numpy as np
import numpy.testing as npt

from scipy import integrate
from scipy import stats
from scipy.special import betainc
from common_tests import (check_normalization, check_moment, check_mean_expect,
        check_var_expect, check_skew_expect, check_kurt_expect,
        check_entropy, check_private_entropy, NUMPY_BELOW_1_7,
        check_edge_support, check_named_args, check_random_state_property,
        check_meth_dtype, check_ppf_dtype)

from scipy.stats._distr_params import distcont

"""
Test all continuous distributions.

Parameters were chosen for those distributions that pass the
Kolmogorov-Smirnov test.  This provides safe parameters for each
distributions so that we can perform further testing of class methods.

These tests currently check only/mostly for serious errors and exceptions,
not for numerically exact results.
"""

## Note that you need to add new distributions you want tested
## to _distr_params

DECIMAL = 5  # specify the precision of the tests  # increased from 0 to 5

## Last four of these fail all around. Need to be checked
distcont_extra = [
    ['betaprime', (100, 86)],
    ['fatiguelife', (5,)],
    ['mielke', (4.6420495492121487, 0.59707419545516938)],
    ['invweibull', (0.58847112119264788,)],
    # burr: sample mean test fails still for c<1
    ['burr', (0.94839838075366045, 4.3820284068855795)],
    # genextreme: sample mean test, sf-logsf test fail
    ['genextreme', (3.3184017469423535,)],
]


# for testing only specific functions
# distcont = [
##    ['fatiguelife', (29,)],   #correction numargs = 1
##    ['loggamma', (0.41411931826052117,)]]

# for testing ticket:767
# distcont = [
##    ['genextreme', (3.3184017469423535,)],
##    ['genextreme', (0.01,)],
##    ['genextreme', (0.00001,)],
##    ['genextreme', (0.0,)],
##    ['genextreme', (-0.01,)]
##    ]

# distcont = [['gumbel_l', ()],
##            ['gumbel_r', ()],
##            ['norm', ()]
##            ]

# distcont = [['norm', ()]]

distmissing = ['wald', 'gausshyper', 'genexpon', 'rv_continuous',
    'loglaplace', 'rdist', 'semicircular', 'invweibull', 'ksone',
    'cosine', 'kstwobign', 'truncnorm', 'mielke', 'recipinvgauss', 'levy',
    'johnsonsu', 'levy_l', 'powernorm', 'wrapcauchy',
    'johnsonsb', 'truncexpon', 'invgauss', 'invgamma',
    'powerlognorm']

distmiss = [[dist,args] for dist,args in distcont if dist in distmissing]
distslow = ['rdist', 'gausshyper', 'recipinvgauss', 'ksone', 'genexpon',
            'vonmises', 'vonmises_line', 'mielke', 'semicircular',
            'cosine', 'invweibull', 'powerlognorm', 'johnsonsu', 'kstwobign']
# distslow are sorted by speed (very slow to slow)


fails_cmplx = set(['alpha', 'beta', 'betaprime', 'chi', 'chi2', 'dgamma', 
    'dweibull', 'erlang', 'expon', 'exponnorm', 'exponpow', 'exponweib', 'f',
    'fatiguelife', 'foldnorm', 'frechet_l', 'frechet_r', 'gamma', 'gausshyper',
    'genexpon', 'genextreme', 'gengamma', 'genlogistic', 'gennorm', 'genpareto',
    'gilbrat', 'gompertz', 'halfcauchy', 'halfgennorm', 'halflogistic',
    'halfnorm', 'invgamma', 'invgauss', 'johnsonsb', 'johnsonsu', 'ksone',
    'kstwobign', 'levy_l', 'loggamma', 'logistic', 'lognorm', 'lomax',
    'maxwell', 'nakagami', 'ncf', 'nct', 'ncx2', 'norm', 'pearson3',
    'powerlognorm', 'powernorm', 'rayleigh', 'recipinvgauss', 'rice', 't',
    'truncexpon', 'truncnorm', 'tukeylambda', 'vonmises', 'vonmises_line',
    'wald', 'weibull_min'])


# NB: not needed anymore?
def _silence_fp_errors(func):
    # warning: don't apply to test_ functions as is, then those will be skipped
    def wrap(*a, **kw):
        olderr = np.seterr(all='ignore')
        try:
            return func(*a, **kw)
        finally:
            np.seterr(**olderr)
    wrap.__name__ = func.__name__
    return wrap


def test_cont_basic():
    # this test skips slow distributions
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=integrate.IntegrationWarning)
        for distname, arg in distcont[:]:
            if distname in distslow:
                continue
            if distname is 'levy_stable':
                continue
            distfn = getattr(stats, distname)
            np.random.seed(765456)
            sn = 500
            rvs = distfn.rvs(size=sn, *arg)
            sm = rvs.mean()
            sv = rvs.var()
            m, v = distfn.stats(*arg)

            yield check_sample_meanvar_, distfn, arg, m, v, sm, sv, sn, \
                   distname + 'sample mean test'
            yield check_cdf_ppf, distfn, arg, distname
            yield check_sf_isf, distfn, arg, distname
            yield check_pdf, distfn, arg, distname
            yield check_pdf_logpdf, distfn, arg, distname
            yield check_cdf_logcdf, distfn, arg, distname
            yield check_sf_logsf, distfn, arg, distname
            if distname in distmissing:
                alpha = 0.01
                yield check_distribution_rvs, distname, arg, alpha, rvs

            locscale_defaults = (0, 1)
            meths = [distfn.pdf, distfn.logpdf, distfn.cdf, distfn.logcdf,
                     distfn.logsf]
            # make sure arguments are within support
            spec_x = {'frechet_l': -0.5, 'weibull_max': -0.5, 'levy_l': -0.5,
                      'pareto': 1.5, 'tukeylambda': 0.3}
            x = spec_x.get(distname, 0.5)
            yield check_named_args, distfn, x, arg, locscale_defaults, meths
            yield check_random_state_property, distfn, arg

            # Entropy
            skp = npt.dec.skipif
            yield check_entropy, distfn, arg, distname

            if distfn.numargs == 0:
                yield skp(NUMPY_BELOW_1_7)(check_vecentropy), distfn, arg
            if distfn.__class__._entropy != stats.rv_continuous._entropy:
                yield check_private_entropy, distfn, arg, stats.rv_continuous

            yield check_edge_support, distfn, arg

            yield check_meth_dtype, distfn, arg, meths
            yield check_ppf_dtype, distfn, arg

            knf = npt.dec.knownfailureif
            yield knf(distname == 'truncnorm')(check_ppf_private), distfn, \
                      arg, distname


@npt.dec.slow
def test_cont_basic_slow():
    # same as above for slow distributions
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=integrate.IntegrationWarning)
        for distname, arg in distcont[:]:
            if distname not in distslow:
                continue
            if distname is 'levy_stable':
                continue
            distfn = getattr(stats, distname)
            np.random.seed(765456)
            sn = 500
            rvs = distfn.rvs(size=sn,*arg)
            sm = rvs.mean()
            sv = rvs.var()
            m, v = distfn.stats(*arg)
            yield check_sample_meanvar_, distfn, arg, m, v, sm, sv, sn, \
                  distname + 'sample mean test'
            yield check_cdf_ppf, distfn, arg, distname
            yield check_sf_isf, distfn, arg, distname
            yield check_pdf, distfn, arg, distname
            yield check_pdf_logpdf, distfn, arg, distname
            yield check_cdf_logcdf, distfn, arg, distname
            yield check_sf_logsf, distfn, arg, distname
            # yield check_oth, distfn, arg # is still missing
            if distname in distmissing:
                alpha = 0.01
                yield check_distribution_rvs, distname, arg, alpha, rvs

            locscale_defaults = (0, 1)
            meths = [distfn.pdf, distfn.logpdf, distfn.cdf, distfn.logcdf,
                     distfn.logsf]
            # make sure arguments are within support
            x = 0.5
            if distname == 'invweibull':
                arg = (1,)
            elif distname == 'ksone':
                arg = (3,)
            yield check_named_args, distfn, x, arg, locscale_defaults, meths
            yield check_random_state_property, distfn, arg

            # Entropy
            skp = npt.dec.skipif
            ks_cond = distname in ['ksone', 'kstwobign']
            yield skp(ks_cond)(check_entropy), distfn, arg, distname

            if distfn.numargs == 0:
                yield skp(NUMPY_BELOW_1_7)(check_vecentropy), distfn, arg
            if distfn.__class__._entropy != stats.rv_continuous._entropy:
                yield check_private_entropy, distfn, arg, stats.rv_continuous

            yield check_edge_support, distfn, arg

            yield check_meth_dtype, distfn, arg, meths
            yield check_ppf_dtype, distfn, arg


@npt.dec.slow
def test_moments():
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=integrate.IntegrationWarning)
        knf = npt.dec.knownfailureif
        fail_normalization = set(['vonmises', 'ksone'])
        fail_higher = set(['vonmises', 'ksone', 'ncf'])
        for distname, arg in distcont[:]:
            if distname is 'levy_stable':
                continue
            distfn = getattr(stats, distname)
            m, v, s, k = distfn.stats(*arg, moments='mvsk')
            cond1, cond2 = distname in fail_normalization, distname in fail_higher
            msg = distname + ' fails moments'
            yield knf(cond1, msg)(check_normalization), distfn, arg, distname
            yield knf(cond2, msg)(check_mean_expect), distfn, arg, m, distname
            yield knf(cond2, msg)(check_var_expect), distfn, arg, m, v, distname
            yield knf(cond2, msg)(check_skew_expect), distfn, arg, m, v, s, \
                  distname
            yield knf(cond2, msg)(check_kurt_expect), distfn, arg, m, v, k, \
                  distname
            yield check_loc_scale, distfn, arg, m, v, distname
            yield check_moment, distfn, arg, m, v, distname


def check_sample_meanvar_(distfn, arg, m, v, sm, sv, sn, msg):
    # this did not work, skipped silently by nose
    if np.isfinite(m):
        check_sample_mean(sm, sv, sn, m)
    if np.isfinite(v):
        check_sample_var(sv, sn, v)


def check_sample_mean(sm,v,n, popmean):
    # from stats.stats.ttest_1samp(a, popmean):
    # Calculates the t-obtained for the independent samples T-test on ONE group
    # of scores a, given a population mean.
    #
    # Returns: t-value, two-tailed prob
    df = n-1
    svar = ((n-1)*v) / float(df)    # looks redundant
    t = (sm-popmean) / np.sqrt(svar*(1.0/n))
    prob = betainc(0.5*df, 0.5, df/(df + t*t))

    # return t,prob
    npt.assert_(prob > 0.01, 'mean fail, t,prob = %f, %f, m, sm=%f,%f' %
            (t, prob, popmean, sm))


def check_sample_var(sv,n, popvar):
    # two-sided chisquare test for sample variance equal to hypothesized variance
    df = n-1
    chi2 = (n-1)*popvar/float(popvar)
    pval = stats.distributions.chi2.sf(chi2, df) * 2
    npt.assert_(pval > 0.01, 'var fail, t, pval = %f, %f, v, sv=%f, %f' %
                (chi2, pval, popvar, sv))


def check_cdf_ppf(distfn,arg,msg):
    values = [0.001, 0.5, 0.999]
    npt.assert_almost_equal(distfn.cdf(distfn.ppf(values, *arg), *arg),
                            values, decimal=DECIMAL, err_msg=msg +
                            ' - cdf-ppf roundtrip')


def check_sf_isf(distfn,arg,msg):
    npt.assert_almost_equal(distfn.sf(distfn.isf([0.1,0.5,0.9], *arg), *arg),
                            [0.1,0.5,0.9], decimal=DECIMAL, err_msg=msg +
                            ' - sf-isf roundtrip')
    npt.assert_almost_equal(distfn.cdf([0.1,0.9], *arg),
                            1.0-distfn.sf([0.1,0.9], *arg),
                            decimal=DECIMAL, err_msg=msg +
                            ' - cdf-sf relationship')


def check_pdf(distfn, arg, msg):
    # compares pdf at median with numerical derivative of cdf
    median = distfn.ppf(0.5, *arg)
    eps = 1e-6
    pdfv = distfn.pdf(median, *arg)
    if (pdfv < 1e-4) or (pdfv > 1e4):
        # avoid checking a case where pdf is close to zero or huge (singularity)
        median = median + 0.1
        pdfv = distfn.pdf(median, *arg)
    cdfdiff = (distfn.cdf(median + eps, *arg) -
               distfn.cdf(median - eps, *arg))/eps/2.0
    # replace with better diff and better test (more points),
    # actually, this works pretty well
    npt.assert_almost_equal(pdfv, cdfdiff,
                decimal=DECIMAL, err_msg=msg + ' - cdf-pdf relationship')


def check_pdf_logpdf(distfn, args, msg):
    # compares pdf at several points with the log of the pdf
    points = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
    vals = distfn.ppf(points, *args)
    pdf = distfn.pdf(vals, *args)
    logpdf = distfn.logpdf(vals, *args)
    pdf = pdf[pdf != 0]
    logpdf = logpdf[np.isfinite(logpdf)]
    npt.assert_almost_equal(np.log(pdf), logpdf, decimal=7, err_msg=msg + " - logpdf-log(pdf) relationship")


def check_sf_logsf(distfn, args, msg):
    # compares sf at several points with the log of the sf
    points = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
    vals = distfn.ppf(points, *args)
    sf = distfn.sf(vals, *args)
    logsf = distfn.logsf(vals, *args)
    sf = sf[sf != 0]
    logsf = logsf[np.isfinite(logsf)]
    npt.assert_almost_equal(np.log(sf), logsf, decimal=7, err_msg=msg + " - logsf-log(sf) relationship")


def check_cdf_logcdf(distfn, args, msg):
    # compares cdf at several points with the log of the cdf
    points = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
    vals = distfn.ppf(points, *args)
    cdf = distfn.cdf(vals, *args)
    logcdf = distfn.logcdf(vals, *args)
    cdf = cdf[cdf != 0]
    logcdf = logcdf[np.isfinite(logcdf)]
    npt.assert_almost_equal(np.log(cdf), logcdf, decimal=7, err_msg=msg + " - logcdf-log(cdf) relationship")


def check_distribution_rvs(dist, args, alpha, rvs):
    # test from scipy.stats.tests
    # this version reuses existing random variables
    D,pval = stats.kstest(rvs, dist, args=args, N=1000)
    if (pval < alpha):
        D,pval = stats.kstest(dist,'',args=args, N=1000)
        npt.assert_(pval > alpha, "D = " + str(D) + "; pval = " + str(pval) +
               "; alpha = " + str(alpha) + "\nargs = " + str(args))


def check_vecentropy(distfn, args):
    npt.assert_equal(distfn.vecentropy(*args), distfn._entropy(*args))


@npt.dec.skipif(NUMPY_BELOW_1_7)
def check_loc_scale(distfn, arg, m, v, msg):
    loc, scale = 10.0, 10.0
    mt, vt = distfn.stats(loc=loc, scale=scale, *arg)
    npt.assert_allclose(m*scale + loc, mt)
    npt.assert_allclose(v*scale*scale, vt)


def check_ppf_private(distfn, arg, msg):
    #fails by design for truncnorm self.nb not defined
    ppfs = distfn._ppf(np.array([0.1, 0.5, 0.9]), *arg)
    npt.assert_(not np.any(np.isnan(ppfs)), msg + 'ppf private is nan')


if __name__ == "__main__":
    npt.run_module_suite()
