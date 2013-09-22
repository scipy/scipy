from __future__ import division, print_function, absolute_import

import warnings
import inspect

import numpy as np
import numpy.testing as npt

from scipy import stats
from common_tests import (check_normalization, check_moment, check_mean_expect,
        check_var_expect, check_skew_expect, check_kurt_expect)

"""
Test all continuous distributions.

Parameters were chosen for those distributions that pass the
Kolmogorov-Smirnov test.  This provides safe parameters for each
distributions so that we can perform further testing of class methods.

These tests currently check only/mostly for serious errors and exceptions,
not for numerically exact results.
"""

DECIMAL = 5  # specify the precision of the tests  # increased from 0 to 5

distcont = [
    ['alpha', (3.5704770516650459,)],
    ['anglit', ()],
    ['arcsine', ()],
    ['beta', (2.3098496451481823, 0.62687954300963677)],
    ['betaprime', (5, 6)],
    ['bradford', (0.29891359763170633,)],
    ['burr', (10.5, 4.3)],
    ['cauchy', ()],
    ['chi', (78,)],
    ['chi2', (55,)],
    ['cosine', ()],
    ['dgamma', (1.1023326088288166,)],
    ['dweibull', (2.0685080649914673,)],
    ['erlang', (10,)],
    ['expon', ()],
    ['exponpow', (2.697119160358469,)],
    ['exponweib', (2.8923945291034436, 1.9505288745913174)],
    ['f', (29, 18)],
    ['fatiguelife', (29,)],   # correction numargs = 1
    ['fisk', (3.0857548622253179,)],
    ['foldcauchy', (4.7164673455831894,)],
    ['foldnorm', (1.9521253373555869,)],
    ['frechet_l', (3.6279911255583239,)],
    ['frechet_r', (1.8928171603534227,)],
    ['gamma', (1.9932305483800778,)],
    ['gausshyper', (13.763771604130699, 3.1189636648681431,
                    2.5145980350183019, 5.1811649903971615)],  # veryslow
    ['genexpon', (9.1325976465418908, 16.231956600590632, 3.2819552690843983)],
    ['genextreme', (-0.1,)],
    ['gengamma', (4.4162385429431925, 3.1193091679242761)],
    ['genhalflogistic', (0.77274727809929322,)],
    ['genlogistic', (0.41192440799679475,)],
    ['genpareto', (0.1,)],   # use case with finite moments
    ['gilbrat', ()],
    ['gompertz', (0.94743713075105251,)],
    ['gumbel_l', ()],
    ['gumbel_r', ()],
    ['halfcauchy', ()],
    ['halflogistic', ()],
    ['halfnorm', ()],
    ['hypsecant', ()],
    ['invgamma', (4.0668996136993067,)],
    ['invgauss', (0.14546264555347513,)],
    ['invweibull', (10.58,)],
    ['johnsonsb', (4.3172675099141058, 3.1837781130785063)],
    ['johnsonsu', (2.554395574161155, 2.2482281679651965)],
    ['ksone', (1000,)],  # replace 22 by 100 to avoid failing range, ticket 956
    ['kstwobign', ()],
    ['laplace', ()],
    ['levy', ()],
    ['levy_l', ()],
#    ['levy_stable', (0.35667405469844993,
#                     -0.67450531578494011)], #NotImplementedError
    #           rvs not tested
    ['loggamma', (0.41411931826052117,)],
    ['logistic', ()],
    ['loglaplace', (3.2505926592051435,)],
    ['lognorm', (0.95368226960575331,)],
    ['lomax', (1.8771398388773268,)],
    ['maxwell', ()],
    ['mielke', (10.4, 3.6)],
    ['nakagami', (4.9673794866666237,)],
    ['ncf', (27, 27, 0.41578441799226107)],
    ['nct', (14, 0.24045031331198066)],
    ['ncx2', (21, 1.0560465975116415)],
    ['norm', ()],
    ['pareto', (2.621716532144454,)],
    ['pearson3', (0.1,)],
    ['powerlaw', (1.6591133289905851,)],
    ['powerlognorm', (2.1413923530064087, 0.44639540782048337)],
    ['powernorm', (4.4453652254590779,)],
    ['rayleigh', ()],
    ['rdist', (0.9,)],   # feels also slow
    ['recipinvgauss', (0.63004267809369119,)],
    ['reciprocal', (0.0062309367010521255, 1.0062309367010522)],
    ['rice', (0.7749725210111873,)],
    ['semicircular', ()],
    ['t', (2.7433514990818093,)],
    ['triang', (0.15785029824528218,)],
    ['truncexpon', (4.6907725456810478,)],
    ['truncnorm', (-1.0978730080013919, 2.7306754109031979)],
    ['truncnorm', (0.1, 2.)],
    ['tukeylambda', (3.1321477856738267,)],
    ['uniform', ()],
    ['vonmises', (3.9939042581071398,)],
    ['vonmises_line', (3.9939042581071398,)],
    ['wald', ()],
    ['weibull_max', (2.8687961709100187,)],
    ['weibull_min', (1.7866166930421596,)],
    ['wrapcauchy', (0.031071279018614728,)]]

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
    'johnsonsb', 'truncexpon', 'rice', 'invgauss', 'invgamma',
    'powerlognorm']

distmiss = [[dist,args] for dist,args in distcont if dist in distmissing]
distslow = ['rdist', 'gausshyper', 'recipinvgauss', 'ksone', 'genexpon',
            'vonmises', 'vonmises_line', 'mielke', 'semicircular',
            'cosine', 'invweibull', 'powerlognorm', 'johnsonsu', 'kstwobign']
# distslow are sorted by speed (very slow to slow)


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
    for distname, arg in distcont[:]:
        if distname in distslow:
            continue
        distfn = getattr(stats, distname)
        np.random.seed(765456)
        sn = 500
        rvs = distfn.rvs(size=sn, *arg)
        sm = rvs.mean()
        sv = rvs.var()
        skurt = stats.kurtosis(rvs)
        sskew = stats.skew(rvs)
        m, v = distfn.stats(*arg)

        yield check_sample_meanvar_, distfn, arg, m, v, sm, sv, sn, distname + \
              'sample mean test'
        # the sample skew kurtosis test has known failures, not very good distance measure
        # yield check_sample_skew_kurt, distfn, arg, sskew, skurt, distname
        yield check_moment, distfn, arg, m, v, distname
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

        # this asserts the vectorization of _entropy w/ no shape parameters
        # NB: broken for older versions of numpy
        if distfn.numargs == 0:
            if np.__version__ > '1.7':
                yield check_vecentropy, distfn, arg
        # compare a generic _entropy w/ distribution-specific implementation,
        # if available
        if distfn.__class__._entropy != stats.rv_continuous._entropy:
            yield check_private_entropy, distfn, arg
        yield check_edge_support, distfn, arg


@npt.dec.slow
def test_cont_basic_slow():
    # same as above for slow distributions
    for distname, arg in distcont[:]:
        if distname not in distslow:
            continue
        distfn = getattr(stats, distname)
        np.random.seed(765456)
        sn = 500
        rvs = distfn.rvs(size=sn,*arg)
        sm = rvs.mean()
        sv = rvs.var()
        skurt = stats.kurtosis(rvs)
        sskew = stats.skew(rvs)
        m, v = distfn.stats(*arg)
        yield check_sample_meanvar_, distfn, arg, m, v, sm, sv, sn, distname + \
              'sample mean test'
        # the sample skew kurtosis test has known failures, not very good distance measure
        # yield check_sample_skew_kurt, distfn, arg, sskew, skurt, distname
        # vonmises and ksone are not supposed to fully work
        yield check_moment, distfn, arg, m, v, distname
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

        # this asserts the vectorization of _entropy w/ no shape parameters
        # NB: broken for older versions of numpy
        if distfn.numargs == 0:
            if np.__version__ > '1.7':
                yield check_vecentropy, distfn, arg
        # compare a generic _entropy w/ distribution-specific implementation,
        # if available
        if distfn.__class__._entropy != stats.rv_continuous._entropy:
            yield check_private_entropy, distfn, arg
        yield check_edge_support, distfn, arg


@npt.dec.slow
def test_moments():
    knf = npt.dec.knownfailureif
    fail_normalization = set(['vonmises', 'ksone'])
    fail_higher = set(['vonmises', 'ksone', 'ncf'])
    for distname, arg in distcont[:]:
        distfn = getattr(stats, distname)
        m, v, s, k = distfn.stats(*arg, moments='mvsk')
        cond1, cond2 = distname in fail_normalization, distname in fail_higher
        msg = distname + ' fails moments'
        yield knf(cond1, msg)(check_normalization), distfn, arg, distname
        yield knf(cond2, msg)(check_mean_expect), distfn, arg, m, distname
        yield knf(cond2, msg)(check_var_expect), distfn, arg, m, v, distname
        yield knf(cond2, msg)(check_skew_expect), distfn, arg, m, v, s, distname
        yield knf(cond2, msg)(check_kurt_expect), distfn, arg, m, v, k, distname


def check_sample_meanvar_(distfn, arg, m, v, sm, sv, sn, msg):
    # this did not work, skipped silently by nose
    if not np.isinf(m):
        check_sample_mean(sm, sv, sn, m)
    if not np.isinf(v):
        check_sample_var(sv, sn, v)


def check_sample_mean(sm,v,n, popmean):
    """
from stats.stats.ttest_1samp(a, popmean):
Calculates the t-obtained for the independent samples T-test on ONE group
of scores a, given a population mean.

Returns: t-value, two-tailed prob
"""
##    a = asarray(a)
##    x = np.mean(a)
##    v = np.var(a, ddof=1)
##    n = len(a)
    df = n-1
    svar = ((n-1)*v) / float(df)    # looks redundant
    t = (sm-popmean)/np.sqrt(svar*(1.0/n))
    prob = stats.betai(0.5*df,0.5,df/(df+t*t))

    # return t,prob
    npt.assert_(prob > 0.01, 'mean fail, t,prob = %f, %f, m,sm=%f,%f' % (t,prob,popmean,sm))


def check_sample_var(sv,n, popvar):
    '''
two-sided chisquare test for sample variance equal to hypothesized variance
    '''
    df = n-1
    chi2 = (n-1)*popvar/float(popvar)
    pval = stats.chisqprob(chi2,df)*2
    npt.assert_(pval > 0.01, 'var fail, t,pval = %f, %f, v,sv=%f,%f' % (chi2,pval,popvar,sv))


def check_sample_skew_kurt(distfn, arg, ss, sk, msg):
    skew,kurt = distfn.stats(moments='sk',*arg)
##    skew = distfn.stats(moment='s',*arg)[()]
##    kurt = distfn.stats(moment='k',*arg)[()]
    check_sample_meanvar(sk, kurt, msg + 'sample kurtosis test')
    check_sample_meanvar(ss, skew, msg + 'sample skew test')


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
    npt.assert_equal( distfn.vecentropy(*args), distfn._entropy(*args) )


def check_named_args(distfn, x, shape_args, defaults, meths):
    """Check calling w/ named arguments."""

    # check consistency of shapes, numargs and _parse signature
    signature = inspect.getargspec(distfn._parse_args)
    npt.assert_(signature.varargs is None)
    npt.assert_(signature.keywords is None)
    npt.assert_(signature.defaults == defaults)

    shape_argnames = signature.args[1:-len(defaults)]  # self, a, b, loc=0, scale=1
    if distfn.shapes:
        shapes_ = distfn.shapes.replace(',',' ').split()
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
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                npt.assert_equal(distfn.moment(3, *a, **k),
                                 distfn.moment(3, *shape_args))

    # unknown arguments should not go through:
    k.update({'kaboom': 42})
    npt.assert_raises(TypeError, distfn.cdf, x, **k)


def check_edge_support(distfn, args):
    """Make sure the x=self.a and self.b are handled correctly."""
    x = [distfn.a, distfn.b]
    npt.assert_equal(distfn.cdf(x, *args), [0.0, 1.0])
    npt.assert_equal(distfn.logcdf(x, *args), [-np.inf, 0.0])

    npt.assert_equal(distfn.sf(x, *args), [1.0, 0.0])
    npt.assert_equal(distfn.logsf(x, *args), [0.0, -np.inf])

    npt.assert_equal(distfn.ppf([0.0, 1.0], *args), x)
    npt.assert_equal(distfn.isf([0.0, 1.0], *args), x[::-1])
    # pdf(x=[a, b], *args) depends on the distribution


def check_private_entropy(distfn, args):
    # compare a generic _entropy with the distribution-specific implementation
    npt.assert_allclose(distfn._entropy(*args),
                        stats.rv_continuous._entropy(distfn, *args))




def check_edge_support(distfn, args):
    """Make sure the x=self.a and self.b are handled correctly."""
    x = [distfn.a, distfn.b]
    npt.assert_equal(distfn.cdf(x, *args), [0.0, 1.0])
    npt.assert_equal(distfn.logcdf(x, *args), [-np.inf, 0.0])

    npt.assert_equal(distfn.sf(x, *args), [1.0, 0.0])
    npt.assert_equal(distfn.logsf(x, *args), [0.0, -np.inf])

    npt.assert_equal(distfn.ppf([0.0, 1.0], *args), x)
    npt.assert_equal(distfn.isf([0.0, 1.0], *args), x[::-1])
    # pdf(x=[a, b], *args) depends on the distribution


if __name__ == "__main__":
    npt.run_module_suite()
