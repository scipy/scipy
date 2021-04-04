import os

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy import stats

from .test_continuous_basic import distcont

# this is not a proper statistical test for convergence, but only
# verifies that the estimate and true values don't differ by too much

fit_sizes = [1000, 5000]  # sample sizes to try

thresh_percent = 0.25  # percent of true parameters for fail cut-off
thresh_min = 0.75  # minimum difference estimate - true to fail test

mle_failing_fits = [
        'gausshyper',
        'genexpon',
        'gengamma',
        'kappa4',
        'ksone',
        'kstwo',
        'ncf',
        'ncx2',
        'truncexpon',
        'tukeylambda',
        'vonmises',
        'levy_stable',
        'trapezoid',
]

# The MLE fit method of these distributions doesn't perform well when all
# parameters are fit, so test them with the location fixed at 0.
mle_use_floc0 = [
    'burr',
    'chi',
    'chi2',
    'mielke',
    'pearson3',
    'genhalflogistic',
    'rdist',
    'pareto',
    'powerlognorm',
    'wrapcauchy',
]

mm_failing_fits = ['alpha', 'betaprime', 'burr', 'burr12', 'cauchy', 'chi',
                   'chi2', 'crystalball', 'dgamma', 'dweibull', 'f',
                   'fatiguelife', 'fisk', 'foldcauchy', 'genextreme',
                   'gengamma', 'gennorm', 'genpareto', 'halfcauchy',
                   'invgamma', 'invweibull', 'johnsonsu',
                   'kappa3', 'ksone', 'kstwo', 'levy', 'levy_l',
                   'levy_stable', 'loglaplace', 'lomax', 'mielke', 'nakagami',
                   'ncf', 'nct', 'ncx2', 'pareto', 'powerlognorm', 'powernorm',
                   'skewcauchy', 't',
                   'trapezoid', 'triang', 'tukeylambda']

# not sure if these fail, but they caused my patience to fail
mm_slow_fits = ['argus', 'exponpow', 'exponweib', 'gausshyper', 'genexpon',
                'genhalflogistic', 'halfgennorm', 'gompertz', 'johnsonsb',
                'kappa4', 'kstwobign', 'recipinvgauss', 'skewnorm',
                'truncexpon', 'vonmises', 'vonmises_line']

failing_fits = {"MM": mm_failing_fits + mm_slow_fits, "MLE": mle_failing_fits}

# Don't run the fit test on these:
skip_fit = [
    'erlang',  # Subclass of gamma, generates a warning.
]


def cases_test_cont_fit():
    # this tests the closeness of the estimated parameters to the true
    # parameters with fit method of continuous distributions
    # Note: is slow, some distributions don't converge with sample
    # size <= 10000
    for distname, arg in distcont:
        if distname not in skip_fit:
            yield distname, arg


@pytest.mark.slow
@pytest.mark.parametrize('distname,arg', cases_test_cont_fit())
@pytest.mark.parametrize('method', ["MLE", 'MM'])
def test_cont_fit(distname, arg, method):
    if distname in failing_fits[method]:
        # Skip failing fits unless overridden
        try:
            xfail = not int(os.environ['SCIPY_XFAIL'])
        except Exception:
            xfail = True
        if xfail:
            msg = "Fitting %s doesn't work reliably yet" % distname
            msg += (" [Set environment variable SCIPY_XFAIL=1 to run this"
                    " test nevertheless.]")
            pytest.xfail(msg)

    distfn = getattr(stats, distname)

    truearg = np.hstack([arg, [0.0, 1.0]])
    diffthreshold = np.max(np.vstack([truearg*thresh_percent,
                                      np.full(distfn.numargs+2, thresh_min)]),
                           0)

    for fit_size in fit_sizes:
        # Note that if a fit succeeds, the other fit_sizes are skipped
        np.random.seed(1234)

        with np.errstate(all='ignore'):
            rvs = distfn.rvs(size=fit_size, *arg)
            if method == 'MLE' and distfn.name in mle_use_floc0:
                kwds = {'floc': 0}
            else:
                kwds = {}
            # start with default values
            est = distfn.fit(rvs, method=method, **kwds)
            if method == 'MLE':
                # Trivial test of the use of CensoredData.  The fit() method
                # will check that data contains no actual censored data, and
                # do a regular uncensored fit.
                data1 = stats.CensoredData(rvs)
                est1 = distfn.fit(data1, **kwds)
                msg = ('Different results fitting uncensored data wrapped as'
                       f' CensoredData: {distfn.name}: est={est} est1={est1}')
                assert_allclose(est1, est, rtol=1e-10, err_msg=msg)
                # Convert the first `nic` values in rvs to interval-censored
                # values. The interval is small, so est2 should be close to
                # est.
                nic = 20
                intervals = np.column_stack((rvs, rvs))
                intervals[:nic, 0] *= 0.975
                intervals[:nic, 1] *= 1.025
                intervals.sort(axis=1)
                data2 = stats.CensoredData(intervals=intervals)
                est2 = distfn.fit(data2, **kwds)
                msg = ('Different results fitting interval-censored'
                       f' data: {distfn.name}: est={est} est2={est2}')
                assert_allclose(est2, est, rtol=0.05, err_msg=msg)

        diff = est - truearg

        # threshold for location
        diffthreshold[-2] = np.max([np.abs(rvs.mean())*thresh_percent,
                                    thresh_min])

        if np.any(np.isnan(est)):
            raise AssertionError('nan returned in fit')
        else:
            if np.all(np.abs(diff) <= diffthreshold):
                break
    else:
        txt = 'parameter: %s\n' % str(truearg)
        txt += 'estimated: %s\n' % str(est)
        txt += 'diff     : %s\n' % str(diff)
        raise AssertionError('fit not very good in %s\n' % distfn.name + txt)


def _check_loc_scale_mle_fit(name, data, desired, atol=None):
    d = getattr(stats, name)
    actual = d.fit(data)[-2:]
    assert_allclose(actual, desired, atol=atol,
                    err_msg='poor mle fit of (loc, scale) in %s' % name)


def test_non_default_loc_scale_mle_fit():
    data = np.array([1.01, 1.78, 1.78, 1.78, 1.88, 1.88, 1.88, 2.00])
    _check_loc_scale_mle_fit('uniform', data, [1.01, 0.99], 1e-3)
    _check_loc_scale_mle_fit('expon', data, [1.01, 0.73875], 1e-3)


def test_expon_fit():
    """gh-6167"""
    data = [0, 0, 0, 0, 2, 2, 2, 2]
    phat = stats.expon.fit(data, floc=0)
    assert_allclose(phat, [0, 1.0], atol=1e-3)
