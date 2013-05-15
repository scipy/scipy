from __future__ import division, print_function, absolute_import

import os

import numpy as np
from numpy.testing import dec

from scipy import stats

from test_continuous_basic import distcont

# this is not a proper statistical test for convergence, but only
# verifies that the estimate and true values don't differ by too much

fit_sizes = [1000, 5000]  # sample sizes to try
thresh_percent = 0.25  # percent of true parameters for fail cut-off
thresh_min = 0.75  # minimum difference estimate - true to fail test

failing_fits = [
        'burr',
        'chi',
        'chi2',
        'gausshyper',
        'genexpon',
        'gengamma',
        'ksone',
        'mielke',
        'ncf',
        'ncx2',
        'pearson3',
        'powerlognorm',
        'truncexpon',
        'tukeylambda',
        'vonmises',
        'wrapcauchy',
]

# Don't run the fit test on these:
skip_fit = [
    'erlang',  # Subclass of gamma, generates a warning.
]


@dec.slow
def test_cont_fit():
    # this tests the closeness of the estimated parameters to the true
    # parameters with fit method of continuous distributions
    # Note: is slow, some distributions don't converge with sample size <= 10000

    for distname, arg in distcont:
        if distname not in skip_fit:
            yield check_cont_fit, distname,arg


def check_cont_fit(distname,arg):
    if distname in failing_fits:
        # Skip failing fits unless overridden
        xfail = True
        try:
            xfail = not int(os.environ['SCIPY_XFAIL'])
        except:
            pass
        if xfail:
            msg = "Fitting %s doesn't work reliably yet" % distname
            msg += " [Set environment variable SCIPY_XFAIL=1 to run this test nevertheless.]"
            dec.knownfailureif(True, msg)(lambda: None)()

    distfn = getattr(stats, distname)

    truearg = np.hstack([arg,[0.0,1.0]])
    diffthreshold = np.max(np.vstack([truearg*thresh_percent,
                                      np.ones(distfn.numargs+2)*thresh_min]),0)

    for fit_size in fit_sizes:
        # Note that if a fit succeeds, the other fit_sizes are skipped
        np.random.seed(1234)

        with np.errstate(all='ignore'):
            rvs = distfn.rvs(size=fit_size, *arg)
            est = distfn.fit(rvs)  # start with default values

        diff = est - truearg

        # threshold for location
        diffthreshold[-2] = np.max([np.abs(rvs.mean())*thresh_percent,thresh_min])

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


if __name__ == "__main__":
    np.testing.run_module_suite()
