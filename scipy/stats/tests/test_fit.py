import os
import numpy as np
import numpy.testing as npt
from numpy.testing import assert_allclose
import pytest
from scipy import stats
from scipy.optimize import differential_evolution

from .test_continuous_basic import distcont
from scipy.stats._distn_infrastructure import FitError
from scipy.stats._distr_params import distdiscrete


# this is not a proper statistical test for convergence, but only
# verifies that the estimate and true values don't differ by too much

fit_sizes = [1000, 5000, 10000]  # sample sizes to try

thresh_percent = 0.25  # percent of true parameters for fail cut-off
thresh_min = 0.75  # minimum difference estimate - true to fail test

mle_failing_fits = [
        'burr',
        'chi2',
        'gausshyper',
        'genexpon',
        'gengamma',
        'kappa4',
        'ksone',
        'kstwo',
        'mielke',
        'ncf',
        'ncx2',
        'pearson3',
        'powerlognorm',
        'truncexpon',
        'tukeylambda',
        'vonmises',
        'levy_stable',
        'trapezoid',
        'truncweibull_min',
        'studentized_range',
]

mm_failing_fits = ['alpha', 'betaprime', 'burr', 'burr12', 'cauchy', 'chi',
                   'chi2', 'crystalball', 'dgamma', 'dweibull', 'f',
                   'fatiguelife', 'fisk', 'foldcauchy', 'genextreme',
                   'gengamma', 'genhyperbolic', 'gennorm', 'genpareto',
                   'halfcauchy', 'invgamma', 'invweibull', 'johnsonsu',
                   'kappa3', 'ksone', 'kstwo', 'levy', 'levy_l',
                   'levy_stable', 'loglaplace', 'lomax', 'mielke', 'nakagami',
                   'ncf', 'nct', 'ncx2', 'pareto', 'powerlognorm', 'powernorm',
                   'skewcauchy', 't', 'trapezoid', 'triang',
                   'truncweibull_min', 'tukeylambda', 'studentized_range']

# not sure if these fail, but they caused my patience to fail
mm_slow_fits = ['argus', 'exponpow', 'exponweib', 'gausshyper', 'genexpon',
                'genhalflogistic', 'halfgennorm', 'gompertz', 'johnsonsb',
                'kappa4', 'kstwobign', 'recipinvgauss', 'skewnorm',
                'truncexpon', 'vonmises', 'vonmises_line']

failing_fits = {"MM": mm_failing_fits + mm_slow_fits, "MLE": mle_failing_fits}

# Don't run the fit test on these:
skip_fit = [
    'erlang',  # Subclass of gamma, generates a warning.
    'genhyperbolic',  # too slow
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
@pytest.mark.parametrize('method', ["MLE", "MM"])
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
            est = distfn.fit(rvs, method=method)  # start with default values

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


def test_fit_error():
    data = np.concatenate([np.zeros(29), np.ones(21)])
    message = "Optimization converged to parameters that are..."
    with pytest.raises(FitError, match=message), \
            pytest.warns(RuntimeWarning):
        stats.beta.fit(data)


@pytest.mark.parametrize("dist, params",
                         [(stats.norm, (0.5, 2.5)),  # type: ignore[attr-defined] # noqa
                          (stats.binom, (10, 0.3, 2))])  # type: ignore[attr-defined] # noqa
def test_nnlf_and_related_methods(dist, params):
    rng = np.random.default_rng(983459824)

    if hasattr(dist, 'pdf'):
        logpxf = dist.logpdf
    else:
        logpxf = dist.logpmf

    x = dist.rvs(*params, size=100, random_state=rng)
    ref = -logpxf(x, *params).sum()
    res1 = dist.nnlf(params, x)
    res2 = dist._penalized_nnlf(params, x)
    assert_allclose(res1, ref)
    assert_allclose(res2, ref)


def cases_test_fit():
    # These three fail default test; check separately
    skip_basic_fit = {'argus', 'foldnorm', 'truncweibull_min'}
    # status of 'studentized_range', 'ksone', 'kstwo' unknown; all others pass
    slow_basic_fit = {'burr12', 'johnsonsb', 'bradford', 'fisk', 'mielke',
                      'exponpow', 'rdist', 'norminvgauss', 'betaprime',
                      'powerlaw', 'pareto', 'johnsonsu', 'loglaplace',
                      'wrapcauchy', 'weibull_max', 'arcsine', 'binom', 'rice',
                      'uniform', 'f', 'invweibull', 'genpareto', 'weibull_min',
                      'nbinom', 'kappa3', 'lognorm', 'halfgennorm', 'pearson3',
                      'alpha', 't', 'crystalball', 'fatiguelife', 'nakagami',
                      'kstwobign', 'gompertz', 'dweibull', 'lomax', 'invgauss',
                      'recipinvgauss', 'chi', 'foldcauchy', 'powernorm',
                      'gennorm', 'skewnorm', 'randint', 'genextreme'}
    xslow_basic_fit = {'studentized_range', 'ksone', 'kstwo', 'levy_stable',
                       'nchypergeom_fisher', 'nchypergeom_wallenius',
                       'gausshyper', 'genexpon', 'gengamma', 'genhyperbolic',
                       'geninvgauss', 'tukeylambda', 'skellam', 'ncx2',
                       'hypergeom', 'nhypergeom', 'zipfian', 'ncf',
                       'truncnorm', 'powerlognorm', 'beta',
                       'loguniform', 'reciprocal', 'trapezoid', 'nct',
                       'kappa4', 'betabinom', 'exponweib', 'genhalflogistic',
                       'burr', 'triang'}

    for dist in dict(distdiscrete + distcont):
        if dist in skip_basic_fit or not isinstance(dist, str):
            reason = "tested separately"
            yield pytest.param(dist, marks=pytest.mark.skip(reason=reason))
        elif dist in slow_basic_fit:
            reason = "too slow (>= 0.25s)"
            yield pytest.param(dist, marks=pytest.mark.slow(reason=reason))
        elif dist in xslow_basic_fit:
            reason = "too slow (>= 1.0s)"
            yield pytest.param(dist, marks=pytest.mark.xslow(reason=reason))
        else:
            yield dist


def cases_test_fitstart():
    for distname, shapes in dict(distcont).items():
        if not isinstance(distname, str) or distname in {'studentized_range'}:
            continue
        yield distname, shapes


@pytest.mark.parametrize('distname, shapes', cases_test_fitstart())
def test_fitstart(distname, shapes):
    dist = getattr(stats, distname)
    rng = np.random.default_rng(216342614)
    data = rng.random(10)

    with np.errstate(invalid='ignore', divide='ignore'):  # irrelevant to test
        guess = dist._fitstart(data)

    assert dist._argcheck(*guess[:-2])


def assert_nllf_less_or_close(dist, data, params1, params0, rtol=1e-7, atol=0):
    nllf1 = dist.nnlf(params1, data)
    nllf0 = dist.nnlf(params0, data)
    if not (nllf1 < nllf0):
        np.testing.assert_allclose(nllf1, nllf0, rtol=rtol, atol=atol)


class TestFit:
    dist = stats.binom  # type: ignore[attr-defined]
    seed = 654634816187
    rng = np.random.default_rng(seed)
    data = stats.binom.rvs(5, 0.5, size=100, random_state=rng)  # type: ignore[attr-defined] # noqa
    shape_bounds_a = [(1, 10), (0, 1)]
    shape_bounds_d = {'n': (1, 10), 'p': (0, 1)}
    atol = 5e-2
    rtol = 1e-2
    tols = {'atol': atol, 'rtol': rtol}

    def opt(self, *args, **kwds):
        return differential_evolution(*args, seed=0, **kwds)

    def test_dist_iv(self):
        message = "`dist` must be an instance of..."
        with pytest.raises(ValueError, match=message):
            stats.fit(10, self.data, self.shape_bounds_a)

    def test_data_iv(self):
        message = "`data` must be exactly one-dimensional."
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, [[1, 2, 3]], self.shape_bounds_a)

        message = "All elements of `data` must be finite numbers."
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, [1, 2, 3, np.nan], self.shape_bounds_a)
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, [1, 2, 3, np.inf], self.shape_bounds_a)
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, ['1', '2', '3'], self.shape_bounds_a)

    def test_bounds_iv(self):
        message = "Bounds provided for the following unrecognized..."
        shape_bounds = {'n': (1, 10), 'p': (0, 1), '1': (0, 10)}
        with pytest.warns(RuntimeWarning, match=message):
            stats.fit(self.dist, self.data, shape_bounds)

        message = "Each element of a `bounds` sequence must be a tuple..."
        shape_bounds = [(1, 10, 3), (0, 1)]
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, self.data, shape_bounds)

        message = "Each element of `bounds` must be a tuple specifying..."
        shape_bounds = [(1, 10, 3), (0, 1, 0.5)]
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, self.data, shape_bounds)
        shape_bounds = [1, 0]
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, self.data, shape_bounds)

        message = "A `bounds` sequence must contain at least 2 elements..."
        shape_bounds = [(1, 10)]
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, self.data, shape_bounds)

        message = "A `bounds` sequence may not contain more than 3 elements..."
        bounds = [(1, 10), (1, 10), (1, 10), (1, 10)]
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, self.data, bounds)

        message = "There are no values for `p` on the interval..."
        shape_bounds = {'n': (1, 10), 'p': (1, 0)}
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, self.data, shape_bounds)

        message = "There are no values for `n` on the interval..."
        shape_bounds = [(10, 1), (0, 1)]
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, self.data, shape_bounds)

        message = "There are no integer values for `n` on the interval..."
        shape_bounds = [(1.4, 1.6), (0, 1)]
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, self.data, shape_bounds)

        message = "The intersection of user-provided bounds for `n`"
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, self.data)
        shape_bounds = [(-np.inf, np.inf), (0, 1)]
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, self.data, shape_bounds)

    def test_guess_iv(self):
        message = "Guesses provided for the following unrecognized..."
        guess = {'n': 1, 'p': 0.5, '1': 255}
        with pytest.warns(RuntimeWarning, match=message):
            stats.fit(self.dist, self.data, self.shape_bounds_d, guess=guess)

        message = "Each element of `guess` must be a scalar..."
        guess = {'n': 1, 'p': 'hi'}
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, self.data, self.shape_bounds_d, guess=guess)
        guess = [1, 'f']
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, self.data, self.shape_bounds_d, guess=guess)
        guess = [[1, 2]]
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, self.data, self.shape_bounds_d, guess=guess)

        message = "A `guess` sequence must contain at least 2..."
        guess = [1]
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, self.data, self.shape_bounds_d, guess=guess)

        message = "A `guess` sequence may not contain more than 3..."
        guess = [1, 2, 3, 4]
        with pytest.raises(ValueError, match=message):
            stats.fit(self.dist, self.data, self.shape_bounds_d, guess=guess)

        message = "Guess for parameter `n` rounded..."
        guess = {'n': 4.5, 'p': -0.5}
        with pytest.warns(RuntimeWarning, match=message):
            stats.fit(self.dist, self.data, self.shape_bounds_d, guess=guess)

        message = "Guess for parameter `loc` rounded..."
        guess = [5, 0.5, 0.5]
        with pytest.warns(RuntimeWarning, match=message):
            stats.fit(self.dist, self.data, self.shape_bounds_d, guess=guess)

        message = "Guess for parameter `p` clipped..."
        guess = {'n': 5, 'p': -0.5}
        with pytest.warns(RuntimeWarning, match=message):
            stats.fit(self.dist, self.data, self.shape_bounds_d, guess=guess)

        message = "Guess for parameter `loc` clipped..."
        guess = [5, 0.5, 1]
        with pytest.warns(RuntimeWarning, match=message):
            stats.fit(self.dist, self.data, self.shape_bounds_d, guess=guess)

    @pytest.mark.parametrize("dist_name", cases_test_fit())
    def test_basic_fit(self, dist_name):

        N = 5000
        dist_data = dict(distcont + distdiscrete)
        rng = np.random.default_rng(self.seed)
        dist = getattr(stats, dist_name)
        shapes = np.array(dist_data[dist_name])
        bounds = np.empty((len(shapes) + 2, 2), dtype=np.float64)
        bounds[:-2, 0] = shapes/10.**np.sign(shapes)
        bounds[:-2, 1] = shapes*10.**np.sign(shapes)
        bounds[-2] = (0, 10)
        bounds[-1] = (0, 10)
        loc = rng.uniform(*bounds[-2])
        scale = rng.uniform(*bounds[-1])
        ref = list(dist_data[dist_name]) + [loc, scale]

        if getattr(dist, 'pmf', False):
            ref = ref[:-1]
            ref[-1] = np.floor(loc)
            data = dist.rvs(*ref, size=N, random_state=rng)
            bounds = bounds[:-1]
        if getattr(dist, 'pdf', False):
            data = dist.rvs(*ref, size=N, random_state=rng)

        with npt.suppress_warnings() as sup:
            sup.filter(RuntimeWarning, "overflow encountered")
            res = stats.fit(dist, data, bounds, optimizer=self.opt)

        assert_nllf_less_or_close(dist, data, res.params, ref, **self.tols)

    def test_argus(self):
        # Can't guarantee that all distributions will fit all data with
        # arbitrary bounds. This distribution just happens to fail above.
        # Try something slightly different.
        N = 1000
        rng = np.random.default_rng(self.seed)
        dist = stats.argus
        shapes = (1., 2., 3.)
        data = dist.rvs(*shapes, size=N, random_state=rng)
        shape_bounds = {'chi': (0.1, 10), 'loc': (0.1, 10), 'scale': (0.1, 10)}
        res = stats.fit(dist, data, shape_bounds, optimizer=self.opt)

        assert_nllf_less_or_close(dist, data, res.params, shapes, **self.tols)

    def test_foldnorm(self):
        # Can't guarantee that all distributions will fit all data with
        # arbitrary bounds. This distribution just happens to fail above.
        # Try something slightly different.
        N = 1000
        rng = np.random.default_rng(self.seed)
        dist = stats.foldnorm
        shapes = (1.952125337355587, 2., 3.)
        data = dist.rvs(*shapes, size=N, random_state=rng)
        shape_bounds = {'c': (0.1, 10), 'loc': (0.1, 10), 'scale': (0.1, 10)}
        res = stats.fit(dist, data, shape_bounds, optimizer=self.opt)

        assert_nllf_less_or_close(dist, data, res.params, shapes, **self.tols)

    def test_truncweibull_min(self):
        # Can't guarantee that all distributions will fit all data with
        # arbitrary bounds. This distribution just happens to fail above.
        # Try something slightly different.
        N = 1000
        rng = np.random.default_rng(self.seed)
        dist = stats.truncweibull_min
        shapes = (2.5, 0.25, 1.75, 2., 3.)
        data = dist.rvs(*shapes, size=N, random_state=rng)
        shape_bounds = [(0.1, 10)]*5
        res = stats.fit(dist, data, shape_bounds, optimizer=self.opt)

        assert_nllf_less_or_close(dist, data, res.params, shapes, **self.tols)

    def test_missing_shape_bounds(self):
        # some distributions have a small domain w.r.t. a parameter, e.g.
        # $p \in [0, 1]$ for binomial distribution
        # User does not need to provide these because the intersection of the
        # user's bounds (none) and the distribution's domain is finite
        N = 1000
        rng = np.random.default_rng(self.seed)

        dist = stats.binom
        n, p, loc = 10, 0.65, 0
        data = dist.rvs(n, p, loc=loc, size=N, random_state=rng)
        shape_bounds = {'n': np.array([0, 20])}  # check arrays are OK, too
        res = stats.fit(dist, data, shape_bounds, optimizer=self.opt)
        assert_allclose(res.params, (n, p, loc), **self.tols)

        dist = stats.bernoulli
        p, loc = 0.314159, 0
        data = dist.rvs(p, loc=loc, size=N, random_state=rng)
        res = stats.fit(dist, data, optimizer=self.opt)
        assert_allclose(res.params, (p, loc), **self.tols)

    def test_fit_only_loc_scale(self):
        # fit only loc
        N = 5000
        rng = np.random.default_rng(self.seed)

        dist = stats.norm
        loc, scale = 1.5, 1
        data = dist.rvs(loc=loc, size=N, random_state=rng)
        loc_bounds = (0, 5)
        bounds = {'loc': loc_bounds}
        res = stats.fit(dist, data, bounds, optimizer=self.opt)
        assert_allclose(res.params, (loc, scale), **self.tols)

        # fit only scale
        loc, scale = 0, 2.5
        data = dist.rvs(scale=scale, size=N, random_state=rng)
        scale_bounds = (0, 5)
        bounds = {'scale': scale_bounds}
        res = stats.fit(dist, data, bounds, optimizer=self.opt)
        assert_allclose(res.params, (loc, scale), **self.tols)

        # fit only loc and scale
        dist = stats.norm
        loc, scale = 1.5, 2.5
        data = dist.rvs(loc=loc, scale=scale, size=N, random_state=rng)
        bounds = {'loc': loc_bounds, 'scale': scale_bounds}
        res = stats.fit(dist, data, bounds, optimizer=self.opt)
        assert_allclose(res.params, (loc, scale), **self.tols)

    def test_everything_fixed(self):
        N = 5000
        rng = np.random.default_rng(self.seed)

        dist = stats.norm
        loc, scale = 1.5, 2.5
        data = dist.rvs(loc=loc, scale=scale, size=N, random_state=rng)

        # loc, scale fixed to 0, 1 by default
        res = stats.fit(dist, data)
        assert_allclose(res.params, (0, 1), **self.tols)

        # loc, scale explicitly fixed
        bounds = {'loc': (loc, loc), 'scale': (scale, scale)}
        res = stats.fit(dist, data, bounds)
        assert_allclose(res.params, (loc, scale), **self.tols)

        # `n` gets fixed during polishing
        dist = stats.binom
        n, p, loc = 10, 0.65, 0
        data = dist.rvs(n, p, loc=loc, size=N, random_state=rng)
        shape_bounds = {'n': (0, 20), 'p': (0.65, 0.65)}
        res = stats.fit(dist, data, shape_bounds, optimizer=self.opt)
        assert_allclose(res.params, (n, p, loc), **self.tols)

    def test_failure(self):
        N = 5000
        rng = np.random.default_rng(self.seed)

        dist = stats.nbinom
        shapes = (5, 0.5)
        data = dist.rvs(*shapes, size=N, random_state=rng)

        assert data.min() == 0
        # With lower bounds on location at 0.5, likelihood is zero
        bounds = [(0, 30), (0, 1), (0.5, 10)]
        res = stats.fit(dist, data, bounds)
        message = "Optimization converged to parameter values that are"
        assert res.message.startswith(message)
        assert res.success is False

    @pytest.mark.xslow
    def test_guess(self):
        # Test that guess helps DE find the desired solution
        N = 2000
        rng = np.random.default_rng(self.seed)
        dist = stats.nhypergeom
        params = (20, 7, 12, 0)
        bounds = [(2, 200), (0.7, 70), (1.2, 120), (0, 10)]

        data = dist.rvs(*params, size=N, random_state=rng)

        res = stats.fit(dist, data, bounds, optimizer=self.opt)
        assert not np.allclose(res.params, params, **self.tols)

        res = stats.fit(dist, data, bounds, guess=params, optimizer=self.opt)
        assert_allclose(res.params, params, **self.tols)
