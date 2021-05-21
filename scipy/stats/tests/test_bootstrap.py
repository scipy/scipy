import numpy as np
import pytest
from scipy.stats import bootstrap
from numpy.testing import assert_allclose
from scipy import stats
from .. import _bootstrap as _bootstrap
from scipy._lib._util import rng_integers


def test_bootstrap_iv():
    message = "`data` must be a sequence of samples."
    with pytest.raises(ValueError, match=message):
        bootstrap(1, np.mean)

    message = "`data` must contain at least one sample."
    with pytest.raises(ValueError, match=message):
        bootstrap(tuple(), np.mean)

    message = "each sample in `data` must contain two or more observations..."
    with pytest.raises(ValueError, match=message):
        bootstrap(([1, 2, 3], [1]), np.mean)

    message = "`axis` must be an integer."
    with pytest.raises(ValueError, match=message):
        bootstrap(([1, 2, 3],), np.mean, axis=1.5)

    message = "could not convert string to float"
    with pytest.raises(ValueError, match=message):
        bootstrap(([1, 2, 3],), np.mean, confidence_level='ni')

    message = "`n_resamples` must be a positive integer."
    with pytest.raises(ValueError, match=message):
        bootstrap(([1, 2, 3],), np.mean, n_resamples=-1000)

    message = "`n_resamples` must be a positive integer."
    with pytest.raises(ValueError, match=message):
        bootstrap(([1, 2, 3],), np.mean, n_resamples=1000.5)

    message = "`method` must be in"
    with pytest.raises(ValueError, match=message):
        bootstrap(([1, 2, 3],), np.mean, method='ekki')

    message = "`method = 'BCa' is only available for one-sample statistics"

    def statistic(x, y, axis):
        mean1 = np.mean(x, axis)
        mean2 = np.mean(y, axis)
        return mean1 - mean2

    with pytest.raises(ValueError, match=message):
        bootstrap(([.1, .2, .3], [.1, .2, .3]), statistic, method='BCa')

    message = "'herring' cannot be used to seed a"
    with pytest.raises(ValueError, match=message):
        bootstrap(([1, 2, 3],), np.mean, random_state='herring')


@pytest.mark.parametrize("method", ['basic', 'percentile', 'BCa'])
def test_bootstrap_against_theory(method):
    # based on https://www.statology.org/confidence-intervals-python/
    data = stats.norm.rvs(loc=5, scale=2, size=5000, random_state=0)
    alpha = 0.95
    dist = stats.t(df=len(data)-1, loc=np.mean(data), scale=stats.sem(data))
    expected_interval = dist.interval(alpha=alpha)
    expected_se = dist.std()

    res = bootstrap((data,), np.mean, n_resamples=5000,
                    confidence_level=alpha, method=method,
                    random_state=0)
    assert_allclose(res.confidence_interval, expected_interval, rtol=5e-4)
    assert_allclose(res.standard_error, expected_se, atol=3e-4)


tests_R = {"basic": (23.77, 79.12),
           "percentile": (28.86, 84.21),
           "BCa": (32.31, 91.43)}


@pytest.mark.parametrize("method, expected", tests_R.items())
def test_bootstrap_against_R(method, expected):
    # Compare against R's "boot" library
    # library(boot)

    # stat <- function (x, a) {
    #     mean(x[a])
    # }

    # x <- c(10, 12, 12.5, 12.5, 13.9, 15, 21, 22,
    #        23, 34, 50, 81, 89, 121, 134, 213)

    # # Use a large value so we get a few significant digits for the CI.
    # n = 1000000
    # bootresult = boot(x, stat, n)
    # result <- boot.ci(bootresult)
    # print(result)
    x = np.array([10, 12, 12.5, 12.5, 13.9, 15, 21, 22,
                  23, 34, 50, 81, 89, 121, 134, 213])
    res = bootstrap((x,), np.mean, n_resamples=1000000, method=method,
                    random_state=0)
    assert_allclose(res.confidence_interval, expected, rtol=0.005)


tests_against_itself_1samp = {"basic": 1780,
                              "percentile": 1784,
                              "BCa": 1784}


@pytest.mark.xfail_on_32bit("Uses too much memory")
@pytest.mark.parametrize("method, expected",
                         tests_against_itself_1samp.items())
def test_bootstrap_against_itself_1samp(method, expected):
    # The expected values in this test were generated using bootstrap
    # to check for unintended changes in behavior. The test also makes sure
    # that bootstrap works with multi-sample statistics and that the
    # `axis` argument works as expected / function is vectorized.
    np.random.seed(0)

    n = 100  # size of sample
    n_resamples = 1000  # number of bootstrap resamples used to form each CI
    confidence_level = 0.9

    # The true mean is 5
    dist = stats.norm(loc=5, scale=1)
    stat_true = dist.mean()

    # Do the same thing 1000 times. (The code is fully vectorized.)
    n_replications = 2000
    data = dist.rvs(size=(n_replications, n))
    res = bootstrap((data,),
                    statistic=np.mean,
                    confidence_level=confidence_level,
                    n_resamples=n_resamples,
                    method=method,
                    axis=-1)
    ci = res.confidence_interval

    # ci contains vectors of lower and upper confidence interval bounds
    ci_contains_true = np.sum((ci[0] < stat_true) & (stat_true < ci[1]))
    assert ci_contains_true == expected

    # ci_contains_true is not inconsistent with confidence_level
    pvalue = stats.binomtest(ci_contains_true, n_replications,
                             confidence_level).pvalue
    assert pvalue > 0.1


tests_against_itself_2samp = {"basic": 888,
                              "percentile": 886}


@pytest.mark.xfail_on_32bit("Uses too much memory")
@pytest.mark.parametrize("method, expected",
                         tests_against_itself_2samp.items())
def test_bootstrap_against_itself_2samp(method, expected):
    # The expected values in this test were generated using bootstrap
    # to check for unintended changes in behavior. The test also makes sure
    # that bootstrap works with multi-sample statistics and that the
    # `axis` argument works as expected / function is vectorized.
    np.random.seed(0)

    n1 = 100  # size of sample 1
    n2 = 120  # size of sample 2
    n_resamples = 1000  # number of bootstrap resamples used to form each CI
    confidence_level = 0.9

    # The statistic we're interested in is the difference in means
    def my_stat(data1, data2, axis=-1):
        mean1 = np.mean(data1, axis=axis)
        mean2 = np.mean(data2, axis=axis)
        return mean1 - mean2

    # The true difference in the means is -0.1
    dist1 = stats.norm(loc=0, scale=1)
    dist2 = stats.norm(loc=0.1, scale=1)
    stat_true = dist1.mean() - dist2.mean()

    # Do the same thing 1000 times. (The code is fully vectorized.)
    n_replications = 1000
    data1 = dist1.rvs(size=(n_replications, n1))
    data2 = dist2.rvs(size=(n_replications, n2))
    res = bootstrap((data1, data2),
                    statistic=my_stat,
                    confidence_level=confidence_level,
                    n_resamples=n_resamples,
                    method=method,
                    axis=-1)
    ci = res.confidence_interval

    # ci contains vectors of lower and upper confidence interval bounds
    ci_contains_true = np.sum((ci[0] < stat_true) & (stat_true < ci[1]))
    assert ci_contains_true == expected

    # ci_contains_true is not inconsistent with confidence_level
    pvalue = stats.binomtest(ci_contains_true, n_replications,
                             confidence_level).pvalue
    assert pvalue > 0.1


def test_jackknife_resample():
    shape = 3, 4, 5, 6
    np.random.seed(0)
    x = np.random.rand(*shape)
    y = _bootstrap._jackknife_resample(x)

    for i in range(shape[-1]):
        # each resample is indexed along second to last axis
        # (last axis is the one the statistic will be taken over / consumed)
        slc = y[..., i, :]
        expected = np.delete(x, i, axis=-1)

        assert np.array_equal(slc, expected)


@pytest.mark.parametrize("rng_name", ["RandomState", "default_rng"])
def test_bootstrap_resample(rng_name):
    rng = getattr(np.random, rng_name, None)
    if rng is None:
        pytest.skip(f"{rng_name} not available.")
    rng1 = rng(0)
    rng2 = rng(0)

    n_resamples = 10
    shape = 3, 4, 5, 6

    np.random.seed(0)
    x = np.random.rand(*shape)
    y = _bootstrap._bootstrap_resample(x, n_resamples, random_state=rng1)

    for i in range(n_resamples):
        # each resample is indexed along second to last axis
        # (last axis is the one the statistic will be taken over / consumed)
        slc = y[..., i, :]

        js = rng_integers(rng2, 0, shape[-1], shape[-1])
        expected = x[..., js]

        assert np.array_equal(slc, expected)


@pytest.mark.parametrize("score", [0, 0.5, 1])
@pytest.mark.parametrize("axis", [0, 1, 2])
def test_percentile_of_score(score, axis):
    shape = 10, 20, 30
    np.random.seed(0)
    x = np.random.rand(*shape)
    p = _bootstrap._percentile_of_score(x, score, axis=-1)

    def vectorized_pos(a, score, axis):
        return np.apply_along_axis(stats.percentileofscore, axis, a, score)

    p2 = vectorized_pos(x, score, axis=-1)/100

    assert_allclose(p, p2, 1e-15)


def test_percentile_along_axis():
    # the difference between _percentile_along_axis and np.percentile is that
    # np.percentile gets _all_ the qs for each axis slice, whereas
    # _percentile_along_axis gets the q corresponding with each axis slice

    shape = 10, 20
    np.random.seed(0)
    x = np.random.rand(*shape)
    q = np.random.rand(*shape[:-1]) * 100
    y = _bootstrap._percentile_along_axis(x, q)

    for i in range(shape[0]):
        res = y[i]
        expected = np.percentile(x[i], q[i], axis=-1)
        assert_allclose(res, expected, 1e-15)
