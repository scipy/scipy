import pytest
import numpy as np
from numpy.testing import assert_equal, assert_allclose
from scipy import stats
from scipy.stats import _survival


def _kaplan_meier_reference(times, censored):
    # This is a very straightforward implementation of the Kaplan-Meier
    # estimator that does almost everything differently from the implementation
    # in stats.ecdf.

    # Begin by sorting the raw data. Note that the order of death and loss
    # at a given time matters: death happens first. See [2] page 461:
    # "These conventions may be paraphrased by saying that deaths recorded as
    # of an age t are treated as if they occurred slightly before t, and losses
    # recorded as of an age t are treated as occurring slightly after t."
    # We implement this by sorting the data first by time, then by `censored`,
    # (which is 0 when there is a death and 1 when there is only a loss).
    dtype = [('time', float), ('censored', int)]
    data = np.array([(t, d) for t, d in zip(times, censored)], dtype=dtype)
    data = np.sort(data, order=('time', 'censored'))
    times = data['time']
    died = np.logical_not(data['censored'])

    m = times.size
    n = np.arange(m, 0, -1)  # number at risk
    sf = np.cumprod((n - died) / n)

    # Find the indices of the *last* occurrence of unique times. The
    # corresponding entries of `times` and `sf` are what we want.
    _, indices = np.unique(times[::-1], return_index=True)
    ref_times = times[-indices - 1]
    ref_sf = sf[-indices - 1]
    return ref_times, ref_sf


class TestSurvival:

    @staticmethod
    def get_random_sample(rng, n_unique):
        # generate random sample
        unique_times = rng.random(n_unique)
        # convert to `np.int32` to resolve `np.repeat` failure in 32-bit CI
        repeats = rng.integers(1, 4, n_unique).astype(np.int32)
        times = rng.permuted(np.repeat(unique_times, repeats))
        censored = rng.random(size=times.size) > rng.random()
        sample = stats.CensoredData.right_censored(times, censored)
        return sample, times, censored

    def test_input_validation(self):
        message = '`sample` must be a one-dimensional sequence.'
        with pytest.raises(ValueError, match=message):
            stats.ecdf([[1]])
        with pytest.raises(ValueError, match=message):
            stats.ecdf(1)

        message = '`sample` must not contain nan'
        with pytest.raises(ValueError, match=message):
            stats.ecdf([np.nan])

        message = 'Currently, only uncensored and right-censored data...'
        with pytest.raises(NotImplementedError, match=message):
            stats.ecdf(stats.CensoredData.left_censored([1], censored=[True]))

        message = 'method` must be one of...'
        res = stats.ecdf([1, 2, 3])
        with pytest.raises(ValueError, match=message):
            res.cdf.confidence_interval(method='ekki-ekki')
        with pytest.raises(ValueError, match=message):
            res.sf.confidence_interval(method='shrubbery')

        message = 'confidence_level` must be a scalar between 0 and 1'
        with pytest.raises(ValueError, match=message):
            res.cdf.confidence_interval(-1)
        with pytest.raises(ValueError, match=message):
            res.sf.confidence_interval([0.5, 0.6])

    def test_edge_cases(self):
        res = stats.ecdf([])
        assert_equal(res.x, [])
        assert_equal(res.cdf.points, [])

        res = stats.ecdf([1])
        assert_equal(res.x, [1])
        assert_equal(res.cdf.points, [1])

    def test_unique(self):
        # Example with unique observations; `stats.ecdf` ref. [1] page 80
        sample = [6.23, 5.58, 7.06, 6.42, 5.20]
        res = stats.ecdf(sample)
        ref_x = np.sort(np.unique(sample))
        ref_cdf = np.arange(1, 6) / 5
        ref_sf = 1 - ref_cdf
        assert_equal(res.x, ref_x)
        assert_equal(res.cdf.points, ref_cdf)
        assert_equal(res.sf.points, ref_sf)

    def test_nonunique(self):
        # Example with non-unique observations; `stats.ecdf` ref. [1] page 82
        sample = [0, 2, 1, 2, 3, 4]
        res = stats.ecdf(sample)
        ref_x = np.sort(np.unique(sample))
        ref_cdf = np.array([1/6, 2/6, 4/6, 5/6, 1])
        ref_sf = 1 - ref_cdf
        assert_equal(res.x, ref_x)
        assert_equal(res.cdf.points, ref_cdf)
        assert_equal(res.sf.points, ref_sf)

    def test_call_methods(self):
        # Test CDF and SF call methods
        rng = np.random.default_rng(1162729143302572461)
        sample, _, _ = self.get_random_sample(rng, 15)
        res = stats.ecdf(sample)
        x = res.x
        xr = res.x + np.diff(res.x, append=x[-1]+1)/2  # right shifted points

        assert_equal(res.cdf(x), res.cdf.points)
        assert_equal(res.cdf(xr), res.cdf.points)
        assert_equal(res.cdf(x[0]-1), 0)  # CDF starts at 0
        assert_equal(res.cdf([-np.inf, np.inf]), [0, 1])

        assert_equal(res.sf(x), res.sf.points)
        assert_equal(res.sf(xr), res.sf.points)
        assert_equal(res.sf(x[0]-1), 1)  # CDF starts at 1
        assert_equal(res.sf([-np.inf, np.inf]), [1, 0])

    # ref. [1] page 91
    t1 = [37, 43, 47, 56, 60, 62, 71, 77, 80, 81]  # times
    d1 = [0, 0, 1, 1, 0, 0, 0, 1, 1, 1]  # 1 means deaths (not censored)
    r1 = [1, 1, 0.875, 0.75, 0.75, 0.75, 0.75, 0.5, 0.25, 0]  # reference SF

    # https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_survival/BS704_Survival5.html  # noqa
    t2 = [8, 12, 26, 14, 21, 27, 8, 32, 20, 40]
    d2 = [1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
    r2 = [0.9, 0.788, 0.675, 0.675, 0.54, 0.405, 0.27, 0.27, 0.27]
    t3 = [33, 28, 41, 48, 48, 25, 37, 48, 25, 43]
    d3 = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
    r3 = [1, 0.875, 0.75, 0.75, 0.6, 0.6, 0.6]

    # https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_survival/bs704_survival4.html  # noqa
    t4 = [24, 3, 11, 19, 24, 13, 14, 2, 18, 17,
          24, 21, 12, 1, 10, 23, 6, 5, 9, 17]
    d4 = [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1]
    r4 = [0.95, 0.95, 0.897, 0.844, 0.844, 0.844, 0.844, 0.844, 0.844,
          0.844, 0.76, 0.676, 0.676, 0.676, 0.676, 0.507, 0.507]

    # https://www.real-statistics.com/survival-analysis/kaplan-meier-procedure/confidence-interval-for-the-survival-function/  # noqa
    t5 = [3, 5, 8, 10, 5, 5, 8, 12, 15, 14, 2, 11, 10, 9, 12, 5, 8, 11]
    d5 = [1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1]
    r5 = [0.944, 0.889, 0.722, 0.542, 0.542, 0.542, 0.361, 0.181, 0.181, 0.181]

    @pytest.mark.parametrize("case", [(t1, d1, r1), (t2, d2, r2), (t3, d3, r3),
                                      (t4, d4, r4), (t5, d5, r5)])
    def test_right_censored_against_examples(self, case):
        # test `ecdf` against other implementations on example problems
        times, died, ref = case
        sample = stats.CensoredData.right_censored(times, np.logical_not(died))
        res = stats.ecdf(sample)
        assert_allclose(res.sf.points, ref, atol=1e-3)
        assert_equal(res.x, np.sort(np.unique(times)))

        # test reference implementation against other implementations
        res = _kaplan_meier_reference(times, np.logical_not(died))
        assert_equal(res[0], np.sort(np.unique(times)))
        assert_allclose(res[1], ref, atol=1e-3)

    @pytest.mark.parametrize('seed', [182746786639392128, 737379171436494115,
                                      576033618403180168, 308115465002673650])
    def test_right_censored_against_reference_implementation(self, seed):
        # test `ecdf` against reference implementation on random problems
        rng = np.random.default_rng(seed)
        n_unique = rng.integers(10, 100)
        sample, times, censored = self.get_random_sample(rng, n_unique)
        res = stats.ecdf(sample)
        ref = _kaplan_meier_reference(times, censored)
        assert_allclose(res.x, ref[0])
        assert_allclose(res.sf.points, ref[1])

        # If all observations are uncensored, the KM estimate should match
        # the usual estimate for uncensored data
        sample = stats.CensoredData(uncensored=times)
        res = _survival._ecdf_right_censored(sample)  # force Kaplan-Meier
        ref = stats.ecdf(times)
        assert_equal(res[0], ref.x)
        assert_allclose(res[1], ref.cdf.points, rtol=1e-14)
        assert_allclose(res[2], ref.sf.points, rtol=1e-14)

    def test_right_censored_ci(self):
        # test "greenwood" confidence interval against example 5 (URL above).
        times, died = self.t4, self.d4
        sample = stats.CensoredData.right_censored(times, np.logical_not(died))
        res = stats.ecdf(sample)
        ref_allowance = [0.096, 0.096, 0.135, 0.162, 0.162, 0.162, 0.162,
                         0.162, 0.162, 0.162, 0.214, 0.246, 0.246, 0.246,
                         0.246, 0.341, 0.341]

        sf_ci = res.sf.confidence_interval()
        cdf_ci = res.cdf.confidence_interval()
        allowance = res.sf.points - sf_ci.low

        assert_allclose(allowance, ref_allowance, atol=1e-3)
        assert_allclose(sf_ci.low, np.clip(res.sf.points - allowance, 0, 1))
        assert_allclose(sf_ci.high, np.clip(res.sf.points + allowance, 0, 1))
        assert_allclose(cdf_ci.low, np.clip(res.cdf.points - allowance, 0, 1))
        assert_allclose(cdf_ci.high, np.clip(res.cdf.points + allowance, 0, 1))

    def test_right_censored_exponential_greenwood_ci(self):
        # test "exponential greenwood" confidence interval against example 5
        times, died = self.t5, self.d5
        sample = stats.CensoredData.right_censored(times, np.logical_not(died))
        res = stats.ecdf(sample)
        lower = np.array([0.66639, 0.624174, 0.456179, 0.287822, 0.287822,
                          0.287822, 0.128489, 0.030957, 0.030957, 0.030957])
        upper = np.array([0.991983, 0.970995, 0.87378, 0.739467, 0.739467,
                          0.739467, 0.603133, 0.430365, 0.430365, 0.430365])

        sf_ci = res.sf.confidence_interval(method='exponential greenwood')
        cdf_ci = res.cdf.confidence_interval(method='exponential greenwood')

        assert_allclose(sf_ci.low, lower, atol=1e-5)
        assert_allclose(sf_ci.high, upper, atol=1e-5)
        assert_allclose(cdf_ci.low, 1-upper, atol=1e-5)
        assert_allclose(cdf_ci.high, 1-lower, atol=1e-5)

    def test_right_censored_ci_nans(self):
        # test `ecdf` confidence interval (with NaNs) against Matlab
        times, died = self.t1, self.d1
        sample = stats.CensoredData.right_censored(times, np.logical_not(died))
        res = stats.ecdf(sample)

        # Reference values generated with Matlab
        # format long
        # t = [37 43 47 56 60 62 71 77 80 81];
        # d = [0 0 1 1 0 0 0 1 1 1];
        # censored = ~d1;
        # [f, x, flo, fup] = ecdf(t, 'Censoring', censored, 'Function',
        #                        'survivor', 'Alpha', 0.05);
        x = [37, 47, 56, 77, 80, 81]
        flo = [np.nan, 0.64582769623, 0.449943020228, 0.05270146407, 0, np.nan]
        fup = [np.nan, 1.0, 1.0, 0.947298535929289, 0.662388873768210, np.nan]
        i = np.searchsorted(res.x, x)

        message = "The variance estimate is undefined at some observations"
        with pytest.warns(RuntimeWarning, match=message):
            ci = res.sf.confidence_interval()

        # Matlab gives NaN as the first element of the CIs. It makes some
        # sense, but it's not what the formula gives, so skip it.
        assert_allclose(ci.low[i][1:], flo[1:])
        assert_allclose(ci.high[i][1:], fup[1:])

        # [f, x, flo, fup] = ecdf(t, 'Censoring', censored, 'Alpha', 0.05);
        flo = [np.nan, 0, 0, 0.052701464070711, 0.337611126231790, np.nan]
        fup = [np.nan, 0.35417230377, 0.5500569798, 0.9472985359, 1.0, np.nan]
        i = np.searchsorted(res.x, x)

        with pytest.warns(RuntimeWarning, match=message):
            ci = res.cdf.confidence_interval()

        assert_allclose(ci.low[i][1:], flo[1:])
        assert_allclose(ci.high[i][1:], fup[1:])
