import math
import re
import pytest
import numpy as np
from numpy.testing import assert_, assert_allclose
from scipy._lib._array_api import make_xp_test_case, eager_warns
from scipy._lib._array_api_no_0d import xp_assert_close, xp_assert_equal, xp_assert_less
from scipy.stats._axis_nan_policy import (SmallSampleWarning, too_small_1d_not_omit,
                                          too_small_1d_omit, too_small_nd_omit)
from scipy import stats


@make_xp_test_case(stats.circmean, stats.circvar, stats.circstd)
class TestCircFuncs:
    # In gh-5747, the R package `circular` was used to calculate reference
    # values for the circular variance, e.g.:
    # library(circular)
    # options(digits=16)
    # x = c(0, 2*pi/3, 5*pi/3)
    # var.circular(x)
    @pytest.mark.parametrize("test_func,expected",
                             [(stats.circmean, 0.167690146),
                              (stats.circvar, 0.006455174000787767),
                              (stats.circstd, 6.520702116)])
    def test_circfuncs(self, test_func, expected, xp):
        x = xp.asarray([355., 5., 2., 359., 10., 350.])
        xp_assert_close(test_func(x, high=360), xp.asarray(expected))

    def test_circfuncs_small(self, xp):
        # Default tolerances won't work here because the reference values
        # are approximations. Ensure all array types work in float64 to
        # avoid needing separate float32 and float64 tolerances.
        x = xp.asarray([20, 21, 22, 18, 19, 20.5, 19.2], dtype=xp.float64)
        M1 = xp.mean(x)
        M2 = stats.circmean(x, high=360)
        xp_assert_close(M2, M1, rtol=1e-5)

        V1 = xp.var(x*xp.pi/180, correction=0)
        # for small variations, circvar is approximately half the
        # linear variance
        V1 = V1 / 2.
        V2 = stats.circvar(x, high=360)
        xp_assert_close(V2, V1, rtol=1e-4)

        S1 = xp.std(x, correction=0)
        S2 = stats.circstd(x, high=360)
        xp_assert_close(S2, S1, rtol=1e-4)

    @pytest.mark.parametrize("test_func, numpy_func",
                             [(stats.circmean, np.mean),
                              (stats.circvar, np.var),
                              (stats.circstd, np.std)])
    def test_circfuncs_close(self, test_func, numpy_func, xp):
        # circfuncs should handle very similar inputs (gh-12740)
        x = np.asarray([0.12675364631578953] * 10 + [0.12675365920187928] * 100)
        circstat = test_func(xp.asarray(x))
        normal = xp.asarray(numpy_func(x))
        xp_assert_close(circstat, normal, atol=2e-8)

    @pytest.mark.parametrize('circfunc', [stats.circmean,
                                          stats.circvar,
                                          stats.circstd])
    def test_circmean_axis(self, xp, circfunc):
        x = xp.asarray([[355, 5, 2, 359, 10, 350],
                        [351, 7, 4, 352, 9, 349],
                        [357, 9, 8, 358, 4, 356.]])
        res = circfunc(x, high=360)
        ref = circfunc(xp.reshape(x, (-1,)), high=360)
        xp_assert_close(res, xp.asarray(ref))

        res = circfunc(x, high=360, axis=1)
        ref = [circfunc(x[i, :], high=360) for i in range(x.shape[0])]
        xp_assert_close(res, xp.stack(ref))

        res = circfunc(x, high=360, axis=0)
        ref = [circfunc(x[:, i], high=360) for i in range(x.shape[1])]
        xp_assert_close(res, xp.stack(ref))

    @pytest.mark.parametrize("test_func,expected",
                             [(stats.circmean, 0.167690146),
                              (stats.circvar, 0.006455174270186603),
                              (stats.circstd, 6.520702116)])
    def test_circfuncs_array_like(self, test_func, expected, xp):
        x = xp.asarray([355, 5, 2, 359, 10, 350.])
        xp_assert_close(test_func(x, high=360), xp.asarray(expected))

    @pytest.mark.parametrize("test_func", [stats.circmean, stats.circvar,
                                           stats.circstd])
    def test_empty(self, test_func, xp):
        dtype = xp.float64
        x = xp.asarray([], dtype=dtype)
        with eager_warns(SmallSampleWarning, match=too_small_1d_not_omit, xp=xp):
            res = test_func(x)

        xp_assert_equal(res, xp.asarray(xp.nan, dtype=dtype))

    @pytest.mark.parametrize("test_func", [stats.circmean, stats.circvar,
                                           stats.circstd])
    def test_nan_propagate(self, test_func, xp):
        x = xp.asarray([355, 5, 2, 359, 10, 350, np.nan])
        xp_assert_equal(test_func(x, high=360), xp.asarray(xp.nan))

    @pytest.mark.parametrize("test_func,expected",
                             [(stats.circmean,
                               {None: np.nan, 0: 355.66582264, 1: 0.28725053}),
                              (stats.circvar,
                               {None: np.nan,
                                0: 0.002570671054089924,
                                1: 0.005545914017677123}),
                              (stats.circstd,
                               {None: np.nan, 0: 4.11093193, 1: 6.04265394})])
    def test_nan_propagate_array(self, test_func, expected, xp):
        x = xp.asarray([[355, 5, 2, 359, 10, 350, 1],
                        [351, 7, 4, 352, 9, 349, np.nan],
                        [1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]])
        for axis in expected.keys():
            out = test_func(x, high=360, axis=axis)
            if axis is None:
                xp_assert_equal(out, xp.asarray(xp.nan))
            else:
                xp_assert_close(out[0], xp.asarray(expected[axis]))
                xp_assert_equal(out[1:], xp.full_like(out[1:], xp.nan))

    def test_circmean_scalar(self, xp):
        x = xp.asarray(1.)[()]
        M1 = x
        M2 = stats.circmean(x)
        xp_assert_close(M2, M1, rtol=1e-5)

    def test_circmean_range(self, xp):
        # regression test for gh-6420: circmean(..., high, low) must be
        # between `high` and `low`
        m = stats.circmean(xp.arange(0, 2, 0.1), xp.pi, -xp.pi)
        xp_assert_less(m, xp.asarray(xp.pi))
        xp_assert_less(-m, xp.asarray(xp.pi))

    def test_circfuncs_uint8(self, xp):
        # regression test for gh-7255: overflow when working with
        # numpy uint8 data type
        x = xp.asarray([150, 10], dtype=xp.uint8)
        xp_assert_close(stats.circmean(x, high=180), xp.asarray(170.0))
        xp_assert_close(stats.circvar(x, high=180), xp.asarray(0.2339555554617))
        xp_assert_close(stats.circstd(x, high=180), xp.asarray(20.91551378))

    def test_circstd_zero(self, xp):
        # circstd() of a single number should return positive zero.
        y = stats.circstd(xp.asarray([0]))
        assert math.copysign(1.0, y) == 1.0

    def test_circmean_accuracy_tiny_input(self, xp):
        # For tiny x such that sin(x) == x and cos(x) == 1.0 numerically,
        # circmean(x) should return x because atan2(sin(x), cos(x)) == x.
        # This test verifies this.
        #
        # The purpose of this test is not to show that circmean() is
        # accurate in the last digit for certain input, because this is
        # neither guaranteed not particularly useful.  Rather, it is a
        # "white-box" sanity check that no undue loss of precision is
        # introduced by conversion between (high - low) and (2 * pi).

        x = xp.linspace(1e-9, 6e-9, 50)
        assert xp.all(xp.sin(x) == x) and xp.all(xp.cos(x) == 1.0)

        m = (x * (2 * xp.pi) / (2 * xp.pi)) != x
        assert xp.any(m)
        x = x[m]

        y = stats.circmean(x[:, None], axis=1)
        assert xp.all(y == x)

    def test_circmean_accuracy_huge_input(self, xp):
        # White-box test that circmean() does not introduce undue loss of
        # numerical accuracy by eagerly rotating the input.  This is detected
        # by supplying a huge input x such that (x - low) == x numerically.
        x = xp.asarray(1e17, dtype=xp.float64)
        y = math.atan2(xp.sin(x), xp.cos(x))  # -2.6584887370946806
        expected = xp.asarray(y, dtype=xp.float64)
        actual = stats.circmean(x, high=xp.pi, low=-xp.pi)
        xp_assert_close(actual, expected, rtol=1e-15, atol=0.0)


class TestCircFuncsNanPolicy:
    # `nan_policy` is implemented by the `_axis_nan_policy` decorator, which is
    # not yet array-API compatible. When it is array-API compatible, the generic
    # tests run on every function will be much stronger than these, so these
    # will not be necessary. So I don't see a need to make these array-API compatible;
    # when the time comes, they can just be removed.
    @pytest.mark.parametrize("test_func,expected",
                             [(stats.circmean,
                               {None: 359.4178026893944,
                                0: np.array([353.0, 6.0, 3.0, 355.5, 9.5,
                                             349.5]),
                                1: np.array([0.16769015, 358.66510252])}),
                              (stats.circvar,
                               {None: 0.008396678483192477,
                                0: np.array([1.9997969, 0.4999873, 0.4999873,
                                             6.1230956, 0.1249992, 0.1249992]
                                            )*(np.pi/180)**2,
                                1: np.array([0.006455174270186603,
                                             0.01016767581393285])}),
                              (stats.circstd,
                               {None: 7.440570778057074,
                                0: np.array([2.00020313, 1.00002539, 1.00002539,
                                             3.50108929, 0.50000317,
                                             0.50000317]),
                                1: np.array([6.52070212, 8.19138093])})])
    def test_nan_omit_array(self, test_func, expected):
        x = np.array([[355, 5, 2, 359, 10, 350, np.nan],
                      [351, 7, 4, 352, 9, 349, np.nan],
                      [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]])
        for axis in expected.keys():
            if axis is None:
                out = test_func(x, high=360, nan_policy='omit', axis=axis)
                assert_allclose(out, expected[axis], rtol=1e-7)
            else:
                with pytest.warns(SmallSampleWarning, match=too_small_nd_omit):
                    out = test_func(x, high=360, nan_policy='omit', axis=axis)
                    assert_allclose(out[:-1], expected[axis], rtol=1e-7)
                    assert_(np.isnan(out[-1]))

    @pytest.mark.parametrize("test_func,expected",
                             [(stats.circmean, 0.167690146),
                              (stats.circvar, 0.006455174270186603),
                              (stats.circstd, 6.520702116)])
    def test_nan_omit(self, test_func, expected):
        x = [355, 5, 2, 359, 10, 350, np.nan]
        assert_allclose(test_func(x, high=360, nan_policy='omit'),
                        expected, rtol=1e-7)

    @pytest.mark.parametrize("test_func", [stats.circmean, stats.circvar,
                                           stats.circstd])
    def test_nan_omit_all(self, test_func):
        x = [np.nan, np.nan, np.nan, np.nan, np.nan]
        with pytest.warns(SmallSampleWarning, match=too_small_1d_omit):
            assert_(np.isnan(test_func(x, nan_policy='omit')))

    @pytest.mark.parametrize("test_func", [stats.circmean, stats.circvar,
                                           stats.circstd])
    def test_nan_omit_all_axis(self, test_func):
        with pytest.warns(SmallSampleWarning, match=too_small_nd_omit):
            x = np.array([[np.nan, np.nan, np.nan, np.nan, np.nan],
                          [np.nan, np.nan, np.nan, np.nan, np.nan]])
            out = test_func(x, nan_policy='omit', axis=1)
            assert_(np.isnan(out).all())
            assert_(len(out) == 2)

    @pytest.mark.parametrize("x",
                             [[355, 5, 2, 359, 10, 350, np.nan],
                              np.array([[355, 5, 2, 359, 10, 350, np.nan],
                                        [351, 7, 4, 352, np.nan, 9, 349]])])
    @pytest.mark.parametrize("test_func", [stats.circmean, stats.circvar,
                                           stats.circstd])
    def test_nan_raise(self, test_func, x):
        with pytest.raises(ValueError):
            test_func(x, high=360, nan_policy='raise')

    @pytest.mark.parametrize("x",
                             [[355, 5, 2, 359, 10, 350, np.nan],
                              np.array([[355, 5, 2, 359, 10, 350, np.nan],
                                        [351, 7, 4, 352, np.nan, 9, 349]])])
    @pytest.mark.parametrize("test_func", [stats.circmean, stats.circvar,
                                           stats.circstd])
    def test_bad_nan_policy(self, test_func, x):
        with pytest.raises(ValueError):
            test_func(x, high=360, nan_policy='foobar')



@make_xp_test_case(stats.directional_stats)
class TestDirectionalStats:
    # Reference implementations are not available
    def test_directional_stats_correctness(self, xp):
        # Data from Fisher: Dispersion on a sphere, 1953 and
        # Mardia and Jupp, Directional Statistics.
        decl = -np.deg2rad(np.array([343.2, 62., 36.9, 27., 359.,
                                     5.7, 50.4, 357.6, 44.]))
        incl = -np.deg2rad(np.array([66.1, 68.7, 70.1, 82.1, 79.5,
                                     73., 69.3, 58.8, 51.4]))
        data = np.stack((np.cos(incl) * np.cos(decl),
                         np.cos(incl) * np.sin(decl),
                         np.sin(incl)),
                        axis=1)

        decl = xp.asarray(decl.tolist())
        incl = xp.asarray(incl.tolist())
        data = xp.asarray(data.tolist())

        dirstats = stats.directional_stats(data)
        directional_mean = dirstats.mean_direction

        reference_mean = xp.asarray([0.2984, -0.1346, -0.9449])
        xp_assert_close(directional_mean, reference_mean, atol=1e-4)

    @pytest.mark.parametrize('angles, ref', [
        ([-np.pi/2, np.pi/2], 1.),
        ([0, 2 * np.pi], 0.)
    ])
    def test_directional_stats_2d_special_cases(self, angles, ref, xp):
        angles = xp.asarray(angles)
        ref = xp.asarray(ref)
        data = xp.stack([xp.cos(angles), xp.sin(angles)], axis=1)
        res = 1 - stats.directional_stats(data).mean_resultant_length
        xp_assert_close(res, ref)

    def test_directional_stats_2d(self, xp):
        # Test that for circular data directional_stats
        # yields the same result as circmean/circvar
        rng = np.random.default_rng(0xec9a6899d5a2830e0d1af479dbe1fd0c)
        testdata = xp.asarray(2 * xp.pi * rng.random((1000, )))
        testdata_vector = xp.stack((xp.cos(testdata),
                                    xp.sin(testdata)),
                                   axis=1)
        dirstats = stats.directional_stats(testdata_vector)
        directional_mean = dirstats.mean_direction
        directional_mean_angle = xp.atan2(directional_mean[1], directional_mean[0])
        directional_mean_angle = directional_mean_angle % (2 * xp.pi)
        circmean = stats.circmean(testdata)
        xp_assert_close(directional_mean_angle, circmean)

        directional_var = 1. - dirstats.mean_resultant_length
        circular_var = stats.circvar(testdata)
        xp_assert_close(directional_var, circular_var)

    def test_directional_mean_higher_dim(self, xp):
        # test that directional_stats works for higher dimensions
        # here a 4D array is reduced over axis = 2
        data = xp.asarray([[0.8660254, 0.5, 0.],
                           [0.8660254, -0.5, 0.]])
        full_array = xp.asarray(xp.tile(data, (2, 2, 2, 1)))
        expected = xp.asarray([[[1., 0., 0.],
                                [1., 0., 0.]],
                               [[1., 0., 0.],
                                [1., 0., 0.]]])
        dirstats = stats.directional_stats(full_array, axis=2)
        xp_assert_close(dirstats.mean_direction, expected)

    @pytest.mark.skip_xp_backends(np_only=True, reason='checking array-like input')
    def test_directional_stats_list_ndarray_input(self, xp):
        # test that list and numpy array inputs yield same results
        data = [[0.8660254, 0.5, 0.], [0.8660254, -0.5, 0]]
        data_array = xp.asarray(data, dtype=xp.float64)
        ref = stats.directional_stats(data)
        res = stats.directional_stats(data_array)
        xp_assert_close(res.mean_direction,
                        xp.asarray(ref.mean_direction))
        xp_assert_close(res.mean_resultant_length,
                        xp.asarray(res.mean_resultant_length))

    def test_directional_stats_1d_error(self, xp):
        # test that one-dimensional data raises ValueError
        data = xp.ones((5, ))
        message = (r"samples must at least be two-dimensional. "
                   r"Instead samples has shape: (5,)")
        with pytest.raises(ValueError, match=re.escape(message)):
            stats.directional_stats(data)

    @pytest.mark.parametrize("dtype", ["float32", "float64"])
    def test_directional_stats_normalize(self, dtype, xp):
        # test that directional stats calculations yield same results
        # for unnormalized input with normalize=True and normalized
        # input with normalize=False
        data = np.array([[0.8660254, 0.5, 0.],
                         [1.7320508, -1., 0.]], dtype=dtype)
        res = stats.directional_stats(xp.asarray(data), normalize=True)
        normalized_data = data / np.linalg.norm(data, axis=-1,
                                                keepdims=True)
        ref = stats.directional_stats(normalized_data, normalize=False)
        xp_assert_close(res.mean_direction,
                        xp.asarray(ref.mean_direction))
        xp_assert_close(res.mean_resultant_length,
                        xp.asarray(ref.mean_resultant_length))
