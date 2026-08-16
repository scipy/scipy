import math
from math import pi
import re
import pytest
import numpy as np
from numpy.testing import assert_, assert_allclose
from scipy._lib._array_api import make_xp_test_case, eager_warns
from scipy._lib._array_api_no_0d import xp_assert_close, xp_assert_equal, xp_assert_less
from scipy.stats._axis_nan_policy import (SmallSampleWarning, too_small_1d_not_omit,
                                          too_small_1d_omit, too_small_nd_omit)
from scipy import stats

lazy_xp_modules = [stats]


@make_xp_test_case(stats.circmedian)
class TestCircMedian:
    # Example data from references; e.g. `example_3_2` is [3] Example #2
    example_2 = [43, 45, 52, 61, 75, 88, 88, 279, 375]
    example_3_1 = [-9 * pi / 16, -9 * pi / 16, 0, 9 * pi / 16, 9 * pi / 16]
    example_3_2 = [-pi / 8, pi / 8, pi, pi]
    example_3_3 = [-3 * pi / 8, 0, 2 * pi / 3]
    example_3_4 = [0, 2 * pi / 3, -2 * pi / 3]
    example_3_5 = [0, pi / 6, 2 * pi / 3, 5 * pi / 6, -2 * pi / 3, -3 * pi / 6]
    # [7] is Otieno & Cook "A More Efficient Way Of Obtaining A Unique Median..."
    # [8] is Otieno & Cook "Measures of preferred direction for environmental..."
    # In the paper, 144 is listed instead of 145, but 145 seems to be used in analysis
    example_7_1 = [104, 110, 117, 121, 127, 130, 136, 145, 152, 178, 184, 192, 200, 316]
    example_8_1 = [2, 9, 18, 24, 30, 35, 35, 39, 39, 44, 44, 49, 56, 70, 76, 76, 81, 86,
                   91, 112, 121, 127, 133, 134, 138, 147, 152, 157, 166, 171, 177, 187,
                   206, 210, 211, 215, 238, 246, 269, 270, 285, 292, 305, 315, 325, 328,
                   329, 343, 354, 359]
    example_8_2 = [85, 135, 135, 140, 145, 150, 150, 150,
                   160, 285, 200, 210, 220, 225, 270]
    example_8_3 = [67, 66, 74, 61, 58, 60, 100, 89, 171, 166, 98, 60, 197, 98, 86, 123,
                   165, 133, 101, 105, 71, 84, 75, 98, 83, 71, 74, 91, 38, 200, 56]
    _bisect = dict(convention='bisecting')
    _deg_bisect = dict(high=360, convention='bisecting')

    @pytest.mark.parametrize("sample, ref, kwargs", [
        (example_2, 52., _deg_bisect),
        (example_3_1, 0.0, {}),
        (example_3_1, pi, _bisect),
        (example_3_2, pi, {}),
        (example_3_3, 0.0, {}),
        (example_3_3, float(stats.circmean([-3*pi/8, 0, 2*pi/3 + pi])), _bisect),
        (example_3_4, math.nan, {}),
        (example_3_4, math.nan, _bisect),
        (example_3_5, math.nan, {}),
        (example_3_5, math.nan, _bisect),
        pytest.param(example_7_1, 46.5, _deg_bisect, marks=pytest.mark.xfail),
        pytest.param(example_8_1, 136.75, _deg_bisect, marks=pytest.mark.xfail),
        (example_8_2, 150., _deg_bisect),
        (example_8_3, 86., _deg_bisect),
    ])
    def test_against_references(self, sample, ref, kwargs, xp):
        # Test against examples from references
        res = stats.circmedian(xp.asarray(sample), **kwargs)
        xp_assert_close(res, xp.asarray(ref))

    @pytest.mark.parametrize("sample, ref", [
        (example_2, 52.00000000000002), (example_7_1, 136.99714451257873),
        (example_8_1, 45.6661962795516), (example_8_2, 150.00000000000003),
        (example_8_3, 86.0),
    ])
    def test_against_pycircstat2_method_deviation(self, sample, ref, xp):
        # test against pycircstat2 method='deviation', e.g.
        # from pycircstat2.descriptive import circ_median
        # d = np.deg2rad(example_5_1)
        # np.rad2deg(circ_median(alpha=d, method='deviation'))
        res = stats.circmedian(xp.asarray(sample), high=360)
        xp_assert_close(res, xp.asarray(ref))

    @pytest.mark.parametrize("sample, ref", [
        (example_2, 52.00000000000002),
        pytest.param(example_7_1, math.nan, marks=pytest.mark.xfail),
        pytest.param(example_8_1, 52.86336626974574, marks=pytest.mark.xfail),
        pytest.param(example_8_2, 152.495231278463, marks=pytest.mark.xfail),
        (example_8_3, 86.0),
    ])
    def test_against_pycircstat2_method_count(self, sample, ref, xp):
        # test against pycircstat2 method='count', e.g.
        # from pycircstat2.descriptive import circ_median
        # d = np.deg2rad(example_5_1)
        # np.rad2deg(circ_median(alpha=d, method='count'))
        res = stats.circmedian(xp.asarray(sample), high=360, convention='bisecting')
        xp_assert_close(res, xp.asarray(ref))

    def test_input_validation(self):
        message = "`convention` must be either 'arc-distance' or 'bisecting'."
        with pytest.raises(ValueError, match=message):
            stats.circmedian([1, 2, 3], convention='ekki-ekki')


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
