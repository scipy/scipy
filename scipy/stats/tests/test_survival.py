import pytest
import numpy as np
from numpy.testing import assert_equal
from scipy import stats

class TestSurvival:
    def test_input_validation(self):
        message = '`sample` must be a one-dimensional sequence.'
        with pytest.raises(ValueError, match=message):
            stats.ecdf([[1]])
        with pytest.raises(ValueError, match=message):
            stats.ecdf(1)

        message = '`sample` must not contain nan'
        with pytest.raises(ValueError, match=message):
            stats.ecdf([np.nan])

        message = 'Currently, only uncensored data is supported.'
        with pytest.raises(NotImplementedError, match=message):
            stats.ecdf(stats.CensoredData.left_censored([1], censored=[True]))

    def test_edge_cases(self):
        res = stats.ecdf([])
        assert_equal(res.x, [])
        assert_equal(res.cdf, [])

        res = stats.ecdf([1])
        assert_equal(res.x, [1])
        assert_equal(res.cdf, [1])

    def test_unique(self):
        # Example with unique observations; `stats.ecdf` ref. [1] page 80
        sample = [6.23, 5.58, 7.06, 6.42, 5.20]
        res = stats.ecdf(sample)
        ref_x = np.sort(np.unique(sample))
        ref_cdf = np.arange(1, 6) / 5
        ref_sf = 1 - ref_cdf
        assert_equal(res.x, ref_x)
        assert_equal(res.cdf, ref_cdf)
        assert_equal(res.sf, ref_sf)

    def test_nonunique(self):
        # Example with non-unique observations; `stats.ecdf` ref. [1] page 82
        sample = [0, 2, 1, 2, 3, 4]
        res = stats.ecdf(sample)
        ref_x = np.sort(np.unique(sample))
        ref_cdf = np.array([1/6, 2/6, 4/6, 5/6, 1])
        ref_sf = 1 - ref_cdf
        assert_equal(res.x, ref_x)
        assert_equal(res.cdf, ref_cdf)
        assert_equal(res.sf, ref_sf)
