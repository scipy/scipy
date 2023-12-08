import numpy as np
from numpy.testing import assert_equal, assert_array_equal, assert_array_almost_equal, assert_allclose
import pytest
from pytest import raises as assert_raises

from scipy.interpolate import InverseDistanceWeightedNDInterpolator


class TestInverseDistanceWeightedNDInterpolator:
    def _get_sample_data(self):
        points = np.array([[1., 1], [3, 3], [5, 5], [7, 7], [9, 9]])
        values = np.array([1., 3, 5, 7, 9])
        return points, values

    def test_initialization_with_valid_data(self):
        points, values = self._get_sample_data()
        interpolator = InverseDistanceWeightedNDInterpolator(points, values)
        assert interpolator is not None, "Failed to initialize with valid data"
        
    def test_initialization_parameters(self):
        points, values = self._get_sample_data()
        interpolator = InverseDistanceWeightedNDInterpolator(points, values, rescale=True, fill_value=0, local=False)
        assert interpolator._local is False, "Initialization parameter 'local' not set correctly"
        assert interpolator._fill_value == 0, "Initialization parameter 'fill_value' not set correctly"

    def test_fill_value_usage(self):
        points, values = self._get_sample_data()
        fill_value = -1
        interpolator = InverseDistanceWeightedNDInterpolator(points, values, fill_value=fill_value)
        query_points = np.array([[2, 2]])
        # Outside the convex hull of the data points
        interpolated_values = interpolator(query_points, distance_upper_bound=1e-7)
        assert_array_equal(interpolated_values, fill_value, err_msg="Incorrect usage of fill_value")

    def test_local_global_flag(self):
        points, values = self._get_sample_data()
        local_interpolator = InverseDistanceWeightedNDInterpolator(points, values, local=True)
        global_interpolator = InverseDistanceWeightedNDInterpolator(points, values, local=False)
        assert local_interpolator._local is True, "Local flag not set correctly"
        assert global_interpolator._local is False, "Global flag not set correctly"

    def test_input_dimensionality(self):
        points = np.random.rand(10, 3)
        values = np.random.rand(10)
        interpolator = InverseDistanceWeightedNDInterpolator(points, values)
        query_points = np.random.rand(5, 3)  # Matching dimensionality
        ip = interpolator(query_points)
        assert_equal(query_points.shape[0], len(ip))



    def test_list_input(self):
        points, values = self._get_sample_data()
        query_points = np.array([[1, 1], [2, 2]])

        list_interpolator = InverseDistanceWeightedNDInterpolator(points, values.tolist())
        norm_interpolator = InverseDistanceWeightedNDInterpolator(points, values)
        assert_allclose(list_interpolator(query_points.tolist()), norm_interpolator(query_points))

    def test_basic_functionality(self):
        points, values = self.get_sample_data()
        interpolator = InverseDistanceWeightedNDInterpolator(points, values)
        query_points = np.array([[1, 1], [2, 2]])
        expected_values = np.array([1, 2])
        interpolated_values = interpolator(query_points, k=2)  # k=2, it will consider [1,1] and [3,3]
        assert_array_almost_equal(interpolated_values, expected_values, err_msg="Basic functionality failed")


    def test_basic_functionality_power_1(self):
        points, values = self.get_sample_data()
        interpolator = InverseDistanceWeightedNDInterpolator(points, values)
        query_points = np.array([[1, 1], [2, 2]])
        expected_values = np.array([1, 2])
        interpolated_values = interpolator(query_points, k=2, p=3)
        np.testing.assert_array_almost_equal(interpolated_values, expected_values, err_msg="Basic functionality failed")
