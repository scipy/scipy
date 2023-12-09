import numpy as np
from numpy.testing import (assert_equal, assert_array_equal,
                             assert_array_almost_equal,
                           assert_allclose, assert_array_less)
import pytest
from pytest import raises as assert_raises

from scipy.interpolate import InverseDistanceWeightedNDInterpolator


class TestInverseDistanceWeightedNDInterpolator:
    def __init__(self):
        pass

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
        interpolator = InverseDistanceWeightedNDInterpolator(points, values,
                                                             rescale=True, fill_value=0,
                                                             local=False)
        assert interpolator._local is False, ("Initialization parameter 'local' "
                                              "not set correctly")
        assert interpolator._fill_value == 0, ("Initialization parameter 'fill_value' "
                                               "not set correctly")

    def test_local_global_flag(self):
        points, values = self._get_sample_data()
        local_interpolator = InverseDistanceWeightedNDInterpolator(points, values,
                                                                   local=True)
        global_interpolator = InverseDistanceWeightedNDInterpolator(points, values,
                                                                    local=False)
        assert local_interpolator._local is True, "Local flag not set correctly"
        assert global_interpolator._local is False, "Global flag not set correctly"

    def test_fill_value_usage(self):
        points, values = self._get_sample_data()
        fill_value = -1
        interpolator = InverseDistanceWeightedNDInterpolator(points, values,
                                                             fill_value=fill_value)
        query_points = np.array([[2, 2]])
        # Outside the convex hull of the data points
        interpolated_values = interpolator(query_points, distance_upper_bound=1e-7)
        assert_array_equal(interpolated_values, fill_value,
                           err_msg="Incorrect usage of fill_value")

    def test_input_dimensionality(self):
        points = np.random.rand(10, 3)
        values = np.random.rand(10)
        interpolator = InverseDistanceWeightedNDInterpolator(points, values)
        query_points = np.random.rand(5, 3)  # Matching dimensionality
        ip = interpolator(query_points)
        assert_equal(query_points.shape[0], len(ip))

    def test_input_dimensionality_leading_dim_stack(self):
        points = np.random.rand(10, 3)
        values = np.random.rand(10)
        interpolator = InverseDistanceWeightedNDInterpolator(points, values)
        query_points = np.random.rand(10, 10, 3)  # Matching dimensionality
        ip = interpolator(query_points)
        assert_equal(query_points.shape[0:2], ip.shape[0:2])

    def test_list_input(self):
        points, values = self._get_sample_data()
        query_points = np.array([[1, 1], [2, 2]])

        list_interpolator = InverseDistanceWeightedNDInterpolator(points,
                                                                  values.tolist())
        norm_interpolator = InverseDistanceWeightedNDInterpolator(points, values)
        assert_allclose(list_interpolator(query_points.tolist()),
                        norm_interpolator(query_points))

    def test_basic_functionality(self):
        points, values = self._get_sample_data()
        interpolator = InverseDistanceWeightedNDInterpolator(points, values)
        query_points = np.array([[1, 1], [2, 2]])
        expected_values = np.array([1, 2])
        interpolated_values = interpolator(query_points,
                                           k=2)  # k=2, it will consider [1,1] and [3,3]
        assert_array_almost_equal(interpolated_values, expected_values,
                                  err_msg="Basic functionality failed")

    def test_basic_functionality_power(self):
        points, values = self._get_sample_data()
        interpolator = InverseDistanceWeightedNDInterpolator(points, values)
        query_points = np.array([[1, 1], [2, 2]])
        expected_values = np.array([1, 2])
        interpolated_values = interpolator(query_points, k=2, p=3)
        assert_array_almost_equal(interpolated_values, expected_values,
                                  err_msg="Basic functionality failed")

    def test_basic_functionality_power_values(self):
        points, values = self._get_sample_data()
        interpolator = InverseDistanceWeightedNDInterpolator(points, values)
        query_points = np.array([[2, 2]])
        interpolated_values_p2 = interpolator(query_points, k=3, p=2)
        interpolated_values_p1 = interpolator(query_points, k=3, p=1)
        assert_array_less(interpolated_values_p2, interpolated_values_p1,
                          err_msg="Basic functionality failed")

    def test_local_global(self):
        points, values = self._get_sample_data()
        interpolator_local = InverseDistanceWeightedNDInterpolator(points, values,
                                                                   local=True)
        interpolator_global = InverseDistanceWeightedNDInterpolator(points, values,
                                                                    local=False)
        query_point = np.array([[0.5, 0.5]])
        result_all_k = interpolator_local(query_point, k=len(points))
        result_global = interpolator_global(query_point)
        assert_allclose(result_all_k, result_global,
                        err_msg="Local/global comparison failed")

    def test_large_number_of_neighbors(self):
        points = np.c_[np.arange(100), np.arange(100)]
        values = np.concatenate([np.ones(50) * 5.0, np.ones(50) * 10.0])
        interpolator = InverseDistanceWeightedNDInterpolator(points, values, local=True)
        query_points = np.array([[0.5, 0.5], [99.5, 99.5]])
        expected_values = np.array([5, 10])
        interpolated_values = interpolator(query_points, k=50)
        assert_array_almost_equal(interpolated_values, expected_values,
                                  err_msg="local functionality failed")
        interpolated_values = interpolator(query_points, k=51)
        with assert_raises(AssertionError):
            assert_array_almost_equal(interpolated_values, expected_values)

    def test_basic_functionality_complex_values(self):
        points = np.array([[1., 1], [3, 3], [5, 5], [7, 7], [9, 9]])
        values = np.array([1., 3, 5, 7, 9]) + 1e-14j
        interpolator = InverseDistanceWeightedNDInterpolator(points, values)
        query_points = np.array([[1, 1], [2, 2]])
        expected_values = np.array([1, 2])
        interpolated_values = interpolator(query_points,
                                           k=2)  # k=2, it will consider [1,1] and [3,3]
        assert_array_almost_equal(interpolated_values, expected_values,
                                  err_msg="Basic complex functionality failed")

    def test_with_different_rescale_options(self):
        points = np.array([[0, 0], [100, 100], [200, 200]])
        values = np.array([0, 1, 2])
        interpolator_rescaled = InverseDistanceWeightedNDInterpolator(points, values,
                                                                      rescale=True)
        interpolator_non_rescaled = InverseDistanceWeightedNDInterpolator(points,
                                                                          values,
                                                                          rescale=False)
        query_points = np.array([[50, 50], [150, 150]])
        interpolated_values_rescaled = interpolator_rescaled(query_points)
        interpolated_values_non_rescaled = interpolator_non_rescaled(query_points)
        assert_array_almost_equal(interpolated_values_rescaled,
                                  interpolated_values_non_rescaled,
                                  err_msg="rescale functionality failed")

    def test_all_points_at_zero_distance_local(self):
        points = np.array([[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]])
        values = np.array([1, 5, 9, 32, 80])
        interpolator = InverseDistanceWeightedNDInterpolator(points, values)
        query_points = np.array([[0, 0]])
        interpolated_values = interpolator(query_points, k=len(points))
        expected_value = values.mean()  # or some other expected behavior
        assert_array_almost_equal(interpolated_values, expected_value,
                                  err_msg="Failed when all points are at zero distance")

    def test_all_points_at_zero_distance_global(self):
        points = np.array([[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]])
        values = np.array([1, 5, 9, 32, 80])
        interpolator = InverseDistanceWeightedNDInterpolator(points, values,
                                                             local=False)
        query_points = np.array([[0, 0]])
        interpolated_values = interpolator(query_points)
        expected_value = values.mean()  # or some other expected behavior
        assert_array_almost_equal(interpolated_values, expected_value,
                                  err_msg="Failed when all points are at zero distance")

    def test_weighting_func_functionality(self):
        points, values = self._get_sample_data()

        def exp_weight_func(dists, p, decorrelation_length = 1):
            return np.exp(-(dists / decorrelation_length) ** p)

        query_points = np.array([[0.5, 0.5]])
        interpolator = InverseDistanceWeightedNDInterpolator(points, values)
        interpolated_values_normal = interpolator(query_points, weight_func=None)
        interpolated_values_exp = interpolator(query_points,
                                               weight_func=exp_weight_func)

        with assert_raises(AssertionError):
            assert_array_almost_equal(interpolated_values_normal,
                                      interpolated_values_exp)

    def test_idw_kdtree_options(self):
        points, values = self._get_sample_data()
        idw_normal = InverseDistanceWeightedNDInterpolator(points, values)

        opts = {'balanced_tree': False, 'compact_nodes': False}
        idw_options = InverseDistanceWeightedNDInterpolator(points, values,
                                                            tree_options=opts)

        query_point = np.array([[1.5, 1.5]])
        assert_allclose(idw_options(query_point), idw_normal(query_point))
