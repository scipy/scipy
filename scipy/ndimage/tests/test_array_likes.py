"""
Smoke test "array_like" inputs.
"""
import numpy as np
from scipy import ndimage
from scipy._lib._array_api import xp_assert_close


def _assert_same_result(result, expected):
    if isinstance(result, tuple):
        assert isinstance(expected, tuple)
        assert len(result) == len(expected)
        for actual_item, expected_item in zip(result, expected):
            _assert_same_result(actual_item, expected_item)
        return

    if isinstance(result, list):
        assert isinstance(expected, list)
        assert len(result) == len(expected)
        for actual_item, expected_item in zip(result, expected):
            _assert_same_result(actual_item, expected_item)
        return

    if isinstance(result, dict):
        assert result.keys() == expected.keys()
        for key in result:
            _assert_same_result(result[key], expected[key])
        return

    if isinstance(result, slice):
        assert result == expected
        return

    if result is None or expected is None:
        assert result is expected
        return

    xp_assert_close(result, expected, atol=1e-14)


def _generic_filter1d(line, output):
    output[...] = line


def _derivative(input, axis, output, mode, cval):
    return ndimage.sobel(input, axis=axis, output=output, mode=mode, cval=cval)


def _derivative2(input, axis, output, mode, cval):
    return ndimage.correlate1d(
        input,
        [1, -2, 1],
        axis=axis,
        output=output,
        mode=mode,
        cval=cval,
    )


def _vectorized_mean(window, *, axis):
    return np.mean(window, axis=axis)


def test_correlate1d_accepts_lists():
    input_list = [1, 2, 3, 4]
    input_array = np.asarray(input_list)
    weights_list = [1, 0, -1]
    weights_array = np.asarray(weights_list)
    result = ndimage.correlate1d(input_list, weights_list)
    expected = ndimage.correlate1d(input_array, weights_array)
    _assert_same_result(result, expected)


def test_convolve1d_accepts_lists():
    input_list = [1, 2, 3, 4]
    input_array = np.asarray(input_list)
    weights_list = [1, 0, -1]
    weights_array = np.asarray(weights_list)
    result = ndimage.convolve1d(input_list, weights_list)
    expected = ndimage.convolve1d(input_array, weights_array)
    _assert_same_result(result, expected)


def test_gaussian_filter1d_accepts_lists():
    input_list = [1, 2, 3, 4]
    input_array = np.asarray(input_list)
    sigma = 1
    result = ndimage.gaussian_filter1d(input_list, sigma)
    expected = ndimage.gaussian_filter1d(input_array, sigma)
    _assert_same_result(result, expected)


def test_gaussian_filter_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    sigma = 1
    result = ndimage.gaussian_filter(input_list, sigma)
    expected = ndimage.gaussian_filter(input_array, sigma)
    _assert_same_result(result, expected)


def test_prewitt_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    result = ndimage.prewitt(input_list)
    expected = ndimage.prewitt(input_array)
    _assert_same_result(result, expected)


def test_sobel_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    result = ndimage.sobel(input_list)
    expected = ndimage.sobel(input_array)
    _assert_same_result(result, expected)


def test_generic_laplace_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    derivative2 = _derivative2
    result = ndimage.generic_laplace(input_list, derivative2)
    expected = ndimage.generic_laplace(input_array, derivative2)
    _assert_same_result(result, expected)


def test_laplace_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    result = ndimage.laplace(input_list)
    expected = ndimage.laplace(input_array)
    _assert_same_result(result, expected)


def test_gaussian_laplace_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    sigma = 1
    result = ndimage.gaussian_laplace(input_list, sigma)
    expected = ndimage.gaussian_laplace(input_array, sigma)
    _assert_same_result(result, expected)


def test_generic_gradient_magnitude_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    derivative = _derivative
    result = ndimage.generic_gradient_magnitude(input_list, derivative)
    expected = ndimage.generic_gradient_magnitude(input_array, derivative)
    _assert_same_result(result, expected)


def test_gaussian_gradient_magnitude_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    sigma = 1
    result = ndimage.gaussian_gradient_magnitude(input_list, sigma)
    expected = ndimage.gaussian_gradient_magnitude(input_array, sigma)
    _assert_same_result(result, expected)


def test_correlate_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    weights_list = [[1, 0, -1], [1, 0, -1], [1, 0, -1]]
    weights_array = np.asarray(weights_list)
    result = ndimage.correlate(input_list, weights_list)
    expected = ndimage.correlate(input_array, weights_array)
    _assert_same_result(result, expected)


def test_convolve_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    weights_list = [[1, 0, -1], [1, 0, -1], [1, 0, -1]]
    weights_array = np.asarray(weights_list)
    result = ndimage.convolve(input_list, weights_list)
    expected = ndimage.convolve(input_array, weights_array)
    _assert_same_result(result, expected)


def test_uniform_filter1d_accepts_lists():
    input_list = [1, 2, 3, 4]
    input_array = np.asarray(input_list)
    size = 2
    result = ndimage.uniform_filter1d(input_list, size)
    expected = ndimage.uniform_filter1d(input_array, size)
    _assert_same_result(result, expected)


def test_uniform_filter_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    size = 2
    result = ndimage.uniform_filter(input_list, size)
    expected = ndimage.uniform_filter(input_array, size)
    _assert_same_result(result, expected)


def test_minimum_filter1d_accepts_lists():
    input_list = [1, 2, 3, 4]
    input_array = np.asarray(input_list)
    size = 2
    result = ndimage.minimum_filter1d(input_list, size)
    expected = ndimage.minimum_filter1d(input_array, size)
    _assert_same_result(result, expected)


def test_maximum_filter1d_accepts_lists():
    input_list = [1, 2, 3, 4]
    input_array = np.asarray(input_list)
    size = 2
    result = ndimage.maximum_filter1d(input_list, size)
    expected = ndimage.maximum_filter1d(input_array, size)
    _assert_same_result(result, expected)


def test_minimum_filter_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    size = 2
    result = ndimage.minimum_filter(input_list, size)
    expected = ndimage.minimum_filter(input_array, size)
    _assert_same_result(result, expected)


def test_maximum_filter_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    size = 2
    result = ndimage.maximum_filter(input_list, size)
    expected = ndimage.maximum_filter(input_array, size)
    _assert_same_result(result, expected)


def test_rank_filter_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    rank = 1
    size = 2
    result = ndimage.rank_filter(input_list, rank, size)
    expected = ndimage.rank_filter(input_array, rank, size)
    _assert_same_result(result, expected)


def test_median_filter_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    size = 2
    result = ndimage.median_filter(input_list, size)
    expected = ndimage.median_filter(input_array, size)
    _assert_same_result(result, expected)


def test_percentile_filter_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    percentile = 50
    size = 2
    result = ndimage.percentile_filter(input_list, percentile, size)
    expected = ndimage.percentile_filter(input_array, percentile, size)
    _assert_same_result(result, expected)


def test_generic_filter1d_accepts_lists():
    input_list = [1, 2, 3, 4]
    input_array = np.asarray(input_list)
    function = _generic_filter1d
    filter_size = 1
    result = ndimage.generic_filter1d(input_list, function, filter_size)
    expected = ndimage.generic_filter1d(input_array, function, filter_size)
    _assert_same_result(result, expected)


def test_generic_filter_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    function = np.mean
    size = 2
    result = ndimage.generic_filter(input_list, function, size)
    expected = ndimage.generic_filter(input_array, function, size)
    _assert_same_result(result, expected)


def test_vectorized_filter_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    function = _vectorized_mean
    size = 2
    result = ndimage.vectorized_filter(input_list, function, size=size)
    expected = ndimage.vectorized_filter(input_array, function, size=size)
    _assert_same_result(result, expected)


def test_fourier_gaussian_accepts_lists():
    input_list = [[1.0, 2.0], [3.0, 4.0]]
    input_array = np.asarray(input_list)
    sigma = 1
    result = ndimage.fourier_gaussian(input_list, sigma)
    expected = ndimage.fourier_gaussian(input_array, sigma)
    _assert_same_result(result, expected)


def test_fourier_uniform_accepts_lists():
    input_list = [[1.0, 2.0], [3.0, 4.0]]
    input_array = np.asarray(input_list)
    size = 1
    result = ndimage.fourier_uniform(input_list, size)
    expected = ndimage.fourier_uniform(input_array, size)
    _assert_same_result(result, expected)


def test_fourier_ellipsoid_accepts_lists():
    input_list = [[1.0, 2.0], [3.0, 4.0]]
    input_array = np.asarray(input_list)
    size = 1
    result = ndimage.fourier_ellipsoid(input_list, size)
    expected = ndimage.fourier_ellipsoid(input_array, size)
    _assert_same_result(result, expected)


def test_fourier_shift_accepts_lists():
    input_list = [[1.0, 2.0], [3.0, 4.0]]
    input_array = np.asarray(input_list)
    shift = 1
    result = ndimage.fourier_shift(input_list, shift)
    expected = ndimage.fourier_shift(input_array, shift)
    _assert_same_result(result, expected)


def test_spline_filter1d_accepts_lists():
    input_list = [1, 2, 3, 4]
    input_array = np.asarray(input_list)
    result = ndimage.spline_filter1d(input_list)
    expected = ndimage.spline_filter1d(input_array)
    _assert_same_result(result, expected)


def test_spline_filter_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    result = ndimage.spline_filter(input_list)
    expected = ndimage.spline_filter(input_array)
    _assert_same_result(result, expected)


def test_geometric_transform_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    result = ndimage.geometric_transform(input_list, lambda x: x)
    expected = ndimage.geometric_transform(input_array, lambda x: x)
    _assert_same_result(result, expected)


def test_map_coordinates_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    coordinates_list = [[0.0, 1.0], [0.0, 1.0]]
    coordinates_array = np.asarray(coordinates_list)
    result = ndimage.map_coordinates(input_list, coordinates_list)
    expected = ndimage.map_coordinates(input_array, coordinates_array)
    _assert_same_result(result, expected)


def test_affine_transform_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    matrix_list = [[1, 0], [0, 1]]
    matrix_array = np.asarray(matrix_list)
    result = ndimage.affine_transform(input_list, matrix_list)
    expected = ndimage.affine_transform(input_array, matrix_array)
    _assert_same_result(result, expected)


def test_shift_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    shift = 1
    result = ndimage.shift(input_list, shift)
    expected = ndimage.shift(input_array, shift)
    _assert_same_result(result, expected)


def test_zoom_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    zoom = 2
    result = ndimage.zoom(input_list, zoom)
    expected = ndimage.zoom(input_array, zoom)
    _assert_same_result(result, expected)


def test_rotate_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    angle = 0
    result = ndimage.rotate(input_list, angle)
    expected = ndimage.rotate(input_array, angle)
    _assert_same_result(result, expected)


def test_label_accepts_lists():
    input_list = [[0, 1, 0], [1, 1, 0], [0, 0, 1]]
    input_array = np.asarray(input_list)
    result = ndimage.label(input_list)
    expected = ndimage.label(input_array)
    _assert_same_result(result, expected)


def test_find_objects_accepts_lists():
    input_list = [[1, 1, 0], [0, 2, 2], [0, 0, 0]]
    input_array = np.asarray(input_list)
    result = ndimage.find_objects(input_list)
    expected = ndimage.find_objects(input_array)
    _assert_same_result(result, expected)


def test_labeled_comprehension_accepts_lists():
    input_list = [1, 2, 3, 4]
    input_array = np.asarray(input_list)
    labels_list = [1, 1, 2, 2]
    labels_array = np.asarray(labels_list)
    index_list = [1, 2]
    index_array = np.asarray(index_list)
    func = np.sum
    out_dtype = float
    default = 0
    result = ndimage.labeled_comprehension(
        input_list,
        labels_list,
        index_list,
        func,
        out_dtype,
        default,
    )
    expected = ndimage.labeled_comprehension(
        input_array,
        labels_array,
        index_array,
        func,
        out_dtype,
        default,
    )
    _assert_same_result(result, expected)


def test_sum_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6]]
    input_array = np.asarray(input_list)
    result = ndimage.sum(input_list)
    expected = ndimage.sum(input_array)
    _assert_same_result(result, expected)


def test_mean_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6]]
    input_array = np.asarray(input_list)
    result = ndimage.mean(input_list)
    expected = ndimage.mean(input_array)
    _assert_same_result(result, expected)


def test_variance_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6]]
    input_array = np.asarray(input_list)
    result = ndimage.variance(input_list)
    expected = ndimage.variance(input_array)
    _assert_same_result(result, expected)


def test_standard_deviation_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6]]
    input_array = np.asarray(input_list)
    result = ndimage.standard_deviation(input_list)
    expected = ndimage.standard_deviation(input_array)
    _assert_same_result(result, expected)


def test_minimum_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6]]
    input_array = np.asarray(input_list)
    result = ndimage.minimum(input_list)
    expected = ndimage.minimum(input_array)
    _assert_same_result(result, expected)


def test_maximum_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6]]
    input_array = np.asarray(input_list)
    result = ndimage.maximum(input_list)
    expected = ndimage.maximum(input_array)
    _assert_same_result(result, expected)


def test_median_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6]]
    input_array = np.asarray(input_list)
    result = ndimage.median(input_list)
    expected = ndimage.median(input_array)
    _assert_same_result(result, expected)


def test_minimum_position_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6]]
    input_array = np.asarray(input_list)
    result = ndimage.minimum_position(input_list)
    expected = ndimage.minimum_position(input_array)
    _assert_same_result(result, expected)


def test_maximum_position_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6]]
    input_array = np.asarray(input_list)
    result = ndimage.maximum_position(input_list)
    expected = ndimage.maximum_position(input_array)
    _assert_same_result(result, expected)


def test_extrema_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6]]
    input_array = np.asarray(input_list)
    result = ndimage.extrema(input_list)
    expected = ndimage.extrema(input_array)
    _assert_same_result(result, expected)


def test_center_of_mass_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6]]
    input_array = np.asarray(input_list)
    result = ndimage.center_of_mass(input_list)
    expected = ndimage.center_of_mass(input_array)
    _assert_same_result(result, expected)


def test_histogram_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6]]
    input_array = np.asarray(input_list)
    min = 1
    max = 6
    bins = 3
    result = ndimage.histogram(input_list, min, max, bins)
    expected = ndimage.histogram(input_array, min, max, bins)
    _assert_same_result(result, expected)


def test_watershed_ift_accepts_list_markers():
    input_list = [
        [5, 4, 5],
        [4, 1, 4],
        [5, 4, 5],
    ]
    input_array = np.asarray(input_list, dtype=np.uint8)
    markers_list = [[1, 0, 0], [0, 0, 0], [0, 0, 2]]
    markers_array = np.asarray(markers_list)
    result = ndimage.watershed_ift(input_array, markers_list)
    expected = ndimage.watershed_ift(input_array, markers_array)
    _assert_same_result(result, expected)


def test_sum_labels_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6]]
    input_array = np.asarray(input_list)
    result = ndimage.sum_labels(input_list)
    expected = ndimage.sum_labels(input_array)
    _assert_same_result(result, expected)


def test_value_indices_accepts_lists():
    arr_list = [[1, 0, 1], [2, 2, 0]]
    arr_array = np.asarray(arr_list)
    result = ndimage.value_indices(arr_list)
    expected = ndimage.value_indices(arr_array)
    _assert_same_result(result, expected)


def test_iterate_structure_accepts_lists():
    structure_list = [[0, 1, 0], [1, 1, 1], [0, 1, 0]]
    structure_array = np.asarray(structure_list)
    iterations = 2
    result = ndimage.iterate_structure(structure_list, iterations)
    expected = ndimage.iterate_structure(structure_array, iterations)
    _assert_same_result(result, expected)


def test_generate_binary_structure_accepts_lists():
    rank = 2
    connectivity = 1
    result = ndimage.generate_binary_structure(rank, connectivity)
    expected = ndimage.generate_binary_structure(rank, connectivity)
    _assert_same_result(result, expected)


def test_binary_erosion_accepts_lists():
    input_list = [[0, 1, 0], [1, 1, 0], [0, 0, 1]]
    input_array = np.asarray(input_list)
    result = ndimage.binary_erosion(input_list)
    expected = ndimage.binary_erosion(input_array)
    _assert_same_result(result, expected)


def test_binary_dilation_accepts_lists():
    input_list = [[0, 1, 0], [1, 1, 0], [0, 0, 1]]
    input_array = np.asarray(input_list)
    result = ndimage.binary_dilation(input_list)
    expected = ndimage.binary_dilation(input_array)
    _assert_same_result(result, expected)


def test_binary_opening_accepts_lists():
    input_list = [[0, 1, 0], [1, 1, 0], [0, 0, 1]]
    input_array = np.asarray(input_list)
    result = ndimage.binary_opening(input_list)
    expected = ndimage.binary_opening(input_array)
    _assert_same_result(result, expected)


def test_binary_closing_accepts_lists():
    input_list = [[0, 1, 0], [1, 0, 1], [0, 1, 0]]
    input_array = np.asarray(input_list)
    result = ndimage.binary_closing(input_list)
    expected = ndimage.binary_closing(input_array)
    _assert_same_result(result, expected)


def test_binary_hit_or_miss_accepts_lists():
    input_list = [[0, 1, 0], [1, 1, 0], [0, 0, 1]]
    input_array = np.asarray(input_list)
    result = ndimage.binary_hit_or_miss(input_list)
    expected = ndimage.binary_hit_or_miss(input_array)
    _assert_same_result(result, expected)


def test_binary_propagation_accepts_lists():
    input_list = [[0, 1, 0], [1, 1, 0], [0, 0, 1]]
    input_array = np.asarray(input_list)
    result = ndimage.binary_propagation(input_list)
    expected = ndimage.binary_propagation(input_array)
    _assert_same_result(result, expected)


def test_binary_fill_holes_accepts_lists():
    input_list = [[1, 1, 1], [1, 0, 1], [1, 1, 1]]
    input_array = np.asarray(input_list)
    result = ndimage.binary_fill_holes(input_list)
    expected = ndimage.binary_fill_holes(input_array)
    _assert_same_result(result, expected)


def test_grey_erosion_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    size = 2
    result = ndimage.grey_erosion(input_list, size=size)
    expected = ndimage.grey_erosion(input_array, size=size)
    _assert_same_result(result, expected)


def test_grey_dilation_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    size = 2
    result = ndimage.grey_dilation(input_list, size=size)
    expected = ndimage.grey_dilation(input_array, size=size)
    _assert_same_result(result, expected)


def test_grey_opening_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    size = 2
    result = ndimage.grey_opening(input_list, size=size)
    expected = ndimage.grey_opening(input_array, size=size)
    _assert_same_result(result, expected)


def test_grey_closing_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    size = 2
    result = ndimage.grey_closing(input_list, size=size)
    expected = ndimage.grey_closing(input_array, size=size)
    _assert_same_result(result, expected)


def test_morphological_gradient_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    size = 2
    result = ndimage.morphological_gradient(input_list, size=size)
    expected = ndimage.morphological_gradient(input_array, size=size)
    _assert_same_result(result, expected)


def test_morphological_laplace_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    size = 2
    result = ndimage.morphological_laplace(input_list, size=size)
    expected = ndimage.morphological_laplace(input_array, size=size)
    _assert_same_result(result, expected)


def test_white_tophat_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    size = 2
    result = ndimage.white_tophat(input_list, size=size)
    expected = ndimage.white_tophat(input_array, size=size)
    _assert_same_result(result, expected)


def test_black_tophat_accepts_lists():
    input_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    input_array = np.asarray(input_list)
    size = 2
    result = ndimage.black_tophat(input_list, size=size)
    expected = ndimage.black_tophat(input_array, size=size)
    _assert_same_result(result, expected)


def test_distance_transform_bf_accepts_lists():
    input_list = [[0, 1, 0], [1, 1, 0], [0, 0, 1]]
    input_array = np.asarray(input_list)
    result = ndimage.distance_transform_bf(input_list)
    expected = ndimage.distance_transform_bf(input_array)
    _assert_same_result(result, expected)


def test_distance_transform_cdt_accepts_lists():
    input_list = [[0, 1, 0], [1, 1, 0], [0, 0, 1]]
    input_array = np.asarray(input_list)
    result = ndimage.distance_transform_cdt(input_list)
    expected = ndimage.distance_transform_cdt(input_array)
    _assert_same_result(result, expected)


def test_distance_transform_edt_accepts_lists():
    input_list = [[0, 1, 0], [1, 1, 0], [0, 0, 1]]
    input_array = np.asarray(input_list)
    result = ndimage.distance_transform_edt(input_list)
    expected = ndimage.distance_transform_edt(input_array)
    _assert_same_result(result, expected)
