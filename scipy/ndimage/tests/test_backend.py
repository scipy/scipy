import numpy as np
import scipy.ndimage
from scipy.ndimage import set_backend
from scipy.ndimage.tests import mock_backend

from numpy.testing import assert_equal
import pytest


fnames = ('affine_transform', 'binary_closing', 'binary_dilation',
          'binary_erosion', 'binary_fill_holes', 'binary_hit_or_miss',
          'binary_opening', 'binary_propagation', 'black_tophat',
          'center_of_mass', 'convolve', 'convolve1d', 'correlate',
          'correlate1d', 'distance_transform_bf', 'distance_transform_cdt',
          'distance_transform_edt', 'extrema', 'find_objects',
          'fourier_ellipsoid', 'fourier_gaussian', 'fourier_shift',
          'fourier_uniform', 'gaussian_filter', 'gaussian_filter1d',
          'gaussian_gradient_magnitude', 'gaussian_laplace',
          'generate_binary_structure', 'generic_filter', 'generic_filter1d',
          'generic_gradient_magnitude', 'generic_laplace',
          'geometric_transform', 'grey_closing', 'grey_dilation',
          'grey_erosion', 'grey_opening', 'histogram', 'iterate_structure',
          'label', 'labeled_comprehension', 'laplace', 'map_coordinates',
          'maximum', 'maximum_filter', 'maximum_filter1d', 'maximum_position',
          'mean', 'median', 'median_filter', 'minimum', 'minimum_filter',
          'minimum_filter1d', 'minimum_position', 'morphological_gradient',
          'morphological_laplace', 'percentile_filter', 'prewitt',
          'rank_filter', 'rotate', 'shift', 'sobel', 'spline_filter',
          'spline_filter1d', 'standard_deviation', 'sum_labels',
          'uniform_filter', 'uniform_filter1d', 'variance', 'watershed_ift',
          'white_tophat', 'zoom')


funcs = [getattr(scipy.ndimage, fname) for fname in fnames]
mocks = [getattr(mock_backend, fname) for fname in fnames]


@pytest.mark.parametrize("func, mock", zip(funcs, mocks))
def test_backend_call(func, mock):
    """
    Checks fake backend dispatch.
    """
    x = np.array([[0, 0], [1, 1]])

    with set_backend(mock_backend, only=True):
        mock.number_calls = 0
        y = func(x)
        assert_equal(y, mock.return_value)
        assert_equal(mock.number_calls, 1)
