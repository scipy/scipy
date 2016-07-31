from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_allclose

from scipy import ndimage
from scipy.ndimage import _ctest
from scipy.ndimage import _ctest_oldapi
from scipy.ndimage import _cytest

MODULES = [_ctest, _ctest_oldapi, _cytest]


def test_generic_filter():
    def filter2d(footprint_elements, weights):
        return (weights*footprint_elements).sum()

    im = np.ones((20, 20))
    im[:10,:10] = 0
    footprint = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])
    footprint_size = np.count_nonzero(footprint)
    weights = np.ones(footprint_size)/footprint_size
    for mod in MODULES:
        res = ndimage.generic_filter(im, mod.filter2d(weights),
                                     footprint=footprint)
        std = ndimage.generic_filter(im, filter2d, footprint=footprint,
                                     extra_arguments=(weights,))
        assert_allclose(res, std, err_msg="{} failed".format(mod.__name__))


def test_generic_filter1d():
    def filter1d(input_line, output_line, filter_size):
        for i in range(output_line.size):
            output_line[i] = 0
            for j in range(filter_size):
                output_line[i] += input_line[i+j]
        output_line /= filter_size

    im = np.tile(np.hstack((np.zeros(10), np.ones(10))), (10, 1))
    filter_size = 3
    for mod in MODULES:
        res = ndimage.generic_filter1d(im, mod.filter1d(filter_size),
                                       filter_size)
        std = ndimage.generic_filter1d(im, filter1d, filter_size,
                                       extra_arguments=(filter_size,))
        assert_allclose(res, std, err_msg="{} failed".format(mod.__name__))


def test_geometric_transform():
    def transform(output_coordinates, shift):
        return output_coordinates[0] - shift, output_coordinates[1] - shift

    im = np.arange(12).reshape(4, 3).astype(np.float64)
    shift = 0.5
    for mod in MODULES:
        res = ndimage.geometric_transform(im, _ctest.transform(shift))
        std = ndimage.geometric_transform(im, transform, extra_arguments=(shift,))
        assert_allclose(res, std, err_msg="{} failed".format(mod.__name__))
