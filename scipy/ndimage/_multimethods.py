import functools
import numpy as np
from scipy._lib.uarray import Dispatchable, all_of_type, create_multimethod
from scipy.ndimage import _api


mark_non_coercible = lambda x: Dispatchable(x, np.ndarray, coercible=False)
create_ndimage = functools.partial(create_multimethod,
                                   domain="numpy.scipy.ndimage")


def _get_docs(func):
    """
    Decorator to take the docstring from original
    function and assign to the multimethod.
    """
    func.__doc__ = getattr(_api, func.__name__).__doc__
    return func


############################################
""" filters multimethods """
############################################


def _input_weights_axis_output_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for signatures involving required ``input``
    and `weights` arrays with ``axis`` and ``output`` kwargs.
    """
    def self_method(input, weights, axis=-1, output=None, *args, **kwargs):
        return (dispatchables[0], dispatchables[1], axis,
                dispatchables[2]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_input_weights_axis_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def correlate1d(input, weights, axis=-1, output=None,
                mode="reflect", cval=0.0, origin=0):
    return input, weights, mark_non_coercible(output)


@create_ndimage(_input_weights_axis_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def convolve1d(input, weights, axis=-1, output=None,
               mode="reflect", cval=0.0, origin=0):
    return input, weights, mark_non_coercible(output)


def _gaussian_filter1d_arg_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for gaussian_filter1d.
    """
    def self_method(input, sigma, axis=-1, order=0,
                    output=None, *args, **kwargs):
        return (dispatchables[0], sigma, axis, order,
                dispatchables[1]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_gaussian_filter1d_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def gaussian_filter1d(input, sigma, axis=-1, order=0, output=None,
                      mode="reflect", cval=0.0, truncate=4.0):
    return input, mark_non_coercible(output)


def _gaussian_filter_arg_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for gaussian_filter.
    """
    def self_method(input, sigma, order=0, output=None, *args, **kwargs):
        return (dispatchables[0], sigma,
                order, dispatchables[1]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_gaussian_filter_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def gaussian_filter(input, sigma, order=0, output=None,
                    mode="reflect", cval=0.0, truncate=4.0):
    return input, mark_non_coercible(output)


def _input_axis_output_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for signatures involving required
    ``input`` array with ``axis`` and ``output`` kwargs.
    """
    def self_method(input, axis=-1, output=None, *args, **kwargs):
        return (dispatchables[0], axis, dispatchables[1]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_input_axis_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def prewitt(input, axis=-1, output=None, mode="reflect", cval=0.0):
    return input, mark_non_coercible(output)


@create_ndimage(_input_axis_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def sobel(input, axis=-1, output=None, mode="reflect", cval=0.0):
    return input, mark_non_coercible(output)


def _input_arg2_output_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for signatures involving required ``input``
    array and a non-dispatchable second argument with an ``output`` kwarg.
    """
    def self_method(input, arg2, output=None, *args, **kwargs):
        return (dispatchables[0], arg2, dispatchables[1]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_input_arg2_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def generic_laplace(input, derivative2, output=None, mode="reflect",
                    cval=0.0, extra_arguments=(), extra_keywords=None):
    return input, mark_non_coercible(output)


def _input_output_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for signatures involving
    required ``input`` array and an ``output`` kwarg.
    """
    def self_method(input, output=None, *args, **kwargs):
        return (dispatchables[0], dispatchables[1]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_input_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def laplace(input, output=None, mode="reflect", cval=0.0):
    return input, mark_non_coercible(output)


@create_ndimage(_input_arg2_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def gaussian_laplace(input, sigma, output=None, mode="reflect",
                     cval=0.0, **kwargs):
    return input, mark_non_coercible(output)


@create_ndimage(_input_arg2_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def generic_gradient_magnitude(input, derivative, output=None,
                               mode="reflect", cval=0.0,
                               extra_arguments=(), extra_keywords=None):
    return input, mark_non_coercible(output)


@create_ndimage(_input_arg2_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def gaussian_gradient_magnitude(input, sigma, output=None,
                                mode="reflect", cval=0.0, **kwargs):
    return input, mark_non_coercible(output)


def _input_weights_output_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for signatures involving required
    ``input`` and `weights` array an ``output`` kwarg.
    """
    def self_method(input, weights, output=None, *args, **kwargs):
        return (dispatchables[0], dispatchables[1],
                dispatchables[2]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_input_weights_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def correlate(input, weights, output=None, mode='reflect', cval=0.0,
              origin=0):
    return input, weights, mark_non_coercible(output)


@create_ndimage(_input_weights_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def convolve(input, weights, output=None, mode='reflect', cval=0.0,
             origin=0):
    return input, weights, mark_non_coercible(output)


def _input_size_axis_output_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for signatures involving required ``input``
    array with ``size``, ``axis`` and ``output`` kwargs.
    """
    def self_method(input, size, axis=-1, output=None, *args, **kwargs):
        return (dispatchables[0], size, axis, dispatchables[1]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_input_size_axis_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def uniform_filter1d(input, size, axis=-1, output=None,
                     mode="reflect", cval=0.0, origin=0):
    return input, mark_non_coercible(output)


def _uniform_filter_arg_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for uniform_filter.
    """
    def self_method(input, size=3, output=None, *args, **kwargs):
        return (dispatchables[0], size, dispatchables[1]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_uniform_filter_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def uniform_filter(input, size=3, output=None, mode="reflect",
                   cval=0.0, origin=0):
    return input, mark_non_coercible(output)


@create_ndimage(_input_size_axis_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def minimum_filter1d(input, size, axis=-1, output=None,
                     mode="reflect", cval=0.0, origin=0):
    return input, mark_non_coercible(output)


@create_ndimage(_input_size_axis_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def maximum_filter1d(input, size, axis=-1, output=None,
                     mode="reflect", cval=0.0, origin=0):
    return input, mark_non_coercible(output)


def _input_size_footprint_output_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for signatures involving required ``input``
    with ``size``, ``footprint`` and an ``output`` kwarg.
    """
    def self_method(input, size=None, footprint=None,
                    output=None, *args, **kwargs):
      return (dispatchables[0], size, dispatchables[1],
              dispatchables[2]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_input_size_footprint_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def minimum_filter(input, size=None, footprint=None, output=None,
                   mode="reflect", cval=0.0, origin=0):
    return input, footprint, mark_non_coercible(output)


@create_ndimage(_input_size_footprint_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def maximum_filter(input, size=None, footprint=None, output=None,
                   mode="reflect", cval=0.0, origin=0):
    return input, footprint, mark_non_coercible(output)


@create_ndimage(_input_size_footprint_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def median_filter(input, size=None, footprint=None, output=None,
                  mode="reflect", cval=0.0, origin=0):
    return input, footprint, mark_non_coercible(output)


def _input_arg2_size_footprint_output_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for signatures involving required ``input``
    and a non-dispatchable second argument with ``size``, ``footprint``
    and ``output`` kwarg.
    """
    def self_method(input, arg2, size=None, footprint=None,
                    output=None, *args, **kwargs):
      return (dispatchables[0], arg2, size,
              dispatchables[1], dispatchables[2]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_input_arg2_size_footprint_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def rank_filter(input, rank, size=None, footprint=None, output=None,
                mode="reflect", cval=0.0, origin=0):
    return input, footprint, mark_non_coercible(output)


@create_ndimage(_input_arg2_size_footprint_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def percentile_filter(input, percentile, size=None, footprint=None,
                      output=None, mode="reflect", cval=0.0, origin=0):
    return input, footprint, mark_non_coercible(output)


def _generic_filter1d_arg_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for generic_filter1d
    """
    def self_method(input, function, filter_size, axis=-1,
                    output=None, *args, **kwargs):
      return (dispatchables[0], arg2, arg3,
              axis, dispatchables[1]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_generic_filter1d_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def generic_filter1d(input, function, filter_size, axis=-1,
                     output=None, mode="reflect", cval=0.0, origin=0,
                     extra_arguments=(), extra_keywords=None):
    return input, mark_non_coercible(output)


@create_ndimage(_input_arg2_size_footprint_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def generic_filter(input, function, size=None, footprint=None,
                   output=None, mode="reflect", cval=0.0, origin=0,
                   extra_arguments=(), extra_keywords=None):
    return input, footprint, mark_non_coercible(output)


############################################
""" fourier multimethods """
############################################

def _fourier_multimethods_arg_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for fourier multimethods.
    """
    def self_method(input, arg2, n=-1, axis=-1, output=None, *args, **kwargs):
      return (dispatchables[0], arg2, n, axis,
              dispatchables[1]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_fourier_multimethods_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def fourier_gaussian(input, sigma, n=-1, axis=-1, output=None):
    return input, mark_non_coercible(output)


@create_ndimage(_fourier_multimethods_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def fourier_uniform(input, size, n=-1, axis=-1, output=None):
    return input, mark_non_coercible(output)


@create_ndimage(_fourier_multimethods_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def fourier_ellipsoid(input, size, n=-1, axis=-1, output=None):
    return input, mark_non_coercible(output)


@create_ndimage(_fourier_multimethods_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def fourier_shift(input, shift, n=-1, axis=-1, output=None):
    return input, mark_non_coercible(output)


############################################
""" interpolation multimethods """
############################################


def _spline_filter1d_arg_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for spline_filter1d.
    """
    def self_method(input, order=3, axis=-1, output=np.float64,
                    *args, **kwargs):
      return (dispatchables[0], order, axis,
              dispatchables[1]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_spline_filter1d_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def spline_filter1d(input, order=3, axis=-1, output=np.float64,
                    mode='mirror'):
    return input, mark_non_coercible(output)


def _input_order_output_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for signatures involving required ``input``
    with ``order`` and ``output`` kwargs.
    """
    def self_method(input, order=3, output=np.float64, *args, **kwargs):
      return (dispatchables[0], order, dispatchables[1]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_input_order_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def spline_filter(input, order=3, output=np.float64, mode='mirror'):
    return input, mark_non_coercible(output)


def _geometric_transform_arg_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for geometric_transform.
    """
    def self_method(input, mapping, output_shape=None,
                    output=None, *args, **kwargs):
      return (dispatchables[0], mapping,
              output_shape, dispatchables[1]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_geometric_transform_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def geometric_transform(input, mapping, output_shape=None,
                        output=None, order=3,
                        mode='constant', cval=0.0, prefilter=True,
                        extra_arguments=(), extra_keywords={}):
    return input, mark_non_coercible(output)


def _map_coordinates_arg_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for map_coordinates.
    """
    def self_method(input, coordinates, output=None, *args, **kwargs):
      return (dispatchables[0], coordinates, dispatchables[1]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_map_coordinates_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def map_coordinates(input, coordinates, output=None, order=3,
                    mode='constant', cval=0.0, prefilter=True):
    return input, mark_non_coercible(output)


def _affine_transform_arg_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for affine_transform.
    """
    def self_method(input, matrix, offset=0.0, output_shape=None,
                    output=None, *args, **kwargs):
      return (dispatchables[0], dispatchables[1], offset,
              output_shape, dispatchables[2]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_affine_transform_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def affine_transform(input, matrix, offset=0.0, output_shape=None,
                     output=None, order=3, mode='constant',
                     cval=0.0, prefilter=True):
    return (input, matrix, mark_non_coercible(output))


@create_ndimage(_input_arg2_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def shift(input, shift, output=None, order=3, mode='constant', cval=0.0,
          prefilter=True):
    return input, mark_non_coercible(output)


@create_ndimage(_input_arg2_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def zoom(input, zoom, output=None, order=3, mode='constant', cval=0.0,
         prefilter=True, *, grid_mode=False):
    return input, mark_non_coercible(output)


def _rotate_arg_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for rotate.
    """
    def self_method(input, angle, axes=(1,0), reshape=True,
                    output=None, *args, **kwargs):
      return (dispatchables[0], angle, axes, reshape,
              dispatchables[2]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_rotate_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def rotate(input, angle, axes=(1, 0), reshape=True, output=None, order=3,
           mode='constant', cval=0.0, prefilter=True):
    return input, mark_non_coercible(output)


############################################
""" measurement multimethods """
############################################


def _input_structure_output_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for signatures involving required
    ``input`` with ``structure`` and ``output`` kwargs.
    """
    def self_method(input, structure=None, output=None, *args, **kwargs):
      return (dispatchables[0], dispatchables[1],
              dispatchables[2]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_input_structure_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def label(input, structure=None, output=None):
    return input, structure, mark_non_coercible(output)


def _input_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for signatures involving
    the required ``input`` as the only dispatchable.
    """
    def self_method(input, *args, **kwargs):
      return (dispatchables[0],) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_input_replacer)
@all_of_type(np.ndarray)
@_get_docs
def find_objects(input, max_label=0):
    return (input,)


def _input_labels_index_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for signatures involving required
    ``input`` with ``labels`` and ``index`` kwargs.
    """
    def self_method(input, labels=None, index=None, *args, **kwargs):
      return (dispatchables[0], dispatchables[1],
              dispatchables[2]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_input_labels_index_replacer)
@all_of_type(np.ndarray)
@_get_docs
def labeled_comprehension(input, labels, index, func, out_dtype,
                          default, pass_positions=False):
    return input, labels, index


@create_ndimage(_input_labels_index_replacer)
@all_of_type(np.ndarray)
@_get_docs
def sum_labels(input, labels=None, index=None):
    return input, labels, index


# alias for sum_labels; kept for backward compatibility
@_get_docs
def sum(input, labels=None, index=None):
    return sum_labels(input, labels=labels, index=index)


@create_ndimage(_input_labels_index_replacer)
@all_of_type(np.ndarray)
@_get_docs
def mean(input, labels=None, index=None):
    return input, labels, index


@create_ndimage(_input_labels_index_replacer)
@all_of_type(np.ndarray)
@_get_docs
def variance(input, labels=None, index=None):
    return input, labels, index


@create_ndimage(_input_labels_index_replacer)
@all_of_type(np.ndarray)
@_get_docs
def standard_deviation(input, labels=None, index=None):
    return input, labels, index


@create_ndimage(_input_labels_index_replacer)
@all_of_type(np.ndarray)
@_get_docs
def minimum(input, labels=None, index=None):
    return input, labels, index


@create_ndimage(_input_labels_index_replacer)
@all_of_type(np.ndarray)
@_get_docs
def maximum(input, labels=None, index=None):
    return input, labels, index


@create_ndimage(_input_labels_index_replacer)
@all_of_type(np.ndarray)
@_get_docs
def median(input, labels=None, index=None):
    return input, labels, index


@create_ndimage(_input_labels_index_replacer)
@all_of_type(np.ndarray)
@_get_docs
def minimum_position(input, labels=None, index=None):
    return input, labels, index


@create_ndimage(_input_labels_index_replacer)
@all_of_type(np.ndarray)
@_get_docs
def maximum_position(input, labels=None, index=None):
    return input, labels, index


@create_ndimage(_input_labels_index_replacer)
@all_of_type(np.ndarray)
@_get_docs
def extrema(input, labels=None, index=None):
    return input, labels, index


@create_ndimage(_input_labels_index_replacer)
@all_of_type(np.ndarray)
@_get_docs
def center_of_mass(input, labels=None, index=None):
    return input, labels, index


def _histogram_arg_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for histogram.
    """
    def self_method(input, min, max, bins, labels=None, index=None,
                    *args, **kwargs):
        return (dispatchables[0], min, max, bins,
                dispatchables[1], dispatchables[2]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_histogram_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def histogram(input, min, max, bins, labels=None, index=None):
    return input, labels, index


def _watershed_ift_arg_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for watershed_ift.
    """
    def self_method(input, markers, structure=None, output=None,
                    *args, **kwargs):
        return (dispatchables[0], dispatchables[1],
                dispatchables[2], dispatchables[3]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_watershed_ift_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def watershed_ift(input, markers, structure=None, output=None):
    return input, markers, structure, mark_non_coercible(output)


############################################
""" morphology multimethods """
############################################


@create_ndimage(_input_replacer)
@all_of_type(np.ndarray)
@_get_docs
def iterate_structure(structure, iterations, origin=None):
    return (structure,)


def _identity_arg_replacer(args, kwargs, arrays):
    """
    uarray argument replacer when dispatchables are empty.
    """
    return args, kwargs


@create_ndimage(_identity_arg_replacer)
def generate_binary_structure(rank, connectivity):
    return


def _binary_erosion_dilation_arg_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for binary_erosion and binary_dilation.
    """
    def self_method(input, structure=None, iterations=1, mask=None,
                    output=None, *args, **kwargs):
        return (dispatchables[0], dispatchables[1], iterations, mask,
                dispatchables[2], dispatchables[3]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_binary_erosion_dilation_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def binary_erosion(input, structure=None, iterations=1, mask=None,
                   output=None, border_value=0, origin=0, brute_force=False):
    return input, structure, mask, mark_non_coercible(output)


@create_ndimage(_binary_erosion_dilation_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def binary_dilation(input, structure=None, iterations=1, mask=None,
                    output=None, border_value=0, origin=0,
                    brute_force=False):

    return input, structure, mask, mark_non_coercible(output)


def _binary_opening_closing_arg_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for binary_opening and binary_closing.
    """
    def self_method(input, structure=None, iterations=1,
                    output=None, *args, **kwargs):
        return (dispatchables[0], dispatchables[1],
                iterations, dispatchables[2]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_binary_opening_closing_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def binary_opening(input, structure=None, iterations=1, output=None,
                   origin=0, mask=None, border_value=0, brute_force=False):
    return input, structure, mask, mark_non_coercible(output)


@create_ndimage(_binary_opening_closing_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def binary_closing(input, structure=None, iterations=1, output=None,
                   origin=0, mask=None, border_value=0, brute_force=False):
    return input, structure, mask, mark_non_coercible(output)


def _binary_hit_or_miss_arg_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for binary_hit_or_miss.
    """
    def self_method(input, structure1=None, structure2=None,
                    output=None, *args, **kwargs):
        return (dispatchables[0], dispatchables[1],
                dispatchables[2], dispatchables[3]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_binary_hit_or_miss_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def binary_hit_or_miss(input, structure1=None, structure2=None,
                       output=None, origin1=0, origin2=None):
    return input, structure1, structure2, mark_non_coercible(output)


def _binary_propagation_arg_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer for binary_propagation.
    """
    def self_method(input, structure=None, mask=None,
                    output=None, *args, **kwargs):
        return (dispatchables[0], dispatchables[1],
                dispatchables[2], dispatchables[3]) + args, kwargs
    return self_method(*args, **kwargs)


@create_ndimage(_binary_propagation_arg_replacer)
@all_of_type(np.ndarray)
@_get_docs
def binary_propagation(input, structure=None, mask=None,
                       output=None, border_value=0, origin=0):
    return input, structure, mask, mark_non_coercible(output)


@create_ndimage(_input_structure_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def binary_fill_holes(input, structure=None, output=None, origin=0):
    return input, structure, mark_non_coercible(output)


def _input_size_footprint_structure_output_replacer(args, kwargs,
                                                    dispatchables):
    """
    uarray argument replacer for signatures involving required ``input``
    with ``size``, ``footprint``, ``structure`` and ``output`` kwargs.
    """
    def self_method(input, size=None, footprint=None, structure=None,
                    output=None, *args, **kwargs):
        return (dispatchables[0], size, dispatchables[1],
                dispatchables[2], dispatchables[3]) + args, kwargs

    return self_method(*args, **kwargs)


@create_ndimage(_input_size_footprint_structure_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def grey_erosion(input, size=None, footprint=None, structure=None,
                 output=None, mode="reflect", cval=0.0, origin=0):
    return input, footprint, structure, mark_non_coercible(output)


@create_ndimage(_input_size_footprint_structure_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def grey_dilation(input, size=None, footprint=None, structure=None,
                  output=None, mode="reflect", cval=0.0, origin=0):
    return input, footprint, structure, mark_non_coercible(output)


@create_ndimage(_input_size_footprint_structure_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def grey_opening(input, size=None, footprint=None, structure=None,
                 output=None, mode="reflect", cval=0.0, origin=0):
    return input, footprint, structure, mark_non_coercible(output)


@create_ndimage(_input_size_footprint_structure_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def grey_closing(input, size=None, footprint=None, structure=None,
                 output=None, mode="reflect", cval=0.0, origin=0):
    return input, footprint, structure, mark_non_coercible(output)


@create_ndimage(_input_size_footprint_structure_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def morphological_gradient(input, size=None, footprint=None, structure=None,
                           output=None, mode="reflect", cval=0.0, origin=0):
    return input, footprint, structure, mark_non_coercible(output)


@create_ndimage(_input_size_footprint_structure_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def morphological_laplace(input, size=None, footprint=None, structure=None,
                          output=None, mode="reflect", cval=0.0, origin=0):
    return input, footprint, structure, mark_non_coercible(output)


@create_ndimage(_input_size_footprint_structure_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def white_tophat(input, size=None, footprint=None, structure=None,
                 output=None, mode="reflect", cval=0.0, origin=0):
    return input, footprint, structure, mark_non_coercible(output)


@create_ndimage(_input_size_footprint_structure_output_replacer)
@all_of_type(np.ndarray)
@_get_docs
def black_tophat(input, size=None, footprint=None, structure=None,
                 output=None, mode="reflect", cval=0.0, origin=0):
    return input, footprint, structure, mark_non_coercible(output)


@create_ndimage(_input_replacer)
@all_of_type(np.ndarray)
@_get_docs
def distance_transform_bf(input, metric="euclidean", sampling=None,
                          return_distances=True, return_indices=False,
                          distances=None, indices=None):
    return (input,)


@create_ndimage(_input_replacer)
@all_of_type(np.ndarray)
@_get_docs
def distance_transform_cdt(input, metric='chessboard', return_distances=True,
                           return_indices=False, distances=None,
                           indices=None):
    return (input,)


@create_ndimage(_input_replacer)
@all_of_type(np.ndarray)
@_get_docs
def distance_transform_edt(input, sampling=None, return_distances=True,
                           return_indices=False, distances=None,
                           indices=None):
    return (input,)
