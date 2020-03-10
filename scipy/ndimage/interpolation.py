# Copyright (C) 2003-2005 Peter J. Verveer
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above
#    copyright notice, this list of conditions and the following
#    disclaimer in the documentation and/or other materials provided
#    with the distribution.
#
# 3. The name of the author may not be used to endorse or promote
#    products derived from this software without specific prior
#    written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
# OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import itertools
import numpy
import warnings

from . import _ni_support
from . import _nd_image
from ._ni_docstrings import docdict
from scipy._lib import doccer

# Change the default 'reflect' to 'constant' via modifying a copy of docdict
docdict_copy = docdict.copy()
del docdict
docdict_copy['mode'] = docdict_copy['mode'].replace("Default is 'reflect'",
                                                    "Default is 'constant'")

docfiller = doccer.filldoc(docdict_copy)

__all__ = ['spline_filter1d', 'spline_filter', 'geometric_transform',
           'map_coordinates', 'affine_transform', 'shift', 'zoom', 'rotate']


@docfiller
def spline_filter1d(input, order=3, axis=-1, output=numpy.float64,
                    mode='mirror'):
    """
    Calculate a 1-D spline filter along the given axis.

    The lines of the array along the given axis are filtered by a
    spline filter. The order of the spline must be >= 2 and <= 5.

    Parameters
    ----------
    %(input)s
    order : int, optional
        The order of the spline, default is 3.
    axis : int, optional
        The axis along which the spline filter is applied. Default is the last
        axis.
    output : ndarray or dtype, optional
        The array in which to place the output, or the dtype of the returned
        array. Default is ``numpy.float64``.
    %(mode)s

    Returns
    -------
    spline_filter1d : ndarray
        The filtered input.

    Notes
    -----
    All functions in `ndimage.interpolation` do spline interpolation of
    the input image. If using B-splines of `order > 1`, the input image
    values have to be converted to B-spline coefficients first, which is
    done by applying this 1-D filter sequentially along all
    axes of the input. All functions that require B-spline coefficients
    will automatically filter their inputs, a behavior controllable with
    the `prefilter` keyword argument. For functions that accept a `mode`
    parameter, the result will only be correct if it matches the `mode`
    used when filtering.
    """
    if order < 0 or order > 5:
        raise RuntimeError('spline order not supported')
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    output = _ni_support._get_output(output, input)
    if order in [0, 1]:
        output[...] = numpy.array(input)
    else:
        mode = _ni_support._extend_mode_to_code(mode)
        axis = _ni_support._check_axis(axis, input.ndim)
        _nd_image.spline_filter1d(input, order, axis, output, mode)
    return output


def spline_filter(input, order=3, output=numpy.float64, mode='mirror'):
    """
    Multidimensional spline filter.

    For more details, see `spline_filter1d`.

    See Also
    --------
    spline_filter1d

    Notes
    -----
    The multidimensional filter is implemented as a sequence of
    1-D spline filters. The intermediate arrays are stored
    in the same data type as the output. Therefore, for output types
    with a limited precision, the results may be imprecise because
    intermediate results may be stored with insufficient precision.

    """
    if order < 2 or order > 5:
        raise RuntimeError('spline order not supported')
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    output = _ni_support._get_output(output, input)
    if order not in [0, 1] and input.ndim > 0:
        for axis in range(input.ndim):
            spline_filter1d(input, order, axis, output=output, mode=mode)
            input = output
    else:
        output[...] = input[...]
    return output


@docfiller
def geometric_transform(input, mapping, output_shape=None,
                        output=None, order=3,
                        mode='constant', cval=0.0, prefilter=True,
                        extra_arguments=(), extra_keywords={}):
    """
    Apply an arbitrary geometric transform.

    The given mapping function is used to find, for each point in the
    output, the corresponding coordinates in the input. The value of the
    input at those coordinates is determined by spline interpolation of
    the requested order.

    Parameters
    ----------
    %(input)s
    mapping : {callable, scipy.LowLevelCallable}
        A callable object that accepts a tuple of length equal to the output
        array rank, and returns the corresponding input coordinates as a tuple
        of length equal to the input array rank.
    output_shape : tuple of ints, optional
        Shape tuple.
    %(output)s
    order : int, optional
        The order of the spline interpolation, default is 3.
        The order has to be in the range 0-5.
    %(mode)s
    %(cval)s
    %(prefilter)s
    extra_arguments : tuple, optional
        Extra arguments passed to `mapping`.
    extra_keywords : dict, optional
        Extra keywords passed to `mapping`.

    Returns
    -------
    output : ndarray
        The filtered input.

    See Also
    --------
    map_coordinates, affine_transform, spline_filter1d


    Notes
    -----
    This function also accepts low-level callback functions with one
    the following signatures and wrapped in `scipy.LowLevelCallable`:

    .. code:: c

       int mapping(npy_intp *output_coordinates, double *input_coordinates,
                   int output_rank, int input_rank, void *user_data)
       int mapping(intptr_t *output_coordinates, double *input_coordinates,
                   int output_rank, int input_rank, void *user_data)

    The calling function iterates over the elements of the output array,
    calling the callback function at each element. The coordinates of the
    current output element are passed through ``output_coordinates``. The
    callback function must return the coordinates at which the input must
    be interpolated in ``input_coordinates``. The rank of the input and
    output arrays are given by ``input_rank`` and ``output_rank``
    respectively. ``user_data`` is the data pointer provided
    to `scipy.LowLevelCallable` as-is.

    The callback function must return an integer error status that is zero
    if something went wrong and one otherwise. If an error occurs, you should
    normally set the Python error status with an informative message
    before returning, otherwise a default error message is set by the
    calling function.

    In addition, some other low-level function pointer specifications
    are accepted, but these are for backward compatibility only and should
    not be used in new code.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.ndimage import geometric_transform
    >>> a = np.arange(12.).reshape((4, 3))
    >>> def shift_func(output_coords):
    ...     return (output_coords[0] - 0.5, output_coords[1] - 0.5)
    ...
    >>> geometric_transform(a, shift_func)
    array([[ 0.   ,  0.   ,  0.   ],
           [ 0.   ,  1.362,  2.738],
           [ 0.   ,  4.812,  6.187],
           [ 0.   ,  8.263,  9.637]])

    >>> b = [1, 2, 3, 4, 5]
    >>> def shift_func(output_coords):
    ...     return (output_coords[0] - 3,)
    ...
    >>> geometric_transform(b, shift_func, mode='constant')
    array([0, 0, 0, 1, 2])
    >>> geometric_transform(b, shift_func, mode='nearest')
    array([1, 1, 1, 1, 2])
    >>> geometric_transform(b, shift_func, mode='reflect')
    array([3, 2, 1, 1, 2])
    >>> geometric_transform(b, shift_func, mode='wrap')
    array([2, 3, 4, 1, 2])

    """
    if order < 0 or order > 5:
        raise RuntimeError('spline order not supported')
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    if output_shape is None:
        output_shape = input.shape
    if input.ndim < 1 or len(output_shape) < 1:
        raise RuntimeError('input and output rank must be > 0')
    mode = _ni_support._extend_mode_to_code(mode)
    if prefilter and order > 1:
        filtered = spline_filter(input, order, output=numpy.float64)
    else:
        filtered = input
    output = _ni_support._get_output(output, input, shape=output_shape)
    _nd_image.geometric_transform(filtered, mapping, None, None, None, output,
                                  order, mode, cval, extra_arguments,
                                  extra_keywords)
    return output


@docfiller
def map_coordinates(input, coordinates, output=None, order=3,
                    mode='constant', cval=0.0, prefilter=True):
    """
    Map the input array to new coordinates by interpolation.

    The array of coordinates is used to find, for each point in the output,
    the corresponding coordinates in the input. The value of the input at
    those coordinates is determined by spline interpolation of the
    requested order.

    The shape of the output is derived from that of the coordinate
    array by dropping the first axis. The values of the array along
    the first axis are the coordinates in the input array at which the
    output value is found.

    Parameters
    ----------
    %(input)s
    coordinates : array_like
        The coordinates at which `input` is evaluated.
    %(output)s
    order : int, optional
        The order of the spline interpolation, default is 3.
        The order has to be in the range 0-5.
    %(mode)s
    %(cval)s
    %(prefilter)s

    Returns
    -------
    map_coordinates : ndarray
        The result of transforming the input. The shape of the output is
        derived from that of `coordinates` by dropping the first axis.

    See Also
    --------
    spline_filter, geometric_transform, scipy.interpolate

    Examples
    --------
    >>> from scipy import ndimage
    >>> a = np.arange(12.).reshape((4, 3))
    >>> a
    array([[  0.,   1.,   2.],
           [  3.,   4.,   5.],
           [  6.,   7.,   8.],
           [  9.,  10.,  11.]])
    >>> ndimage.map_coordinates(a, [[0.5, 2], [0.5, 1]], order=1)
    array([ 2.,  7.])

    Above, the interpolated value of a[0.5, 0.5] gives output[0], while
    a[2, 1] is output[1].

    >>> inds = np.array([[0.5, 2], [0.5, 4]])
    >>> ndimage.map_coordinates(a, inds, order=1, cval=-33.3)
    array([  2. , -33.3])
    >>> ndimage.map_coordinates(a, inds, order=1, mode='nearest')
    array([ 2.,  8.])
    >>> ndimage.map_coordinates(a, inds, order=1, cval=0, output=bool)
    array([ True, False], dtype=bool)

    """
    if order < 0 or order > 5:
        raise RuntimeError('spline order not supported')
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    coordinates = numpy.asarray(coordinates)
    if numpy.iscomplexobj(coordinates):
        raise TypeError('Complex type not supported')
    output_shape = coordinates.shape[1:]
    if input.ndim < 1 or len(output_shape) < 1:
        raise RuntimeError('input and output rank must be > 0')
    if coordinates.shape[0] != input.ndim:
        raise RuntimeError('invalid shape for coordinate array')
    mode = _ni_support._extend_mode_to_code(mode)
    if prefilter and order > 1:
        filtered = spline_filter(input, order, output=numpy.float64)
    else:
        filtered = input
    output = _ni_support._get_output(output, input,
                                     shape=output_shape)
    _nd_image.geometric_transform(filtered, None, coordinates, None, None,
                                  output, order, mode, cval, None, None)
    return output


@docfiller
def affine_transform(input, matrix, offset=0.0, output_shape=None,
                     output=None, order=3,
                     mode='constant', cval=0.0, prefilter=True):
    """
    Apply an affine transformation.

    Given an output image pixel index vector ``o``, the pixel value
    is determined from the input image at position
    ``np.dot(matrix, o) + offset``.

    This does 'pull' (or 'backward') resampling, transforming the output space
    to the input to locate data. Affine transformations are often described in
    the 'push' (or 'forward') direction, transforming input to output. If you
    have a matrix for the 'push' transformation, use its inverse
    (:func:`numpy.linalg.inv`) in this function.

    Parameters
    ----------
    %(input)s
    matrix : ndarray
        The inverse coordinate transformation matrix, mapping output
        coordinates to input coordinates. If ``ndim`` is the number of
        dimensions of ``input``, the given matrix must have one of the
        following shapes:

            - ``(ndim, ndim)``: the linear transformation matrix for each
              output coordinate.
            - ``(ndim,)``: assume that the 2-D transformation matrix is
              diagonal, with the diagonal specified by the given value. A more
              efficient algorithm is then used that exploits the separability
              of the problem.
            - ``(ndim + 1, ndim + 1)``: assume that the transformation is
              specified using homogeneous coordinates [1]_. In this case, any
              value passed to ``offset`` is ignored.
            - ``(ndim, ndim + 1)``: as above, but the bottom row of a
              homogeneous transformation matrix is always ``[0, 0, ..., 1]``,
              and may be omitted.

    offset : float or sequence, optional
        The offset into the array where the transform is applied. If a float,
        `offset` is the same for each axis. If a sequence, `offset` should
        contain one value for each axis.
    output_shape : tuple of ints, optional
        Shape tuple.
    %(output)s
    order : int, optional
        The order of the spline interpolation, default is 3.
        The order has to be in the range 0-5.
    %(mode)s
    %(cval)s
    %(prefilter)s

    Returns
    -------
    affine_transform : ndarray
        The transformed input.

    Notes
    -----
    The given matrix and offset are used to find for each point in the
    output the corresponding coordinates in the input by an affine
    transformation. The value of the input at those coordinates is
    determined by spline interpolation of the requested order. Points
    outside the boundaries of the input are filled according to the given
    mode.

    .. versionchanged:: 0.18.0
        Previously, the exact interpretation of the affine transformation
        depended on whether the matrix was supplied as a 1-D or a
        2-D array. If a 1-D array was supplied
        to the matrix parameter, the output pixel value at index ``o``
        was determined from the input image at position
        ``matrix * (o + offset)``.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Homogeneous_coordinates
    """
    if order < 0 or order > 5:
        raise RuntimeError('spline order not supported')
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    if output_shape is None:
        output_shape = input.shape
    if input.ndim < 1 or len(output_shape) < 1:
        raise RuntimeError('input and output rank must be > 0')
    mode = _ni_support._extend_mode_to_code(mode)
    if prefilter and order > 1:
        filtered = spline_filter(input, order, output=numpy.float64)
    else:
        filtered = input
    output = _ni_support._get_output(output, input,
                                     shape=output_shape)
    matrix = numpy.asarray(matrix, dtype=numpy.float64)
    if matrix.ndim not in [1, 2] or matrix.shape[0] < 1:
        raise RuntimeError('no proper affine matrix provided')
    if (matrix.ndim == 2 and matrix.shape[1] == input.ndim + 1 and
            (matrix.shape[0] in [input.ndim, input.ndim + 1])):
        if matrix.shape[0] == input.ndim + 1:
            exptd = [0] * input.ndim + [1]
            if not numpy.all(matrix[input.ndim] == exptd):
                msg = ('Expected homogeneous transformation matrix with '
                       'shape %s for image shape %s, but bottom row was '
                       'not equal to %s' % (matrix.shape, input.shape, exptd))
                raise ValueError(msg)
        # assume input is homogeneous coordinate transformation matrix
        offset = matrix[:input.ndim, input.ndim]
        matrix = matrix[:input.ndim, :input.ndim]
    if matrix.shape[0] != input.ndim:
        raise RuntimeError('affine matrix has wrong number of rows')
    if matrix.ndim == 2 and matrix.shape[1] != output.ndim:
        raise RuntimeError('affine matrix has wrong number of columns')
    if not matrix.flags.contiguous:
        matrix = matrix.copy()
    offset = _ni_support._normalize_sequence(offset, input.ndim)
    offset = numpy.asarray(offset, dtype=numpy.float64)
    if offset.ndim != 1 or offset.shape[0] < 1:
        raise RuntimeError('no proper offset provided')
    if not offset.flags.contiguous:
        offset = offset.copy()
    if matrix.ndim == 1:
        warnings.warn(
            "The behavior of affine_transform with a 1-D "
            "array supplied for the matrix parameter has changed in "
            "SciPy 0.18.0."
        )
        _nd_image.zoom_shift(filtered, matrix, offset/matrix, output, order,
                             mode, cval)
    else:
        _nd_image.geometric_transform(filtered, None, None, matrix, offset,
                                      output, order, mode, cval, None, None)
    return output


@docfiller
def shift(input, shift, output=None, order=3, mode='constant', cval=0.0,
          prefilter=True):
    """
    Shift an array.

    The array is shifted using spline interpolation of the requested order.
    Points outside the boundaries of the input are filled according to the
    given mode.

    Parameters
    ----------
    %(input)s
    shift : float or sequence
        The shift along the axes. If a float, `shift` is the same for each
        axis. If a sequence, `shift` should contain one value for each axis.
    %(output)s
    order : int, optional
        The order of the spline interpolation, default is 3.
        The order has to be in the range 0-5.
    %(mode)s
    %(cval)s
    %(prefilter)s

    Returns
    -------
    shift : ndarray
        The shifted input.

    """
    if order < 0 or order > 5:
        raise RuntimeError('spline order not supported')
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    if input.ndim < 1:
        raise RuntimeError('input and output rank must be > 0')
    mode = _ni_support._extend_mode_to_code(mode)
    if prefilter and order > 1:
        filtered = spline_filter(input, order, output=numpy.float64)
    else:
        filtered = input
    output = _ni_support._get_output(output, input)
    shift = _ni_support._normalize_sequence(shift, input.ndim)
    shift = [-ii for ii in shift]
    shift = numpy.asarray(shift, dtype=numpy.float64)
    if not shift.flags.contiguous:
        shift = shift.copy()
    _nd_image.zoom_shift(filtered, None, shift, output, order, mode, cval)
    return output


@docfiller
def zoom(input, zoom, output=None, order=3, mode='constant', cval=0.0,
         prefilter=True):
    """
    Zoom an array.

    The array is zoomed using spline interpolation of the requested order.

    Parameters
    ----------
    %(input)s
    zoom : float or sequence
        The zoom factor along the axes. If a float, `zoom` is the same for each
        axis. If a sequence, `zoom` should contain one value for each axis.
    %(output)s
    order : int, optional
        The order of the spline interpolation, default is 3.
        The order has to be in the range 0-5.
    %(mode)s
    %(cval)s
    %(prefilter)s

    Returns
    -------
    zoom : ndarray
        The zoomed input.

    Examples
    --------
    >>> from scipy import ndimage, misc
    >>> import matplotlib.pyplot as plt

    >>> fig = plt.figure()
    >>> ax1 = fig.add_subplot(121)  # left side
    >>> ax2 = fig.add_subplot(122)  # right side
    >>> ascent = misc.ascent()
    >>> result = ndimage.zoom(ascent, 3.0)
    >>> ax1.imshow(ascent)
    >>> ax2.imshow(result)
    >>> plt.show()

    >>> print(ascent.shape)
    (512, 512)

    >>> print(result.shape)
    (1536, 1536)
    """
    if order < 0 or order > 5:
        raise RuntimeError('spline order not supported')
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError('Complex type not supported')
    if input.ndim < 1:
        raise RuntimeError('input and output rank must be > 0')
    mode = _ni_support._extend_mode_to_code(mode)
    if prefilter and order > 1:
        filtered = spline_filter(input, order, output=numpy.float64)
    else:
        filtered = input
    zoom = _ni_support._normalize_sequence(zoom, input.ndim)
    output_shape = tuple(
            [int(round(ii * jj)) for ii, jj in zip(input.shape, zoom)])

    zoom_div = numpy.array(output_shape, float) - 1
    # Zooming to infinite values is unpredictable, so just choose
    # zoom factor 1 instead
    zoom = numpy.divide(numpy.array(input.shape) - 1, zoom_div,
                        out=numpy.ones_like(input.shape, dtype=numpy.float64),
                        where=zoom_div != 0)

    output = _ni_support._get_output(output, input,
                                     shape=output_shape)
    zoom = numpy.ascontiguousarray(zoom)
    _nd_image.zoom_shift(filtered, zoom, None, output, order, mode, cval)
    return output


@docfiller
def rotate(input, angle, axes=(1, 0), reshape=True, output=None, order=3,
           mode='constant', cval=0.0, prefilter=True):
    """
    Rotate an array.

    The array is rotated in the plane defined by the two axes given by the
    `axes` parameter using spline interpolation of the requested order.

    Parameters
    ----------
    %(input)s
    angle : float
        The rotation angle in degrees.
    axes : tuple of 2 ints, optional
        The two axes that define the plane of rotation. Default is the first
        two axes.
    reshape : bool, optional
        If `reshape` is true, the output shape is adapted so that the input
        array is contained completely in the output. Default is True.
    %(output)s
    order : int, optional
        The order of the spline interpolation, default is 3.
        The order has to be in the range 0-5.
    %(mode)s
    %(cval)s
    %(prefilter)s

    Returns
    -------
    rotate : ndarray
        The rotated input.

    Examples
    --------
    >>> from scipy import ndimage, misc
    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure(figsize=(10, 3))
    >>> ax1, ax2, ax3 = fig.subplots(1, 3)
    >>> img = misc.ascent()
    >>> img_45 = ndimage.rotate(img, 45, reshape=False)
    >>> full_img_45 = ndimage.rotate(img, 45, reshape=True)
    >>> ax1.imshow(img, cmap='gray')
    >>> ax1.set_axis_off()
    >>> ax2.imshow(img_45, cmap='gray')
    >>> ax2.set_axis_off()
    >>> ax3.imshow(full_img_45, cmap='gray')
    >>> ax3.set_axis_off()
    >>> fig.set_tight_layout(True)
    >>> plt.show()
    >>> print(img.shape)
    (512, 512)
    >>> print(img_45.shape)
    (512, 512)
    >>> print(full_img_45.shape)
    (724, 724)

    """
    input_arr = numpy.asarray(input)
    ndim = input_arr.ndim

    if ndim < 2:
        raise ValueError('input array should be at least 2D')

    axes = list(axes)

    if len(axes) != 2:
        raise ValueError('axes should contain exactly two values')

    if not all([float(ax).is_integer() for ax in axes]):
        raise ValueError('axes should contain only integer values')

    if axes[0] < 0:
        axes[0] += ndim
    if axes[1] < 0:
        axes[1] += ndim
    if axes[0] < 0 or axes[1] < 0 or axes[0] >= ndim or axes[1] >= ndim:
        raise ValueError('invalid rotation plane specified')

    axes.sort()

    angle_rad = numpy.deg2rad(angle)
    c, s = numpy.cos(angle_rad), numpy.sin(angle_rad)

    rot_matrix = numpy.array([[c, s],
                              [-s, c]])

    img_shape = numpy.asarray(input_arr.shape)
    in_plane_shape = img_shape[axes]
    if reshape:
        # Compute transformed input bounds
        iy, ix = in_plane_shape
        out_bounds = rot_matrix @ [[0, 0, iy, iy],
                                   [0, ix, 0, ix]]
        # Compute the shape of the transformed input plane
        out_plane_shape = (out_bounds.ptp(axis=1) + 0.5).astype(int)
    else:
        out_plane_shape = img_shape[axes]

    out_center = rot_matrix @ ((out_plane_shape - 1) / 2)
    in_center = (in_plane_shape - 1) / 2
    offset = in_center - out_center

    output_shape = img_shape
    output_shape[axes] = out_plane_shape
    output_shape = tuple(output_shape)

    output = _ni_support._get_output(output, input_arr, shape=output_shape)

    if ndim <= 2:
        affine_transform(input_arr, rot_matrix, offset, output_shape, output,
                         order, mode, cval, prefilter)
    else:
        # If ndim > 2, the rotation is applied over all the planes
        # parallel to axes
        planes_coord = itertools.product(
            *[[slice(None)] if ax in axes else range(img_shape[ax])
              for ax in range(ndim)])

        out_plane_shape = tuple(out_plane_shape)

        for coordinates in planes_coord:
            ia = input_arr[coordinates]
            oa = output[coordinates]
            affine_transform(ia, rot_matrix, offset, out_plane_shape,
                             oa, order, mode, cval, prefilter)

    return output
