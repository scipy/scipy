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

import types
import math
import warnings
import numpy
import _ni_support
import _nd_image

def _extend_mode_to_code(mode):
    mode = _ni_support._extend_mode_to_code(mode)
    return mode

def spline_filter1d(input, order = 3, axis = -1, output = numpy.float64,
                    output_type = None):
    """Calculates a one-dimensional spline filter along the given axis.

    The lines of the array along the given axis are filtered by a
    spline filter. The order of the spline must be >= 2 and <= 5.
    """
    if order < 0 or order > 5:
        raise RuntimeError, 'spline order not supported'
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    output, return_value = _ni_support._get_output(output, input,
                                                    output_type)
    if order in [0, 1]:
        output[...] = numpy.array(input)
    else:
        axis = _ni_support._check_axis(axis, input.ndim)
        _nd_image.spline_filter1d(input, order, axis, output)
    return return_value


def spline_filter(input, order = 3, output = numpy.float64,
                  output_type = None):
    """Multi-dimensional spline filter.

    Note: The multi-dimensional filter is implemented as a sequence of
    one-dimensional spline filters. The intermediate arrays are stored
    in the same data type as the output. Therefore, for output types
    with a limited precision, the results may be imprecise because
    intermediate results may be stored with insufficient precision.
    """
    if order < 2 or order > 5:
        raise RuntimeError, 'spline order not supported'
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    output, return_value = _ni_support._get_output(output, input,
                                                    output_type)
    if order not in [0, 1] and input.ndim > 0:
        for axis in range(input.ndim):
            spline_filter1d(input, order, axis, output = output)
            input = output
    else:
        output[...] = input[...]
    return return_value

def geometric_transform(input, mapping, output_shape = None,
                        output_type = None, output = None, order = 3,
                        mode = 'constant', cval = 0.0, prefilter = True,
                        extra_arguments = (), extra_keywords = {}):
    """Apply an arbritrary geometric transform.

    The given mapping function is used to find, for each point in the
    output, the corresponding coordinates in the input. The value of the
    input at those coordinates is determined by spline interpolation of
    the requested order.

    mapping must be a callable object that accepts a tuple of length
    equal to the output array rank and returns the corresponding input
    coordinates as a tuple of length equal to the input array
    rank. Points outside the boundaries of the input are filled
    according to the given mode ('constant', 'nearest', 'reflect' or
    'wrap'). The output shape can optionally be given. If not given,
    it is equal to the input shape. The parameter prefilter determines
    if the input is pre-filtered before interpolation (necessary for
    spline interpolation of order > 1).  If False it is assumed that
    the input is already filtered. The extra_arguments and
    extra_keywords arguments can be used to provide extra arguments
    and keywords that are passed to the mapping function at each call.

    Example
    -------
    >>> a = arange(12.).reshape((4,3))
    >>> def shift_func(output_coordinates):
    ...     return (output_coordinates[0]-0.5, output_coordinates[1]-0.5)
    ...
    >>> print geometric_transform(a,shift_func)
    array([[ 0.    ,  0.    ,  0.    ],
           [ 0.    ,  1.3625,  2.7375],
           [ 0.    ,  4.8125,  6.1875],
           [ 0.    ,  8.2625,  9.6375]])
    """
    if order < 0 or order > 5:
        raise RuntimeError, 'spline order not supported'
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    if output_shape is None:
        output_shape = input.shape
    if input.ndim < 1 or len(output_shape) < 1:
        raise RuntimeError, 'input and output rank must be > 0'
    mode = _extend_mode_to_code(mode)
    if prefilter and order > 1:
        filtered = spline_filter(input, order, output = numpy.float64)
    else:
        filtered = input
    output, return_value = _ni_support._get_output(output, input,
                                        output_type, shape = output_shape)
    _nd_image.geometric_transform(filtered, mapping, None, None, None,
               output, order, mode, cval, extra_arguments, extra_keywords)
    return return_value


def map_coordinates(input, coordinates, output_type = None, output = None,
                order = 3, mode = 'constant', cval = 0.0, prefilter = True):
    """Apply an arbritrary coordinate transformation.

    The array of coordinates is used to find, for each point in the output,
    the corresponding coordinates in the input. The value of the input at
    those coordinates is determined by spline interpolation of the
    requested order.

    The shape of the output is derived from that of the coordinate
    array by dropping the first axis. The values of the array along
    the first axis are the coordinates in the input array at which the
    output value is found.  For example, if the input has dimensions
    (100,200,3), then the shape of coordinates will be (3,100,200,3),
    where coordinates[:,1,2,3] specify the input coordinate at which
    output[1,2,3] is found.

    Points outside the boundaries of the input are filled according to
    the given mode ('constant', 'nearest', 'reflect' or 'wrap'). The
    parameter prefilter determines if the input is pre-filtered before
    interpolation (necessary for spline interpolation of order >
    1). If False it is assumed that the input is already filtered.

    Example
    -------
    >>> a = arange(12.).reshape((4,3))
    >>> print a
    [[  0.   1.   2.]
     [  3.   4.   5.]
     [  6.   7.   8.]
     [  9.  10.  11.]]
    >>> output = map_coordinates(a,[[0.5, 2], [0.5, 1]],order=1)
    >>> print output
    [ 2. 7.]

    Here, the interpolated value of a[0.5,0.5] gives output[0], while
    a[2,1] is output[1].

    """
    if order < 0 or order > 5:
        raise RuntimeError, 'spline order not supported'
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    coordinates = numpy.asarray(coordinates)
    if numpy.iscomplexobj(coordinates):
        raise TypeError, 'Complex type not supported'
    output_shape = coordinates.shape[1:]
    if input.ndim < 1 or len(output_shape) < 1:
        raise RuntimeError, 'input and output rank must be > 0'
    if coordinates.shape[0] != input.ndim:
        raise RuntimeError, 'invalid shape for coordinate array'
    mode = _extend_mode_to_code(mode)
    if prefilter and order > 1:
        filtered = spline_filter(input, order, output = numpy.float64)
    else:
        filtered = input
    output, return_value = _ni_support._get_output(output, input,
                                        output_type, shape = output_shape)
    _nd_image.geometric_transform(filtered, None, coordinates, None, None,
               output, order, mode, cval, None, None)
    return return_value


def affine_transform(input, matrix, offset = 0.0, output_shape = None,
                     output_type = None, output = None, order = 3,
                     mode = 'constant', cval = 0.0, prefilter = True):
    """Apply an affine transformation.

    The given matrix and offset are used to find for each point in the
    output the corresponding coordinates in the input by an affine
    transformation. The value of the input at those coordinates is
    determined by spline interpolation of the requested order. Points
    outside the boundaries of the input are filled according to the given
    mode. The output shape can optionally be given. If not given it is
    equal to the input shape. The parameter prefilter determines if the
    input is pre-filtered before interpolation, if False it is assumed
    that the input is already filtered.

    The matrix must be two-dimensional or can also be given as a
    one-dimensional sequence or array. In the latter case, it is
    assumed that the matrix is diagonal. A more efficient algorithms
    is then applied that exploits the separability of the problem.
    """
    if order < 0 or order > 5:
        raise RuntimeError, 'spline order not supported'
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    if output_shape is None:
        output_shape = input.shape
    if input.ndim < 1 or len(output_shape) < 1:
        raise RuntimeError, 'input and output rank must be > 0'
    mode = _extend_mode_to_code(mode)
    if prefilter and order > 1:
        filtered = spline_filter(input, order, output = numpy.float64)
    else:
        filtered = input
    output, return_value = _ni_support._get_output(output, input,
                                        output_type, shape = output_shape)
    matrix = numpy.asarray(matrix, dtype = numpy.float64)
    if matrix.ndim not in [1, 2] or matrix.shape[0] < 1:
        raise RuntimeError, 'no proper affine matrix provided'
    if matrix.shape[0] != input.ndim:
        raise RuntimeError, 'affine matrix has wrong number of rows'
    if matrix.ndim == 2 and matrix.shape[1] != output.ndim:
        raise RuntimeError, 'affine matrix has wrong number of columns'
    if not matrix.flags.contiguous:
        matrix = matrix.copy()
    offset = _ni_support._normalize_sequence(offset, input.ndim)
    offset = numpy.asarray(offset, dtype = numpy.float64)
    if offset.ndim != 1 or offset.shape[0] < 1:
        raise RuntimeError, 'no proper offset provided'
    if not offset.flags.contiguous:
        offset = offset.copy()
    if matrix.ndim == 1:
        _nd_image.zoom_shift(filtered, matrix, offset, output, order,
                             mode, cval)
    else:
        _nd_image.geometric_transform(filtered, None, None, matrix, offset,
                            output, order, mode, cval, None, None)
    return return_value


def shift(input, shift, output_type = None, output = None, order = 3,
          mode = 'constant', cval = 0.0, prefilter = True):
    """Shift an array.

    The array is shifted using spline interpolation of the requested
    order. Points outside the boundaries of the input are filled according
    to the given mode. The parameter prefilter determines if the input is
    pre-filtered before interpolation, if False it is assumed that the
    input is already filtered.
    """
    if order < 0 or order > 5:
        raise RuntimeError, 'spline order not supported'
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    if input.ndim < 1:
        raise RuntimeError, 'input and output rank must be > 0'
    mode = _extend_mode_to_code(mode)
    if prefilter and order > 1:
        filtered = spline_filter(input, order, output = numpy.float64)
    else:
        filtered = input
    output, return_value = _ni_support._get_output(output, input,
                                                    output_type)
    shift = _ni_support._normalize_sequence(shift, input.ndim)
    shift = [-ii for ii in shift]
    shift = numpy.asarray(shift, dtype = numpy.float64)
    if not shift.flags.contiguous:
        shift = shift.copy()
    _nd_image.zoom_shift(filtered, None, shift, output, order, mode, cval)
    return return_value


def zoom(input, zoom, output_type = None, output = None, order = 3,
         mode = 'constant', cval = 0.0, prefilter = True):
    """Zoom an array.

    The array is zoomed using spline interpolation of the requested order.
    Points outside the boundaries of the input are filled according to the
    given mode. The parameter prefilter determines if the input is pre-
    filtered before interpolation, if False it is assumed that the input
    is already filtered.
    """
    if order < 0 or order > 5:
        raise RuntimeError, 'spline order not supported'
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    if input.ndim < 1:
        raise RuntimeError, 'input and output rank must be > 0'
    mode = _extend_mode_to_code(mode)
    if prefilter and order > 1:
        filtered = spline_filter(input, order, output = numpy.float64)
    else:
        filtered = input
    zoom = _ni_support._normalize_sequence(zoom, input.ndim)
    output_shape = [int(ii * jj) for ii, jj in zip(input.shape, zoom)]
    zoom = (numpy.array(input.shape)-1)/(numpy.array(output_shape,float)-1)
    output, return_value = _ni_support._get_output(output, input,
                                        output_type, shape = output_shape)
    zoom = numpy.asarray(zoom, dtype = numpy.float64)
    zoom = numpy.ascontiguousarray(zoom)
    _nd_image.zoom_shift(filtered, zoom, None, output, order, mode, cval)
    return return_value

def _minmax(coor, minc, maxc):
    if coor[0] < minc[0]:
        minc[0] = coor[0]
    if coor[0] > maxc[0]:
        maxc[0] = coor[0]
    if coor[1] < minc[1]:
        minc[1] = coor[1]
    if coor[1] > maxc[1]:
        maxc[1] = coor[1]
    return minc, maxc

def rotate(input, angle, axes = (1, 0), reshape = True,
           output_type = None, output = None, order = 3,
           mode = 'constant', cval = 0.0, prefilter = True):
    """Rotate an array.

    The array is rotated in the plane defined by the two axes given by the
    axes parameter using spline interpolation of the requested order. The
    angle is given in degrees. Points outside the boundaries of the input
    are filled according to the given mode. If reshape is true, the output
    shape is adapted so that the input array is contained completely in
    the output. The parameter prefilter determines if the input is pre-
    filtered before interpolation, if False it is assumed that the input
    is already filtered.
    """
    input = numpy.asarray(input)
    axes = list(axes)
    rank = input.ndim
    if axes[0] < 0:
        axes[0] += rank
    if axes[1] < 0:
        axes[1] += rank
    if axes[0] < 0 or axes[1] < 0 or axes[0] > rank or axes[1] > rank:
        raise RuntimeError, 'invalid rotation plane specified'
    if axes[0] > axes[1]:
        axes = axes[1], axes[0]
    angle = numpy.pi / 180 * angle
    m11 = math.cos(angle)
    m12 = math.sin(angle)
    m21 = -math.sin(angle)
    m22 = math.cos(angle)
    matrix = numpy.array([[m11, m12],
                             [m21, m22]], dtype = numpy.float64)
    iy = input.shape[axes[0]]
    ix = input.shape[axes[1]]
    if reshape:
        mtrx = numpy.array([[ m11, -m21],
                               [-m12,  m22]], dtype = numpy.float64)
        minc = [0, 0]
        maxc = [0, 0]
        coor = numpy.dot(mtrx, [0, ix])
        minc, maxc = _minmax(coor, minc, maxc)
        coor = numpy.dot(mtrx, [iy, 0])
        minc, maxc = _minmax(coor, minc, maxc)
        coor = numpy.dot(mtrx, [iy, ix])
        minc, maxc = _minmax(coor, minc, maxc)
        oy = int(maxc[0] - minc[0] + 0.5)
        ox = int(maxc[1] - minc[1] + 0.5)
    else:
        oy = input.shape[axes[0]]
        ox = input.shape[axes[1]]
    offset = numpy.zeros((2,), dtype = numpy.float64)
    offset[0] = float(oy) / 2.0 - 0.5
    offset[1] = float(ox) / 2.0 - 0.5
    offset = numpy.dot(matrix, offset)
    tmp = numpy.zeros((2,), dtype = numpy.float64)
    tmp[0] = float(iy) / 2.0 - 0.5
    tmp[1] = float(ix) / 2.0 - 0.5
    offset = tmp - offset
    output_shape = list(input.shape)
    output_shape[axes[0]] = oy
    output_shape[axes[1]] = ox
    output_shape = tuple(output_shape)
    output, return_value = _ni_support._get_output(output, input,
                                        output_type, shape = output_shape)
    if input.ndim <= 2:
        affine_transform(input, matrix, offset, output_shape, None, output,
                         order, mode, cval, prefilter)
    else:
        coordinates = []
        size = numpy.product(input.shape,axis=0)
        size /= input.shape[axes[0]]
        size /= input.shape[axes[1]]
        for ii in range(input.ndim):
            if ii not in axes:
                coordinates.append(0)
            else:
                coordinates.append(slice(None, None, None))
        iter_axes = range(input.ndim)
        iter_axes.reverse()
        iter_axes.remove(axes[0])
        iter_axes.remove(axes[1])
        os = (output_shape[axes[0]], output_shape[axes[1]])
        for ii in range(size):
            ia = input[tuple(coordinates)]
            oa = output[tuple(coordinates)]
            affine_transform(ia, matrix, offset, os, None, oa, order, mode,
                             cval, prefilter)
            for jj in iter_axes:
                if coordinates[jj] < input.shape[jj] - 1:
                    coordinates[jj] += 1
                    break
                else:
                    coordinates[jj] = 0
    return return_value
