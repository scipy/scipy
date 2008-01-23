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

import math
import numpy
import _ni_support
import _nd_image

_mode_doc = \
"""The mode parameter determines how the array borders are handled,
    where cval is the value when mode is equal to 'constant'.  Other
    modes are 'nearest', 'mirror', 'reflect' and 'wrap'."""

_origin_doc = """

    The origin parameter controls the placement of the filter."""

def moredoc(*args):
    def decorate(f):
        if not f.__doc__: f.__doc__ = ""
        for a in args: f.__doc__ += a
        return f
    return decorate

@moredoc(_mode_doc, _origin_doc)
def correlate1d(input, weights, axis = -1, output = None, mode = "reflect",
                cval = 0.0, origin = 0):
    """Calculate a one-dimensional correlation along the given axis.

    The lines of the array along the given axis are correlated with the
    given weights. The weights parameter must be a one-dimensional sequence
    of numbers.

    """
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    output, return_value = _ni_support._get_output(output, input)
    weights = numpy.asarray(weights, dtype = numpy.float64)
    if weights.ndim != 1 or weights.shape[0] < 1:
        raise RuntimeError, 'no filter weights given'
    if not weights.flags.contiguous:
        weights = weights.copy()
    axis = _ni_support._check_axis(axis, input.ndim)
    if ((len(weights) // 2 + origin < 0) or
        (len(weights) // 2 + origin > len(weights))):
        raise ValueError, 'invalid origin'
    mode = _ni_support._extend_mode_to_code(mode)
    _nd_image.correlate1d(input, weights, axis, output, mode, cval,
                          origin)
    return return_value

@moredoc(_mode_doc, _origin_doc)
def convolve1d(input, weights, axis = -1, output = None, mode = "reflect",
               cval = 0.0, origin = 0):
    """Calculate a one-dimensional convolution along the given axis.

    The lines of the array along the given axis are convolved with the
    given weights. The weights parameter must be a one-dimensional sequence
    of numbers.

    """
    weights = weights[::-1]
    origin = -origin
    if not len(weights) & 1:
        origin -= 1
    return correlate1d(input, weights, axis, output, mode, cval, origin)

@moredoc(_mode_doc, _origin_doc)
def gaussian_filter1d(input, sigma, axis = -1, order = 0, output = None,
                      mode = "reflect", cval = 0.0):
    """One-dimensional Gaussian filter.

    The standard-deviation of the Gaussian filter is given by
    sigma. An order of 0 corresponds to convolution with a Gaussian
    kernel. An order of 1, 2, or 3 corresponds to convolution with the
    first, second or third derivatives of a Gaussian. Higher order
    derivatives are not implemented.

    """
    sd = float(sigma)
    # make the length of the filter equal to 4 times the standard
    # deviations:
    lw = int(4.0 * sd + 0.5)
    weights = [0.0] * (2 * lw + 1)
    weights[lw] = 1.0
    sum = 1.0
    sd = sd * sd
    # calculate the kernel:
    for ii in range(1, lw + 1):
        tmp = math.exp(-0.5 * float(ii * ii) / sd)
        weights[lw + ii] = tmp
        weights[lw - ii] = tmp
        sum += 2.0 * tmp
    for ii in range(2 * lw + 1):
        weights[ii] /= sum
    # implement first, second and third order derivatives:
    if order == 1 : # first derivative
        weights[lw] = 0.0
        for ii in range(1, lw + 1):
            x = float(ii)
            tmp = -x / sd * weights[lw + ii]
            weights[lw + ii] = -tmp
            weights[lw - ii] = tmp
    elif order == 2: # second derivative
        weights[lw] *= -1.0 / sd
        for ii in range(1, lw + 1):
            x = float(ii)
            tmp = (x * x / sd - 1.0) * weights[lw + ii] / sd
            weights[lw + ii] = tmp
            weights[lw - ii] = tmp
    elif order == 3: # third derivative
        weights[lw] = 0.0
        sd2 = sd * sd
        for ii in range(1, lw + 1):
            x = float(ii)
            tmp = (3.0 - x * x / sd) * x * weights[lw + ii] / sd / sd
            weights[lw + ii] = -tmp
            weights[lw - ii] = tmp
    return correlate1d(input, weights, axis, output, mode, cval, 0)

@moredoc(_mode_doc, _origin_doc)
def gaussian_filter(input, sigma, order = 0, output = None,
                  mode = "reflect", cval = 0.0):
    """Multi-dimensional Gaussian filter.

    The standard-deviations of the Gaussian filter are given for each
    axis as a sequence, or as a single number, in which case it is
    equal for all axes. The order of the filter along each axis is
    given as a sequence of integers, or as a single number. An order
    of 0 corresponds to convolution with a Gaussian kernel. An order
    of 1, 2, or 3 corresponds to convolution with the first, second or
    third derivatives of a Gaussian. Higher order derivatives are not
    implemented.'

    Note: The multi-dimensional filter is implemented as a sequence of
    one-dimensional convolution filters. The intermediate arrays are
    stored in the same data type as the output. Therefore, for output
    types with a limited precision, the results may be imprecise
    because intermediate results may be stored with insufficient
    precision.

    """
    input = numpy.asarray(input)
    output, return_value = _ni_support._get_output(output, input)
    orders = _ni_support._normalize_sequence(order, input.ndim)
    sigmas = _ni_support._normalize_sequence(sigma, input.ndim)
    axes = range(input.ndim)
    axes = [(axes[ii], sigmas[ii], orders[ii])
                        for ii in range(len(axes)) if sigmas[ii] > 1e-15]
    if len(axes) > 0:
        for axis, sigma, order in axes:
            gaussian_filter1d(input, sigma, axis, order, output,
                              mode, cval)
            input = output
    else:
        output[...] = input[...]
    return return_value

@moredoc(_mode_doc, _origin_doc)
def prewitt(input, axis = -1, output = None, mode = "reflect", cval = 0.0):
    """Calculate a Prewitt filter.

    """
    input = numpy.asarray(input)
    axis = _ni_support._check_axis(axis, input.ndim)
    output, return_value = _ni_support._get_output(output, input)
    correlate1d(input, [-1, 0, 1], axis, output, mode, cval, 0)
    axes = [ii for ii in range(input.ndim) if ii != axis]
    for ii in axes:
        correlate1d(output, [1, 1, 1], ii, output, mode, cval, 0,)
    return return_value

@moredoc(_mode_doc, _origin_doc)
def sobel(input, axis = -1, output = None, mode = "reflect", cval = 0.0):
    """Calculate a Sobel filter.

    """
    input = numpy.asarray(input)
    axis = _ni_support._check_axis(axis, input.ndim)
    output, return_value = _ni_support._get_output(output, input)
    correlate1d(input, [-1, 0, 1], axis, output, mode, cval, 0)
    axes = [ii for ii in range(input.ndim) if ii != axis]
    for ii in axes:
        correlate1d(output, [1, 2, 1], ii, output, mode, cval, 0)
    return return_value

@moredoc(_mode_doc, _origin_doc)
def generic_laplace(input, derivative2, output = None, mode = "reflect",
                    cval = 0.0, extra_arguments = (), extra_keywords = {}):
    """Calculate a multidimensional laplace filter using the provided
    second derivative function.

    The derivative2 parameter must be a callable with the following
    signature:

    derivative2(input, axis, output, mode, cval,
                *extra_arguments, **extra_keywords)

    The extra_arguments and extra_keywords arguments can be used to pass
    extra arguments and keywords that are passed to derivative2 at each
    call.

    """
    input = numpy.asarray(input)
    output, return_value = _ni_support._get_output(output, input)
    axes = range(input.ndim)
    if len(axes) > 0:
        derivative2(input, axes[0], output, mode, cval,
                    *extra_arguments, **extra_keywords)
        for ii in range(1, len(axes)):
            tmp = derivative2(input, axes[ii], output.dtype, mode, cval,
                              *extra_arguments, **extra_keywords)
            output += tmp
    else:
        output[...] = input[...]
    return return_value

@moredoc(_mode_doc, _origin_doc)
def laplace(input, output = None, mode = "reflect", cval = 0.0):
    """Calculate a multidimensional laplace filter using an estimation
    for the second derivative based on differences.

    """
    def derivative2(input, axis, output, mode, cval):
        return correlate1d(input, [1, -2, 1], axis, output, mode, cval, 0)
    return generic_laplace(input, derivative2, output, mode, cval)

@moredoc(_mode_doc, _origin_doc)
def gaussian_laplace(input, sigma, output = None, mode = "reflect",
                     cval = 0.0):
    """Calculate a multidimensional laplace filter using gaussian
    second derivatives.

    The standard-deviations of the Gaussian filter are given for each
    axis as a sequence, or as a single number, in which case it is
    equal for all axes..

    """
    input = numpy.asarray(input)
    def derivative2(input, axis, output, mode, cval, sigma):
        order = [0] * input.ndim
        order[axis] = 2
        return gaussian_filter(input, sigma, order, output, mode, cval)
    return generic_laplace(input, derivative2, output, mode, cval,
                           extra_arguments = (sigma,))

@moredoc(_mode_doc, _origin_doc)
def generic_gradient_magnitude(input, derivative, output = None,
                mode = "reflect", cval = 0.0,
                extra_arguments = (), extra_keywords = {}):
    """Calculate a gradient magnitude using the provdide function for
    the gradient.

    The derivative parameter must be a callable with the following
    signature:

    derivative(input, axis, output, mode, cval,
               *extra_arguments, **extra_keywords)

    The extra_arguments and extra_keywords arguments can be used to pass
    extra arguments and keywords that are passed to derivative2 at each
    call.

    """
    input = numpy.asarray(input)
    output, return_value = _ni_support._get_output(output, input)
    axes = range(input.ndim)
    if len(axes) > 0:
        derivative(input, axes[0], output, mode, cval,
                   *extra_arguments, **extra_keywords)
        numpy.multiply(output, output, output)
        for ii in range(1, len(axes)):
            tmp = derivative(input, axes[ii], output.dtype, mode, cval,
                             *extra_arguments, **extra_keywords)
            numpy.multiply(tmp, tmp, tmp)
            output += tmp
        numpy.sqrt(output, output)
    else:
        output[...] = input[...]
    return return_value

@moredoc(_mode_doc, _origin_doc)
def gaussian_gradient_magnitude(input, sigma, output = None,
                mode = "reflect", cval = 0.0):
    """Calculate a multidimensional gradient magnitude using gaussian
    derivatives.

    The standard-deviations of the Gaussian filter are given for each
    axis as a sequence, or as a single number, in which case it is
    equal for all axes..

    """
    input = numpy.asarray(input)
    def derivative(input, axis, output, mode, cval, sigma):
        order = [0] * input.ndim
        order[axis] = 1
        return gaussian_filter(input, sigma, order, output, mode, cval)
    return generic_gradient_magnitude(input, derivative, output, mode,
                            cval, extra_arguments = (sigma,))

def _correlate_or_convolve(input, weights, output, mode, cval, origin,
                           convolution):
    input = numpy.asarray(input)
    if numpy.iscomplexobj(int):
        raise TypeError, 'Complex type not supported'
    origins = _ni_support._normalize_sequence(origin, input.ndim)
    weights = numpy.asarray(weights, dtype = numpy.float64)
    wshape = [ii for ii in weights.shape if ii > 0]
    if len(wshape) != input.ndim:
        raise RuntimeError, 'filter weights array has incorrect shape.'
    if convolution:
        weights = weights[tuple([slice(None, None, -1)] * weights.ndim)]
        for ii in range(len(origins)):
            origins[ii] = -origins[ii]
            if not weights.shape[ii] & 1:
                origins[ii] -= 1
    for origin, lenw in zip(origins, wshape):
        if (lenw // 2 + origin < 0) or (lenw // 2 + origin > lenw):
            raise ValueError, 'invalid origin'
    if not weights.flags.contiguous:
        weights = weights.copy()
    output, return_value = _ni_support._get_output(output, input)
    mode = _ni_support._extend_mode_to_code(mode)
    _nd_image.correlate(input, weights, output, mode, cval, origins)
    return return_value

@moredoc(_mode_doc, _origin_doc)
def correlate(input, weights, output = None, mode = 'reflect', cval = 0.0,
              origin = 0):
    """Multi-dimensional correlation.

    The array is correlated with the given kernel.

    """
    return _correlate_or_convolve(input, weights, output, mode, cval,
                                  origin, False)

@moredoc(_mode_doc, _origin_doc)
def convolve(input, weights, output = None, mode = 'reflect', cval = 0.0,
             origin = 0):
    """Multi-dimensional convolution.

    The array is convolved with the given kernel.

    """
    return _correlate_or_convolve(input, weights, output, mode, cval,
                                  origin, True)

@moredoc(_mode_doc, _origin_doc)
def uniform_filter1d(input, size, axis = -1, output = None,
                     mode = "reflect", cval = 0.0, origin = 0):
    """Calculate a one-dimensional uniform filter along the given axis.

    The lines of the array along the given axis are filtered with a
    uniform filter of given size.

    """
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    axis = _ni_support._check_axis(axis, input.ndim)
    if size < 1:
        raise RuntimeError, 'incorrect filter size'
    output, return_value = _ni_support._get_output(output, input)
    if (size // 2 + origin < 0) or (size // 2 + origin > size):
        raise ValueError, 'invalid origin'
    mode = _ni_support._extend_mode_to_code(mode)
    _nd_image.uniform_filter1d(input, size, axis, output, mode, cval,
                               origin)
    return return_value

@moredoc(_mode_doc, _origin_doc)
def uniform_filter(input, size = 3, output = None, mode = "reflect",
                   cval = 0.0, origin = 0):
    """Multi-dimensional uniform filter.

    The sizes of the uniform filter are given for each axis as a
    sequence, or as a single number, in which case the size is equal
    for all axes.

    The multi-dimensional filter is implemented as a sequence of
    one-dimensional uniform filters. The intermediate arrays are stored
    in the same data type as the output. Therefore, for output types
    with a limited precision, the results may be imprecise because
    intermediate results may be stored with insufficient precision.

    """
    input = numpy.asarray(input)
    output, return_value = _ni_support._get_output(output, input)
    sizes = _ni_support._normalize_sequence(size, input.ndim)
    origins = _ni_support._normalize_sequence(origin, input.ndim)
    axes = range(input.ndim)
    axes = [(axes[ii], sizes[ii], origins[ii])
                           for ii in range(len(axes)) if sizes[ii] > 1]
    if len(axes) > 0:
        for axis, size, origin in axes:
            uniform_filter1d(input, int(size), axis, output, mode,
                             cval, origin)
            input = output
    else:
        output[...] = input[...]
    return return_value

@moredoc(_mode_doc, _origin_doc)
def minimum_filter1d(input, size, axis = -1, output = None,
                     mode = "reflect", cval = 0.0, origin = 0):
    """Calculate a one-dimensional minimum filter along the given axis.

    The lines of the array along the given axis are filtered with a
    minimum filter of given size.

    """
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    axis = _ni_support._check_axis(axis, input.ndim)
    if size < 1:
        raise RuntimeError, 'incorrect filter size'
    output, return_value = _ni_support._get_output(output, input)
    if (size // 2 + origin < 0) or (size // 2 + origin > size):
        raise ValueError, 'invalid origin'
    mode = _ni_support._extend_mode_to_code(mode)
    _nd_image.min_or_max_filter1d(input, size, axis, output, mode, cval,
                                  origin, 1)
    return return_value

@moredoc(_mode_doc, _origin_doc)
def maximum_filter1d(input, size, axis = -1, output = None,
                     mode = "reflect", cval = 0.0, origin = 0):
    """Calculate a one-dimensional maximum filter along the given axis.

    The lines of the array along the given axis are filtered with a
    maximum filter of given size.

    """
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    axis = _ni_support._check_axis(axis, input.ndim)
    if size < 1:
        raise RuntimeError, 'incorrect filter size'
    output, return_value = _ni_support._get_output(output, input)
    if (size // 2 + origin < 0) or (size // 2 + origin > size):
        raise ValueError, 'invalid origin'
    mode = _ni_support._extend_mode_to_code(mode)
    _nd_image.min_or_max_filter1d(input, size, axis, output, mode, cval,
                                  origin, 0)
    return return_value

def _min_or_max_filter(input, size, footprint, structure, output, mode,
                       cval, origin, minimum):
    if structure is None:
        if footprint is None:
            if size is None:
                raise RuntimeError, "no footprint provided"
            separable= True
        else:
            footprint = numpy.asarray(footprint)
            footprint = footprint.astype(bool)
            if numpy.alltrue(numpy.ravel(footprint),axis=0):
                size = footprint.shape
                footprint = None
                separable = True
            else:
                separable = False
    else:
        structure = numpy.asarray(structure, dtype = numpy.float64)
        separable = False
        if footprint is None:
            footprint = numpy.ones(structure.shape, bool)
        else:
            footprint = numpy.asarray(footprint)
            footprint = footprint.astype(bool)
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    output, return_value = _ni_support._get_output(output, input)
    origins = _ni_support._normalize_sequence(origin, input.ndim)
    if separable:
        sizes = _ni_support._normalize_sequence(size, input.ndim)
        axes = range(input.ndim)
        axes = [(axes[ii], sizes[ii], origins[ii])
                               for ii in range(len(axes)) if sizes[ii] > 1]
        if minimum:
            filter = minimum_filter1d
        else:
            filter = maximum_filter1d
        if len(axes) > 0:
            for axis, size, origin in axes:
                filter(input, int(size), axis, output, mode, cval, origin)
                input = output
        else:
            output[...] = input[...]
    else:
        fshape = [ii for ii in footprint.shape if ii > 0]
        if len(fshape) != input.ndim:
            raise RuntimeError, 'footprint array has incorrect shape.'
        for origin, lenf in zip(origins, fshape):
            if (lenf // 2 + origin < 0) or (lenf // 2 + origin > lenf):
                raise ValueError, 'invalid origin'
        if not footprint.flags.contiguous:
            footprint = footprint.copy()
        if structure is not None:
            if len(structure.shape) != input.ndim:
                raise RuntimeError, 'structure array has incorrect shape'
            if not structure.flags.contiguous:
                structure = structure.copy()
        mode = _ni_support._extend_mode_to_code(mode)
        _nd_image.min_or_max_filter(input, footprint, structure, output,
                                    mode, cval, origins, minimum)
    return return_value

@moredoc(_mode_doc, _origin_doc)
def minimum_filter(input, size = None, footprint = None, output = None,
      mode = "reflect", cval = 0.0, origin = 0):
    """Calculates a multi-dimensional minimum filter.

    Either a size or a footprint with the filter must be
    provided. An output array can optionally be provided.

    """
    return _min_or_max_filter(input, size, footprint, None, output, mode,
                              cval, origin, 1)

@moredoc(_mode_doc, _origin_doc)
def maximum_filter(input, size = None, footprint = None, output = None,
      mode = "reflect", cval = 0.0, origin = 0):
    """Calculates a multi-dimensional maximum filter.

    Either a size or a footprint with the filter must be
    provided. An output array can optionally be provided.

    """
    return _min_or_max_filter(input, size, footprint, None, output, mode,
                              cval, origin, 0)


def _rank_filter(input, rank, size = None, footprint = None, output = None,
     mode = "reflect", cval = 0.0, origin = 0, operation = 'rank'):
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    origins = _ni_support._normalize_sequence(origin, input.ndim)
    if footprint is None:
        if size is None:
            raise RuntimeError, "no footprint or filter size provided"
        sizes = _ni_support._normalize_sequence(size, input.ndim)
        footprint = numpy.ones(sizes, dtype = bool)
    else:
        footprint = numpy.asarray(footprint, dtype = bool)
    fshape = [ii for ii in footprint.shape if ii > 0]
    if len(fshape) != input.ndim:
        raise RuntimeError, 'filter footprint array has incorrect shape.'
    for origin, lenf in zip(origins, fshape):
        if (lenf // 2 + origin < 0) or (lenf // 2 + origin > lenf):
            raise ValueError, 'invalid origin'
    if not footprint.flags.contiguous:
        footprint = footprint.copy()
    filter_size = numpy.where(footprint, 1, 0).sum()
    if operation == 'median':
        rank = filter_size // 2
    elif operation == 'percentile':
        percentile = rank
        if percentile < 0.0:
            percentile += 100.0
        if percentile < 0 or percentile > 100:
            raise RuntimeError, 'invalid percentile'
        if percentile == 100.0:
            rank = filter_size - 1
        else:
            rank = int(float(filter_size) * percentile / 100.0)
    if rank < 0:
        rank += filter_size
    if rank < 0  or rank >= filter_size:
        raise RuntimeError, 'rank not within filter footprint size'
    if rank == 0:
        return minimum_filter(input, None, footprint, output, mode, cval,
                              origin)
    elif rank == filter_size - 1:
        return maximum_filter(input, None, footprint, output, mode, cval,
                              origin)
    else:
        output, return_value = _ni_support._get_output(output, input)
        mode = _ni_support._extend_mode_to_code(mode)
        _nd_image.rank_filter(input, rank, footprint, output, mode, cval,
                              origins)
        return return_value

@moredoc(_mode_doc, _origin_doc)
def rank_filter(input, rank, size = None, footprint = None, output = None,
      mode = "reflect", cval = 0.0, origin = 0):
    """Calculates a multi-dimensional rank filter.

    The rank parameter may be less then zero, i.e., rank = -1
    indicates the larges element. Either a size or a footprint with
    the filter must be provided. An output array can optionally be
    provided.

    """
    return _rank_filter(input, rank, size, footprint, output, mode, cval,
                        origin, 'rank')

@moredoc(_mode_doc, _origin_doc)
def median_filter(input, size = None, footprint = None, output = None,
      mode = "reflect", cval = 0.0, origin = 0):
    """Calculates a multi-dimensional median filter.

    Either a size or a footprint with the filter must be provided. An
    output array can optionally be provided.

    """
    return _rank_filter(input, 0, size, footprint, output, mode, cval,
                        origin, 'median')

@moredoc(_mode_doc, _origin_doc)
def percentile_filter(input, percentile, size = None, footprint = None,
                 output = None, mode = "reflect", cval = 0.0, origin = 0):
    """Calculates a multi-dimensional percentile filter.

    The percentile parameter may be less then zero, i.e., percentile =
    -20 equals percentile = 80. Either a size or a footprint with the
    filter must be provided. An output array can optionally be
    provided.

    """
    return _rank_filter(input, percentile, size, footprint, output, mode,
                                   cval, origin, 'percentile')

@moredoc(_mode_doc, _origin_doc)
def generic_filter1d(input, function, filter_size, axis = -1,
                 output = None, mode = "reflect", cval = 0.0, origin = 0,
                 extra_arguments = (), extra_keywords = {}):
    """Calculate a one-dimensional filter along the given axis.

    The function iterates over the lines of the array, calling the
    given function at each line. The arguments of the line are the
    input line, and the output line. The input and output lines are 1D
    double arrays.  The input line is extended appropiately according
    to the filter size and origin. The output line must be modified
    in-place with the result.  The extra_arguments and extra_keywords
    arguments can be used to pass extra arguments and keywords that
    are passed to the function at each call.

    """
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    output, return_value = _ni_support._get_output(output, input)
    if filter_size < 1:
        raise RuntimeError, 'invalid filter size'
    axis = _ni_support._check_axis(axis, input.ndim)
    if ((filter_size // 2 + origin < 0) or
        (filter_size // 2 + origin > filter_size)):
        raise ValueError, 'invalid origin'
    mode = _ni_support._extend_mode_to_code(mode)
    _nd_image.generic_filter1d(input, function, filter_size, axis, output,
                      mode, cval, origin, extra_arguments, extra_keywords)
    return return_value

@moredoc(_mode_doc, _origin_doc)
def generic_filter(input, function, size = None, footprint = None,
                   output = None, mode = "reflect", cval = 0.0, origin = 0,
                   extra_arguments = (), extra_keywords = {}):
    """Calculates a multi-dimensional filter using the given function.

    At each element the provided function is called. The input values
    within the filter footprint at that element are passed to the function
    as a 1D array of double values.

    Either a size or a footprint with the filter must be provided. An
    output array can optionally be provided. The extra_arguments and
    extra_keywords arguments can be used to pass extra arguments and
    keywords that are passed to the function at each call.

    """
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    origins = _ni_support._normalize_sequence(origin, input.ndim)
    if footprint is None:
        if size is None:
            raise RuntimeError, "no footprint or filter size provided"
        sizes = _ni_support._normalize_sequence(size, input.ndim)
        footprint = numpy.ones(size, dtype = bool)
    else:
        footprint = numpy.asarray(footprint)
        footprint = footprint.astype(bool)
    fshape = [ii for ii in footprint.shape if ii > 0]
    if len(fshape) != input.ndim:
        raise RuntimeError, 'filter footprint array has incorrect shape.'
    for origin, lenf in zip(origins, fshape):
        if (lenf // 2 + origin < 0) or (lenf // 2 + origin > lenf):
            raise ValueError, 'invalid origin'
    if not footprint.flags.contiguous:
        footprint = footprint.copy()
    output, return_value = _ni_support._get_output(output, input)
    mode = _ni_support._extend_mode_to_code(mode)
    _nd_image.generic_filter(input, function, footprint, output, mode,
                         cval, origins, extra_arguments, extra_keywords)
    return return_value
