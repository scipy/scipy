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

import numpy
import _ni_support
import _nd_image
import filters
import types


def _center_is_true(structure, origin):
    structure = numpy.array(structure)
    coor = tuple([oo + ss // 2 for ss, oo in zip(structure.shape,
                                                 origin)])
    return bool(structure[coor])

def iterate_structure(structure, iterations, origin = None):
    """Iterate a structure by dilating it with itself.

    If origin is None, only the iterated structure is returned. If
    not, a tuple of the iterated structure and the modified origin is
    returned.
    """
    structure = numpy.asarray(structure)
    if iterations < 2:
        return structure.copy()
    ni = iterations - 1
    shape = [ii + ni * (ii - 1) for ii in structure.shape]
    pos = [ni * (structure.shape[ii] / 2) for ii in range(len(shape))]
    slc = [slice(pos[ii], pos[ii] + structure.shape[ii], None)
           for ii in range(len(shape))]
    out = numpy.zeros(shape, bool)
    out[slc] = structure != 0
    out = binary_dilation(out, structure, iterations = ni)
    if origin is None:
        return out
    else:
        origin = _ni_support._normalize_sequence(origin, structure.ndim)
        origin = [iterations * o for o in origin]
        return out, origin

def generate_binary_structure(rank, connectivity):
    """Generate a binary structure for binary morphological operations.

    The inputs are the rank of the array to which the structure will
    be applied and the square of the connectivity of the structure.
    """
    if connectivity < 1:
        connectivity = 1
    if rank < 1:
        if connectivity < 1:
            return numpy.array(0, dtype = bool)
        else:
            return numpy.array(1, dtype = bool)
    output = numpy.zeros([3] * rank, bool)
    output = numpy.fabs(numpy.indices([3] * rank) - 1)
    output = numpy.add.reduce(output, 0)
    return numpy.asarray(output <= connectivity, dtype = bool)


def _binary_erosion(input, structure, iterations, mask, output,
                    border_value, origin, invert, brute_force):
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    if structure is None:
        structure = generate_binary_structure(input.ndim, 1)
    else:
        structure = numpy.asarray(structure)
        structure = structure.astype(bool)
    if structure.ndim != input.ndim:
        raise RuntimeError, 'structure rank must equal input rank'
    if not structure.flags.contiguous:
        structure = structure.copy()
    if numpy.product(structure.shape,axis=0) < 1:
        raise RuntimeError, 'structure must not be empty'
    if mask is not None:
        mask = numpy.asarray(mask)
        if mask.shape != input.shape:
            raise RuntimeError, 'mask and input must have equal sizes'
    origin = _ni_support._normalize_sequence(origin, input.ndim)
    cit = _center_is_true(structure, origin)
    if isinstance(output, numpy.ndarray):
        if numpy.iscomplexobj(output):
            raise TypeError, 'Complex output type not supported'
    else:
        output = bool
    output, return_value = _ni_support._get_output(output, input)


    if iterations == 1:
        _nd_image.binary_erosion(input, structure, mask, output,
                                     border_value, origin, invert, cit, 0)
        return return_value
    elif cit and not brute_force:
        changed, coordinate_list = _nd_image.binary_erosion(input,
             structure, mask, output, border_value, origin, invert, cit, 1)
        structure = structure[tuple([slice(None, None, -1)] *
                                    structure.ndim)]
        for ii in range(len(origin)):
            origin[ii] = -origin[ii]
            if not structure.shape[ii] & 1:
                origin[ii] -= 1
        if mask is not None:
            msk = numpy.asarray(mask)
            msk = mask.astype(numpy.int8)
            if msk is mask:
                msk = mask.copy()
            mask = msk
        if not structure.flags.contiguous:
            structure = structure.copy()
        _nd_image.binary_erosion2(output, structure, mask, iterations - 1,
                                  origin, invert, coordinate_list)
        return return_value
    else:
        tmp_in = numpy.zeros(input.shape, bool)
        if return_value is None:
            tmp_out = output
        else:
            tmp_out = numpy.zeros(input.shape, bool)
        if not iterations & 1:
            tmp_in, tmp_out = tmp_out, tmp_in
        changed = _nd_image.binary_erosion(input, structure, mask,
                            tmp_out, border_value, origin, invert, cit, 0)
        ii = 1
        while (ii < iterations) or (iterations < 1) and changed:
            tmp_in, tmp_out = tmp_out, tmp_in
            changed = _nd_image.binary_erosion(tmp_in, structure, mask,
                            tmp_out, border_value, origin, invert, cit, 0)
            ii += 1
        if return_value is not None:
            return tmp_out


def binary_erosion(input, structure = None, iterations = 1, mask = None,
        output = None, border_value = 0, origin = 0, brute_force = False):
    """Multi-dimensional binary erosion with the given structure.

    An output array can optionally be provided. The origin parameter
    controls the placement of the filter. If no structuring element is
    provided an element is generated with a squared connectivity equal
    to one. The border_value parameter gives the value of the array
    outside the border. The erosion operation is repeated iterations
    times. If iterations is less than 1, the erosion is repeated until
    the result does not change anymore. If a mask is given, only those
    elements with a true value at the corresponding mask element are
    modified at each iteration.
    """
    return _binary_erosion(input, structure, iterations, mask,
                           output, border_value, origin, 0, brute_force)

def binary_dilation(input, structure = None, iterations = 1, mask = None,
        output = None, border_value = 0, origin = 0, brute_force = False):
    """Multi-dimensional binary dilation with the given structure.

    An output array can optionally be provided. The origin parameter
    controls the placement of the filter. If no structuring element is
    provided an element is generated with a squared connectivity equal
    to one. The dilation operation is repeated iterations times.  If
    iterations is less than 1, the dilation is repeated until the
    result does not change anymore.  If a mask is given, only those
    elements with a true value at the corresponding mask element are
    modified at each iteration.
    """
    input = numpy.asarray(input)
    if structure is None:
        structure = generate_binary_structure(input.ndim, 1)
    origin = _ni_support._normalize_sequence(origin, input.ndim)
    structure = numpy.asarray(structure)
    structure = structure[tuple([slice(None, None, -1)] *
                                structure.ndim)]
    for ii in range(len(origin)):
        origin[ii] = -origin[ii]
        if not structure.shape[ii] & 1:
            origin[ii] -= 1
    return _binary_erosion(input, structure, iterations, mask,
                           output, border_value, origin, 1, brute_force)


def binary_opening(input, structure = None, iterations = 1, output = None,
                   origin = 0):
    """Multi-dimensional binary opening with the given structure.

    An output array can optionally be provided. The origin parameter
    controls the placement of the filter. If no structuring element is
    provided an element is generated with a squared connectivity equal
    to one. The iterations parameter gives the number of times the
    erosions and then the dilations are done.
    """
    input = numpy.asarray(input)
    if structure is None:
        rank = input.ndim
        structure = generate_binary_structure(rank, 1)
    tmp = binary_erosion(input, structure, iterations, None, None, 0,
                         origin)
    return binary_dilation(tmp, structure, iterations, None, output, 0,
                           origin)


def binary_closing(input, structure = None, iterations = 1, output = None,
                   origin = 0):
    """Multi-dimensional binary closing with the given structure.

    An output array can optionally be provided. The origin parameter
    controls the placement of the filter. If no structuring element is
    provided an element is generated with a squared connectivity equal
    to one. The iterations parameter gives the number of times the
    dilations and then the erosions are done.
    """
    input = numpy.asarray(input)
    if structure is None:
        rank = input.ndim
        structure = generate_binary_structure(rank, 1)
    tmp = binary_dilation(input, structure, iterations, None, None, 0,
                          origin)
    return binary_erosion(tmp, structure, iterations, None, output, 0,
                          origin)


def binary_hit_or_miss(input, structure1 = None, structure2 = None,
                       output = None, origin1 = 0, origin2 = None):
    """Multi-dimensional binary hit-or-miss transform.

    An output array can optionally be provided. The origin parameters
    controls the placement of the structuring elements. If the first
    structuring element is not given one is generated with a squared
    connectivity equal to one. If the second structuring element is
    not provided, it set equal to the inverse of the first structuring
    element. If the origin for the second structure is equal to None
    it is set equal to the origin of the first.
    """
    input = numpy.asarray(input)
    if structure1 is None:
        structure1 = generate_binary_structure(input.ndim, 1)
    if structure2 is None:
        structure2 = numpy.logical_not(structure1)
    origin1 = _ni_support._normalize_sequence(origin1, input.ndim)
    if origin2 is None:
        origin2 = origin1
    else:
        origin2 = _ni_support._normalize_sequence(origin2, input.ndim)

    tmp1 = _binary_erosion(input, structure1, 1, None, None, 0, origin1,
                           0, False)
    inplace = isinstance(output, numpy.ndarray)
    result = _binary_erosion(input, structure2, 1, None, output, 0,
                             origin2, 1, False)
    if inplace:
        numpy.logical_not(output, output)
        numpy.logical_and(tmp1, output, output)
    else:
        numpy.logical_not(result, result)
        return numpy.logical_and(tmp1, result)

def binary_propagation(input, structure = None, mask = None,
                       output = None, border_value = 0, origin = 0):
    """Multi-dimensional binary propagation with the given structure.

    An output array can optionally be provided. The origin parameter
    controls the placement of the filter. If no structuring element is
    provided an element is generated with a squared connectivity equal
    to one. If a mask is given, only those elements with a true value at
    the corresponding mask element are.

    This function is functionally equivalent to calling binary_dilation
    with the number of iterations less then one: iterative dilation until
    the result does not change anymore.
    """
    return binary_dilation(input, structure, -1, mask, output,
                           border_value, origin)

def binary_fill_holes(input, structure = None, output = None, origin = 0):
    """Fill the holes in binary objects.

    An output array can optionally be provided. The origin parameter
    controls the placement of the filter. If no structuring element is
    provided an element is generated with a squared connectivity equal
    to one.
    """
    mask = numpy.logical_not(input)
    tmp = numpy.zeros(mask.shape, bool)
    inplace = isinstance(output, numpy.ndarray)
    if inplace:
        binary_dilation(tmp, structure, -1, mask, output, 1, origin)
        numpy.logical_not(output, output)
    else:
        output = binary_dilation(tmp, structure, -1, mask, None, 1,
                                 origin)
        numpy.logical_not(output, output)
        return output

def grey_erosion(input,  size = None, footprint = None, structure = None,
                 output = None, mode = "reflect", cval = 0.0, origin = 0):
    """Calculate a grey values erosion.

    Either a size or a footprint, or the structure must be provided. An
    output array can optionally be provided. The origin parameter
    controls the placement of the filter. The mode parameter
    determines how the array borders are handled, where cval is the
    value when mode is equal to 'constant'.
    """
    return filters._min_or_max_filter(input, size, footprint, structure,
                                      output, mode, cval, origin, 1)


def grey_dilation(input,  size = None, footprint = None, structure = None,
                 output = None, mode = "reflect", cval = 0.0, origin = 0):
    """Calculate a grey values dilation.

    Either a size or a footprint, or the structure must be
    provided. An output array can optionally be provided. The origin
    parameter controls the placement of the filter. The mode parameter
    determines how the array borders are handled, where cval is the
    value when mode is equal to 'constant'.
    """
    if structure is not None:
        structure = numpy.asarray(structure)
        structure = structure[tuple([slice(None, None, -1)] *
                                    structure.ndim)]
    if footprint is not None:
        footprint = numpy.asarray(footprint)
        footprint = footprint[tuple([slice(None, None, -1)] *
                                    footprint.ndim)]
    input = numpy.asarray(input)
    origin = _ni_support._normalize_sequence(origin, input.ndim)
    for ii in range(len(origin)):
        origin[ii] = -origin[ii]
        if footprint is not None:
            sz = footprint.shape[ii]
        else:
            sz = size[ii]
        if not sz & 1:
            origin[ii] -= 1
    return filters._min_or_max_filter(input, size, footprint, structure,
                                      output, mode, cval, origin, 0)


def grey_opening(input, size = None, footprint = None, structure = None,
                 output = None, mode = "reflect", cval = 0.0, origin = 0):
    """Multi-dimensional grey valued opening.

    Either a size or a footprint, or the structure must be provided. An
    output array can optionally be provided. The origin parameter
    controls the placement of the filter. The mode parameter
    determines how the array borders are handled, where cval is the
    value when mode is equal to 'constant'.
    """
    tmp = grey_erosion(input, size, footprint, structure, None, mode,
                       cval, origin)
    return grey_dilation(tmp, size, footprint, structure, output, mode,
                         cval, origin)


def grey_closing(input, size = None, footprint = None, structure = None,
                 output = None, mode = "reflect", cval = 0.0, origin = 0):
    """Multi-dimensional grey valued closing.

    Either a size or a footprint, or the structure must be provided. An
    output array can optionally be provided. The origin parameter
    controls the placement of the filter. The mode parameter
    determines how the array borders are handled, where cval is the
    value when mode is equal to 'constant'.
    """
    tmp = grey_dilation(input, size, footprint, structure, None, mode,
                        cval, origin)
    return grey_erosion(tmp, size, footprint, structure, output, mode,
                        cval, origin)


def morphological_gradient(input, size = None, footprint = None,
                        structure = None, output = None, mode = "reflect",
                        cval = 0.0, origin = 0):
    """Multi-dimensional morphological gradient.

    Either a size or a footprint, or the structure must be provided. An
    output array can optionally be provided. The origin parameter
    controls the placement of the filter. The mode parameter
    determines how the array borders are handled, where cval is the
    value when mode is equal to 'constant'.
    """
    tmp = grey_dilation(input, size, footprint, structure, None, mode,
                        cval, origin)
    if isinstance(output, numpy.ndarray):
        grey_erosion(input, size, footprint, structure, output, mode,
                     cval, origin)
        return numpy.subtract(tmp, output, output)
    else:
        return (tmp - grey_erosion(input, size, footprint, structure,
                                   None, mode, cval, origin))


def morphological_laplace(input, size = None, footprint = None,
                          structure = None, output = None,
                          mode = "reflect", cval = 0.0, origin = 0):
    """Multi-dimensional morphological laplace.

    Either a size or a footprint, or the structure must be provided. An
    output array can optionally be provided. The origin parameter
    controls the placement of the filter. The mode parameter
    determines how the array borders are handled, where cval is the
    value when mode is equal to 'constant'.
    """
    tmp1 = grey_dilation(input, size, footprint, structure, None, mode,
                         cval, origin)
    if isinstance(output, numpy.ndarray):
        grey_erosion(input, size, footprint, structure, output, mode,
                     cval, origin)
        numpy.add(tmp1, output, output)
        del tmp1
        numpy.subtract(output, input, output)
        return numpy.subtract(output, input, output)
    else:
        tmp2 = grey_erosion(input, size, footprint, structure, None, mode,
                            cval, origin)
        numpy.add(tmp1, tmp2, tmp2)
        del tmp1
        numpy.subtract(tmp2, input, tmp2)
        numpy.subtract(tmp2, input, tmp2)
        return tmp2


def white_tophat(input, size = None, footprint = None, structure = None,
                 output = None, mode = "reflect", cval = 0.0, origin = 0):
    """Multi-dimensional white tophat filter.

    Either a size or a footprint, or the structure must be provided. An
    output array can optionally be provided. The origin parameter
    controls the placement of the filter. The mode parameter
    determines how the array borders are handled, where cval is the
    value when mode is equal to 'constant'.
    """
    tmp = grey_erosion(input, size, footprint, structure, None, mode,
                       cval, origin)
    if isinstance(output, numpy.ndarray):
        grey_dilation(tmp, size, footprint, structure, output, mode, cval,
                      origin)
        del tmp
        return numpy.subtract(input, output, output)
    else:
        tmp = grey_dilation(tmp, size, footprint, structure, None, mode,
                            cval, origin)
        return input - tmp


def black_tophat(input, size = None, footprint = None,
                 structure = None, output = None, mode = "reflect",
                 cval = 0.0, origin = 0):
    """Multi-dimensional black tophat filter.

    Either a size or a footprint, or the structure must be provided. An
    output array can optionally be provided. The origin parameter
    controls the placement of the filter. The mode parameter
    determines how the array borders are handled, where cval is the
    value when mode is equal to 'constant'.
    """
    tmp = grey_dilation(input, size, footprint, structure, None, mode,
                        cval, origin)
    if isinstance(output, numpy.ndarray):
        grey_erosion(tmp, size, footprint, structure, output, mode, cval,
                      origin)
        del tmp
        return numpy.subtract(output, input, output)
    else:
        tmp = grey_erosion(tmp, size, footprint, structure, None, mode,
                           cval, origin)
        return tmp - input


def distance_transform_bf(input, metric = "euclidean", sampling = None,
                          return_distances = True, return_indices = False,
                          distances = None, indices = None):
    """Distance transform function by a brute force algorithm.

    This function calculates the distance transform of the input, by
    replacing each background element (zero values), with its
    shortest distance to the foreground (any element non-zero). Three
    types of distance metric are supported: 'euclidean', 'city_block'
    and 'chessboard'.

    In addition to the distance transform, the feature transform can
    be calculated. In this case the index of the closest background
    element is returned along the first axis of the result.

    The return_distances, and return_indices flags can be used to
    indicate if the distance transform, the feature transform, or both
    must be returned.

    Optionally the sampling along each axis can be given by the
    sampling parameter which should be a sequence of length equal to
    the input rank, or a single number in which the sampling is assumed
    to be equal along all axes. This parameter is only used in the
    case of the euclidean distance transform.

    This function employs a slow brute force algorithm, see also the
    function distance_transform_cdt for more efficient city_block and
    chessboard algorithms.

    the distances and indices arguments can be used to give optional
    output arrays that must be of the correct size and type (float64
    and int32).
    """
    if (not return_distances) and (not return_indices):
        msg = 'at least one of distances/indices must be specified'
        raise RuntimeError, msg
    tmp1 = numpy.asarray(input) != 0
    struct = generate_binary_structure(tmp1.ndim, tmp1.ndim)
    tmp2 = binary_dilation(tmp1, struct)
    tmp2 = numpy.logical_xor(tmp1, tmp2)
    tmp1 = tmp1.astype(numpy.int8) - tmp2.astype(numpy.int8)
    del tmp2
    metric = metric.lower()
    if metric == 'euclidean':
        metric = 1
    elif metric == 'cityblock':
        metric = 2
    elif metric == 'chessboard':
        metric = 3
    else:
        raise RuntimeError, 'distance metric not supported'
    if sampling is not None:
        sampling = _ni_support._normalize_sequence(sampling, tmp1.ndim)
        sampling = numpy.asarray(sampling, dtype = numpy.float64)
        if not sampling.flags.contiguous:
            sampling = sampling.copy()
    if return_indices:
        ft = numpy.zeros(tmp1.shape, dtype = numpy.int32)
    else:
        ft = None
    if return_distances:
        if distances is None:
            if metric == 1:
                dt = numpy.zeros(tmp1.shape, dtype = numpy.float64)
            else:
                dt = numpy.zeros(tmp1.shape, dtype = numpy.uint32)
        else:
            if distances.shape != tmp1.shape:
                raise RuntimeError, 'distances array has wrong shape'
            if metric == 1:
                if distances.dtype.type != numpy.float64:
                    raise RuntimeError, 'distances array must be float64'
            else:
                if distances.dtype.type != numpy.uint32:
                    raise RuntimeError, 'distances array must be uint32'
            dt = distances
    else:
        dt = None
    _nd_image.distance_transform_bf(tmp1, metric, sampling, dt, ft)
    if return_indices:
        if isinstance(indices, numpy.ndarray):
            if indices.dtype.type != numpy.int32:
                raise RuntimeError, 'indices must of int32 type'
            if indices.shape != (tmp1.ndim,) + tmp1.shape:
                raise RuntimeError, 'indices has wrong shape'
            tmp2 = indices
        else:
            tmp2 = numpy.indices(tmp1.shape, dtype = numpy.int32)
        ft = numpy.ravel(ft)
        for ii in range(tmp2.shape[0]):
            rtmp = numpy.ravel(tmp2[ii, ...])[ft]
            rtmp.shape = tmp1.shape
            tmp2[ii, ...] = rtmp
        ft = tmp2
    # construct and return the result
    result = []
    if return_distances and not isinstance(distances, numpy.ndarray):
        result.append(dt)
    if return_indices and not isinstance(indices, numpy.ndarray):
        result.append(ft)
    if len(result) == 2:
        return tuple(result)
    elif len(result) == 1:
        return result[0]
    else:
        return None

def distance_transform_cdt(input, structure = 'chessboard',
                        return_distances = True, return_indices = False,
                        distances = None, indices = None):
    """Distance transform for chamfer type of transforms.

    The structure determines the type of chamfering that is done. If
    the structure is equal to 'cityblock' a structure is generated
    using generate_binary_structure with a squared distance equal to
    1. If the structure is equal to 'chessboard', a structure is
    generated using generate_binary_structure with a squared distance
    equal to the rank of the array. These choices correspond to the
    common interpretations of the cityblock and the chessboard
    distance metrics in two dimensions.

    In addition to the distance transform, the feature transform can
    be calculated. In this case the index of the closest background
    element is returned along the first axis of the result.

    The return_distances, and return_indices flags can be used to
    indicate if the distance transform, the feature transform, or both
    must be returned.

    The distances and indices arguments can be used to give optional
    output arrays that must be of the correct size and type (both int32).
    """
    if (not return_distances) and (not return_indices):
        msg = 'at least one of distances/indices must be specified'
        raise RuntimeError, msg
    ft_inplace = isinstance(indices, numpy.ndarray)
    dt_inplace = isinstance(distances, numpy.ndarray)
    input = numpy.asarray(input)
    if structure == 'cityblock':
        rank = input.ndim
        structure = generate_binary_structure(rank, 1)
    elif structure == 'chessboard':
        rank = input.ndim
        structure = generate_binary_structure(rank, rank)
    else:
        try:
            structure = numpy.asarray(structure)
        except:
            raise RuntimeError, 'invalid structure provided'
        for s in structure.shape:
            if s != 3:
                raise RuntimeError, 'structure sizes must be equal to 3'
    if not structure.flags.contiguous:
        structure = structure.copy()
    if dt_inplace:
        if distances.dtype.type != numpy.int32:
            raise RuntimeError, 'distances must be of int32 type'
        if distances.shape != input.shape:
            raise RuntimeError, 'distances has wrong shape'
        dt = distances
        dt[...] = numpy.where(input, -1, 0).astype(numpy.int32)
    else:
        dt = numpy.where(input, -1, 0).astype(numpy.int32)
    rank = dt.ndim
    if return_indices:
        sz = numpy.product(dt.shape,axis=0)
        ft = numpy.arange(sz, dtype = numpy.int32)
        ft.shape = dt.shape
    else:
        ft = None
    _nd_image.distance_transform_op(structure, dt, ft)
    dt = dt[tuple([slice(None, None, -1)] * rank)]
    if return_indices:
        ft = ft[tuple([slice(None, None, -1)] * rank)]
    _nd_image.distance_transform_op(structure, dt, ft)
    dt = dt[tuple([slice(None, None, -1)] * rank)]
    if return_indices:
        ft = ft[tuple([slice(None, None, -1)] * rank)]
        ft = numpy.ravel(ft)
        if ft_inplace:
            if indices.dtype.type != numpy.int32:
                raise RuntimeError, 'indices must of int32 type'
            if indices.shape != (dt.ndim,) + dt.shape:
                raise RuntimeError, 'indices has wrong shape'
            tmp = indices
        else:
            tmp = numpy.indices(dt.shape, dtype = numpy.int32)
        for ii in range(tmp.shape[0]):
            rtmp = numpy.ravel(tmp[ii, ...])[ft]
            rtmp.shape = dt.shape
            tmp[ii, ...] = rtmp
        ft = tmp

    # construct and return the result
    result = []
    if return_distances and not dt_inplace:
        result.append(dt)
    if return_indices and not ft_inplace:
        result.append(ft)
    if len(result) == 2:
        return tuple(result)
    elif len(result) == 1:
        return result[0]
    else:
        return None


def distance_transform_edt(input, sampling = None,
                        return_distances = True, return_indices = False,
                        distances = None, indices = None):
    """Exact euclidean distance transform.

    In addition to the distance transform, the feature transform can
    be calculated. In this case the index of the closest background
    element is returned along the first axis of the result.

    The return_distances, and return_indices flags can be used to
    indicate if the distance transform, the feature transform, or both
    must be returned.

    Optionally the sampling along each axis can be given by the
    sampling parameter which should be a sequence of length equal to
    the input rank, or a single number in which the sampling is assumed
    to be equal along all axes.

    the distances and indices arguments can be used to give optional
    output arrays that must be of the correct size and type (float64
    and int32).
    """
    if (not return_distances) and (not return_indices):
        msg = 'at least one of distances/indices must be specified'
        raise RuntimeError, msg
    ft_inplace = isinstance(indices, numpy.ndarray)
    dt_inplace = isinstance(distances, numpy.ndarray)
    # calculate the feature transform
    input = numpy.where(input, 1, 0).astype(numpy.int8)
    if sampling is not None:
        sampling = _ni_support._normalize_sequence(sampling, input.ndim)
        sampling = numpy.asarray(sampling, dtype = numpy.float64)
        if not sampling.flags.contiguous:
            sampling = sampling.copy()
    if ft_inplace:
        ft = indices
        if ft.shape != (input.ndim,) + input.shape:
            raise RuntimeError, 'indices has wrong shape'
        if ft.dtype.type != numpy.int32:
            raise RuntimeError, 'indices must be of int32 type'
    else:
        ft = numpy.zeros((input.ndim,) + input.shape,
                            dtype = numpy.int32)
    _nd_image.euclidean_feature_transform(input, sampling, ft)
    # if requested, calculate the distance transform
    if return_distances:
        dt = ft - numpy.indices(input.shape, dtype = ft.dtype)
        dt = dt.astype(numpy.float64)
        if sampling is not None:
            for ii in range(len(sampling)):
                dt[ii, ...] *= sampling[ii]
        numpy.multiply(dt, dt, dt)
        if dt_inplace:
            dt = numpy.add.reduce(dt, axis = 0)
            if distances.shape != dt.shape:
                raise RuntimeError, 'indices has wrong shape'
            if distances.dtype.type != numpy.float64:
                raise RuntimeError, 'indices must be of float64 type'
            numpy.sqrt(dt, distances)
            del dt
        else:
            dt = numpy.add.reduce(dt, axis = 0)
            dt = numpy.sqrt(dt)
    # construct and return the result
    result = []
    if return_distances and not dt_inplace:
        result.append(dt)
    if return_indices and not ft_inplace:
        result.append(ft)
    if len(result) == 2:
        return tuple(result)
    elif len(result) == 1:
        return result[0]
    else:
        return None
