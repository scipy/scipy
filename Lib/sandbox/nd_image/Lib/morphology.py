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

import numarray
import _ni_support
import _nd_image
import filters

def _center_is_true(structure, origin):
    structure = numarray.array(structure)
    coor = tuple([oo + ss // 2 for ss, oo in zip(structure.shape,
                                                 origin)])
    return bool(structure[coor])

def iterate_structure(structure, iterations, origin = None):
    """Iterate a structure by dilating it with itself.

    If origin is None, only the iterated structure is returned. If
    not, a tuple of the iterated structure and the modified origin is
    returned.
    """
    structure = numarray.asarray(structure)
    if iterations < 2:
        return structure.copy()
    ni = iterations - 1
    shape = [ii + ni * (ii - 1) for ii in structure.shape]
    pos = [ni * (structure.shape[ii] / 2) for ii in range(len(shape))]
    slc = [slice(pos[ii], pos[ii] + structure.shape[ii], None)
           for ii in range(len(shape))]
    out = numarray.zeros(shape, numarray.Bool)
    out[slc] = structure != 0
    out = binary_dilation(out, structure, iterations = ni)
    if origin is None:
        return out
    else:
        origin = _ni_support._normalize_sequence(origin, structure.rank)
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
            return numarray.array(0, type = numarray.Bool)
        else:
            return numarray.array(1, type = numarray.Bool)
    output = numarray.zeros([3] * rank, numarray.Bool)
    output = numarray.abs(numarray.indices([3] * rank) - 1)
    output = numarray.add.reduce(output, 0)
    return numarray.asarray(output <= connectivity, type = numarray.Bool)


def _binary_erosion(input, structure, iterations, mask, output,
                    border_value, origin, invert, brute_force):
    input = numarray.asarray(input)
    if isinstance(input.type(), numarray.ComplexType):
        raise TypeError, 'Complex type not supported'
    if structure is None:
        structure = generate_binary_structure(input.rank, 1)
    else:
        structure = numarray.asarray(structure, type = numarray.Bool)
    if structure.rank != input.rank:
        raise RuntimeError, 'structure rank must equal input rank'
    if not structure.iscontiguous():
        structure = structure.copy()
    if structure.nelements() < 1:
        raise RuntimeError, 'structure must not be empty'
    if mask is not None:
        mask = numarray.asarray(mask)
        if mask.shape != input.shape:
            raise RuntimeError, 'mask and input must have equal sizes'
    origin = _ni_support._normalize_sequence(origin, input.rank)
    cit = _center_is_true(structure, origin)
    if isinstance(output, numarray.NumArray):
        if isinstance(output.type(), numarray.ComplexType):
            raise TypeError, 'Complex output type not supported'
    else:
        output = numarray.Bool
    output, return_value = _ni_support._get_output(output, input)
    if iterations == 1:
        _nd_image.binary_erosion(input, structure, mask, output,
                                     border_value, origin, invert, cit, 0)
        return return_value
    elif cit and not brute_force:
        changed, coordinate_list = _nd_image.binary_erosion(input,
             structure, mask, output, border_value, origin, invert, cit, 1)
        structure = structure[tuple([slice(None, None, -1)] *
                                    structure.rank)]
        for ii in range(len(origin)):
            origin[ii] = -origin[ii]
            if not structure.shape[ii] & 1:
                origin[ii] -= 1
        if mask != None:
            msk = numarray.asarray(mask, type = numarray.Int8)
            if msk is mask:
                msk = mask.copy()
            mask = msk
        if not structure.iscontiguous():
            structure = structure.copy()
        _nd_image.binary_erosion2(output, structure, mask, iterations - 1,
                                  origin, invert, coordinate_list)
        return return_value
    else:
        tmp_in = numarray.zeros(input.shape, numarray.Bool)
        if return_value == None:
            tmp_out = output
        else:
            tmp_out = numarray.zeros(input.shape, numarray.Bool)
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
        if return_value != None:
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
    input = numarray.asarray(input)
    if structure == None:
        structure = generate_binary_structure(input.rank, 1)
    origin = _ni_support._normalize_sequence(origin, input.rank)
    structure = numarray.asarray(structure)
    structure = structure[tuple([slice(None, None, -1)] *
                                structure.rank)]
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
    input = numarray.asarray(input)
    if structure is None:
        rank = input.rank
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
    input = numarray.asarray(input)
    if structure is None:
        rank = input.rank
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
    input = numarray.asarray(input)
    if structure1 is None:
        structure1 = generate_binary_structure(input.rank, 1)
    if structure2 is None:
        structure2 = numarray.logical_not(structure1)
    origin1 = _ni_support._normalize_sequence(origin1, input.rank)
    if origin2 is None:
        origin2 = origin1
    else:
        origin2 = _ni_support._normalize_sequence(origin2, input.rank)

    tmp1 = _binary_erosion(input, structure1, 1, None, None, 0, origin1,
                           0, False)
    inplace = isinstance(output, numarray.NumArray)
    result = _binary_erosion(input, structure2, 1, None, output, 0,
                             origin2, 1, False)
    if inplace:
        numarray.logical_not(output, output)
        numarray.logical_and(tmp1, output, output)
    else:
        numarray.logical_not(result, result)
        return numarray.logical_and(tmp1, result)
        
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
    mask = numarray.logical_not(input)
    tmp = numarray.zeros(mask.shape, numarray.Bool)
    inplace = isinstance(output, numarray.NumArray)
    if inplace:
        binary_dilation(tmp, structure, -1, mask, output, 1, origin)
        numarray.logical_not(output, output)
    else:
        output = binary_dilation(tmp, structure, -1, mask, None, 1,
                                 origin)
        numarray.logical_not(output, output)
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
        structure = numarray.asarray(structure)
        structure = structure[tuple([slice(None, None, -1)] *
                                    structure.rank)]
    if footprint is not None:
        footprint = numarray.asarray(footprint)
        footprint = footprint[tuple([slice(None, None, -1)] *
                                    footprint.rank)]
    input = numarray.asarray(input)
    origin = _ni_support._normalize_sequence(origin, input.rank)
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
    if isinstance(output, numarray.NumArray):
        grey_erosion(input, size, footprint, structure, output, mode,
                     cval, origin)
        return numarray.subtract(tmp, output, output)
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
    if isinstance(output, numarray.NumArray):
        grey_erosion(input, size, footprint, structure, output, mode,
                     cval, origin)
        numarray.add(tmp1, output, output)
        del tmp1
        numarray.subtract(output, input, output)
        return numarray.subtract(output, input, output)
    else:
        tmp2 = grey_erosion(input, size, footprint, structure, None, mode,
                            cval, origin)
        numarray.add(tmp1, tmp2, tmp2)
        del tmp1
        numarray.subtract(tmp2, input, tmp2)
        numarray.subtract(tmp2, input, tmp2)
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
    if isinstance(output, numarray.NumArray):
        grey_dilation(tmp, size, footprint, structure, output, mode, cval,
                      origin)
        del tmp
        return numarray.subtract(input, output, output)
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
    if isinstance(output, numarray.NumArray):
        grey_erosion(tmp, size, footprint, structure, output, mode, cval,
                      origin)
        del tmp
        return numarray.subtract(output, input, output)
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
    output arrays that must be of the correct size and type (Float64
    and Int32).
    """
    if (not return_distances) and (not return_indices):
        msg = 'at least one of distances/indices must be specified'
        raise RuntimeError, msg
    tmp1 = numarray.asarray(input) != 0
    struct = generate_binary_structure(tmp1.rank, tmp1.rank)
    tmp2 = binary_dilation(tmp1, struct)
    tmp2 = numarray.logical_xor(tmp1, tmp2)
    tmp1 = tmp1.astype(numarray.Int8) - tmp2.astype(numarray.Int8)
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
    if sampling != None:
        sampling = _ni_support._normalize_sequence(sampling, tmp1.rank)
        sampling = numarray.asarray(sampling, type = numarray.Float64)
        if not sampling.iscontiguous():
            sampling = sampling.copy()
    if return_indices:
        ft = numarray.zeros(tmp1.shape, type = numarray.Int32)
    else:
        ft = None
    if return_distances:
        if distances == None:
            if metric == 1:
                dt = numarray.zeros(tmp1.shape, type = numarray.Float64)
            else:
                dt = numarray.zeros(tmp1.shape, type = numarray.UInt32)
        else:
            if distances.shape != tmp1.shape:
                raise RuntimeError, 'distances array has wrong shape'
            if metric == 1:
                if distances.type() != numarray.Float64:
                    raise RuntimeError, 'distances array must be Float64'
            else:
                if distances.type() != numarray.UInt32:
                    raise RuntimeError, 'distances array must be UInt32'
            dt = distances
    else:
        dt = None
    _nd_image.distance_transform_bf(tmp1, metric, sampling, dt, ft)
    if return_indices:
        if isinstance(indices, numarray.NumArray):
            if indices.type() != numarray.Int32:
                raise RuntimeError, 'indices must of Int32 type'
            if indices.shape != (tmp1.rank,) + tmp1.shape:
                raise RuntimeError, 'indices has wrong shape'
            tmp2 = indices
        else:
            tmp2 = numarray.indices(tmp1.shape, type = numarray.Int32)
        ft = numarray.ravel(ft)
        for ii in range(tmp2.shape[0]):
            rtmp = numarray.ravel(tmp2[ii, ...])[ft]
            rtmp.setshape(tmp1.shape)
            tmp2[ii, ...] = rtmp
        ft = tmp2
    # construct and return the result
    result = []
    if return_distances and not isinstance(distances, numarray.NumArray):
        result.append(dt)
    if return_indices and not isinstance(indices, numarray.NumArray):
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
    output arrays that must be of the correct size and type (both Int32).
    """
    if (not return_distances) and (not return_indices):
        msg = 'at least one of distances/indices must be specified'
        raise RuntimeError, msg    
    ft_inplace = isinstance(indices, numarray.NumArray)
    dt_inplace = isinstance(distances, numarray.NumArray)
    input = numarray.asarray(input)
    if structure == 'cityblock':
        rank = input.rank
        structure = generate_binary_structure(rank, 1)
    elif structure == 'chessboard':
        rank = input.rank
        structure = generate_binary_structure(rank, rank)
    else:
        try:
            structure = numarray.asarray(structure)
        except:
            raise RuntimeError, 'invalid structure provided'
        for s in structure.shape:
            if s != 3:
                raise RuntimeError, 'structure sizes must be equal to 3'
    if not structure.iscontiguous():
        structure = structure.copy()
    if dt_inplace:
        if distances.type() != numarray.Int32:
            raise RuntimeError, 'distances must be of Int32 type'    
        if distances.shape != input.shape:
            raise RuntimeError, 'distances has wrong shape'
        dt = distances
        dt[...] = numarray.where(input, -1, 0).astype(numarray.Int32)
    else:
        dt = numarray.where(input, -1, 0).astype(numarray.Int32)
    rank = dt.rank
    if return_indices:
        sz = dt.nelements()
        ft = numarray.arange(sz, shape=dt.shape, type = numarray.Int32)
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
        ft = numarray.ravel(ft)
        if ft_inplace:
            if indices.type() != numarray.Int32:
                raise RuntimeError, 'indices must of Int32 type'
            if indices.shape != (dt.rank,) + dt.shape:
                raise RuntimeError, 'indices has wrong shape'
            tmp = indices
        else:
            tmp = numarray.indices(dt.shape, type = numarray.Int32)
        for ii in range(tmp.shape[0]):
            rtmp = numarray.ravel(tmp[ii, ...])[ft]
            rtmp.setshape(dt.shape)
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
    output arrays that must be of the correct size and type (Float64
    and Int32).
    """
    if (not return_distances) and (not return_indices):
        msg = 'at least one of distances/indices must be specified'
        raise RuntimeError, msg
    ft_inplace = isinstance(indices, numarray.NumArray)
    dt_inplace = isinstance(distances, numarray.NumArray)
    # calculate the feature transform
    input = numarray.where(input, 1, 0).astype(numarray.Int8)
    if sampling is not None:
        sampling = _ni_support._normalize_sequence(sampling, input.rank)
        sampling = numarray.asarray(sampling, type = numarray.Float64)
        if not sampling.iscontiguous():
            sampling = sampling.copy()
    if ft_inplace:
        ft = indices
        if ft.shape != (input.rank,) + input.shape:
            raise RuntimeError, 'indices has wrong shape'
        if ft.type() != numarray.Int32:
            raise RuntimeError, 'indices must be of Int32 type'
    else:
        ft = numarray.zeros((input.rank,) + input.shape,
                            type = numarray.Int32) 
    _nd_image.euclidean_feature_transform(input, sampling, ft)
    # if requested, calculate the distance transform
    if return_distances:
        dt = ft - numarray.indices(input.shape, type = ft.type())
        dt = dt.astype(numarray.Float64)
        if sampling is not None:
            for ii in range(len(sampling)):
                dt[ii, ...] *= sampling[ii]
        numarray.multiply(dt, dt, dt)
        if dt_inplace:
            dt = numarray.add.reduce(dt, axis = 0)
            if distances.shape != dt.shape:
                raise RuntimeError, 'indices has wrong shape'
            if distances.type() != numarray.Float64:
                raise RuntimeError, 'indices must be of Float64 type'
            numarray.sqrt(dt, distances)
            del dt
        else:
            dt = numarray.add.reduce(dt, axis = 0)
            dt = numarray.sqrt(dt)
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
