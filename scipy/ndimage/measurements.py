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
import numpy
import _ni_support
import _nd_image
import morphology

def label(input, structure = None, output = None):
    """Label an array of objects.

    The structure that defines the object connections must be
    symmetric.  If no structuring element is provided an element is
    generated with a squared connectivity equal to one. This function
    returns a tuple consisting of the array of labels and the number
    of objects found. If an output array is provided only the number of
    objects found is returned.
    """
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    if structure is None:
        structure = morphology.generate_binary_structure(input.ndim, 1)
    structure = numpy.asarray(structure, dtype = bool)
    if structure.ndim != input.ndim:
        raise RuntimeError, 'structure and input must have equal rank'
    for ii in structure.shape:
        if ii != 3:
            raise  RuntimeError, 'structure dimensions must be equal to 3'
    if not structure.flags.contiguous:
        structure = structure.copy()
    if isinstance(output, numpy.ndarray):
        if output.dtype.type != numpy.int32:
            raise RuntimeError, 'output type must be int32'
    else:
        output = numpy.int32
    output, return_value = _ni_support._get_output(output, input)
    max_label = _nd_image.label(input, structure, output)
    if return_value is None:
        return max_label
    else:
        return return_value, max_label

def find_objects(input, max_label = 0):
    """Find objects in a labeled array.

    The input must be an array with labeled objects. A list of slices
    into the array is returned that contain the objects. The list
    represents a sequence of the numbered objects. If a number is
    missing, None is returned instead of a slice. If max_label > 0, it
    gives the largest object number that is searched for, otherwise
    all are returned.
    """
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    if max_label < 1:
        max_label = input.max()
    return _nd_image.find_objects(input, max_label)

def sum(input, labels=None, index=None):
    """Calculate the sum of the values of the array.

    :Parameters:
        labels : array of integers, same shape as input
            Assign labels to the values of the array.

        index : scalar or array
            A single label number or a sequence of label numbers of
            the objects to be measured. If index is None, all
            values are used where 'labels' is larger than zero.

    Examples
    --------

    >>> input =  [0,1,2,3]
    >>> labels = [1,1,2,2]
    >>> sum(input, labels, index=[1,2])
    [1.0, 5.0]

    """
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    if labels is not None:
        labels = numpy.asarray(labels)
        labels = _broadcast(labels, input.shape)

        if labels.shape != input.shape:
            raise RuntimeError, 'input and labels shape are not equal'
    if index is not None:
        T = getattr(index,'dtype',numpy.int32)
        if T not in [numpy.int8, numpy.int16, numpy.int32,
                     numpy.uint8, numpy.uint16, numpy.bool]:
            raise ValueError("Invalid index type")
        index = numpy.asarray(index,dtype=T)
    return _nd_image.statistics(input, labels, index, 0)


def mean(input, labels = None, index = None):
    """Calculate the mean of the values of the array.

    The index parameter is a single label number or a sequence of
    label numbers of the objects to be measured. If index is None, all
    values are used where labels is larger than zero.
    """
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    if labels is not None:
        labels = numpy.asarray(labels)
        labels = _broadcast(labels, input.shape)

        if labels.shape != input.shape:
            raise RuntimeError, 'input and labels shape are not equal'
    return _nd_image.statistics(input, labels, index, 1)


def variance(input, labels = None, index = None):
    """Calculate the variance of the values of the array.

    The index parameter is a single label number or a sequence of
    label numbers of the objects to be measured. If index is None, all
    values are used where labels is larger than zero.
    """
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    if labels is not None:
        labels = numpy.asarray(labels)
        labels = _broadcast(labels, input.shape)

        if labels.shape != input.shape:
            raise RuntimeError, 'input and labels shape are not equal'
    return _nd_image.statistics(input, labels, index, 2)


def standard_deviation(input, labels = None, index = None):
    """Calculate the standard deviation of the values of the array.

    The index parameter is a single label number or a sequence of
    label numbers of the objects to be measured. If index is None, all
    values are used where labels is larger than zero.
    """
    var = variance(input, labels, index)
    if (isinstance(var, types.ListType)):
        return [math.sqrt(x) for x in var]
    else:
        return math.sqrt(var)


def minimum(input, labels = None, index = None):
    """Calculate the minimum of the values of the array.

    The index parameter is a single label number or a sequence of
    label numbers of the objects to be measured. If index is None, all
    values are used where labels is larger than zero.
    """
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    if labels is not None:
        labels = numpy.asarray(labels)
        labels = _broadcast(labels, input.shape)

        if labels.shape != input.shape:
            raise RuntimeError, 'input and labels shape are not equal'
    return _nd_image.statistics(input, labels, index, 3)


def maximum(input, labels=None, index=None):
    """Return the maximum input value.

    The index parameter is a single label number or a sequence of
    label numbers of the objects to be measured. If index is None, all
    values are used where labels is larger than zero.

    """
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    if labels is not None:
        labels = numpy.asarray(labels)
        labels = _broadcast(labels, input.shape)

        if labels.shape != input.shape:
            raise RuntimeError, 'input and labels shape are not equal'
    return _nd_image.statistics(input, labels, index, 4)


def _index_to_position(index, shape):
    """Convert a linear index to a position"""
    if len(shape) > 0:
        pos = []
        stride = numpy.multiply.reduce(shape)
        for size in shape:
            stride = stride // size
            pos.append(index // stride)
            index -= pos[-1] * stride
        return tuple(pos)
    else:
        return 0


def minimum_position(input, labels = None, index = None):
    """Find the position of the minimum of the values of the array.

    The index parameter is a single label number or a sequence of
    label numbers of the objects to be measured. If index is None, all
    values are used where labels is larger than zero.
    """
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    if labels is not None:
        labels = numpy.asarray(labels)
        labels = _broadcast(labels, input.shape)

        if labels.shape != input.shape:
            raise RuntimeError, 'input and labels shape are not equal'
    pos = _nd_image.statistics(input, labels, index, 5)
    if (isinstance(pos, types.ListType)):
        return [_index_to_position(x, input.shape) for x in pos]
    else:
        return _index_to_position(pos, input.shape)


def maximum_position(input, labels = None, index = None):
    """Find the position of the maximum of the values of the array.

    The index parameter is a single label number or a sequence of
    label numbers of the objects to be measured. If index is None, all
    values are used where labels is larger than zero.
    """
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    if labels is not None:
        labels = numpy.asarray(labels)
        labels = _broadcast(labels, input.shape)

        if labels.shape != input.shape:
            raise RuntimeError, 'input and labels shape are not equal'
    pos = _nd_image.statistics(input, labels, index, 6)
    if (isinstance(pos, types.ListType)):
        return [_index_to_position(x, input.shape) for x in pos]
    else:
        return _index_to_position(pos, input.shape)


def extrema(input, labels = None, index = None):
    """Calculate the minimum, the maximum and their positions of the
       values of the array.

    The index parameter is a single label number or a sequence of
    label numbers of the objects to be measured. If index is None, all
    values are used where labels is larger than zero.
    """
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    if labels is not None:
        labels = numpy.asarray(labels)
        labels = _broadcast(labels, input.shape)

        if labels.shape != input.shape:
            raise RuntimeError, 'input and labels shape are not equal'


    min, max, minp, maxp = _nd_image.statistics(input, labels, index, 7)
    if (isinstance(minp, types.ListType)):
        minp = [_index_to_position(x, input.shape) for x in minp]
        maxp = [_index_to_position(x, input.shape) for x in maxp]
    else:
        minp = _index_to_position(minp, input.shape)
        maxp = _index_to_position(maxp, input.shape)
    return min, max, minp, maxp


def center_of_mass(input, labels = None, index = None):
    """Calculate the center of mass of of the array.

    The index parameter is a single label number or a sequence of
    label numbers of the objects to be measured. If index is None, all
    values are used where labels is larger than zero.
    """
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    if labels is not None:
        labels = numpy.asarray(labels)
        labels = _broadcast(labels, input.shape)

        if labels.shape != input.shape:
            raise RuntimeError, 'input and labels shape are not equal'
    return _nd_image.center_of_mass(input, labels, index)


def histogram(input, min, max, bins, labels = None, index = None):
    """Calculate a histogram of of the array.

    The histogram is defined by its minimum and maximum value and the
    number of bins.

    The index parameter is a single label number or a sequence of
    label numbers of the objects to be measured. If index is None, all
    values are used where labels is larger than zero.
    """
    input = numpy.asarray(input)
    if numpy.iscomplexobj(input):
        raise TypeError, 'Complex type not supported'
    if labels is not None:
        labels = numpy.asarray(labels)
        labels = _broadcast(labels, input.shape)

        if labels.shape != input.shape:
            raise RuntimeError, 'input and labels shape are not equal'
    if bins < 1:
        raise RuntimeError, 'number of bins must be >= 1'
    if min >= max:
        raise RuntimeError, 'min must be < max'
    return _nd_image.histogram(input, min, max, bins, labels, index)

def watershed_ift(input, markers, structure = None, output = None):
    """Apply watershed from markers using a iterative forest transform
    algorithm.

    Negative markers are considered background markers which are
    processed after the other markers. A structuring element defining
    the connectivity of the object can be provided. If none is
    provided an element is generated iwth a squared connecitiviy equal
    to one. An output array can optionally be provided.
    """
    input = numpy.asarray(input)
    if input.dtype.type not in [numpy.uint8, numpy.uint16]:
        raise TypeError, 'only 8 and 16 unsigned inputs are supported'
    if structure is None:
        structure = morphology.generate_binary_structure(input.ndim, 1)
    structure = numpy.asarray(structure, dtype = bool)
    if structure.ndim != input.ndim:
        raise RuntimeError, 'structure and input must have equal rank'
    for ii in structure.shape:
        if ii != 3:
            raise  RuntimeError, 'structure dimensions must be equal to 3'
    if not structure.flags.contiguous:
        structure = structure.copy()
    markers = numpy.asarray(markers)
    if input.shape != markers.shape:
        raise RuntimeError, 'input and markers must have equal shape'

    integral_types = [numpy.int0,
                      numpy.int8,
                      numpy.int16,
                      numpy.int32,
                      numpy.int_,
                      numpy.int64,
                      numpy.intc,
                      numpy.intp]

    if markers.dtype.type not in integral_types:
        raise RuntimeError, 'marker should be of integer type'
    if isinstance(output, numpy.ndarray):
        if output.dtype.type not in integral_types:
            raise RuntimeError, 'output should be of integer type'
    else:
        output = markers.dtype
    output, return_value = _ni_support._get_output(output, input)
    _nd_image.watershed_ift(input, markers, structure, output)
    return return_value

def _broadcast(arr, sshape):
    """Return broadcast view of arr, else return None."""
    ashape = arr.shape
    return_value = numpy.zeros(sshape, arr.dtype)
    # Just return arr if they have the same shape
    if sshape == ashape:
        return arr
    srank = len(sshape)
    arank = len(ashape)

    aslices = []
    sslices = []
    for i in range(arank):
        aslices.append(slice(0, ashape[i], 1))

    for i in range(srank):
        sslices.append(slice(0, sshape[i], 1))
    return_value[sslices] = arr[aslices]
    return return_value
