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
import numpy as np
import _ni_support
import _nd_image
import morphology
import time

def label(input, structure = None, output = None):
    """
    Label features in an array.

    Parameters
    ----------

    input : array_like
        An array-like object to be labeled.  Any non-zero values in `input` are
        counted as features and zero values are considered the background.


    structure : array_like, optional
        A structuring element that defines feature connections.

        `structure` must be symmetric.  If no structuring element is provided,
        one is automatically generated with a squared connectivity equal to
        one.

        That is, for a 2D `input` array, the default structuring element is::

            [[0,1,0],
             [1,1,1],
             [0,1,0]]


    output : (None, data-type, array_like), optional
        If `output` is a data type, it specifies the type of the resulting
        labeled feature array

        If `output` is an array-like object, then `output` will be updated
        with the labeled features from this function

    Returns
    -------
    labeled_array : array_like
        An array-like object where each unique feature has a unique value

    num_features : int


    If `output` is None or a data type, this function returns a tuple,
    (`labeled_array`, `num_features`).

    If `output` is an array, then it will be updated with values in
    `labeled_array` and only `num_features` will be returned by this function.


    See Also
    --------
    find_objects : generate a list of slices for the labeled features (or
                   objects); useful for finding features' position or
                   dimensions

    Examples
    --------

    Create an image with some features, then label it using the default
    (cross-shaped) structuring element:

    >>> a = array([[0,0,1,1,0,0],
    ...            [0,0,0,1,0,0],
    ...            [1,1,0,0,1,0],
    ...            [0,0,0,1,0,0]])
    >>> labeled_array, num_features = label(a)

    Each of the 4 features are labeled with a different integer:

    >>> print num_features
    4
    >>> print labeled_array
    array([[0, 0, 1, 1, 0, 0],
           [0, 0, 0, 1, 0, 0],
           [2, 2, 0, 0, 3, 0],
           [0, 0, 0, 4, 0, 0]])

    Generate a structuring element that will consider features connected even
    if they touch diagonally:

    >>> s = generate_binary_structure(2,2)

    or,

    >>> s = [[1,1,1],
             [1,1,1],
             [1,1,1]]

    Label the image using the new structuring element:

    >>> labeled_array, num_features = label(a, structure=s)

    Show the 2 labeled features (note that features 1, 3, and 4 from above are
    now considered a single feature):

    >>> print num_features
    2
    >>> print labeled_array
    array([[0, 0, 1, 1, 0, 0],
           [0, 0, 0, 1, 0, 0],
           [2, 2, 0, 0, 1, 0],
           [0, 0, 0, 1, 0, 0]])

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

def labeled_comprehension(input, labels, index, func, out_dtype, default, pass_positions=False):
    '''Roughly equivalent to [func(input[labels == i]) for i in index].

    Special cases:
      - index a scalar: returns a single value
      - index is None: returns func(inputs[labels > 0])

    func will be called with linear indices as a second argument if
    pass_positions is True.
    '''
    
    as_scalar = numpy.isscalar(index)
    input = numpy.asarray(input)

    if pass_positions:
        positions = numpy.arange(input.size).reshape(input.shape)

    if labels is None:
        if index is not None:
            raise ValueError, "index without defined labels"
        if not pass_positions:
            return func(input.ravel())
        else:
            return func(input.ravel(), positions.ravel())

    try:
        input, labels = numpy.broadcast_arrays(input, labels)
    except ValueError:
        raise ValueError, "input and labels must have the same shape (excepting dimensions with width 1)"

    if index is None:
        if not pass_positions:
            return func(input[labels > 0])
        else:
            return func(input[labels > 0], positions[labels > 0])

    index = numpy.atleast_1d(index)
    if np.any(index.astype(labels.dtype).astype(index.dtype) != index):
        raise ValueError, "Cannot convert index values from <%s> to <%s> (labels' type) without loss of precision"%(index.dtype, labels.dtype)
    index = index.astype(labels.dtype)

    # optimization: find min/max in index, and select those parts of labels, input, and positions
    lo = index.min()
    hi = index.max()
    mask = (labels >= lo) & (labels <= hi)

    # this also ravels the arrays
    labels = labels[mask]
    input = input[mask]
    if pass_positions:
        positions = positions[mask]

    # sort everything by labels
    label_order = labels.argsort()
    labels = labels[label_order]
    input = input[label_order]
    if pass_positions:
        positions = positions[label_order]
    
    index_order = index.argsort()
    sorted_index = index[index_order]

    def do_map(inputs, output):
        '''labels must be sorted'''
        
        nlabels = labels.size
        nidx = sorted_index.size

        # Find boundaries for each stretch of constant labels
        # This could be faster, but we already paid N log N to sort labels.
        lo = numpy.searchsorted(labels, sorted_index, side='left')
        hi = numpy.searchsorted(labels, sorted_index, side='right')
    
        for i, l, h in zip(range(nidx), lo, hi):
            if l == h:
                continue
            idx = sorted_index[i]
            output[i] = func(*[inp[l:h] for inp in inputs])
            
    temp = numpy.empty(index.shape, out_dtype)
    temp[:] = default
    if not pass_positions:
        do_map([input], temp)
    else:
        do_map([input, positions], temp)
    output = numpy.zeros(index.shape, out_dtype)
    output[index_order] = temp

    if as_scalar:
        output = output[0]

    return output

def _stats(input, labels = None, index = None, do_sum2=False):
    '''returns count, sum, and optionally sum^2 by label'''

    def single_group(vals):
        if do_sum2:
            return vals.size, vals.sum(), (vals * vals.conjugate()).sum()
        else:
            return vals.size, vals.sum()
        
    if labels is None:
        return single_group(input)

    # ensure input and labels match sizes
    input, labels = numpy.broadcast_arrays(input, labels)

    if index is None:
        return single_group(input[labels > 0])

    if numpy.isscalar(index):
        return single_group(input[labels == index])

    # remap labels to unique integers if necessary, or if the largest
    # label is larger than the number of values.
    if ((not numpy.issubdtype(labels.dtype, numpy.int)) or 
        (labels.min() < 0) or (labels.max() > labels.size)):
        unique_labels, new_labels = numpy.unique1d(labels, return_inverse=True)

        counts = numpy.bincount(new_labels)
        sums = numpy.bincount(new_labels, weights=input.ravel())
        if do_sum2:
            sums2 = numpy.bincount(new_labels, weights=(input * input.conjugate()).ravel())

        idxs = numpy.searchsorted(unique_labels, index)
        # make all of idxs valid
        idxs[idxs >= unique_labels.size] = 0
        found = (unique_labels[idxs] == index)
    else:
        # labels are an integer type, and there aren't too many, so
        # call bincount directly.
        counts = numpy.bincount(labels.ravel())
        sums = numpy.bincount(labels.ravel(), weights=input.ravel())
        if do_sum2:
            sums2 = numpy.bincount(labels.ravel(), weights=(input * input.conjugate()).ravel())

        # make sure all index values are valid
        idxs = numpy.asanyarray(index, numpy.int).copy()
        found = (idxs >= 0) & (idxs < counts.size)
        idxs[~ found] = 0

    counts = counts[idxs]
    counts[~ found] = 0
    sums = sums[idxs]
    sums[~ found] = 0
    if not do_sum2:
        return (counts, sums)    
    sums2 = sums2[idxs]
    sums2[~ found] = 0
    return (counts, sums, sums2)

        
def sum(input, labels = None, index = None):
    """
    Calculate the sum of the values of the array.

    Parameters
    ----------

    input : array_like
        Values of `input` inside the regions defined by `labels`
        are summed together.

    labels : array of integers, same shape as input
        Assign labels to the values of the array.

    index : scalar or array
        A single label number or a sequence of label numbers of
        the objects to be measured.

    Returns
    -------

    output : list
        A list of the sums of the values of `input` inside the regions
        defined by `labels`.

    See also
    --------

    mean

    Examples
    --------

    >>> input =  [0,1,2,3]
    >>> labels = [1,1,2,2]
    >>> sum(input, labels, index=[1,2])
    [1.0, 5.0]

    """
    count, sum = _stats(input, labels, index)
    return sum

def mean(input, labels = None, index = None):
    """Calculate the mean of the values of an array at labels.

    Labels must be None or an array that can be broadcast to the input.

    Index must be None, a single label or sequence of labels.  If
    None, the mean for all values where label is greater than 0 is
    calculated.
    """

    count, sum = _stats(input, labels, index)
    return sum / numpy.asanyarray(count).astype(numpy.float)

def variance(input, labels = None, index = None):
    """Calculate the variance of the values of an array at labels.

    Labels must be None or an array of the same dimensions as the input.  

    Index must be None, a single label or sequence of labels.  If
    none, all values where label is greater than zero are used.
    """

    count, sum, sum2 = _stats(input, labels, index, do_sum2=True)
    mean = sum / numpy.asanyarray(count).astype(numpy.float)
    mean2 = sum2 / numpy.asanyarray(count).astype(numpy.float)

    return mean2 - (mean * mean.conjugate())

def standard_deviation(input, labels = None, index = None):
    """Calculate the standard deviation of the values of an array at labels.

    Labels must be None or an array of the same dimensions as the input.  

    Index must be None, a single label or sequence of labels.  If
    none, all values where label is greater than zero are used.
    """

    return numpy.sqrt(variance(input, labels, index))

def _select(input, labels = None, index = None, find_min=False, find_max=False, find_min_positions=False, find_max_positions=False):
    '''returns min, max, or both, plus positions if requested'''


    find_positions = find_min_positions or find_max_positions
    positions = None
    if find_positions:
        positions = numpy.arange(input.size).reshape(input.shape)

    def single_group(vals, positions):
        result = []
        if find_min:
            result += [vals.min()]
        if find_min_positions:
            result += [positions[vals == vals.min()][0]]
        if find_max:
            result += [vals.max()]
        if find_max_positions:
            result += [positions[vals == vals.max()][0]]
        return result
        
    if labels is None:
        return single_group(input, positions)

    # ensure input and labels match sizes
    input, labels = numpy.broadcast_arrays(input, labels)

    if index is None:
        mask = (labels > 0)
        masked_positions = None
        if find_positions:
            masked_positions = positions[mask]
        return single_group(input[mask], masked_positions)

    if numpy.isscalar(index):
        mask = (labels == index)
        masked_positions = None
        if find_positions:
            masked_positions = positions[mask]
        return single_group(input[mask], masked_positions)

    order = input.ravel().argsort()
    input = input.ravel()[order]
    labels = labels.ravel()[order]
    if find_positions:
        positions = positions.ravel()[order]

    # remap labels to unique integers if necessary, or if the largest
    # label is larger than the number of values.
    if ((not numpy.issubdtype(labels.dtype, numpy.int)) or 
        (labels.min() < 0) or (labels.max() > labels.size)):
        # remap labels, and indexes
        unique_labels, labels = numpy.unique1d(labels, return_inverse=True)
        idxs = numpy.searchsorted(unique_labels, index)

        # make all of idxs valid
        idxs[idxs >= unique_labels.size] = 0
        found = (unique_labels[idxs] == index)
    else:
        # labels are an integer type, and there aren't too many.
        idxs = numpy.asanyarray(index, numpy.int).copy()
        found = (idxs >= 0) & (idxs <= labels.max())
    
    idxs[~ found] = labels.max() + 1

    result = []
    if find_min:
        mins = numpy.zeros(labels.max() + 2, input.dtype)
        mins[labels[::-1]] = input[::-1]
        result += [mins[idxs]]
    if find_min_positions:
        minpos = numpy.zeros(labels.max() + 2)
        minpos[labels[::-1]] = positions[::-1]
        result += [minpos[idxs]]
    if find_max:
        maxs = numpy.zeros(labels.max() + 2, input.dtype)
        maxs[labels] = input
        result += [maxs[idxs]]
    if find_max_positions:
        maxpos = numpy.zeros(labels.max() + 2)
        maxpos[labels] = positions
        result += [maxpos[idxs]]
    return result

def minimum(input, labels = None, index = None):
    """
    Calculate the minimum of the values of an array over labeled regions.

    Parameters
    ----------

    input: array-like
        Array-like of values. For each region specified by `labels`, the
        minimal values of `input` over the region is computed.

    labels: array-like, optional
        An array-like of integers marking different regions over which the
        minimum value of `input` is to be computed. `labels` must have the
        same shape as `input`. If `labels` is not specified, the minimum
        over the whole array is returned.

    index: array-like, optional
        A list of region labels that are taken into account for computing the
        minima. If index is None, the minimum over all elements where `labels`
        is non-zero is returned.

    Returns
    -------
    output : float or list of floats
        List of minima of `input` over the regions determined by `labels` and
        whose index is in `index`. If `index` or `labels` are not specified, a
        float is returned: the minimal value of `input` if `labels` is None,
        and the minimal value of elements where `labels` is greater than zero
        if `index` is None.

    See also
    --------

    label, maximum, minimum_position, extrema, sum, mean, variance,
    standard_deviation

    Notes
    -----

    The function returns a Python list and not a Numpy array, use
    `np.array` to convert the list to an array.

    Examples
    --------

    >>> a = np.array([[1, 2, 0, 0],
    ...               [5, 3, 0, 4],
    ...               [0, 0, 0, 7],
    ...               [9, 3, 0, 0]])
    >>> labels, labels_nb = ndimage.label(a)
    >>> labels
    array([[1, 1, 0, 0],
           [1, 1, 0, 2],
           [0, 0, 0, 2],
           [3, 3, 0, 0]])
    >>> ndimage.minimum(a, labels=labels, index=np.arange(1, labels_nb + 1))
    [1.0, 4.0, 3.0]
    >>> ndimage.minimum(a)
    0.0
    >>> ndimage.minimum(a, labels=labels)
    1.0

    """
    return _select(input, labels, index, find_min=True)[0]

def maximum(input, labels = None, index = None):
    """
    Calculate the maximum of the values of an array over labeled regions.

    Parameters
    ----------
    input : array_like
        Array-like of values. For each region specified by `labels`, the
        maximal values of `input` over the region is computed.
    labels : array_like, optional
        An array of integers marking different regions over which the
        maximum value of `input` is to be computed. `labels` must have the
        same shape as `input`. If `labels` is not specified, the maximum
        over the whole array is returned.
    index : array_like, optional
        A list of region labels that are taken into account for computing the
        maxima. If index is None, the maximum over all elements where `labels`
        is non-zero is returned.

    Returns
    -------
    output : float or list of floats
        List of maxima of `input` over the regions determined by `labels` and
        whose index is in `index`. If `index` or `labels` are not specified, a
        float is returned: the maximal value of `input` if `labels` is None,
        and the maximal value of elements where `labels` is greater than zero
        if `index` is None.

    See also
    --------
    label, minimum, maximum_position, extrema, sum, mean, variance,
    standard_deviation

    Notes
    -----
    The function returns a Python list and not a Numpy array, use
    `np.array` to convert the list to an array.

    Examples
    --------

    >>> a = np.arange(16).reshape((4,4))
    >>> a
    array([[ 0,  1,  2,  3],
           [ 4,  5,  6,  7],
           [ 8,  9, 10, 11],
           [12, 13, 14, 15]])
    >>> labels = np.zeros_like(a)
    >>> labels[:2,:2] = 1
    >>> labels[2:, 1:3] = 2
    >>> labels
    array([[1, 1, 0, 0],
           [1, 1, 0, 0],
           [0, 2, 2, 0],
           [0, 2, 2, 0]])
    >>> from scipy import ndimage
    >>> ndimage.maximum(a)
    15.0
    >>> ndimage.maximum(a, labels=labels, index=[1,2])
    [5.0, 14.0]
    >>> ndimage.maximum(a, labels=labels)
    14.0

    >>> b = np.array([[1, 2, 0, 0],
                      [5, 3, 0, 4],
                      [0, 0, 0, 7],
                      [9, 3, 0, 0]])
    >>> labels, labels_nb = ndimage.label(b)
    >>> labels
    array([[1, 1, 0, 0],
           [1, 1, 0, 2],
           [0, 0, 0, 2],
           [3, 3, 0, 0]])
    >>> ndimage.maximum(b, labels=labels, index=np.arange(1, labels_nb + 1))
    [5.0, 7.0, 9.0]

    """
    return _select(input, labels, index, find_max=True)[0]

def minimum_position(input, labels = None, index = None):
    """Find the positions of the minimums of the values of an array at labels.

    Labels must be None or an array of the same dimensions as the input.  

    Index must be None, a single label or sequence of labels.  If
    none, all values where label is greater than zero are used.
    """
    
    dims = numpy.array(numpy.asarray(input).shape)
    # see numpy.unravel_index to understand this line.
    dim_prod = numpy.cumprod([1] + list(dims[:0:-1]))[::-1]

    result = _select(input, labels, index, find_min_positions=True)[0]

    if numpy.isscalar(result):
        return tuple((result // dim_prod) % dims)

    return [tuple(v) for v in (result.reshape(-1, 1) // dim_prod) % dims]

def maximum_position(input, labels = None, index = None):
    """Find the positions of the maximums of the values of an array at labels.

    Labels must be None or an array of the same dimensions as the input.  

    Index must be None, a single label or sequence of labels.  If
    none, all values where label is greater than zero are used.
    """
    
    dims = numpy.array(numpy.asarray(input).shape)
    # see numpy.unravel_index to understand this line.
    dim_prod = numpy.cumprod([1] + list(dims[:0:-1]))[::-1]

    result = _select(input, labels, index, find_max_positions=True)[0]

    if numpy.isscalar(result):
        return tuple((result // dim_prod) % dims)

    return [tuple(v) for v in (result.reshape(-1, 1) // dim_prod) % dims]

def extrema(input, labels = None, index = None):
    """Calculate the minimums and maximums of the values of an array
    at labels, along with their positions.

    Labels must be None or an array of the same dimensions as the input.  

    Index must be None, a single label or sequence of labels.  If
    none, all values where label is greater than zero are used.
    
    Returns: minimums, maximums, min_positions, max_positions
    """
    
    dims = numpy.array(numpy.asarray(input).shape)
    # see numpy.unravel_index to understand this line.
    dim_prod = numpy.cumprod([1] + list(dims[:0:-1]))[::-1]

    minimums, min_positions, maximums, max_positions = _select(input, labels, index, 
                                                               find_min=True, find_max=True, 
                                                               find_min_positions=True, find_max_positions=True)


    if numpy.isscalar(minimums):
        return minimums, maximums, tuple((min_positions // dim_prod) % dims), tuple((max_positions // dim_prod) % dims)

    min_positions = [tuple(v) for v in (min_positions.reshape(-1, 1) // dim_prod) % dims]
    max_positions = [tuple(v) for v in (max_positions.reshape(-1, 1) // dim_prod) % dims]

    return minimums, maximums, min_positions, max_positions

def center_of_mass(input, labels = None, index = None):
    """Calculate the center of mass of the values of an array at labels.

    Labels must be None or an array of the same dimensions as the input.  

    Index must be None, a single label or sequence of labels.  If
    none, all values where label is greater than zero are used.
    """

    normalizer = sum(input, labels, index)
    grids = numpy.ogrid[[slice(0, i) for i in input.shape]]

    results = [sum(input * grids[dir].astype(float), labels, index) / normalizer for dir in range(input.ndim)]
    
    if numpy.isscalar(results[0]):
        return tuple(results)

    return [tuple(v) for v in numpy.array(results).T]

def histogram(input, min, max, bins, labels = None, index = None):
    """Calculate the histogram of the values of an array at labels.

    Labels must be None or an array of the same dimensions as the input.  

    The histograms are defined by the minimum and maximum values and the
    number of bins.

    Index must be None, a single label or sequence of labels.  If
    none, all values where label is greater than zero are used.
    """
    
    _bins = numpy.linspace(min, max, bins + 1)

    def _hist(vals):
        return numpy.histogram(vals, _bins)[0]

    return labeled_comprehension(input, labels, index, _hist, object, None, pass_positions=False)

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
