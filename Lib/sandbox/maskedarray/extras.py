"""Masked arrays add-ons.

A collection of utilities for maskedarray

:author: Pierre Gerard-Marchant
:contact: pierregm_at_uga_dot_edu
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

__all__ = ['apply_along_axis', 'atleast_1d', 'atleast_2d', 'atleast_3d',
               'average',
           'vstack', 'hstack', 'dstack', 'row_stack', 'column_stack',
           'compress_rowcols', 'compress_rows', 'compress_cols', 'count_masked', 
           'dot',
           'mask_rowcols', 'mask_rows', 'mask_cols', 'masked_all', 
               'masked_all_like', 'mr_',
           'notmasked_edges', 'notmasked_contiguous',
           'stdu', 'varu',
           ]

from itertools import groupby

import core
#reload(core)
from core import *

import numpy
from numpy import float_
import numpy.core.umath as umath 
import numpy.core.numeric as numeric
from numpy.core.numeric import ndarray
from numpy.core.numeric import array as nxarray
from numpy.core.fromnumeric import asarray as nxasarray

from numpy.lib.index_tricks import concatenator
import numpy.lib.function_base as function_base

#...............................................................................
def issequence(seq):
    """Returns True if the argumnet is a sequence (ndarray, list or tuple)."""
    if isinstance(seq, ndarray):
        return True
    elif isinstance(seq, tuple):
        return True
    elif isinstance(seq, list):
        return True
    return False

def count_masked(arr, axis=None):
    """Counts the number of masked elements along the given axis."""
    m = getmaskarray(arr)
    return m.sum(axis)

def masked_all(shape, dtype=float_):
    """Returns an empty masked array of the given shape and dtype,
    where all the data are masked."""
    a = masked_array(numeric.empty(shape, dtype),
                     mask=numeric.ones(shape, bool_))
    return a

def masked_all_like(arr):
    """Returns an empty masked array of the same shape and dtype as the array `a`,
    where all the data are masked."""
    a = masked_array(numeric.empty_like(arr),
                     mask=numeric.ones(shape, bool_))
    return a

#####--------------------------------------------------------------------------
#---- --- New methods ---
#####-------------------------------------------------------------------------- 
def varu(a, axis=None, dtype=None):
    """a.var(axis=None, dtype=None)
    Returns an unbiased estimate of the variance.
    
    Instead of dividing the sum of squared anomalies (SSA) by n, the number of 
    elements, the SSA is divided by n-1.
        """
    a = asarray(a)
    cnt = a.count(axis=axis)
    anom = a.anom(axis=axis, dtype=dtype)
    anom *= anom
    dvar = anom.sum(axis) / (cnt-1)
    if axis is None:
        return dvar
    return a.__class__(dvar, 
                          mask=mask_or(a._mask.all(axis), (cnt==1)),
                          fill_value=a._fill_value)
            
def stdu(a, axis=None, dtype=None):
    """a.var(axis=None, dtype=None)
    Returns an unbiased estimate of the standard deviation.

    Instead of dividing the sum of squared anomalies (SSA) by n, the number of 
    elements, the SSA is divided by n-1.
        """
    a = asarray(a)
    dvar = a.varu(axis,dtype)
    if axis is None:
        if dvar is masked:
            return masked
        else:
            # Should we use umath.sqrt instead ?
            return sqrt(dvar)
    return a.__class__(sqrt(dvar._data), mask=dvar._mask, 
                          fill_value=a._fill_value)

MaskedArray.stdu = stdu
MaskedArray.varu = varu

#####--------------------------------------------------------------------------
#---- --- Standard functions ---
#####--------------------------------------------------------------------------
class _fromnxfunction:
    """Defines a wrapper to adapt numpy functions to masked arrays."""
    def __init__(self, funcname):
        self._function = funcname
        self.__doc__ = self.getdoc()
    def getdoc(self):
        "Retrieves the __doc__ string from the function."
        return getattr(numpy, self._function).__doc__ +\
            "(The function is applied to both the _data and the mask, if any.)"
    def __call__(self, *args, **params):
        func = getattr(numpy, self._function)
        if len(args)==1:
            x = args[0]
            if isinstance(x,ndarray):
                _d = func(nxasarray(x), **params)
                _m = func(getmaskarray(x), **params)
                return masked_array(_d, mask=_m)
            elif isinstance(x, tuple) or isinstance(x, list):
                _d = func(tuple([nxasarray(a) for a in x]), **params)
                _m = func(tuple([getmaskarray(a) for a in x]), **params)
                return masked_array(_d, mask=_m)
        else:
            arrays = []
            args = list(args)
            while len(args)>0 and issequence(args[0]):
                arrays.append(args.pop(0))
            res = []
            for x in arrays:
                _d = func(nxasarray(x), *args, **params)
                _m = func(getmaskarray(x), *args, **params)
                res.append(masked_array(_d, mask=_m))
            return res
                
atleast_1d = _fromnxfunction('atleast_1d')     
atleast_2d = _fromnxfunction('atleast_2d')  
atleast_3d = _fromnxfunction('atleast_3d')             

vstack = row_stack = _fromnxfunction('vstack')
hstack = _fromnxfunction('hstack')
column_stack = _fromnxfunction('column_stack')
dstack = _fromnxfunction('dstack')

#####--------------------------------------------------------------------------
#---- 
#####--------------------------------------------------------------------------
def apply_along_axis(func1d,axis,arr,*args):
    """ Execute func1d(arr[i],*args) where func1d takes 1-D arrays
        and arr is an N-d array.  i varies so as to apply the function
        along the given axis for each 1-d subarray in arr.
    """
    arr = numeric.asanyarray(arr)
    nd = arr.ndim
    if axis < 0:
        axis += nd
    if (axis >= nd):
        raise ValueError("axis must be less than arr.ndim; axis=%d, rank=%d."
            % (axis,nd))
    ind = [0]*(nd-1)
    i = numeric.zeros(nd,'O')
    indlist = range(nd)
    indlist.remove(axis)
    i[axis] = slice(None,None)
    outshape = numeric.asarray(arr.shape).take(indlist)
    i.put(indlist, ind)
    res = func1d(arr[tuple(i.tolist())],*args)
    #  if res is a number, then we have a smaller output array
    asscalar = numeric.isscalar(res)
    if not asscalar:
        try:
            len(res)
        except TypeError:
            asscalar = True
    # Note: we shouldn't set the dtype of the output from the first result...
    #...so we force the type to object, and build a list of dtypes
    #...we'll just take the largest, to avoid some downcasting
    dtypes = []
    if asscalar:
        dtypes.append(numeric.asarray(res).dtype)
        outarr = zeros(outshape, object_)
        outarr[ind] = res
        Ntot = numeric.product(outshape)
        k = 1
        while k < Ntot:
            # increment the index
            ind[-1] += 1
            n = -1
            while (ind[n] >= outshape[n]) and (n > (1-nd)):
                ind[n-1] += 1
                ind[n] = 0
                n -= 1
            i.put(indlist,ind)
            res = func1d(arr[tuple(i.tolist())],*args)
            outarr[ind] = res
            dtypes.append(asarray(res).dtype)
            k += 1
    else:
        Ntot = numeric.product(outshape)
        holdshape = outshape
        outshape = list(arr.shape)
        outshape[axis] = len(res)
        dtypes.append(asarray(res).dtype)
        outarr = zeros(outshape, object_)
        outarr[tuple(i.tolist())] = res
        k = 1
        while k < Ntot:
            # increment the index
            ind[-1] += 1
            n = -1
            while (ind[n] >= holdshape[n]) and (n > (1-nd)):
                ind[n-1] += 1
                ind[n] = 0
                n -= 1
            i.put(indlist, ind)
            res = func1d(arr[tuple(i.tolist())],*args)
            outarr[tuple(i.tolist())] = res
            dtypes.append(asarray(res).dtype)
            k += 1
    if not hasattr(arr, '_mask'):
        return numeric.asarray(outarr, dtype=max(dtypes))
    else:
        return outarr.astype(max(dtypes))


def average (a, axis=None, weights=None, returned = 0):
    """average(a, axis=None weights=None, returned=False)

    Averages the array over the given axis.  If the axis is None, averages
    over all dimensions of the array.  Equivalent to a.mean(axis)

    If an integer axis is given, this equals:
        a.sum(axis) * 1.0 / size(a, axis)

    If axis is None, this equals:
        a.sum(axis) * 1.0 / a.size

    If weights are given, result is:
        sum(a * weights,axis) / sum(weights,axis),
    where the weights must have a's shape or be 1D with length the
    size of a in the given axis. Integer weights are converted to
    Float.  Not specifying weights is equivalent to specifying
    weights that are all 1.

    If 'returned' is True, return a tuple: the result and the sum of
    the weights or count of values. The shape of these two results
    will be the same.

    Returns masked values instead of ZeroDivisionError if appropriate.
    
    """
    a = asarray(a)
    mask = a.mask
    ash = a.shape
    if ash == ():
        ash = (1,)
    if axis is None:
        if mask is nomask:
            if weights is None:
                n = a.sum(axis=None)
                d = float(a.size)
            else:
                w = filled(weights, 0.0).ravel()
                n = umath.add.reduce(a._data.ravel() * w)
                d = umath.add.reduce(w)
                del w
        else:
            if weights is None:
                n = a.filled(0).sum(axis=None)
                d = umath.add.reduce((-mask).ravel().astype(int_))
            else:
                w = array(filled(weights, 0.0), float, mask=mask).ravel()
                n = add.reduce(a.ravel() * w)
                d = add.reduce(w)
                del w
    else:
        if mask is nomask:
            if weights is None:
                d = ash[axis] * 1.0
                n = add.reduce(a._data, axis, dtype=float_)
            else:
                w = filled(weights, 0.0)
                wsh = w.shape
                if wsh == ():
                    wsh = (1,)
                if wsh == ash:
                    w = numeric.array(w, float_, copy=0)
                    n = add.reduce(a*w, axis)
                    d = add.reduce(w, axis)
                    del w
                elif wsh == (ash[axis],):
                    ni = ash[axis]
                    r = [None]*len(ash)
                    r[axis] = slice(None, None, 1)
                    w = eval ("w["+ repr(tuple(r)) + "] * ones(ash, float)")
                    n = add.reduce(a*w, axis, dtype=float_)
                    d = add.reduce(w, axis, dtype=float_)
                    del w, r
                else:
                    raise ValueError, 'average: weights wrong shape.'
        else:
            if weights is None:
                n = add.reduce(a, axis, dtype=float_)
                d = umath.add.reduce((-mask), axis=axis, dtype=float_)
            else:
                w = filled(weights, 0.0)
                wsh = w.shape
                if wsh == ():
                    wsh = (1,)
                if wsh == ash:
                    w = array(w, dtype=float_, mask=mask, copy=0)
                    n = add.reduce(a*w, axis, dtype=float_)
                    d = add.reduce(w, axis, dtype=float_)
                elif wsh == (ash[axis],):
                    ni = ash[axis]
                    r = [None]*len(ash)
                    r[axis] = slice(None, None, 1)
                    w = eval ("w["+ repr(tuple(r)) + "] * masked_array(ones(ash, float), mask)")
                    n = add.reduce(a*w, axis, dtype=float_)
                    d = add.reduce(w, axis, dtype=float_)
                else:
                    raise ValueError, 'average: weights wrong shape.'
                del w
    if n is masked or d is masked: 
        return masked
    result = n/d
    del n
    
    if isMaskedArray(result):
        if ((axis is None) or (axis==0 and a.ndim == 1)) and \
           (result.mask is nomask):
            result = result._data
        if returned:
            if not isMaskedArray(d):
                d = masked_array(d)
            if isinstance(d, ndarray) and (not d.shape == result.shape):
                d = ones(result.shape, dtype=float_) * d
    if returned:
        return result, d
    else:
        return result
    
#..............................................................................
def compress_rowcols(x, axis=None):
    """Suppresses the rows and/or columns of a 2D array that contains masked values.
    The suppression behavior is selected with the `axis`parameter.
        - If axis is None, rows and columns are suppressed. 
        - If axis is 0, only rows are suppressed. 
        - If axis is 1 or -1, only columns are suppressed.
    Returns a *pure* ndarray.    
    """
    x = asarray(x)
    if x.ndim <> 2:
        raise NotImplementedError, "compress2d works for 2D arrays only."
    m = getmask(x)
    # Nothing is masked: return x
    if m is nomask or not m.any():
        return x._data
    # All is masked: return empty
    if m.all():
        return nxarray([])
    # Builds a list of rows/columns indices
    (idxr, idxc) = (range(len(x)), range(x.shape[1]))
    masked = m.nonzero()
    if not axis:
        for i in function_base.unique(masked[0]):
            idxr.remove(i)
    if axis in [None, 1, -1]:
        for j in function_base.unique(masked[1]):
            idxc.remove(j)
    return x._data[idxr][:,idxc]

def compress_rows(a):
    """Suppresses whole rows of a 2D array that contain masked values."""
    return compress_rowcols(a,0)

def compress_cols(a):
    """Suppresses whole columnss of a 2D array that contain masked values."""
    return compress_rowcols(a,1)

def mask_rowcols(a, axis=None):
    """Masks whole rows and/or columns of a 2D array that contain masked values.
    The masking behavior is selected with the `axis`parameter.
        - If axis is None, rows and columns are suppressed. 
        - If axis is 0, only rows are suppressed. 
        - If axis is 1 or -1, only columns are suppressed.
    Returns a *pure* ndarray.    
    """
    a = asarray(a)
    if a.ndim != 2:
        raise NotImplementedError, "compress2d works for 2D arrays only."
    m = getmask(a)
    # Nothing is masked: return a
    if m is nomask or not m.any():
        return a
    maskedval = m.nonzero()
    if not axis:
        a[function_base.unique(maskedval[0])] = masked
    if axis in [None, 1, -1]:
        a[:,function_base.unique(maskedval[1])] = masked
    return a

def mask_rows(a, axis=None):
    """Masks whole rows of a 2D array that contain masked values."""
    return mask_rowcols(a, 0)

def mask_cols(a, axis=None):
    """Masks whole columns of a 2D array that contain masked values."""
    return mask_rowcols(a, 1)

        
def dot(a,b):
    """Returns the dot product of two 2D masked arrays a and b.
    Like the generic numpy equivalent the product sum is over
    the last dimension of a and the second-to-last dimension of b.
    
    Masked values are propagated: if a masked value appears in a row or column,
    the whole row or column is considered masked.
    
    NB: The first argument is not conjugated.
    """
    #TODO: Works only with 2D arrays. There should be a way to get it to run with higher dimension
    a = mask_rows(a)
    b = mask_cols(b)
    #
    d = numpy.dot(a.filled(0), b.filled(0))
    #
    am = (~getmaskarray(a))
    bm = (~getmaskarray(b))
    m = ~numpy.dot(am,bm)
    return masked_array(d, mask=m)


#####--------------------------------------------------------------------------
#---- --- Concatenation helpers ---
#####--------------------------------------------------------------------------

class mconcatenator(concatenator):
    """Translates slice objects to concatenation along an axis."""

    def __init__(self, axis=0):
        concatenator.__init__(self, axis, matrix=False)

    def __getitem__(self,key):
        if isinstance(key, str):
            raise MAError, "Unavailable for masked array."
        if type(key) is not tuple:
            key = (key,)
        objs = []
        scalars = []
        final_dtypedescr = None
        for k in range(len(key)):
            scalar = False
            if type(key[k]) is slice:
                step = key[k].step
                start = key[k].start
                stop = key[k].stop
                if start is None: 
                    start = 0
                if step is None:
                    step = 1
                if type(step) is type(1j):
                    size = int(abs(step))
                    newobj = function_base.linspace(start, stop, num=size)
                else:
                    newobj = numeric.arange(start, stop, step)
            elif type(key[k]) is str:
                if (key[k] in 'rc'):
                    self.matrix = True
                    self.col = (key[k] == 'c')
                    continue
                try:
                    self.axis = int(key[k])
                    continue
                except (ValueError, TypeError):
                    raise ValueError, "Unknown special directive"
            elif type(key[k]) in numeric.ScalarType:
                newobj = asarray([key[k]])
                scalars.append(k)
                scalar = True
            else:
                newobj = key[k]
            objs.append(newobj)
            if isinstance(newobj, numeric.ndarray) and not scalar:
                if final_dtypedescr is None:
                    final_dtypedescr = newobj.dtype
                elif newobj.dtype > final_dtypedescr:
                    final_dtypedescr = newobj.dtype
        if final_dtypedescr is not None:
            for k in scalars:
                objs[k] = objs[k].astype(final_dtypedescr)
        res = concatenate(tuple(objs),axis=self.axis)
        return self._retval(res)

class mr_class(mconcatenator):
    """Translates slice objects to concatenation along the first axis.

        For example:
        >>> r_[array([1,2,3]), 0, 0, array([4,5,6])]
        array([1, 2, 3, 0, 0, 4, 5, 6])
    """
    def __init__(self):
        mconcatenator.__init__(self, 0)

mr_ = mr_class()

#####--------------------------------------------------------------------------
#---- ---
#####--------------------------------------------------------------------------

def flatnotmasked_edges(a):
    """Finds the indices of the first and last not masked values in a  1D masked array.
    If all values are masked, returns None.
    """
    m = getmask(a)
    if m is nomask or not numpy.any(m):
        return [0,-1]
    unmasked = numeric.flatnonzero(~m)
    if len(unmasked) > 0:
        return unmasked[[0,-1]]
    else:
        return None

def notmasked_edges(a, axis=None):
    """Finds the indices of the first and last not masked values along the given
    axis in a masked array.
    If all values are masked, returns None.
    Otherwise, returns a list of 2 tuples, corresponding to the indices of the
    first and last unmasked values respectively.
    """
    a = asarray(a)
    if axis is None or a.ndim == 1:
        return flatnotmasked_edges(a)
    m = getmask(a)
    idx = array(numpy.indices(a.shape), mask=nxasarray([m]*a.ndim))
    return [tuple([idx[i].min(axis).compressed() for i in range(a.ndim)]),
            tuple([idx[i].max(axis).compressed() for i in range(a.ndim)]),]

def flatnotmasked_contiguous(a):
    """Finds contiguous unmasked data in a flattened masked array.
    Returns a sorted sequence of tuples (size,(start index, end index)).
    """
    m = getmask(a)
    if m is nomask:
        return (a.size, [0,-1])
    unmasked = numeric.flatnonzero(~m)
    if len(unmasked) == 0:
        return None
    result = []
    for k, group in groupby(enumerate(unmasked), lambda (i,x):i-x):
        tmp = numpy.fromiter((g[1] for g in group), int_)
        result.append((tmp.size, tuple(tmp[[0,-1]])))
    result.sort()
    return result

def notmasked_contiguous(a, axis=None):
    """Finds contiguous unmasked data in a masked array along the given axis.
    Returns a sorted sequence of tuples (size,(start index, end index)).
    Note: Only accepts 2D arrays at most.
    """
    a = asarray(a)
    nd = a.ndim
    if nd > 2:
        raise NotImplementedError,"Currently limited to atmost 2D array."
    if axis is None or nd == 1:
        return flatnotmasked_contiguous(a)
    #
    result = []
    #
    other = (axis+1)%2
    idx = [0,0]
    idx[axis] = slice(None,None)
    #
    for i in range(a.shape[other]):
        idx[other] = i
        result.append( flatnotmasked_contiguous(a[idx]) )
    return result
  
