"""Masked arrays add-ons.

A collection of utilities for maskedarray

:author: Pierre Gerard-Marchant
:contact: pierregm_at_uga_dot_edu
:version: $Id: extras.py 38 2006-12-09 23:01:14Z backtopop $
"""
__author__ = "Pierre GF Gerard-Marchant ($Author: backtopop $)"
__version__ = '1.0'
__revision__ = "$Revision: 38 $"
__date__     = '$Date: 2006-12-09 18:01:14 -0500 (Sat, 09 Dec 2006) $'

__all__ = ['apply_along_axis', 'atleast_1d', 'atleast_2d', 'atleast_3d',
           'vstack', 'hstack', 'dstack', 'row_stack', 'column_stack',
           'count_masked', 
           'masked_all', 'masked_all_like', 'mr_',
           'stdu', 'varu',
           ]

import core
reload(core)
from core import *
from core import _arraymethod

import numpy
import numpy.core.numeric as numeric
from numpy.core.numeric import ndarray
from numpy.core.numeric import array as nxarray
from numpy.core.fromnumeric import asarray as nxasarray

from numpy.lib.index_tricks import concatenator
import numpy.lib.function_base as function_base

def issequence(seq):
    """Returns True if the argumnet is a sequence (ndarray, list or tuple."""
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

def masked_all(shape, dtype):
    """Returns an empty masked array of the given shape and dtype,
    where all the data are masked."""
    a = empty(shape, dtype)
    a[:] = masked
    return a

def masked_all_like(arr):
    """Returns an empty masked array of the same shape and dtype as the array `a`,
    where all the data are masked."""
    a = empty_like(arr)
    a[:] = masked
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
                _d = func.__call__(nxasarray(x), **params)
                _m = func.__call__(getmaskarray(x), **params)
                return masked_array(_d, mask=_m)
            elif isinstance(x, tuple):
                _d = func.__call__(tuple([nxasarray(a) for a in x]), **params)
                _m = func.__call__(tuple([getmaskarray(a) for a in x]), **params)
                return masked_array(_d, mask=_m)
        else:
            arrays = []
            args = list(args)
            while len(args)>0 and issequence(args[0]):
                arrays.append(args.pop(0))
            res = []
            for x in arrays:
                _d = func.__call__(nxasarray(x), *args, **params)
                _m = func.__call__(getmaskarray(x), *args, **params)
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
    print dtypes
    if not hasattr(arr, '_mask'):
        return numeric.asarray(outarr, dtype=max(dtypes))
    else:
        return outarr.astype(max(dtypes))

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
