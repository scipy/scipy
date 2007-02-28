"""Time Series add-ons.

A collection of utilities for timeseries

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknox_ca_at_hotmail_dot_com
:version: $Id: tcore.py 2752 2007-02-22 20:50:12Z mattknox_ca $
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author: mattknox_ca $)"
__version__ = '1.0'
__revision__ = "$Revision: 2752 $"
__date__     = '$Date: 2007-02-22 15:50:12 -0500 (Thu, 22 Feb 2007) $'

__all__ = [
    'forward_fill', 'backward_fill', 'interp_masked1d',
    'expmave'
          ]

import numpy as N
from scipy.interpolate import fitpack
import maskedarray as MA
from maskedarray import masked, nomask, getmask, MaskedArray
import numpy.core.numeric as numeric

#####---------------------------------------------------------------------------
#---- --- Functions for filling in masked values in a masked array ---
#####---------------------------------------------------------------------------

def forward_fill(marr, maxgap=None):
    """forward_fill(marr, maxgap=None)

Forward fill masked values in a 1-d array when there are <= maxgap
consecutive masked values. If maxgap is None, then forward fill all
masked values."""

    if numeric.ndim(marr) > 1:
        raise ValueError,"The input array should be 1D only!"

    marr = MA.array(marr, copy=True)
    if getmask(marr) is nomask or marr.size == 0:
        return marr

    currGap = 0

    if maxgap is not None:
        for i in range(1, marr.size):
            if marr._mask[i]:
                currGap += 1
                if currGap <= maxgap and not marr._mask[i-1]:
                    marr._data[i] = marr._data[i-1]
                    marr._mask[i] = False
                elif currGap == maxgap + 1:
                    marr._mask[i-maxgap:i] = True
            else:
                currGap = 0               
    else:
        for i in range(1, marr.size):
            if marr._mask[i] and not marr._mask[i-1]:
                marr._data[i] = marr._data[i-1]
                marr._mask[i] = False
    return marr


def backward_fill(marr, maxgap=None):
    """backward_fill(marr, maxgap=None)

backward fill masked values in a 1-d array when there are <= maxgap
consecutive masked values. If maxgap is None, then backward fill all
masked values."""
    return forward_fill(marr[::-1], maxgap=maxgap)[::-1]
    

def interp_masked1d(marr, kind='linear'):
    """interp_masked1d(marr, king='linear')

Interpolate masked values in marr according to method kind.
kind must be one of 'constant', 'linear', 'cubic', quintic'
"""
    if numeric.ndim(marr) > 1: 
        raise ValueError("array must be 1 dimensional!")
    #
    marr = MA.array(marr, copy=True)
    if getmask(marr) is nomask: 
        return marr
    #
    unmaskedIndices = (~marr._mask).nonzero()[0]
    if unmaskedIndices.size < 2: 
        return marr
    #    
    kind = kind.lower()
    if kind == 'constant': 
        return forward_fill(marr)
    try:
        k = {'linear' : 1,
             'cubic' : 3,
             'quintic' : 5}[kind.lower()]
    except KeyError:
        raise ValueError("Unsupported interpolation type.")
    
    first_unmasked, last_unmasked = MA.extras.flatnotmasked_edges(marr)
    
    vals = marr.data[unmaskedIndices]
    
    tck = fitpack.splrep(unmaskedIndices, vals, k=k)
    
    maskedIndices = marr._mask.nonzero()[0]
    interpIndices = maskedIndices[(maskedIndices > first_unmasked) & \
                                  (maskedIndices < last_unmasked)]
    marr[interpIndices] = fitpack.splev(interpIndices, tck).astype(marr.dtype)
    return marr


#####---------------------------------------------------------------------------
#---- --- Moving average functions                                           ---
#####---------------------------------------------------------------------------

def expmave(data, n, tol=1e-6):
    """calculate exponential moving average of a series.

:Parameters:
    - `data` (ndarray, MaskedArray) : data is a valid ndarray or MaskedArray
      or an instance of a subclass of these types. In particular, TimeSeries
      objects are valid here.
    - `n` (int) : time periods. Where the smoothing factor is 2/(n + 1)
    - `tol` (float, *[1e-6]*) : when `data` contains masked values, this
      parameter will determine what points in the result should be masked.
      Values in the result that would not be "significantly" impacted (as
      determined by this parameter) by the masked values are left unmasked.
"""
    if isinstance(data, MaskedArray):
        ismasked = (data._mask is not nomask)
    else:
        ismasked = False
        
    k = 2./float(n + 1)
    def expmave_sub(a, b):
        return b + k * (a - b)
        
    if ismasked:
        data = data.filled(0)

    result = N.frompyfunc(expmave_sub, 2, 1).accumulate(data).astype(data.dtype)
    if ismasked:
        _unmasked = N.logical_not(MA.getmask(data)).astype(N.float_)
        marker = 1 - N.frompyfunc(expmave_sub, 2, 1).accumulate(_unmasked)
        result[marker > tol] = masked

    return result
