"""

A collection of filters for timeseries

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknox_ca_at_hotmail_dot_com
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import numpy as N
from numpy import bool_, float_
narray = N.array

from scipy.signal import convolve, get_window

import maskedarray as MA
from maskedarray import MaskedArray, nomask, getmask, getmaskarray, masked
marray = MA.array


__all__ = ['expmave'
           'running_window', 'running_mean'           
           ]

#####---------------------------------------------------------------------------
#---- --- Moving average functions ---
#####---------------------------------------------------------------------------
def expmave(data, n, tol=1e-6):
    """Calculates the exponential moving average of a series.

:Parameters:
    data : ndarray
        Data as a valid (subclass of) ndarray or MaskedArray. In particular, 
        TimeSeries objects are valid here.
    n : int 
        Time periods. The smoothing factor is 2/(n + 1)
    tol : float, *[1e-6]*
        Tolerance for the definition of the mask. When data contains masked 
        values, this parameter determinea what points in the result should be masked.
        Values in the result that would not be "significantly" impacted (as 
        determined by this parameter) by the masked values are left unmasked.
"""
    data = marray(data, copy=True, subok=True)
    ismasked = (data._mask is not nomask)
    data._mask = N.zeros(data.shape, bool_)
    _data = data._data
    #
    k = 2./float(n + 1)
    def expmave_sub(a, b):
        return b + k * (a - b)
    #
    data._data.flat = N.frompyfunc(expmave_sub, 2, 1).accumulate(_data)
    if ismasked:
        _unmasked = N.logical_not(data._mask).astype(float_)
        marker = 1. - N.frompyfunc(expmave_sub, 2, 1).accumulate(_unmasked)
        data._mask[marker > tol] = True
    data._mask[0] = True
    #
    return data

def weightmave(data, n):
    data = marray(data, subok=True, copy=True)
    data._mask = N.zeros(data.shape, bool_)
    # Set the data
    _data = data._data
    tmp = N.empty_like(_data)
    tmp[:n] = _data[:n]
    s = 0
    for i in range(n, len(data)):
        s += _data[i] - _data[i-n]
        tmp[i] = n*_data[i] + tmp[i-1] - s
    tmp *= 2./(n*(n+1))
    data._data.flat = tmp
    # Set the mask
    if data._mask is not nomask:
        msk = data._mask.nonzero()[0].repeat(n).reshape(-1,n)
        msk += range(n)
        data._mask[msk.ravel()] = True
    data._mask[:n] = True
    return data


#...............................................................................
def running_window(data, window_type, window_size,
                   centered=False, trailing=False):
    """Applies a running window of type window_type and size window_size on the
data. Returns a (subclass of) MaskedArray.

:Parameters:
    data : ndarray
        Data to process. The array should be at most 2D. On 2D arrays, the
        window
        is applied recursively on each column.
    window_type : string/tuple/float
        Window type (see Notes)
    window_size : integer
        The width of the window.
    centered : boolean, *[False]*
        If both centered and trailing are False, then centered is forced to
        True. If centered, the result at position i uses data points from
        [i-k:i+k+1] in the calculation. The k first and k last data are always
        masked (with k=window_size//2). When data has a missing value at
        position i, the result has missing values in the interval [i-k:i+k+1].
    trailing : boolean, *[False]*
        If trailing is True, the result at position i uses data points from
        [i-window_size:i+1] in the calculation.the first "window_size" data
        points are always masked. When data has a missing value at position i,
        the result has missing values in the interval [i-window_size:i+1].
        
Notes
-----

The recognized window types are: boxcar, triang, blackman, hamming, hanning, 
bartlett, parzen, bohman, blackmanharris, nuttall, barthann, kaiser (needs beta), 
gaussian (needs std), general_gaussian (needs power, width), slepian (needs width).
If the window requires parameters, the window_type argument should be a tuple
with the first argument the string name of the window, and the next arguments 
the needed parameters. If window_type is a floating point number, it is interpreted 
as the beta parameter of the kaiser window.

Note also that only boxcar has been thoroughly tested.
    """

    if not centered and not trailing:
        centered = True
    elif centered and trailing:
        raise ValueError("Cannot specify both centered and trailing")
    
    
    data = marray(data, copy=True, subok=True)
    if data._mask is nomask:
        data._mask = N.zeros(data.shape, bool_)
    window = get_window(window_type, window_size, fftbins=False)
    n = len(data)
    
    if centered: k = window_size//2
    else:        k = 0

    if data.ndim == 1:
        
        data._data.flat = convolve(data._data, window)[k:n+k] / float(window_size)
        data._mask[:] = ((convolve(getmaskarray(data), window) > 0)[k:n+k])
    elif data.ndim == 2:
        for i in range(data.shape[-1]):
            _data = data._data[:,i]
            _data.flat = convolve(_data, window)[k:n+k] / float(window_size)
            data._mask[:,i] = (convolve(data._mask[:,i], window) > 0)[k:n+k]
    else:
        raise ValueError, "Data should be at most 2D"
        
    if centered:
        data._mask[:k] = data._mask[-k:] = True
    else:
        data._mask[:window_size] = True
    return data

def running_mean(data, width,
                 centered=False, trailing=False):
    """Computes the running mean of size width on the data. Returns a
(subclass of) MaskedArray.

:Parameters:
    data : ndarray
        Data to process. The array should be at most 2D. On 2D arrays, the window
        is applied recursively on each column.
    window_size : integer
        The width of the window.
    centered : boolean, *[False]*
        If both centered and trailing are False, then centered is forced to
        True. If centered, the result at position i uses data points from
        [i-k:i+k+1] in the calculation. The k first and k last data are always
        masked (with k=window_size//2). When data has a missing value at
        position i, the result has missing values in the interval [i-k:i+k+1].
    trailing : boolean, *[False]*
        If trailing is True, the result at position i uses data points from
        [i-window_size:i+1] in the calculation.the first "window_size" data
        points are always masked. When data has a missing value at position i,
        the result has missing values in the interval [i-window_size:i+1]."""

    return running_window(data, 'boxcar', width,
                          centered=centered, trailing=trailing)

################################################################################
if __name__ == '__main__':
    from maskedarray.testutils import assert_equal, assert_almost_equal
    from timeseries import time_series, thisday
    #
    data = MA.arange(100)
