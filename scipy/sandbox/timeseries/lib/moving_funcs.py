"""

A collection of moving functions for masked arrays and time series

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknox_ca_at_hotmail_dot_com
:version: $Id: filters.py 2819 2007-03-03 23:00:20Z pierregm $
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author: pierregm $)"
__version__ = '1.0'
__revision__ = "$Revision: 2819 $"
__date__     = '$Date: 2007-03-03 18:00:20 -0500 (Sat, 03 Mar 2007) $'

__all__ = ['mov_sum', 'mov_median', 'mov_min', 'mov_max',
           'mov_average', 'mov_mean', 'mov_average_expw',
           'mov_stddev', 'mov_var', 'mov_covar', 'mov_corr',
           'cmov_average', 'cmov_mean', 'cmov_window'
           ]

import numpy as N
from numpy import bool_, float_
narray = N.array

from scipy.signal import convolve, get_window

import maskedarray as MA
from maskedarray import MaskedArray, nomask, getmask, getmaskarray, masked
marray = MA.array

from timeseries.cseries import MA_mov_stddev, MA_mov_sum, MA_mov_average, \
                               MA_mov_median, MA_mov_min, MA_mov_max

def _process_result_dict(orig_data, result_dict):
    "process the results from the c function"

    rarray = result_dict['array']
    rmask = result_dict['mask']

    # makes a copy of the appropriate type
    data = orig_data.astype(rarray.dtype).copy()
    data.flat = rarray.ravel()
    if not hasattr(data, '__setmask__'):
        data = data.view(MA.MaskedArray)
    data.__setmask__(rmask)
    return data

def _moving_func(data, cfunc, kwargs):

    if data.ndim == 1:
        kwargs['array'] = data

        result_dict = cfunc(**kwargs)
        return _process_result_dict(data, result_dict)

    elif data.ndim == 2:
        for i in range(data.shape[-1]):
            kwargs['array'] = data[:,i]
            result_dict = cfunc(**kwargs)
            
            if i == 0:
                rtype = result_dict['array'].dtype
                result = data.astype(rtype)
                print data.dtype, result.dtype
            
            rmask = result_dict['mask']

            curr_col = marray(result_dict['array'], mask=rmask, copy=False)
            result[:,i] = curr_col

        return result

    else:
        raise ValueError, "Data should be at most 2D"

#...............................................................................
def mov_sum(data, span, dtype=None):
    """Calculates the moving sum of a series.

*Parameters*:
    $$data$$
    $$span$$
    $$dtype$$"""

    kwargs = {'span':span}
    if dtype is not None:
        kwargs['dtype'] = dtype
 
    return _moving_func(data, MA_mov_sum, kwargs)
#...............................................................................
def mov_median(data, span, dtype=None):
    """Calculates the moving median of a series.

*Parameters*:
    $$data$$
    $$span$$
    $$dtype$$"""

    kwargs = {'span':span}
    if dtype is not None:
        kwargs['dtype'] = dtype

    return _moving_func(data, MA_mov_median, kwargs)
#...............................................................................
def mov_min(data, span, dtype=None):
    """Calculates the moving minimum of a series.

*Parameters*:
    $$data$$
    $$span$$
    $$dtype$$"""

    kwargs = {'span':span}
    if dtype is not None:
        kwargs['dtype'] = dtype

    return _moving_func(data, MA_mov_min, kwargs)
#...............................................................................
def mov_max(data, span, dtype=None):
    """Calculates the moving max of a series.

*Parameters*:
    $$data$$
    $$span$$
    $$dtype$$"""

    kwargs = {'span':span}
    if dtype is not None:
        kwargs['dtype'] = dtype

    return _moving_func(data, MA_mov_max, kwargs)
#...............................................................................
def mov_average(data, span, dtype=None):
    """Calculates the moving average of a series.

*Parameters*:
    $$data$$
    $$span$$
    $$dtype$$"""
    
    kwargs = {'span':span}
    if dtype is not None:
        kwargs['dtype'] = dtype

    return _moving_func(data, MA_mov_average, kwargs)
mov_mean = mov_average
#...............................................................................
def _mov_var_stddev(data, span, is_variance, bias, dtype):
    "helper function for mov_var and mov_stddev functions"

    kwargs = {'span':span,
              'is_variance':is_variance,
              'bias':bias}
    if dtype is not None:
        kwargs['dtype'] = dtype

    return _moving_func(data, MA_mov_stddev, kwargs)
#...............................................................................
def mov_var(data, span, bias=False, dtype=None):
    """Calculates the moving variance of a 1-D array.

*Parameters*:
    $$data$$
    $$span$$
    $$bias$$
    $$dtype$$"""
    
    return _mov_var_stddev(data=data, span=span,
                           is_variance=1, bias=int(bias), dtype=dtype)
#...............................................................................
def mov_stddev(data, span, bias=False, dtype=None):
    """Calculates the moving standard deviation of a 1-D array.

*Parameters*:
    $$data$$
    $$span$$
    $$bias$$
    $$dtype$$"""
    
    return _mov_var_stddev(data=data, span=span,
                           is_variance=0, bias=int(bias), dtype=dtype)
#...............................................................................
def mov_covar(x, y, span, bias=False, dtype=None):
    """Calculates the moving covariance of two 1-D arrays.

*Parameters*:
    $$x$$
    $$y$$
    $$span$$
    $$bias$$
    $$dtype$$"""
    
    result = x - mov_average(x, span, dtype=dtype)
    result = result * (y - mov_average(y, span, dtype=dtype))
    
    if bias: denom = span
    else: denom = span - 1
    
    return result/denom
#...............................................................................
def mov_corr(x, y, span, dtype=None):
    """Calculates the moving correlation of two 1-D arrays.

*Parameters*:
    $$x$$
    $$y$$
    $$span$$
    $$dtype$$"""

    result = mov_covar(x, y, span, bias=True, dtype=dtype)
    result = result / mov_stddev(x, span, bias=True, dtype=dtype)
    result = result / mov_stddev(y, span, bias=True, dtype=dtype)
   
    return result
#...............................................................................
def mov_average_expw(data, span, tol=1e-6):
    """Calculates the exponentially weighted moving average of a series.

*Parameters*:
    $$data$$
    span : int 
        Time periods. The smoothing factor is 2/(span + 1)
    tol : float, *[1e-6]*
        Tolerance for the definition of the mask. When data contains masked 
        values, this parameter determinea what points in the result shoud be
        masked. Values in the result that would not be "significantly"
        impacted (as determined by this parameter) by the masked values are
        left unmasked."""

    data = marray(data, copy=True, subok=True)
    ismasked = (data._mask is not nomask)
    data._mask = N.zeros(data.shape, bool_)
    _data = data._data
    #
    k = 2./float(span + 1)
    def expmave_sub(a, b):
        return a + k * (b - a)
    #
    data._data.flat = N.frompyfunc(expmave_sub, 2, 1).accumulate(_data)
    if ismasked:
        _unmasked = N.logical_not(data._mask).astype(float_)
        marker = 1. - N.frompyfunc(expmave_sub, 2, 1).accumulate(_unmasked)
        data._mask[marker > tol] = True
    data._mask[0] = True
    #
    return data
#.............................................................................
def cmov_window(data, span, window_type):
    """Applies a centered moving window of type window_type and size span on
the data.

Returns a (subclass of) MaskedArray. The k first and k last data are always 
masked (with k=span//2). When data has a missing value at position i, the
result has missing values in the interval [i-k:i+k+1].
    
    
*Parameters*:
    data : {ndarray}
        Data to process. The array should be at most 2D. On 2D arrays, the
        window is applied recursively on each column.
    span : {int}
        The width of the window.
    window_type : {string/tuple/float}
        Window type (see Notes)
        
*Notes*:

    The recognized window types are: boxcar, triang, blackman, hamming,
    hanning, bartlett, parzen, bohman, blackmanharris, nuttall, barthann,
    kaiser (needs beta), gaussian (needs std), general_gaussian (needs power,
    width), slepian (needs width). If the window requires parameters, the
    window_type argument should be a tuple with the first argument the string
    name of the window, and the next arguments the needed parameters. If
    window_type is a floating point number, it is interpreted as the beta
    parameter of the kaiser window.

    Note also that only boxcar has been thoroughly tested.
"""

    data = marray(data, copy=True, subok=True)
    if data._mask is nomask:
        data._mask = N.zeros(data.shape, bool_)
    window = get_window(window_type, span, fftbins=False)
    (n, k) = (len(data), span//2)
    #
    if data.ndim == 1:
        data._data.flat = convolve(data._data, window)[k:n+k] / float(span)
        data._mask[:] = ((convolve(getmaskarray(data), window) > 0)[k:n+k])
    elif data.ndim == 2:
        for i in range(data.shape[-1]):
            _data = data._data[:,i]
            _data.flat = convolve(_data, window)[k:n+k] / float(span)
            data._mask[:,i] = (convolve(data._mask[:,i], window) > 0)[k:n+k]
    else:
        raise ValueError, "Data should be at most 2D"
    data._mask[:k] = data._mask[-k:] = True
    return data

def cmov_average(data, span):
    """Computes the centered moving average of size span on the data.
    
*Parameters*:
    data : {ndarray}
        Data to process. The array should be at most 2D. On 2D arrays, the
        window is applied recursively on each column.
    span : {int}
        The width of the window.

*Returns*:    
    A (subclass of) MaskedArray. The k first and k last data are always masked
    (with k=span//2). When data has a missing value at position i, the result
    has missing values in the interval [i-k:i+k+1].
"""
    return cmov_window(data, span, 'boxcar')

cmov_mean = cmov_average

param_doc = {}
param_doc['data'] = \
"""data : {ndarray}
        Data must be an ndarray (or subclass). In particular, note that
        TimeSeries objects are valid here."""

param_doc['x'] = \
"""x : {ndarray}
        First array to be included in the calculation. x must be an ndarray (or
        subclass). In particular, note that TimeSeries objects are valid here."""

param_doc['y'] = \
"""y : {ndarray}
        Second array to be included in the calculation. y must be an ndarray (or
        subclass). In particular, note that TimeSeries objects are valid here."""

param_doc['span'] = \
"""span : {int }
        Time periods to use for each calculation."""

param_doc['bias'] = \
"""bias : {False, True}, optional
        If False, Normalization is by (N-1) where N == span (unbiased
        estimate).  If True then normalization is by N."""

param_doc['dtype'] = \
"""dtype : {numpy data type specification}, optional
        dtype for the result"""

mov_result_doc = \
"""

*Returns*:
    The result is always a masked array (preserves subclass attributes). The
    result at index i uses values from [i-span:i+1], and will be masked for
    the first `span` values. The result will also be masked at i if any of the
    input values in the slice [i-span:i+1] are masked.
"""

_g = globals()

# generate function doc strings
for fn in (x for x in __all__ if x[:4] == 'mov_' and x[4:] != 'mean'):
    fdoc = _g[fn].func_doc
    for prm, dc in param_doc.iteritems():
        fdoc = fdoc.replace('$$'+prm+'$$', dc)
    fdoc += mov_result_doc
    _g[fn].func_doc = fdoc


###############################################################################
if __name__ == '__main__':
    from timeseries import time_series, today
    from maskedarray.testutils import assert_equal, assert_almost_equal
    #
    series = time_series(N.arange(10),start_date=today('D'))
    #
    filtered = mov_sum(series,3)
    assert_equal(filtered, [0,1,3,6,9,12,15,18,21,24])
    assert_equal(filtered._mask, [1,1,0,0,0,0,0,0,0,0])
    assert_equal(filtered._dates, series._dates)
    assert_equal(series, N.arange(10))
    #
    filtered = mov_average(series,3)
    assert_equal(filtered, [0,1,1,2,3,4,5,6,7,8])
    assert_equal(filtered._mask, [1,1,0,0,0,0,0,0,0,0])
    assert_equal(filtered._dates, series._dates)
    assert_equal(series, N.arange(10))
    #
    filtered = mov_average(series._data,3)
    assert_equal(filtered, [0,1,1,2,3,4,5,6,7,8])
    assert_equal(filtered._mask, [1,1,0,0,0,0,0,0,0,0])
    assert_equal(series, N.arange(10))
    
