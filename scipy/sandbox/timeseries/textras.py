"""
Extras functions for time series.

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknox_ca_at_hotmail_dot_com
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'


import numpy
import maskedarray
from maskedarray import masked

import const as _c
from tseries import TimeSeries



__all__ = ['isleapyear', 'count_missing', 'accept_atmost_missing']

#..............................................................................
def isleapyear(year):
    """Returns true if year is a leap year.
    
:Input:
    year : integer / sequence
        A given (list of) year(s).
    """
    year = numpy.asarray(year)
    return numpy.logical_or(year % 400 == 0,
                            numpy.logical_and(year % 4 == 0, year % 100 > 0))

#..............................................................................    
def count_missing(series):
    """Returns the number of missing data per period.
    
    
Notes
----- 
This function is designed to return the actual number of missing values when 
a series has been converted from one frequency to a smaller frequency.
    
For example, converting a 12-month-long daily series to months will yield 
a (12x31) array, with missing values in February, April, June... 
count_missing will discard these extra missing values.
    """
    if not isinstance(series, TimeSeries):
        raise TypeError, "The input data should be a valid TimeSeries object! "\
                         "(got %s instead)" % type(series)
    if series.ndim == 1:
        return len(series) - series.count()
    elif series.ndim != 2:
        raise NotImplementedError
    #
    missing =  series.shape[-1] - series.count(axis=-1)
    period = series.shape[-1]
    freq = series.freq
    if (period == 366) and (freq//_c.FR_ANN == 1):
        # row: years, cols: days
        missing -= ~isleapyear(series.year)
    elif period == 31 and (freq//_c.FR_MTH == 1):
        months = series.months
        # row: months, cols: days
        missing[numpy.array([m in [4,6,9,11] for m in months])] -= 1
        isfeb = (months == 2)
        missing[isfeb] -= 2
        missing[isfeb & ~isleapyear(series.year)] -= 1
    elif period not in (12,7):
        raise NotImplementedError, "Not yet implemented for that frequency..."
    return missing
    
#.............................................................................
def accept_atmost_missing(series, max_missing, strict=False):
    """Masks the rows of the series that contains more than max_missing missing data.
    Returns a new masked series.
    
:Inputs:
    series : TimeSeries
        Input time series.
    max_missing : float
        Number of maximum acceptable missing values per row (if larger than 1),
        or maximum acceptable percentage of missing values (if lower than 1).
    strict : boolean *[False]*
        Whether the
    """
    series = numpy.array(series, copy=True, subok=True)
    if not isinstance(series, TimeSeries):
        raise TypeError, "The input data should be a valid TimeSeries object! "\
                         "(got %s instead)" % type(series)
    # Find the number of missing values ....
    missing = count_missing(series)
    # Transform an acceptable percentage in a number
    if max_missing < 1:
        max_missing = numpy.round(max_missing * series.shape[-1],0)
    #
    series.unshare_mask()
    if strict:
        series[missing > max_missing] = masked
    else:
        series[missing >= max_missing] = masked
    return series
    