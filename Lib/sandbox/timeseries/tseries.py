# pylint: disable-msg=W0201, W0212
"""
Core classes for time/date related arrays.

The `DateArray` class provides a base for the creation of date-based objects,
using days as the base units. This class could be adapted easily to objects
with a smaller temporal resolution (for example, using one hour, one second as the
base unit).

The `TimeSeries` class provides  a base for the definition of time series.
A time series is defined here as the combination of two arrays:
    
    - an array storing the time information (as a `DateArray` instance);
    - an array storing the data (as a `MaskedArray` instance.

These two classes were liberally adapted from `MaskedArray` class.



:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknox_ca_at_hotmail_dot_com
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import numpy
from numpy import ndarray
from numpy.core import bool_, complex_, float_, int_, object_
import numpy.core.fromnumeric as fromnumeric
import numpy.core.numeric as numeric
import numpy.core.umath as umath
from numpy.core.records import recarray
from numpy.core.records import fromarrays as recfromarrays

import maskedarray.core as MA
from maskedarray.core import MaskedArray, MAError, masked, nomask, \
    filled, getmask, getmaskarray, make_mask_none, mask_or, make_mask, \
    masked_array

import tcore as corelib
from tcore import *

import tdates
from tdates import DateError, InsufficientDateError
from tdates import Date, isDate, DateArray, isDateArray, \
    date_array, date_array_fromlist, date_array_fromrange, thisday, today

import cseries



__all__ = [
'TimeSeriesError','TimeSeriesCompatibilityError','TimeSeries','isTimeSeries',
'time_series', 'tsmasked',
'mask_period','mask_inside_period','mask_outside_period','compressed',
'adjust_endpoints','align_series','aligned','convert','group_byperiod',
'tshift','fill_missing_dates', 'stack', 'concatenate_series','empty_like',
'day_of_week','day_of_year','day','month','quarter','year',
'hour','minute','second',  
'tofile','asrecords','flatten',
           ]

#...............................................................................
#
#ufunc_domain = {}
#ufunc_fills = {}

#### --------------------------------------------------------------------------
#--- ... TimeSeriesError class ...
#### --------------------------------------------------------------------------
class TimeSeriesError(Exception):
    "Class for TS related errors."
    def __init__ (self, args=None):
        "Creates an exception."
        Exception.__init__(self)
        self.args = args
    def __str__(self):
        "Calculates the string representation."
        return str(self.args)
    __repr__ = __str__

class TimeSeriesCompatibilityError(TimeSeriesError):
    """Defines the exception raised when series are incompatible."""
    def __init__(self, mode, first, second):
        if mode == 'freq':
            msg = "Incompatible time steps! (%s <> %s)"
        elif mode == 'start_date':
            msg = "Incompatible starting dates! (%s <> %s)"
        elif mode == 'size':
            msg = "Incompatible sizes! (%s <> %s)"
        else:
            msg = "Incompatibility !  (%s <> %s)"
        msg = msg % (first, second)
        TimeSeriesError.__init__(self, msg)

#def _compatibilitycheck(a, b):
def _timeseriescompat(a, b):
    """Checks the date compatibility of two TimeSeries object.
    Returns True if everything's fine, or raises an exception."""
    if not (hasattr(a,'freq') and hasattr(b, 'freq')):
        return True
    if a.freq != b.freq:
        raise TimeSeriesCompatibilityError('freq', a.freq, b.freq)
    elif a.start_date != b.start_date:
        raise TimeSeriesCompatibilityError('start_date', 
                                           a.start_date, b.start_date)
    else:
        step_diff = a._dates.get_steps() != b._dates.get_steps()
        if (step_diff is True) or (hasattr(step_diff, "any") and step_diff.any()):
            raise TimeSeriesCompatibilityError('time_steps', 
                                               a._dates.get_steps(), b._dates.get_steps())
        elif a.shape != b.shape:
            raise TimeSeriesCompatibilityError('size', "1: %s" % str(a.shape), 
                                                   "2: %s" % str(b.shape))
    return True
    

def _timeseriescompat_multiple(*series):
    """Checks the date compatibility of multiple TimeSeries objects.
    Returns True if everything's fine, or raises an exception. Unlike
    the binary version, all items must be TimeSeries objects."""
    
    freqs, start_dates, steps, shapes = \
        zip(*[(ser.freq, ser.start_date,
               (ser._dates.get_steps() != series[0]._dates.get_steps()).any(),
               ser.shape)  for ser in series])
                
    if len(set(freqs)) > 1:
        errItems = tuple(set(freqs))
        raise TimeSeriesCompatibilityError('freq', errItems[0], errItems[1])
    
    if len(set(start_dates)) > 1:
        errItems = tuple(set(start_dates))
        raise TimeSeriesCompatibilityError('start_dates', errItems[0], errItems[1])

            
    if max(steps) == True:
        bad_index = [x for x, val in enumerate(steps) if val][0]
        raise TimeSeriesCompatibilityError('time_steps', 
                series[0]._dates.get_steps(), series[bad_index]._dates.get_steps())        

    if len(set(shapes)) > 1:
        errItems = tuple(set(shapes))
        raise TimeSeriesCompatibilityError('size', "1: %s" % str(errItems[0].shape), 
                                                   "2: %s" % str(errItems[1].shape))

    return True


def _datadatescompat(data,dates):
    """Checks the compatibility of dates and data at the creation of a TimeSeries.
    Returns True if everything's fine, raises an exception otherwise."""    
    # If there's only 1 element, the date is a Date object, which has no size...
    tsize = numeric.size(dates)
    dsize = data.size
    # Only one data
    if dsize == tsize:
        return True
    elif data.ndim > 1:
        #dsize = numeric.asarray(data.shape)[:-1].prod()
        dsize = data.shape[0]
        if dsize == tsize:
            return True 
    elif data.ndim == 0 and tsize <= 1:
        return True   
    raise TimeSeriesCompatibilityError('size', "data: %s" % dsize, 
                                               "dates: %s" % tsize)

def _getdatalength(data):
    "Estimates the length of a series (size/nb of variables)."
    if numeric.ndim(data) >= 2:
        return numeric.asarray(numeric.shape(data))[:-1].prod()
    else:
        return numeric.size(data)
    
def _compare_frequencies(*series):
    """Compares the frequencies of a sequence of series. 
    Returns the common frequency, or raises an exception if series have different 
    frequencies."""
    unique_freqs = numpy.unique([x.freqstr for x in series])
    try:
        common_freq = unique_freqs.item()
    except ValueError:
        raise TimeSeriesError, \
            "All series must have same frequency!"
    return common_freq    

##### --------------------------------------------------------------------------
##--- ... Time Series ...
##### --------------------------------------------------------------------------
class _tsmathmethod(object):
    """Defines a wrapper for arithmetic array methods (add, mul...).
When called, returns a new TimeSeries object, with the new series the result of
the method applied on the original series.
The `_dates` part remains unchanged.
    """
    def __init__ (self, methodname):
        self._name = methodname
    #
    def __get__(self, obj, objtype=None):
        "Gets the calling object."
        self.obj = obj
        return self
    #
    def __call__ (self, other, *args):
        "Execute the call behavior."
        instance = self.obj
        if isinstance(other, TimeSeries):
            assert(_timeseriescompat(instance, other))
        func = getattr(super(TimeSeries, instance), self._name)
        result = func(other, *args).view(type(instance))
        result._dates = instance._dates
        return result

class _tsarraymethod(object):
    """Defines a wrapper for basic array methods.
When called, returns a new TimeSeries object, with the new series the result of
the method applied on the original series.
If `ondates` is True, the same operation is performed on the `_dates`.
If `ondates` is False, the `_dates` part remains unchanged.
    """
    def __init__ (self, methodname, ondates=False):
        """abfunc(fillx, filly) must be defined.
           abinop(x, filly) = x for all x to enable reduce.
        """
        self._name = methodname
        self._ondates = ondates
    #
    def __get__(self, obj, objtype=None):
        self.obj = obj
        return self
    #
    def __call__ (self, *args):
        "Execute the call behavior."
        _name = self._name
        instance = self.obj
        func_series = getattr(super(TimeSeries, instance), _name)
        result = func_series(*args)
        if self._ondates:
            result._dates = getattr(instance._dates, _name)(*args)
        else:
            result._dates = instance._dates
        return result

class _tsaxismethod(object):
    """Defines a wrapper for array methods working on an axis (mean...).
When called, returns a ndarray, as the result of the method applied on the series.
    """
    def __init__ (self, methodname):
        """abfunc(fillx, filly) must be defined.
           abinop(x, filly) = x for all x to enable reduce.
        """
        self._name = methodname
    #
    def __get__(self, obj, objtype=None):
        self.obj = obj
        return self
    #
    def __call__ (self, *args, **params):
        "Execute the call behavior."
        (_dates, _series) = (self.obj._dates, self.obj._series)
        func = getattr(_series, self._name)
        result = func(*args, **params)
        if _dates.size == _series.size:
            return result
        else:
            try:
                axis = params.get('axis', args[0])
                if axis in [-1, _series.ndim-1]:
                    result = TimeSeries(result, dates=_dates)
            except IndexError:
                pass
            return result

class TimeSeries(MaskedArray, object):     
    """Base class for the definition of time series.
A time series is here defined as the combination of three arrays:
    
    - `series` : *[ndarray]*
        Data part
    - `mask` : *[ndarray]*
        Mask part
    - `dates` : *[DateArray]*
        Date part
    
The combination of `series` and `dates` is the `data` part.
    """
    options = None
    _defaultobserved = None
    _genattributes = ['fill_value', 'observed']
    def __new__(cls, data, dates=None, mask=nomask, 
                freq=None, observed=None, start_date=None, 
                dtype=None, copy=False, fill_value=None,
                keep_mask=True, small_mask=True, hard_mask=False, **options):
        maparms = dict(copy=copy, dtype=dtype, fill_value=fill_value,
                       keep_mask=keep_mask, small_mask=small_mask, 
                       hard_mask=hard_mask,)
        # Get the data ...............................
        _data = MaskedArray(data, mask=mask, **maparms).view(cls)
        # Get the frequency ..........................
        freq = corelib.check_freq(freq)
        # Get the dates ..............................
        if dates is None:
            newdates = getattr(data, '_dates', None)
        else:
            newdates = dates
        if newdates is not None:
            if not hasattr(newdates, 'freq'):
                newdates = date_array(dlist=dates, freq=freq)
            if freq is not None and newdates.freq != freq:
                newdates = newdates.tofreq(freq)
        else:
            dshape = _data.shape
            if len(dshape) > 0:
                newdates = date_array(start_date=start_date, length=dshape[0],
                                      freq=freq)
            else:
                newdates = date_array([], freq=freq)          
        # Get observed ...............................
        observed = getattr(data, 'observed', corelib.fmtObserv(observed))
        
        if _data is masked:
            assert(numeric.size(newdates)==1)
            return _data.view(cls)
        assert(_datadatescompat(_data,newdates))
        _data._dates = newdates
        _data._defaultdates = _data._dates
        _data.observed = observed
        return _data
    #............................................
    def __array_finalize__(self,obj):
        MaskedArray.__array_finalize__(self, obj)        
        self._dates = getattr(obj, '_dates', [])
        self.observed = getattr(obj, 'observed', None)
        return
    #..................................
    def __array_wrap__(self, obj, context=None):
        result = super(TimeSeries, self).__array_wrap__(obj, context)
        result._dates = self._dates
        return result
    #............................................
    def _get_series(self):
        if self._mask.ndim == 0 and self._mask:
            return masked
        return self.view(MaskedArray)
    _series = property(fget=_get_series)
    #............................................
    def __checkindex(self, indx):
        "Checks the validity of an index."
        if isinstance(indx, int):
            return (indx, indx)
        elif isinstance(indx, str):
            indx = self._dates.date_to_index(Date(self._dates.freq, string=indx))
            return (indx, indx)
        elif isDate(indx):
            indx = self._dates.date_to_index(indx)
            return (indx, indx)
        elif isDateArray(indx):
            if indx.size == 1:
                indx = self._dates.date_to_index(indx[0])
                return (indx,indx)
            else:
                d2i = self._dates.date_to_index
                tmp = numpy.fromiter((d2i(i) for i in indx),int_)
                return (tmp,tmp)
        elif isinstance(indx,slice):
            slice_start = self.__checkindex(indx.start)[0]
            slice_stop = self.__checkindex(indx.stop)[0]
            indx = slice(slice_start, slice_stop, indx.step)
            return (indx,indx)
        elif isinstance(indx, tuple):
            if len(indx) > self.shape:
                raise IndexError, "Too many indices"
            if self._dates.size == self.size:
                return (indx, indx)
            return (indx,indx[0])
#            elif len(indx)==2:
#                return (indx,indx[0])
#            return (indx,indx[:-1])
        elif isTimeSeries(indx):
            indx = indx._series
        if getmask(indx) is not nomask:
            msg = "Masked arrays must be filled before they can be used as indices!"
            raise IndexError, msg
        return (indx,indx)

    def __getitem__(self, indx):
        """x.__getitem__(y) <==> x[y]
Returns the item described by i. Not a copy as in previous versions.
        """
        (sindx, dindx) = self.__checkindex(indx)
        newdata = numeric.array(self._series[sindx], copy=False, subok=True)
        newdate = self._dates[dindx]
        m = self._mask
        singlepoint = (len(numeric.shape(newdate))==0)
        if singlepoint:
            newdate = DateArray(newdate)
            if newdata is masked:
                newdata = tsmasked
                newdata._dates = newdate
                return newdata
            elif self.ndim > 1:
                # CHECK: use reshape, or set shape ?
                newshape = (list((1,)) + list(newdata.shape))
                newdata.shape = newshape
        newdata = newdata.view(type(self))    
        newdata._dates = newdate
        return newdata
# CHECK : The implementation below should work, but does not. Why ?
#        newdata = numeric.array(self._data[sindx], copy=False)
#        newdates = self._dates[dindx]
#        if self._mask is not nomask:
#            newmask = self._mask.copy()[sindx]
#        else:
#            newmask = nomask
#        singlepoint = (len(numeric.shape(newdates))==0)
#        if singlepoint:
#            if newmask.ndim == 0 and newmask:
#                output = tsmasked
#                output._dates = newdates
#                return output
#            if self.ndim > 1:
#                # CHECK: use reshape, or set shape ?
#                newdata = newdata.reshape((list((1,)) + list(newdata.shape)))
#                if newmask is not nomask:
#                    newmask.shape = newdata.shape
#        newdata = newdata.view(type(self))    
#        newdata._dates = newdates
#        newdata._mask = newmask
#        return newdata



    #........................
    def __setitem__(self, indx, value):
        """x.__setitem__(i, y) <==> x[i]=y
Sets item described by index. If value is masked, masks those locations.
        """
        if self is masked:
            raise MAError, 'Cannot alter the masked element.'
        (sindx, dindx) = self.__checkindex(indx)
        #....        
        super(TimeSeries, self).__setitem__(sindx, value)    
    #........................
    def __getslice__(self, i, j):
        "Gets slice described by i, j"
        (si,di) = self.__checkindex(i)
        (sj,dj) = self.__checkindex(j)
        result = super(TimeSeries, self).__getitem__(slice(si,sj))
        result._dates = self._dates[di:dj]
        return result
    #....
    def __setslice__(self, i, j, value):
        "Gets item described by i. Not a copy as in previous versions."
        (si,_) = self.__checkindex(i)
        (sj,_) = self.__checkindex(j)
        #....
        if isinstance(value, TimeSeries):
            assert(_timeseriescompat(self[si:sj], value))
        super(TimeSeries, self).__setitem__(slice(si,sj), value)
    #......................................................
    def __str__(self):
        """Returns a string representation of self (w/o the dates...)"""
        return str(self._series)
    def __repr__(self):
        """Calculates the repr representation, using masked for fill if
           it is enabled. Otherwise fill with fill value.
        """
        desc = """\
timeseries(data  =
 %(data)s,
           dates = 
 %(time)s, 
           freq  = %(freq)s)
"""
        desc_short = """\
timeseries(data  = %(data)s,
           dates = %(time)s,
           freq  = %(freq)s)
"""
        if numeric.size(self._dates) > 2 and self.isvalid():
            timestr = "[%s ... %s]" % (str(self._dates[0]),str(self._dates[-1]))
        else:
            timestr = str(self.dates)

        if self.ndim <= 1:
            return desc_short % {'data': str(self._series),
                                 'time': timestr,
                                 'freq': self.freqstr, }
        return desc % {'data': str(self._series),
                       'time': timestr,
                       'freq': self.freqstr, }
    #............................................
    __add__ = _tsmathmethod('__add__')
    __radd__ = _tsmathmethod('__add__')
    __sub__ = _tsmathmethod('__sub__')
    __rsub__ = _tsmathmethod('__rsub__')
    __pow__ = _tsmathmethod('__pow__')
    __mul__ = _tsmathmethod('__mul__')
    __rmul__ = _tsmathmethod('__mul__')
    __div__ = _tsmathmethod('__div__')
    __rdiv__ = _tsmathmethod('__rdiv__')
    __truediv__ = _tsmathmethod('__truediv__')
    __rtruediv__ = _tsmathmethod('__rtruediv__')
    __floordiv__ = _tsmathmethod('__floordiv__')
    __rfloordiv__ = _tsmathmethod('__rfloordiv__')
    __eq__ = _tsmathmethod('__eq__')
    __ne__ = _tsmathmethod('__ne__')
    __lt__ = _tsmathmethod('__lt__')
    __le__ = _tsmathmethod('__le__')
    __gt__ = _tsmathmethod('__gt__')
    __ge__ = _tsmathmethod('__ge__')
    
    reshape = _tsarraymethod('reshape', ondates=True)
    copy = _tsarraymethod('copy', ondates=True)
    compress = _tsarraymethod('compress', ondates=True)
    ravel = _tsarraymethod('ravel', ondates=True)
    cumsum = _tsarraymethod('cumsum',ondates=False)
    cumprod = _tsarraymethod('cumprod',ondates=False)
    anom = _tsarraymethod('anom',ondates=False)
    
    sum = _tsaxismethod('sum')
    prod = _tsaxismethod('prod')
    mean = _tsaxismethod('mean')
    var = _tsaxismethod('var')
    varu = _tsaxismethod('varu')
    std = _tsaxismethod('std')
    stdu = _tsaxismethod('stdu')
    all = _tsaxismethod('all')
    any = _tsaxismethod('any')

    
#    def nonzero(self):
#        """Returns a tuple of ndarrays, one for each dimension of the array,
#    containing the indices of the non-zero elements in that dimension."""
#        return self._series.nonzero()
    
#    filled = _tsarraymethod('filled', ondates=False)
    
    #............................................ 
    def ids (self):
        """Return the ids of the data, dates and mask areas"""
        return (id(self._series), id(self.dates),)   
    #------------------------------------------------------
    @property
    def series(self):
        "Returns the series."
        return self._series
    @property
    def dates(self):
        """Returns the dates"""
        return self._dates
    @property
    def freq(self):
        """Returns the corresponding frequency (as an integer)."""
        return self._dates.freq
    @property
    def freqstr(self):
        """Returns the corresponding frequency (as a string)."""
        return self._dates.freqstr
    @property
    def day(self):          
        "Returns the day of month for each date in self._dates."
        return self._dates.day
    @property
    def day_of_week(self):  
        "Returns the day of week for each date in self._dates."
        return self._dates.day_of_week
    @property
    def day_of_year(self):  
        "Returns the day of year for each date in self._dates."
        return self._dates.day_of_year
    @property
    def month(self):        
        "Returns the month for each date in self._dates."
        return self._dates.month
    @property
    def quarter(self):   
        "Returns the quarter for each date in self._dates."   
        return self._dates.quarter
    @property
    def year(self):         
        "Returns the year for each date in self._dates."
        return self._dates.year
    @property
    def second(self):    
        "Returns the seconds for each date in self._dates."  
        return self._dates.second
    @property
    def minute(self):     
        "Returns the minutes for each date in self._dates."  
        return self._dates.minute
    @property
    def hour(self):         
        "Returns the hour for each date in self._dates."
        return self._dates.hour
    @property
    def week(self):
        "Returns the week for each date in self._dates."
        return self._dates.week

    days = day
    weekdays = day_of_week
    yeardays = day_of_year
    months = month
    quarters = quarter
    years = year
    seconds = second
    minutes = minute
    hours = hour
    weeks = week

    @property
    def start_date(self):
        """Returns the first date of the series."""
        if self._dates.size != 0:
            return self._dates[0]
        else:
            return None
    @property
    def end_date(self):
        """Returns the last date of the series."""
        if self._dates.size != 0:
            return self._dates[-1]
        else:
            return None
    
    def isvalid(self):
        """Returns whether the series has no duplicate/missing dates."""
        return self._dates.isvalid()
    
    def has_missing_dates(self):
        """Returns whether there's a date gap in the series."""
        return self._dates.has_missing_dates()
    
    def isfull(self):
        """Returns whether there's no date gap in the series."""
        return self._dates.isfull()
    
    def has_duplicated_dates(self):
        """Returns whether there are duplicated dates in the series."""
        return self._dates.has_duplicated_dates()
    
    def date_to_index(self, date):
        "Returns the index corresponding to a given date, as an integer."
        return self._dates.date_to_index(date) 
    #.....................................................
    def asfreq(self, freq=None):
        "Converts the dates to another frequency."
        if freq is None:
            return self
        return TimeSeries(self._series, dates=self._dates.asfreq(freq))
    
    def convert(self, freq, func='auto', position='END'):
        "Converts the dates to another frequency, and adapt the data."
        return convert(self, freq, func=func, position=position)
    #.....................................................
    def transpose(self, *axes):
        """ a.transpose(*axes)

    Returns a view of 'a' with axes transposed. If no axes are given,
    or None is passed, switches the order of the axes. For a 2-d
    array, this is the usual matrix transpose. If axes are given,
    they describe how the axes are permuted.

        """
        if self._dates.size == self.size:
            result = super(TimeSeries, self).transpose(*axes)
            result._dates = self._dates.transpose(*axes)
        else:
            errmsg = "Operation not permitted on multi-variable series"
            if (len(axes)==0) or axes[0] != 0:
                raise TimeSeriesError, errmsg
            else:
                result = super(TimeSeries, self).transpose(*axes)
                result._dates = self._dates
        return result
    #......................................................
    def copy_attributes(self, oldseries, exclude=[]):
        "Copies the attributes from oldseries if they are not in the exclude list."
        attrlist = type(self)._genattributes
        if not isinstance(oldseries, TimeSeries):
            msg = "Series should be a valid TimeSeries object! (got <%s> instead)"
            raise TimeSeriesError, msg % type(oldseries)
        for attr in attrlist:
            if not attr in exclude:
                setattr(self, attr, getattr(oldseries, attr))
    
        
def _attrib_dict(series, exclude=[]):
    """this function is used for passing through attributes of one
time series to a new one being created"""
    result = {'fill_value':series.fill_value,
              'observed':series.observed}
    return dict(filter(lambda x: x[0] not in exclude, result.iteritems()))
  
        
##### --------------------------------------------------------------------------
##--- ... Additional methods ...
##### --------------------------------------------------------------------------

#.......................................


class _tsblockedmethods(object):
    """Defines a wrapper for array methods that should be temporarily disabled.
    """
    def __init__ (self, methodname):
        """abfunc(fillx, filly) must be defined.
           abinop(x, filly) = x for all x to enable reduce.
        """
        self._name = methodname
    #
    def __get__(self, obj, objtype=None):
        self.obj = obj
        return self
    #
    def __call__ (self, *args, **params):
        raise NotImplementedError
#TimeSeries.transpose = _tsarraymethod('transpose', ondates=True)
TimeSeries.swapaxes = _tsarraymethod('swapaxes', ondates=True)


#####---------------------------------------------------------------------------
#---- --- Definition of functions from the corresponding methods ---
#####---------------------------------------------------------------------------
class _frommethod(object):
    """Defines functions from existing MaskedArray methods.
:ivar _methodname (String): Name of the method to transform.
    """
    def __init__(self, methodname):
        self._methodname = methodname
        self.__doc__ = self.getdoc()
    def getdoc(self):
        "Returns the doc of the function (from the doc of the method)."
        try:
            return getattr(TimeSeries, self._methodname).__doc__
        except:
            return "???"    
    #
    def __call__ (self, caller, *args, **params):
        if hasattr(caller, self._methodname):
            method = getattr(caller, self._methodname)
            # If method is not callable, it's a property, and don't call it
            if hasattr(method, '__call__'):
                return method.__call__(*args, **params)
            return method
        method = getattr(fromnumeric.asarray(caller), self._methodname)
        try:
            return method(*args, **params)
        except SystemError:
            return getattr(numpy,self._methodname).__call__(caller, *args, **params)
#............................
day_of_week = _frommethod('day_of_week')
day_of_year = _frommethod('day_of_year')
year = _frommethod('year')
quarter = _frommethod('quarter')
month = _frommethod('month')
day = _frommethod('day')
hour = _frommethod('hour')
minute = _frommethod('minute')
second = _frommethod('second')
#
##### ---------------------------------------------------------------------------
#---- ... Additional methods ...
##### ---------------------------------------------------------------------------
def tofile(self, output, sep='\t', format_dates=None):
    """Writes the TimeSeries to a file.

:Parameters:
    - `output` (String) : Name or handle of the output file.
    - `sep` (String) : Column separator *['\t']*.
    - `format` (String) : Data format *['%s']*.
    """
    if not hasattr(output, 'writeline'):
        ofile = open(output,'w')
    else:
        ofile = output
    if format_dates is None:
        format_dates = self.dates[0].default_fmtstr()
    oformat = "%%s%s%s" % (sep,format_dates)
    for (_dates,_data) in numpy.broadcast(self._dates.ravel().asstrings(), 
                                          filled(self)):
        ofile.write('%s\n' % sep.join([oformat % (_dates, _data) ]))
    ofile.close()
TimeSeries.tofile = tofile

#............................................
def asrecords(series):
    """Returns the masked time series as a recarray.
Fields are `_dates`, `_data` and _`mask`.
        """
    desctype = [('_dates',int_), ('_series',series.dtype), ('_mask', bool_)]
    flat = series.ravel()
    _dates = numeric.asarray(flat._dates)
    if flat.size > 0:
        return recfromarrays([_dates, flat._data, getmaskarray(flat)],
                             dtype=desctype,
                             shape = (flat.size,),  
                             )
    else:
        return recfromarrays([[], [], []], dtype=desctype, 
                             shape = (flat.size,),  
                             )
TimeSeries.asrecords = asrecords

def flatten(series):
    """Flattens a (multi-) time series to 1D series."""
    shp_ini = series.shape
    # Already flat time series....
    if len(shp_ini) == 1:
        return series
    # Folded single time series ..
    newdates = series._dates.ravel()
    if series._dates.size == series._series.size:
        newshape = (series._series.size,)
    else:
        newshape = (numeric.asarray(shp_ini[:-1]).prod(), shp_ini[-1])
    newseries = series._series.reshape(newshape)
    return time_series(newseries, newdates)
TimeSeries.flatten = flatten



#####---------------------------------------------------------------------------
#---- --- Archiving ---
#####---------------------------------------------------------------------------
def _tsreconstruct(baseclass, datesclass, baseshape, basetype, fill_value):
    """Internal function that builds a new TimeSeries from the information stored
in a pickle."""
#    raise NotImplementedError,"Please use timeseries.archive/unarchive instead."""
    _series = ndarray.__new__(ndarray, baseshape, basetype)
    _dates = ndarray.__new__(datesclass, baseshape, int_)
    _mask = ndarray.__new__(ndarray, baseshape, bool_)
    return baseclass.__new__(baseclass, _series, dates=_dates, mask=_mask, 
                             dtype=basetype, fill_value=fill_value)
#    
def _tsgetstate(a):
    "Returns the internal state of the TimeSeries, for pickling purposes."
#    raise NotImplementedError,"Please use timeseries.archive/unarchive instead."""
    records = a.asrecords()
    state = (1,
             a.shape, 
             a.dtype,
             a.freq,
             records.flags.fnc,
             a.fill_value,
             records
             )
    return state
#    
def _tssetstate(a, state):
    """Restores the internal state of the TimeSeries, for pickling purposes.
`state` is typically the output of the ``__getstate__`` output, and is a 5-tuple:

    - class name
    - a tuple giving the shape of the data
    - a typecode for the data
    - a binary string for the data
    - a binary string for the mask.
    """
    (ver, shp, typ, frq, isf, flv, rec) = state
    a.fill_value = flv
    a._dates = a._dates.__class__(rec['_dates'], freq=frq)
    (a._dates).__tostr = None
    _data = rec['_series'].view(typ)
    _mask = rec['_mask'].view(MA.MaskType)
    a._series = masked_array(_data, mask=_mask)
#    a._data.shape = shp
#    a._dates.shape = shp
#    a._mask = rec['_mask'].view(MA.MaskType)
#    a._mask.shape = shp
#        
def _tsreduce(a):
    """Returns a 3-tuple for pickling a MaskedArray."""
    return (_tsreconstruct,
            (a.__class__, a.dates.__class__, (0,), 'b', -9999),
            a.__getstate__())
#    
TimeSeries.__getstate__ = _tsgetstate
TimeSeries.__setstate__ = _tssetstate
TimeSeries.__reduce__ = _tsreduce
#TimeSeries.__dump__ = dump
#TimeSeries.__dumps__ = dumps


##### -------------------------------------------------------------------------
#---- --- TimeSeries creator ---
##### -------------------------------------------------------------------------
def time_series(data, dates=None, freq=None, observed=None,
                start_date=None, end_date=None, length=None, include_last=True,
                mask=nomask,
                dtype=None, copy=False, fill_value=None,
                keep_mask=True, small_mask=True, hard_mask=False):
    """Creates a TimeSeries object
    
:Parameters:
    `dates` : ndarray
        Array of dates.
    `data` : 
        Array of data.
    """
    data = numeric.array(data, copy=False, subok=True)
    if dates is None:
        dshape = data.shape
        if len(dshape) > 0:
            dates = date_array(start_date=start_date, end_date=end_date,
                               length=dshape[0], include_last=include_last, 
                               freq=freq) 
        else:
            dates = date_array([], freq=freq)   
    elif not isinstance(dates, DateArray):
        dates = date_array(dlist=dates, freq=freq)
    return TimeSeries(data=data, dates=dates, mask=mask, observed=observed,
                      copy=copy, dtype=dtype, fill_value=fill_value,
                      keep_mask=keep_mask, small_mask=small_mask, 
                      hard_mask=hard_mask,)


def isTimeSeries(series):
    "Returns whether the series is a valid TimeSeries object."
    return isinstance(series, TimeSeries)

tsmasked = TimeSeries(masked,dates=Date('D',1))

##### --------------------------------------------------------------------------
#---- ... Additional functions ...
##### --------------------------------------------------------------------------
def mask_period(data, start_date=None, end_date=None, 
                inside=True, include_edges=True, inplace=True):
    """Returns x as an array masked where dates fall outside the selection period,
as well as where data are initially missing (masked).

:Parameters:
    `data` : Timeseries
        Data to process
    `start_date` : Date *[None]*
        Starting date. If None, uses the first date.
    `end_date` : Date *[None]*
        Ending date. If None, uses the last date.
    `inside` : Boolean *[True]*
        Whether the dates inside the range should be masked. If not, masks outside.
    `include_edges` : Boolean *[True]*
        Whether the starting and ending dates should be masked.
    `inplace` : Boolean *[True]*
        Whether the data mask should be modified in place. If not, returns a new
        TimeSeries.
"""
    if not isTimeSeries(data):
        raise ValueError,"Data should be a valid TimeSeries!"
    # Check the starting date ..............
    if start_date is None:
        start_date = data._dates[0]
    elif isinstance(start_date, str):
        start_date = Date(data.freq, string=start_date)
    elif not isinstance(start_date, Date):
        raise DateError,"Starting date should be a valid Date object!"
    start_date = max(start_date, data.dates[0])
    # Check the ending date ................
    if end_date is None:
        end_date = data._dates[-1]
    elif isinstance(end_date, str):
        end_date = Date(data.freq, string=end_date)
    elif not isinstance(end_date, Date):
        raise DateError,"Starting date should be a valid Date object!"
    end_date = min(end_date, data.dates[-1])
    # Constructs the selection mask .........
    if inside:
        if include_edges:
            selection = (data.dates >= start_date) & (data.dates <= end_date)
        else:
            selection = (data.dates > start_date) & (data.dates < end_date)
    else:
        if include_edges:
            selection = (data.dates <= start_date) | (data.dates >= end_date)
        else:
            selection = (data.dates < start_date) | (data.dates > end_date)
    # Process the data:
    if inplace:
        if data._mask is nomask:
            data._mask = selection
        else:
            data._mask += selection
    else:
        return TimeSeries(data, mask=selection, keep_mask=True)
    return data

def mask_inside_period(data, start_date=None, end_date=None, 
                       include_edges=True, inplace=True):
    """Masks values falling inside a given range of dates."""
    return mask_period(data, start_date=start_date, end_date=end_date, 
                       inside=True, include_edges=include_edges, inplace=inplace)
def mask_outside_period(data, start_date=None, end_date=None, 
                       include_edges=True, inplace=True):
    """Masks values falling outside a given range of dates."""
    return mask_period(data, start_date=start_date, end_date=end_date, 
                       inside=False, include_edges=include_edges, inplace=inplace)
    
#...............................................................................
def compressed(series):
    """Suppresses missing values from a time series."""
    if series._mask is nomask:
        return series
    if series.ndim == 1:
        keeper = ~(series._mask)
    elif series.ndim == 2:
        # Both dates and data are 2D: ravel first
        if series._dates.ndim == 2:
            series = series.ravel()
            keeper = ~(series._mask)
        # a 2D series: suppress the rows (dates are in columns)
        else:
            keeper = ~(series._mask.any(-1))
    else:
        raise NotImplementedError
    return series[keeper]
TimeSeries.compressed = compressed
#...............................................................................
def adjust_endpoints(a, start_date=None, end_date=None):
    """Returns a TimeSeries going from `start_date` to `end_date`.
    If `start_date` and `end_date` both fall into the initial range of dates, 
    the new series is NOT a copy.
    """
    # Series validity tests .....................
    if not isinstance(a, TimeSeries):
        raise TypeError,"Argument should be a valid TimeSeries object!"
    if a.freq == 'U':
        raise TimeSeriesError, \
            "Cannot adjust a series with 'Undefined' frequency."
    if not a.dates.isvalid():
        raise TimeSeriesError, \
            "Cannot adjust a series with missing or duplicated dates."
    # Flatten the series if needed ..............
    a = a.flatten()
    shp_flat = a.shape
    # Dates validity checks .,...................
    msg = "%s should be a valid Date object! (got %s instead)"
    if a.dates.size >= 1:
        (dstart, dend) = a.dates[[0,-1]]
    else:
        (dstart, dend) = (None, None)
    # Skip the empty series case
    if dstart is None and (start_date is None or end_date is None):
        raise TimeSeriesError, "Both start_date and end_date must be specified"+\
                               " to adjust endpoints of a zero length series!" 
    #....
    if start_date is None: 
        start_date = dstart
        start_lag = 0
    else:
        if not isinstance(start_date, Date):
            raise TypeError, msg % ('start_date', type(start_date))
        if dstart is not None:
            start_lag = start_date - dstart
        else:
            start_lag = start_date
    #....          
    if end_date is None: 
        end_date = dend
        end_lag = 0
    else:
        if not isinstance(end_date, Date):
            raise TypeError, msg % ('end_date', type(end_date))
        if dend is not None:
            end_lag = end_date - dend    
        else:
            end_lag = end_date
    # Check if the new range is included in the old one
    if start_lag >= 0:
        if end_lag == 0:
            return a[start_lag:]
        elif end_lag < 0:
            return a[start_lag:end_lag]
    # Create a new series .......................
    newdates = date_array(start_date=start_date, end_date=end_date)
    
    newshape = list(shp_flat)
    newshape[0] = len(newdates)
    newshape = tuple(newshape)
    
    newseries = numeric.empty(newshape, dtype=a.dtype).view(type(a))
    newseries.__setmask__(numeric.ones(newseries.shape, dtype=bool_))
    newseries._dates = newdates
    if dstart is not None:
        start_date = max(start_date, dstart)
        end_date = min(end_date, dend) + 1
        newseries[start_date:end_date] = a[start_date:end_date]
    newseries.copy_attributes(a)
    return newseries
#....................................................................
def align_series(*series, **kwargs):
    """Aligns several TimeSeries, so that their starting and ending dates match.
    Series are resized and filled with mased values accordingly.
    
    The function accepts two extras parameters:
    - `start_date` forces the series to start at that given date,
    - `end_date` forces the series to end at that given date.
    By default, `start_date` and `end_date` are set to the smallest and largest
    dates respectively.    
    """
    if len(series) < 2:
        return series  
    unique_freqs = numpy.unique([x.freqstr for x in series])
    common_freq = _compare_frequencies(*series)
    valid_states = [x.isvalid() for x in series]
    if not numpy.all(valid_states):
        raise TimeSeriesError, \
            "Cannot adjust a series with missing or duplicated dates."
    
    start_date = kwargs.pop('start_date', 
                            min([x.start_date for x in series 
                                     if x.start_date is not None]))
    if isinstance(start_date,str):
        start_date = Date(common_freq, string=start_date)
    end_date = kwargs.pop('end_date', 
                          max([x.end_date for x in series
                                   if x.end_date is not None]))
    if isinstance(end_date,str):
        end_date = Date(common_freq, string=end_date)
    
    return [adjust_endpoints(x, start_date, end_date) for x in series]
aligned = align_series
#....................................................................
def convert(series, freq, func='auto', position='END'):
    """Converts a series to a frequency.
       
    When converting to a lower frequency, func is a function that acts
    on a 1-d array and returns a scalar or 1-d array. func should handle
    masked values appropriately. If func is "auto", then an
    appropriate function is determined based on the observed attribute
    of the series. If func is None, then a 2D array is returned, where each
    column represents the values appropriately grouped into the new frequency.
    interp and position will be ignored in this case.
        
    When converting to a higher frequency, position is 'START' or 'END'
    and determines where the data point is in each period (eg. if going
    from monthly to daily, and position is 'END', then each data point is
    placed at the end of the month).
    """
    if not isinstance(series,TimeSeries):
        raise TypeError, "The argument should be a valid TimeSeries!"
    if not series.isvalid():
        raise TimeSeriesError, \
            "Cannot adjust a series with missing or duplicated dates."
    
    if position.upper() not in ('END','START'): 
        raise ValueError("invalid value for position argument: (%s)",str(position))
    
    toFreq = corelib.check_freq(freq)
    fromFreq = series.freq
    start_date = series._dates[0]
    
    if fromFreq == toFreq:
        return series.copy()
    
    if series.size == 0:
        return TimeSeries(series, freq=toFreq, 
                          start_date=start_date.asfreq(toFreq))
    if func == 'auto':
        func = corelib.obs_dict[series.observed]

    tempData = series._series.filled()
    tempMask = getmaskarray(series)

    cRetVal = cseries.convert(tempData, fromFreq, toFreq, position, 
                              int(start_date), tempMask)
    _values = cRetVal['values']
    _mask = cRetVal['mask']
    _startindex = cRetVal['startindex']
    start_date = Date(freq=toFreq, value=_startindex)
        
    tempData = masked_array(_values, mask=_mask)

    if tempData.ndim == 2 and func is not None:
        tempData = MA.apply_along_axis(func, -1, tempData)

    newseries = tempData.view(type(series))
    newseries._dates = date_array(start_date=start_date, length=len(newseries),
                                  freq=toFreq)
    newseries.copy_attributes(series)
    return newseries

def group_byperiod(series, freq, func='auto', position='END'):
    """Converts a series to a frequency, without any processing.
    """
    return convert(series, freq, func=None, position='END')

TimeSeries.convert = convert
TimeSeries.group_byperiod = group_byperiod

#...............................................................................
def tshift(series, nper, copy=True):
    """Returns a series of the same size as `series`, with the same
start_date and end_date, but values shifted by `nper`. This is useful
for doing things like calculating a percentage change.
Eg. pct_change = 100 * (series/tshift(series, -1, copy=False) - 1)
Note: By default the data is copied, but if you are using the result in
a way that is going to make a copy anyway (like the above example) then
you may want to bypass copying the data.

:Parameters:
    - `series` (TimeSeries) : TimeSeries object to shift
    - `nper` (int) : number of periods to shift. Negative numbers
      shift values to the right, positive to the left
    - `copy` (boolean, *[True]*) : copies the data if True, returns
      a view if False.

:Example:
>>> series = tseries.time_series([0,1,2,3], start_date=tdates.Date(freq='A', year=2005))
>>> series
timeseries(data  = [0 1 2 3],
           dates = [2005 ... 2008],
           freq  = A)
>>> tshift(series, -1)
timeseries(data  = [-- 0 1 2],
           dates = [2005 ... 2008],
           freq  = A)
"""
    #Backup series attributes
    options = _attrib_dict(series)
    newdata = masked_array(numeric.empty(series.shape, dtype=series.dtype), 
                           mask=True)
    if copy:
        inidata = series._series.copy()
    else:
        inidata = series._series
    if nper < 0:
        nper = max(-len(series), nper)
        newdata[-nper:] = inidata[:nper]
    elif nper > 0:
        nper = min(len(series), nper)
        newdata[:-nper] = inidata[nper:]
    else:
        newdata = inidata
    newseries = newdata.view(type(series))
    newseries._dates = series._dates
    newseries.copy_attributes(series)
    return newseries
TimeSeries.tshift = tshift
#...............................................................................
def fill_missing_dates(data, dates=None, freq=None,fill_value=None):
    """Finds and fills the missing dates in a time series.
The data corresponding to the initially missing dates are masked, or filled to 
`fill_value`.
    
:Parameters:
    `data`
        Initial array of data.
    `dates` 
        Initial array of dates. 
    `freq` : float *[None]*
        New date resolutions. If *None*, the initial resolution is used instead.
    `fill_value` : float *[None]*
        Default value for missing data. If None, the data are just masked.
    """
    (freq, freqstr) = corelib.check_freq(freq), corelib.check_freqstr(freq)
    if freqstr == 'U':
        raise ValueError,\
              "Unable to define a proper date resolution (found %s)." % freqstr
    if dates is None:
        if not isTimeSeries(data):
            raise InsufficientDateError
        dates = data._dates
        freq = dates.freq
        datad = data._series._data
        datam = data._series._mask
#        if fill_value is None:
#            fill_value = data._fill_value
    elif not isinstance(dates, DateArray):
        dates = DateArray(dates, freq)
        if isinstance(data, MaskedArray):
            datad = data._data
            datam = data._mask
        else:
            datad = data
            datam = nomask
    dflat = dates.asfreq(freq).ravel()
    n = len(dflat)
    if not dflat.has_missing_dates():
        return time_series(data, dflat)
    # ...and now, fill it ! ......
    (tstart, tend) = dflat[[0,-1]]
    newdates = date_array(start_date=tstart, end_date=tend, include_last=True)
    nsize = newdates.size
    #.............................
    # Get the steps between consecutive data. 
    delta = dflat.get_steps()-1
    gap = delta.nonzero()
    slcid = numpy.r_[[0,], numpy.arange(1,n)[gap], [n,]]
    oldslc = numpy.array([slice(i,e) 
                          for (i,e) in numpy.broadcast(slcid[:-1],slcid[1:])])
    addidx = delta[gap].astype(int_).cumsum()
    newslc = numpy.r_[[oldslc[0]], 
                      [slice(i+d,e+d) for (i,e,d) in \
                           numpy.broadcast(slcid[1:-1],slcid[2:],addidx)] 
                     ]
    #.............................
    # Just a quick check
    vdflat = numeric.asarray(dflat)
    vnewdates = numeric.asarray(newdates)
    for (osl,nsl) in zip(oldslc,newslc):
        assert numpy.equal(vdflat[osl],vnewdates[nsl]).all(),\
            "Slicing mishap ! Please check %s (old) and %s (new)" % (osl,nsl)
    #.............................
    data = MA.asarray(data)
    newdatad = numeric.empty(nsize, data.dtype)
    newdatam = numeric.ones(nsize, bool_)
    #....
    if datam is nomask:
        for (new,old) in zip(newslc,oldslc):
            newdatad[new] = datad[old]
            newdatam[new] = False
    else:
        for (new,old) in zip(newslc,oldslc):
            newdatad[new] = datad[old]
            newdatam[new] = datam[old]
    newdata = MA.masked_array(newdatad, mask=newdatam, fill_value=fill_value)    
    # Get new shape ..............
    if data.ndim == 1:
        nshp = (newdates.size,)
    else:
        nshp = tuple([-1,] + list(data.shape[1:]))
    return time_series(newdata.reshape(nshp), newdates)
#...............................................................................
def stack(*series):
    """performs a column_stack on the data from each series, and the
resulting series has the same dates as each individual series. All series
must be date compatible.
    
:Parameters:
    `*series` : the series to be stacked
"""
    _timeseriescompat_multiple(*series)
    return time_series(MA.column_stack(series), series[0]._dates,
                       **_attrib_dict(series[0]))
#...............................................................................
def concatenate_series(series, keep_gap=True):
    """Concatenates a sequence of series, by chronological order.
    Overlapping data are processed in a FIFO basis: the data from the first series
    of the sequence will be overwritten by the data of the second series, and so forth.
    If keep_gap is true, any gap between consecutive, non overlapping series are
    kept: the corresponding data are masked.
    """
    common_f = _compare_frequencies(*series)
    start_date = min([s.start_date for s in series if s.start_date is not None])
    end_date =   max([s.end_date for s in series if s.end_date is not None])
    newdtype = max([s.dtype for s in series])
    whichone = numeric.zeros((end_date-start_date+1), dtype=int_)
    newseries = time_series(numeric.empty((end_date-start_date+1), dtype=newdtype), 
                            dates=date_array(start_date, end_date, freq=common_f),
                            mask=True)
    newdata = newseries._data
    newmask = newseries._mask
    for (k,s) in enumerate(series):
        start = s.start_date - start_date
        end = start + len(s)
        whichone[start:end] = k+1
        newdata[start:end] = s._data
        if s._mask is nomask:
            newmask[start:end] = False
        else:
            newmask[start:end] = s._mask
    keeper = whichone.astype(bool_)
    if not keep_gap:
        newseries = newseries[keeper]
    else:
        newdata[~keeper] = 0
    return newseries
#...............................................................................
def empty_like(series):
    result = N.empty_like(series).view(type(series))
    result._dates = series._dates
    result._mask = series._mask
    return result                           
    
################################################################################
if __name__ == '__main__':
    from maskedarray.testutils import assert_equal, assert_array_equal
    import numpy as N
    
    if 1:
        dlist = ['2007-01-%02i' % i for i in range(1,11)]
        dates = date_array_fromlist(dlist)
        data = masked_array(numeric.arange(10), mask=[1,0,0,0,0]*2, dtype=float_)
        ser1d = time_series(data, dlist)
        
        serfolded = ser1d.reshape((5,2))
        assert_equal(serfolded._dates.shape, (5,2))
        assert_equal(serfolded[0], time_series([0,1],mask=[1,0],
                                               start_date=dates[0]))
        assert_equal(serfolded[:,0], 
                     time_series(ser1d[::2], dates=dates[::2]))
        sertrans = serfolded.transpose()
        assert_equal(sertrans.shape, (2,5))
        
    if 1:
        hodie = today('D')
        ser2d = time_series(N.arange(10).reshape(5,2), start_date=hodie,
                            mask=[[1,1],[0,0],[0,0],[0,0],[0,0]])
        try:
            ser2d_transpose = ser2d.transpose()
        except TimeSeriesError:
            pass
    if 1:
        hodie = today('D')
        ser3d = time_series(N.arange(30).reshape(5,3,2), start_date=hodie,)
        try:
            ser3d_transpose = ser3d.transpose()
        except TimeSeriesError:
            pass
        assert_equal(ser3d.transpose(0,2,1).shape, (5,2,3))
        
    if 1:        
        data = dates
    
        series = time_series(data, dates)
        assert(isinstance(series, TimeSeries))
        assert_equal(series._dates, dates)
        assert_equal(series._data, data)
        assert_equal(series.freqstr, 'D')
        
        series[5] = MA.masked
        
        # ensure that series can be represented by a string after masking a value
        # (there was a bug before that prevented this from working when using a 
        # DateArray for the data)
        strrep = str(series)
