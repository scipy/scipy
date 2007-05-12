"""
The `TimeSeries` class provides  a base for the definition of time series.
A time series is defined here as the combination of two arrays:

    - an array storing the time information (as a `DateArray` instance);
    - an array storing the data (as a `MaskedArray` instance.)

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
from numpy.core.multiarray import dtype
import numpy.core.fromnumeric as fromnumeric
import numpy.core.numeric as numeric
import numpy.core.umath as umath
from numpy.core.records import recarray
from numpy.core.records import fromarrays as recfromarrays

import maskedarray as MA
from maskedarray import MaskedArray, MAError, masked, nomask, \
    filled, getmask, getmaskarray, hsplit, make_mask_none, mask_or, make_mask, \
    masked_array

import tcore as corelib
import const as _c

import tdates
from tdates import DateError, InsufficientDateError
from tdates import Date, isDate, DateArray, isDateArray, \
    date_array, date_array_fromlist, date_array_fromrange, thisday, today, \
    check_freq, check_freq_str

import cseries



__all__ = [
'TimeSeriesError','TimeSeriesCompatibilityError','TimeSeries','isTimeSeries',
'time_series', 'tsmasked',
'mask_period','mask_inside_period','mask_outside_period','compressed',
'adjust_endpoints','align_series','aligned','convert','group_byperiod',
'pct','tshift','fill_missing_dates', 'split', 'stack', 'concatenate_series',
'empty_like',
'day_of_week','day_of_year','day','month','quarter','year',
'hour','minute','second',
'tofile','asrecords','flatten', 'check_observed',
           ]


#####---------------------------------------------------------------------------
#---- --- Observed options ---
#####---------------------------------------------------------------------------
fmtobs_dict = {'UNDEFINED': ['UNDEF','UNDEFINED',None],
               'BEGINNING': ['BEGIN','BEGINNING'],
               'ENDING': ['END','ENDING'],
               'AVERAGED': ['AVERAGE','AVERAGED','MEAN'],
               'SUMMED': ['SUM','SUMMED'],
               'MAXIMUM': ['MAX','MAXIMUM','HIGH'],
               'MINIMUM': ['MIN','MINIMUM','LOW']}

obs_dict = {"UNDEFINED":None,
            "BEGINNING": corelib.first_unmasked_val,
            "ENDING": corelib.last_unmasked_val,
            "AVERAGED": MA.average,
            "SUMMED": MA.sum,
            "MAXIMUM": MA.maximum,
            "MINIMUM": MA.minimum,
            }

alias_obs_dict = {}
for ob, aliases in fmtobs_dict.iteritems():
    for al in aliases:
        alias_obs_dict[al] = obs_dict[ob]
obs_dict.update(alias_obs_dict)
fmtobs_revdict = corelib.reverse_dict(fmtobs_dict)

def fmtObserv(obStr):
    "Converts a possible 'Observed' string into acceptable values."
    if obStr is None:
        return fmtobs_revdict[None]
    elif obStr.upper() in fmtobs_revdict:
        return fmtobs_revdict[obStr.upper()]
    else:
        raise ValueError("Invalid value for observed attribute: %s " % str(obStr))
check_observed = fmtObserv

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
                    result = result.view(type(self.obj))
                    result._dates = _dates
#                    result = TimeSeries(result, dates=_dates)
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
                freq=None, observed=None, start_date=None, length=None,
                dtype=None, copy=False, fill_value=None, subok=True,
                keep_mask=True, small_mask=True, hard_mask=False, **options):
        maparms = dict(copy=copy, dtype=dtype, fill_value=fill_value,subok=subok,
                       keep_mask=keep_mask, small_mask=small_mask,
                       hard_mask=hard_mask,)
        _data = MaskedArray(data, mask=mask, **maparms)
        # Get the frequency ..........................
        freq = check_freq(freq)
        # Get the dates ..............................
        if dates is None:
            newdates = getattr(data, '_dates', None)
        else:
            newdates = dates
        if newdates is not None:
            if not hasattr(newdates, 'freq'):
                newdates = date_array(dlist=dates, freq=freq)
            if freq != _c.FR_UND and newdates.freq != freq:
                newdates = newdates.asfreq(freq)
        else:
            dshape = _data.shape
            if len(dshape) > 0:
                if length is None:
                    length = dshape[0]
                newdates = date_array(start_date=start_date, length=length,
                                      freq=freq)
            else:
                newdates = date_array([], freq=freq)
        # Get observed ...............................
        observed = getattr(data, 'observed', fmtObserv(observed))
        # Get the data ...............................
        if newdates._unsorted is not None:
            _data = _data[newdates._unsorted]
        if not subok or not isinstance(_data,TimeSeries):
            _data = _data.view(cls)
        if _data is masked:
            assert(numeric.size(newdates)==1)
            return _data.view(cls)
        assert(_datadatescompat(_data,newdates))
        _data._dates = newdates
        if _data._dates.size == _data.size and _data.ndim > 1:
            _data._dates.shape = _data.shape
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
        "Returns the series as a regular masked array."
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
timeseries(
 %(data)s,
           dates =
 %(time)s,
           freq  = %(freq)s)
"""
        desc_short = """\
timeseries(%(data)s,
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
        _dates = self._dates
        dsize = _dates.size
        if dsize == 0:
            return None
        elif dsize == 1:
            return _dates[0]
        else:
            return Date(self.freq, _dates.flat[0])
    @property
    def end_date(self):
        """Returns the last date of the series."""
        _dates = self._dates
        dsize = _dates.size
        if dsize == 0:
            return None
        elif dsize == 1:
            return _dates[-1]
        else:
            return Date(self.freq, _dates.flat[-1])

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
    
    def split(self):
        """Split a multiple series into individual columns."""
        if self.ndim == 1:
            return [self]
        else:
            n = self.shape[1]
            arr = hsplit(self, n)[0]
            return [self.__class__(numpy.squeeze(a), 
                                   self._dates, 
                                   **_attrib_dict(self)) for a in arr]        
    
    def filled(self, fill_value=None):
        """Returns an array of the same class as `_data`,
 with masked values filled with `fill_value`.
Subclassing is preserved.

If `fill_value` is None, uses self.fill_value.
        """
        result = self._series.filled(fill_value=fill_value).view(type(self))
        result._dates = self._dates
        result.copy_attributes(self)
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
    #......................................................
    # Pickling
    def __getstate__(self):
        "Returns the internal state of the TimeSeries, for pickling purposes."
    #    raise NotImplementedError,"Please use timeseries.archive/unarchive instead."""
        state = (1,
                 self.shape,
                 self.dtype,
                 self.flags.fnc,
                 self._data.tostring(),
                 getmaskarray(self).tostring(),
                 self._fill_value,
                 self._dates.shape,
                 numeric.asarray(self._dates).tostring(),
                 self.freq,
                 )
        return state
    #
    def __setstate__(self, state):
        """Restores the internal state of the TimeSeries, for pickling purposes.
    `state` is typically the output of the ``__getstate__`` output, and is a 5-tuple:

        - class name
        - a tuple giving the shape of the data
        - a typecode for the data
        - a binary string for the data
        - a binary string for the mask.
        """
        (ver, shp, typ, isf, raw, msk, flv, dsh, dtm, frq) = state
        super(TimeSeries, self).__setstate__((ver, shp, typ, isf, raw, msk, flv))
        self._dates.__setstate__((dsh, dtype(int_), isf, dtm))
        self._dates.freq = frq
#
    def __reduce__(self):
        """Returns a 3-tuple for pickling a MaskedArray."""
        return (_tsreconstruct,
                (self.__class__, self._baseclass,
                 self.shape, self._dates.shape, self.dtype, self._fill_value),
                self.__getstate__())

def _tsreconstruct(genclass, baseclass, baseshape, dateshape, basetype, fill_value):
    """Internal function that builds a new TimeSeries from the information stored
    in a pickle."""
    #    raise NotImplementedError,"Please use timeseries.archive/unarchive instead."""
    _series = ndarray.__new__(baseclass, baseshape, basetype)
    _dates = ndarray.__new__(DateArray, dateshape, int_)
    _mask = ndarray.__new__(ndarray, baseshape, bool_)
    return genclass.__new__(genclass, _series, dates=_dates, mask=_mask,
                            dtype=basetype, fill_value=fill_value)

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

split = _frommethod('split')

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
            if length is None:
                length = dshape[0]
        if len(dshape) > 0:
            dates = date_array(start_date=start_date, end_date=end_date,
                               length=length, include_last=include_last,
                               freq=freq)
        else:
            dates = date_array([], freq=freq)
    elif not isinstance(dates, DateArray):
        dates = date_array(dlist=dates, freq=freq)
    if dates._unsorted is not None:
        idx = dates._unsorted
        data = data[idx]
        if mask is not nomask:
            mask = mask[idx]
        dates._unsorted = None
    return TimeSeries(data=data, dates=dates, mask=mask, 
                      observed=observed, copy=copy, dtype=dtype, 
                      fill_value=fill_value, keep_mask=keep_mask, 
                      small_mask=small_mask, hard_mask=hard_mask,)


def isTimeSeries(series):
    "Returns whether the series is a valid TimeSeries object."
    return isinstance(series, TimeSeries)

tsmasked = TimeSeries(masked,dates=DateArray(Date('D',1)))

##### --------------------------------------------------------------------------
#---- ... Additional functions ...
##### --------------------------------------------------------------------------
def mask_period(data, period=None, start_date=None, end_date=None,
                inside=True, include_edges=True, inplace=False):
    """Returns x as an array masked where dates fall outside the selection period,
as well as where data are initially missing (masked).

:Parameters:
    data : Timeseries
        Data to process
    period : Sequence
        A sequence of (starting date, ending date).
    start_date : string/Date *[None]*
        Starting date. If None, uses the first date of the series.
    end_date : string/Date *[None]*
        Ending date. If None, uses the last date of the series.
    inside : Boolean *[True]*
        Whether the dates inside the range should be masked. If not, masks outside.
    include_edges : Boolean *[True]*
        Whether the starting and ending dates should be masked.
    inplace : Boolean *[True]*
        Whether the data mask should be modified in place. If not, returns a new
        TimeSeries.
    """
    data = masked_array(data, subok=True, copy=not inplace)
    if not isTimeSeries(data):
        raise ValueError,"Data should be a valid TimeSeries!"
    dates = data._dates
    if dates.ndim == 1:
        dates_lims = dates[[0,-1]]
    else:
        dates_lims = dates.ravel()[[0,-1]]
    # Check the period .....................
    if period is not None:
        if isinstance(period, (tuple, list, ndarray)):
            (start_date, end_date) = (period[0], period[-1])
        else:
            (start_date, end_date) = (period, start_date)
    # Check the starting date ..............
    if start_date is None:
        start_date = dates_lims[0]
    elif isinstance(start_date, str):
        start_date = Date(data.freq, string=start_date)
    elif not isinstance(start_date, Date):
        raise DateError,"Starting date should be a valid Date object!"
    # Check the ending date ................
    if end_date is None:
        end_date = dates_lims[-1]
    elif isinstance(end_date, str):
        end_date = Date(data.freq, string=end_date)
    elif not isinstance(end_date, Date):
        raise DateError,"Starting date should be a valid Date object!"
    # Constructs the selection mask .........
    dates = data.dates
    if inside:
        if include_edges:
            selection = (dates >= start_date) & (dates <= end_date)
        else:
            selection = (dates > start_date) & (dates < end_date)
    else:
        if include_edges:
            selection = (dates <= start_date) | (dates >= end_date)
        else:
            selection = (dates < start_date) | (dates > end_date)
    data[selection] = masked
    return data

def mask_inside_period(data, start_date=None, end_date=None,
                       include_edges=True, inplace=False):
    """Masks values falling inside a given range of dates."""
    return mask_period(data, start_date=start_date, end_date=end_date,
                       inside=True, include_edges=include_edges, inplace=inplace)
def mask_outside_period(data, start_date=None, end_date=None,
                       include_edges=True, inplace=False):
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
def _convert1d(series, freq, func='auto', position='END'):
    """Converts a series to a frequency. Private function called by convert

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

    toFreq = check_freq(freq)
    fromFreq = series.freq

    if toFreq == _c.FR_UND:
        raise TimeSeriesError, \
            "Cannot convert a series to UNDEFINED frequency."

    if fromFreq == _c.FR_UND:
        raise TimeSeriesError, \
            "Cannot convert a series with UNDEFINED frequency."

    if not series.isvalid():
        raise TimeSeriesError, \
            "Cannot adjust a series with missing or duplicated dates."

    if position.upper() not in ('END','START'):
        raise ValueError("Invalid value for position argument: (%s). "\
                         "Should be in ['END','START']," % str(position))

    start_date = series._dates[0]

    if series.size == 0:
        return TimeSeries(series, freq=toFreq,
                          start_date=start_date.asfreq(toFreq))
    if func == 'auto':
        func = obs_dict[series.observed]

    tempData = series._series.filled()
    tempMask = getmaskarray(series)

    if (tempData.size // series._dates.size) > 1:
        raise TimeSeriesError("convert works with 1D data only !")

    cRetVal = cseries.TS_convert(tempData, fromFreq, toFreq, position,
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

def convert(series, freq, func='auto', position='END'):
    """Converts a series to a frequency. Private function called by convert

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
    if series.ndim == 1:
        obj = _convert1d(series, freq, func, position)
    elif series.ndim == 2:
        base = _convert1d(series[:,0], freq, func, position)
        obj = MA.column_stack([_convert1d(m,freq,func,position)._series 
                               for m in series.split()]).view(type(series))
        obj._dates = base._dates                        
        if func is None or (func,series.observed) == ('auto','UNDEFINED'):         
            shp = obj.shape
            ncols = base.shape[-1]
            obj.shape = (shp[0], shp[-1]//ncols, ncols)
            obj = numpy.swapaxes(obj,1,2)
    return obj
        

def group_byperiod(series, freq, position='END'):
    """Converts a series to a frequency, without any processing. If the series
    has missing data, it is first filled with masked data. Duplicate values in the
    series will raise an exception.
    """
    if series.has_duplicated_dates():
        raise TimeSeriesError("The input series must not have duplicated dates!")
    elif series.has_missing_dates():
        series = fill_missing_dates(series)
    return convert(series, freq, func=None, position=position)

TimeSeries.convert = convert
TimeSeries.group_byperiod = group_byperiod

#...............................................................................
def tshift(series, nper, copy=True):
    """Returns a series of the same size as `series`, with the same
start_date and end_date, but values shifted by `nper`.

:Parameters:
    - series : (TimeSeries)
        TimeSeries object to shift
    - nper : (int)
        number of periods to shift. Negative numbers shift values to the
        right, positive to the left
    - copy : (boolean, *[True]*)
        copies the data if True, returns a view if False.

:Example:
>>> series = time_series([0,1,2,3], start_date=Date(freq='A', year=2005))
>>> series
timeseries(data  = [0 1 2 3],
           dates = [2005 ... 2008],
           freq  = A-DEC)
>>> tshift(series, -1)
timeseries(data  = [-- 0 1 2],
           dates = [2005 ... 2008],
           freq  = A-DEC)
>>> pct_change = 100 * (series/tshift(series, -1, copy=False) - 1)"""
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
def pct(series, nper=1):
    """Returns the rolling percentage change of the series.

:Parameters:
    - series : (TimeSeries)
        TimeSeries object to to calculate percentage chage for
    - nper : (int)
        number of periods for percentage change

:Example:
>>> series = time_series([2.,1.,2.,3.], start_date=Date(freq='A', year=2005))
>>> pct(series)
timeseries(data  = [-- -50.0 100.0 50.0],
           dates = [2005 ... 2008],
           freq  = A-DEC)
>>> pct(series, 2)
timeseries(data  = [-- -- 0.0 200.0],
           dates = [2005 ... 2008],
           freq  = A-DEC)"""

    newdata = masked_array(numeric.empty(series.shape, dtype=series.dtype),
                           mask=True)
    if nper < newdata.size:
        newdata[nper:] = 100*(series._series[nper:]/series._series[:-nper] - 1)
    newseries = newdata.view(type(series))
    newseries._dates = series._dates
    newseries.copy_attributes(series)
    return newseries
TimeSeries.pct = pct
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

    orig_freq = freq
    freq = check_freq(freq)

    if orig_freq is not None and freq == _c.FR_UND:
        freqstr = check_freq_str(freq)
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
    _data = newdata.reshape(nshp).view(type(data))
    _data._dates = newdates
    return _data
#    return time_series(newdata.reshape(nshp), newdates)
#...............................................................................
def stack(*series):
    """Performs a column_stack on the data from each series, and the
resulting series has the same dates as each individual series. All series
must be date compatible.

:Parameters:
    `*series` : the series to be stacked
"""
    _timeseriescompat_multiple(*series)
    return time_series(MA.column_stack(series), series[0]._dates,
                       **_attrib_dict(series[0]))
#...............................................................................
def concatenate_series(*series, **kwargs):
    """Concatenates a sequence of series, by chronological order.
    Overlapping data are processed in a FIFO basis: the data from the first series
    of the sequence will be overwritten by the data of the second series, and so forth.
    If keep_gap is true, any gap between consecutive, non overlapping series are
    kept: the corresponding data are masked.
    """
    
    keep_gap = kwargs.pop('keep_gap', True)
    if len(kwargs) > 0:
        raise KeyError("unrecognized keyword: %s" % list(kwargs)[0])
    
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
    """Returns an empty series with the same dtype, mask and dates as series."""
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

    if 0:
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
    
    if 0:
        series = time_series(numpy.arange(1,501),
                             start_date=Date('D', string='2007-01-01'))
        mseries = convert(series, 'M')
        aseries = convert(mseries, 'A')
        (freq, func, position) = ('A', None, 'END')
        
        tmp = mseries[:,0].convert('A')
        aseries = MA.concatenate([_convert1d(m,'A')._series for m in mseries.split()],
                                 axis=-1).view(type(series))
        aseries._dates = tmp._dates                                 
        shp = aseries.shape
        aseries.shape = (shp[0], shp[-1]//tmp.shape[-1], tmp.shape[-1])
        numpy.swapaxes(aseries,1,2)
    
    if 1:
        series = time_series(N.arange(124).reshape(62,2), 
                             start_date=Date(freq='d', year=2005, month=7, day=1))
        assert_equal(series.convert('M',sum), [[930,961],[2852,2883]])