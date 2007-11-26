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
from numpy import bool_, complex_, float_, int_, object_
from numpy import dtype
import numpy.core.numeric as numeric
import numpy.core.umath as umath
from numpy.core.records import recarray
from numpy.core.records import fromarrays as recfromarrays

import maskedarray
from maskedarray import MaskedArray, MAError, masked, nomask, \
    filled, getmask, getmaskarray, hsplit, make_mask_none, mask_or, make_mask, \
    masked_array

import const as _c

import dates
from dates import DateError, InsufficientDateError
from dates import Date, isDate, DateArray, isDateArray, \
    date_array, date_array_fromlist, date_array_fromrange, now, \
    check_freq, check_freq_str

import cseries

__all__ = [
'TimeSeriesError','TimeSeriesCompatibilityError','TimeSeries','isTimeSeries',
'time_series', 'tsmasked',
'adjust_endpoints','align_series','align_with','aligned','asrecords',
'compressed','concatenate', 'concatenate_series','convert',
'day_of_week','day_of_year','day',
'empty_like',
'fill_missing_dates','first_unmasked_val','flatten',
'group_byperiod',
'hour',
'last_unmasked_val',
'mask_period','mask_inside_period','mask_outside_period','minute','month',
'pct',
'quarter',
'second','split', 'stack',
'tofile','tshift',
'week',
'year',
]

def _unmasked_val(marray, x):
    "helper function for first_unmasked_val and last_unmasked_val"
    try:
        assert(marray.ndim == 1)
    except AssertionError:
        raise ValueError("array must have ndim == 1")

    idx = maskedarray.extras.flatnotmasked_edges(marray)
    if idx is None:
        return masked
    return marray[idx[x]]

def first_unmasked_val(marray):
    """Retrieve the first unmasked value in a 1d maskedarray.

*Parameters*:
    marray : {MaskedArray}
        marray must be 1 dimensional.

*Returns*:
    val : {singleton of type marray.dtype}
        first unmasked value in marray. If all values in marray are masked,
        the function returns the maskedarray.masked constant
"""
    return _unmasked_val(marray, 0)

def last_unmasked_val(marray):
    """Retrieve the last unmasked value in a 1d maskedarray.

*Parameters*:
    marray : {MaskedArray}
        marray must be 1 dimensional.

*Returns*:
    val : {singleton of type marray.dtype}
        last unmasked value in marray. If all values in marray are masked,
        the function returns the maskedarray.masked constant
"""
    return _unmasked_val(marray, 1)

#### -------------------------------------------------------------------------
#--- ... TimeSeriesError class ...
#### -------------------------------------------------------------------------
class TimeSeriesError(Exception):
    "Class for TS related errors."
    def __init__ (self, value=None):
        "Creates an exception."
        self.value = value
    def __str__(self):
        "Calculates the string representation."
        return str(self.value)
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

def _timeseriescompat(a, b, raise_error=True):
    """Checks the date compatibility of two TimeSeries object.
    Returns True if everything's fine, or raises an exception."""
    if not (hasattr(a,'freq') and hasattr(b, 'freq')):
        return True
    if a.freq != b.freq:
        if raise_error:
            raise TimeSeriesCompatibilityError('freq', a.freq, b.freq)
        else:
            return False
    elif a.start_date != b.start_date:
        if raise_error:
            raise TimeSeriesCompatibilityError('start_date',
                                               a.start_date, b.start_date)
        else:
            return False
    else:
        step_diff = a._dates.get_steps() != b._dates.get_steps()
        if (step_diff is True) or \
           (hasattr(step_diff, "any") and step_diff.any()):
            if raise_error:
                raise TimeSeriesCompatibilityError('time_steps',
                                                   a._dates.get_steps(),
                                                   b._dates.get_steps())
            else:
                return False
        elif a.shape != b.shape:
            if raise_error:
                raise TimeSeriesCompatibilityError(
                        'size', "1: %s" % str(a.shape),
                        "2: %s" % str(b.shape))
            else:
                return False
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
        raise TimeSeriesCompatibilityError('start_dates',
                                           errItems[0], errItems[1])

    if max(steps) == True:
        bad_index = [x for x, val in enumerate(steps) if val][0]
        raise TimeSeriesCompatibilityError('time_steps',
                series[0]._dates.get_steps(),
                series[bad_index]._dates.get_steps())

    if len(set(shapes)) > 1:
        errItems = tuple(set(shapes))
        raise TimeSeriesCompatibilityError('size',
                                           "1: %s" % str(errItems[0].shape),
                                           "2: %s" % str(errItems[1].shape))

    return True

def _datadatescompat(data, dates):
    """Checks the compatibility of dates and data at the creation of a
TimeSeries.

Returns True if everything's fine, raises an exception otherwise.
"""
    # If there's only 1 element, the date is a Date object, which has no
    # size...
    tsize = numeric.size(dates)
    dsize = data.size
    # Only one data
    if dsize == tsize:
        return True
    elif data.ndim > 1:
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
frequencies.
"""
    unique_freqs = numpy.unique([x.freqstr for x in series])
    try:
        common_freq = unique_freqs.item()
    except ValueError:
        raise TimeSeriesError, \
            "All series must have same frequency! (got %s instead)" % \
            unique_freqs
    return common_freq

##### ------------------------------------------------------------------------
##--- ... Time Series ...
##### ------------------------------------------------------------------------
class _tsmathmethod(object):
    """Defines a wrapper for arithmetic array methods (add, mul...).
When called, returns a new TimeSeries object, with the new series the result
of the method applied on the original series. The `_dates` part remains
unchanged.
"""
    def __init__ (self, methodname):
        self._name = methodname

    def __get__(self, obj, objtype=None):
        "Gets the calling object."
        self.obj = obj
        return self

    def __call__ (self, other, *args):
        "Execute the call behavior."
        instance = self.obj
        if isinstance(other, TimeSeries):
            compat = _timeseriescompat(instance, other, raise_error=False)
        else:
            compat = True

        func = getattr(super(TimeSeries, instance), self._name)
        if compat:
            result = func(other, *args).view(type(instance))
            result._dates = instance._dates
        else:
            _result = func(other, *args)
            if hasattr(_result, '_series'):
                result = _result._series
            else:
                result = _result
        return result

class _tsarraymethod(object):
    """Defines a wrapper for basic array methods.
When called, returns a new TimeSeries object, with the new series the result
of the method applied on the original series.
If `ondates` is True, the same operation is performed on the `_dates`.
If `ondates` is False, the `_dates` part remains unchanged.
"""
    def __init__ (self, methodname, ondates=False):
        """abfunc(fillx, filly) must be defined.
           abinop(x, filly) = x for all x to enable reduce.
        """
        self._name = methodname
        self._ondates = ondates

    def __get__(self, obj, objtype=None):
        self.obj = obj
        return self

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

When called, returns a ndarray, as the result of the method applied on the
series.
"""
    def __init__ (self, methodname):
        """abfunc(fillx, filly) must be defined.
           abinop(x, filly) = x for all x to enable reduce.
        """
        self._name = methodname

    def __get__(self, obj, objtype=None):
        self.obj = obj
        return self

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
            except IndexError:
                pass
            return result

class TimeSeries(MaskedArray, object):
    """Base class for the definition of time series.

A time series is here defined as the combination of two arrays:

    series : {MaskedArray}
        Data part
    dates : {DateArray}
        Date part

*Construction*:
    data : {array_like}
        data portion of the array. Any data that is valid for constructing a
        MaskedArray can be used here.
    dates : {DateArray}

*Other Parameters*:
    all other parameters are the same as for MaskedArray. Please see the
    documentation for the MaskedArray class in the maskedarray module
    for details.

*Notes*:
    it is typically recommended to use the `time_series` function for
    construction as it allows greater flexibility and convenience.
"""
    def __new__(cls, data, dates, mask=nomask, dtype=None, copy=False,
                fill_value=None, subok=True, keep_mask=True, hard_mask=False,
                **options):

        maparms = dict(copy=copy, dtype=dtype, fill_value=fill_value, subok=subok,
                       keep_mask=keep_mask, hard_mask=hard_mask)
        _data = MaskedArray(data, mask=mask, **maparms)

        # Get the dates ......................................................
        if not isinstance(dates, (Date, DateArray)):
            raise TypeError("The input dates should be a valid Date or " + \
                            "DateArray object (got %s instead)" % type(dates))

        # Get the data .......................................................
        if not subok or not isinstance(_data,TimeSeries):
            _data = _data.view(cls)
        if _data is masked:
            assert(numeric.size(newdates)==1)
            return _data.view(cls)
        assert(_datadatescompat(_data,dates))
        _data._dates = dates
        if _data._dates.size == _data.size:
            if _data.ndim > 1:
                current_shape = data.shape

                if dates._unsorted is not None:
                    _data.shape = (-1,)
                    _data = _data[dates._unsorted]
                    _data.shape = current_shape
                _data._dates.shape = current_shape
            elif dates._unsorted is not None:
                _data = _data[dates._unsorted]
        return _data
    #.........................................................................
    def __array_finalize__(self,obj):
        MaskedArray.__array_finalize__(self, obj)
        self._dates = getattr(obj, '_dates', DateArray([]))
        return
    #.........................................................................
    def __array_wrap__(self, obj, context=None):
        result = super(TimeSeries, self).__array_wrap__(obj, context)
        result._dates = self._dates
        return result
    #.........................................................................
    def _get_series(self):
        "Returns the series as a regular masked array."
        if self._mask.ndim == 0 and self._mask:
            return masked
        return self.view(MaskedArray)
    _series = property(fget=_get_series)
    #.........................................................................
    def __checkindex(self, indx):
        "Checks the validity of an index."
        if isinstance(indx, int):
            return (indx, indx)
        elif isinstance(indx, str):
            indx = self._dates.date_to_index(
                                    Date(self._dates.freq, string=indx))
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
        elif isTimeSeries(indx):
            indx = indx._series
        if getmask(indx) is not nomask:
            msg = "Masked arrays must be filled before they can be used " + \
                  "as indices!"
            raise IndexError, msg
        return (indx,indx)

    def __getitem__(self, indx):
        """x.__getitem__(y) <==> x[y]
Returns the item described by i. Not a copy.
"""
        (sindx, dindx) = self.__checkindex(indx)
        newdata = numeric.array(self._series[sindx], copy=False, subok=True)
        newdate = self._dates[dindx]
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
    #........................
    def __setitem__(self, indx, value):
        """x.__setitem__(i, y) <==> x[i]=y
Sets item described by index. If value is masked, masks those locations.
"""
        if self is masked:
            raise MAError, 'Cannot alter the masked element.'
        (sindx, _) = self.__checkindex(indx)
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
        """Calculates the repr representation, using masked for fill if it is
enabled. Otherwise fill with fill value.
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
    #.........................................................................
    def ids (self):
        """Return the ids of the data, dates and mask areas"""
        return (id(self._series), id(self.dates),)
    #.........................................................................
    @property
    def series(self):
        """Returns the series."""
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
        """Returns the day of month for each date in self._dates."""
        return self._dates.day
    @property
    def weekday(self):
        """Returns the day of week for each date in self._dates."""
        return self._dates.weekday
    # deprecated alias for weekday
    day_of_week = weekday
    @property
    def day_of_year(self):
        """Returns the day of year for each date in self._dates."""
        return self._dates.day_of_year
    @property
    def month(self):
        """Returns the month for each date in self._dates."""
        return self._dates.month
    @property
    def quarter(self):
        """Returns the quarter for each date in self._dates."""
        return self._dates.quarter
    @property
    def year(self):
        """Returns the year for each date in self._dates."""
        return self._dates.year
    @property
    def second(self):
        """Returns the second for each date in self._dates."""
        return self._dates.second
    @property
    def minute(self):
        """Returns the minute for each date in self._dates."""
        return self._dates.minute
    @property
    def hour(self):
        """Returns the hour for each date in self._dates."""
        return self._dates.hour
    @property
    def week(self):
        """Returns the week for each date in self._dates."""
        return self._dates.week

    days = day
    weekdays = weekday
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
        """Returns the index corresponding to a given date, as an integer."""
        return self._dates.date_to_index(date)
    #.....................................................
    def asfreq(self, freq, relation="END"):
        """Converts the dates portion of the TimeSeries to another frequency.

The resulting TimeSeries will have the same shape and dimensions as the
original series (unlike the `convert` method).

*Parameters*:
    freq : {freq_spec}
    relation : {'END', 'START'} (optional)

*Returns*:
    a new TimeSeries with the .dates DateArray at the specified frequency (the
    .asfreq method of the .dates property will be called). The data in the
    resulting series will be a VIEW of the original series.

*Notes*:
    The parameters are the exact same as for DateArray.asfreq , please see the
    __doc__ string for that method for details on the parameters and how the
    actual conversion is performed.
"""
        if freq is None: return self

        return TimeSeries(self._series,
                          dates=self._dates.asfreq(freq, relation=relation))
    #.....................................................
    def transpose(self, *axes):
        """Returns a view of the series with axes transposed

*Parameters*:
    *axes : {integers}
        the axes to swap

*Returns*:
    a VIEW of the series with axes for both the data and dates transposed

*Notes*:
    If no axes are given, the order of the axes are switches. For a 2-d array,
    this is the usual matrix transpose. If axes are given, they describe how
    the axes are permuted.
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
        """Split a multi-dimensional series into individual columns."""
        if self.ndim == 1:
            return [self]
        else:
            n = self.shape[1]
            arr = hsplit(self, n)[0]
            return [self.__class__(numpy.squeeze(a),
                                   self._dates,
                                   **_attrib_dict(self)) for a in arr]

    def filled(self, fill_value=None):
        """Returns an array of the same class as `_data`,  with masked values
filled with `fill_value`. Subclassing is preserved.

*Parameters*:
    fill_value : {None, singleton of type self.dtype}, optional
        The value to fill in masked values with. If `fill_value` is None, uses
        self.fill_value.
"""
        result = self._series.filled(fill_value=fill_value).view(type(self))
        result._dates = self._dates
        return result
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
    result = {'fill_value':series.fill_value}
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
        method = getattr(numpy.asarray(caller), self._methodname)
        try:
            return method(*args, **params)
        except SystemError:
            return getattr(numpy,self._methodname).__call__(caller, *args, **params)
#............................
weekday = _frommethod('weekday')
# deprecated alias for weekday 
day_of_week = weekday
day_of_year = _frommethod('day_of_year')
week = _frommethod('week')
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

##### -------------------------------------------------------------------------
#---- --- TimeSeries creator ---
##### -------------------------------------------------------------------------
def time_series(data, dates=None, start_date=None, freq=None, mask=nomask,
                dtype=None, copy=False, fill_value=None, keep_mask=True,
                hard_mask=False):
    """Creates a TimeSeries object

*Parameters*:
    data : {array_like}
        data portion of the array. Any data that is valid for constructing a
        MaskedArray can be used here. May also be a TimeSeries object
    dates : {DateArray}, optional
        Date part.
    freq : {freq_spec}, optional
        a valid frequency specification
    start_date : {Date}, optional
        date corresponding to index 0 in the data

*Other Parameters*:
    All other parameters that are accepted by the *array* function in the
    maskedarray module are also accepted by this function.

*Notes*:
    the date portion of the time series must be specified in one of the
    following ways:
        - specify a TimeSeries object for the *data* parameter
        - pass a DateArray for the *dates* parameter
        - specify a start_date (a continuous DateArray will be automatically
          constructed for the dates portion)
        - specify just a frequency (for TimeSeries of size zero)
"""
    maparms = dict(copy=copy, dtype=dtype, fill_value=fill_value, subok=True,
                   keep_mask=keep_mask, hard_mask=hard_mask,)
    data = masked_array(data, mask=mask, **maparms)

    freq = check_freq(freq)

    if dates is None:
        _dates = getattr(data, '_dates', None)
    elif isinstance(dates, (Date, DateArray)):
        _dates = date_array(dates)
    elif isinstance(dates, (tuple, list, ndarray)):
        _dates = date_array(dlist=dates, freq=freq)
    else:
        _dates = date_array([], freq=freq)

    if _dates is not None:
        # Make sure _dates has the proper freqncy
        if (freq != _c.FR_UND) and (_dates.freq != freq):
            _dates = _dates.asfreq(freq)
    else:
        dshape = data.shape
        if len(dshape) > 0:
            length = dshape[0]
            _dates = date_array(start_date=start_date, freq=freq, length=length)
        else:
            _dates = date_array([], freq=freq)

    if _dates._unsorted is not None:
        idx = _dates._unsorted
        data = data[idx]
        _dates._unsorted = None
    return TimeSeries(data=data, dates=_dates, mask=data._mask,
                      copy=copy, dtype=dtype,
                      fill_value=fill_value, keep_mask=keep_mask,
                      hard_mask=hard_mask,)


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
        # 2D series w/ only one date : return a new series ....
        elif series._dates.size == 1:
            result = series._series.compressed().view(type(series))
            result._dates = series.dates
            return result
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
    newseries._update_from(a)
    return newseries
#.....................................................
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

#.....................................................
def align_with(*series):
    """Aligns several TimeSeries to the first of the list, so that their
    starting and ending dates match.
    Series are resized and filled with masked values accordingly.
    """
    if len(series) < 2:
        return series
    dates = series[0]._dates[[0,-1]]
    if len(series) == 2:
        return adjust_endpoints(series[-1], dates[0], dates[-1])
    return [adjust_endpoints(x, dates[0], dates[-1]) for x in series[1:]]


#....................................................................
def _convert1d(series, freq, func, position, *args, **kwargs):
    "helper function for `convert` function"
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
        tempData = maskedarray.apply_along_axis(func, -1, tempData, *args, **kwargs)

    newseries = tempData.view(type(series))
    newseries._dates = date_array(start_date=start_date, length=len(newseries),
                                  freq=toFreq)
    newseries._update_from(series)
    return newseries

def convert(series, freq, func=None, position='END', *args, **kwargs):
    """Converts a series to a frequency. Private function called by convert

*Parameters*:
    series : {TimeSeries}
        the series to convert. Skip this parameter if you are calling this as
        a method of the TimeSeries object instead of the module function.
    freq : {freq_spec}
        Frequency to convert the TimeSeries to. Accepts any valid frequency
        specification (string or integer)
    func : {function} (optional)
        When converting to a lower frequency, func is a function that acts on
        one date's worth of data. func should handle masked values appropriately.
        If func is None, then each data point in the resulting series will a
        group of data points that fall into the date at the lower frequency.

        For example, if converting from monthly to daily and you wanted each
        data point in the resulting series to be the average value for each
        month, you could specify maskedarray.average for the 'func' parameter.
    position : {'END', 'START'} (optional)
        When converting to a higher frequency, position is 'START' or 'END'
        and determines where the data point is in each period. For example, if
        going from monthly to daily, and position is 'END', then each data
        point is placed at the end of the month.
    *args : {extra arguments for func parameter} (optional)
        if a func is specified that requires additional parameters, specify
        them here.
    **kwargs : {extra keyword arguments for func parameter} (optional)
        if a func is specified that requires additional keyword parameters,
        specify them here.
"""
    if series.ndim == 1:
        obj = _convert1d(series, freq, func, position, *args, **kwargs)
    elif series.ndim == 2:
        base = _convert1d(series[:,0], freq, func, position, *args, **kwargs)
        obj = maskedarray.column_stack([_convert1d(m,freq,func,position,
                                          *args, **kwargs)._series
                               for m in series.split()]).view(type(series))
        obj._dates = base._dates
        if func is None:
            shp = obj.shape
            ncols = base.shape[-1]
            obj.shape = (shp[0], shp[-1]//ncols, ncols)
            obj = numpy.swapaxes(obj,1,2)
    else:
        raise ValueError(
            "only series with ndim == 1 or ndim == 2 may be converted")

    return obj

TimeSeries.convert = convert

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

TimeSeries.group_byperiod = group_byperiod

#...............................................................................
def tshift(series, nper, copy=True):
    """Returns a series of the same size as `series`, with the same
start_date and end_date, but values shifted by `nper`.

*Parameters*:
    series : {TimeSeries}
        TimeSeries object to shift
    nper : {int}
        number of periods to shift. Negative numbers shift values to the
        right, positive to the left
    copy : {True, False} (optional)
        copies the data if True, returns a view if False.

*Example*:
>>> series = time_series([0,1,2,3], start_date=Date(freq='A', year=2005))
>>> series
timeseries(data  = [0 1 2 3],
           dates = [2005 ... 2008],
           freq  = A-DEC)
>>> tshift(series, -1)
timeseries(data  = [-- 0 1 2],
           dates = [2005 ... 2008],
           freq  = A-DEC)
>>> pct_change = 100 * (series/tshift(series, -1, copy=False) - 1)
"""
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
    newseries._update_from(series)
    return newseries
TimeSeries.tshift = tshift
#...............................................................................
def pct(series, nper=1):
    """Returns the rolling percentage change of the series.

*Parameters*:
    series : {TimeSeries}
        TimeSeries object to to calculate percentage chage for
    nper : {int}
        number of periods for percentage change

*Example*:
>>> series = time_series([2.,1.,2.,3.], start_date=Date(freq='A', year=2005))
>>> pct(series)
timeseries(data  = [-- -50.0 100.0 50.0],
           dates = [2005 ... 2008],
           freq  = A-DEC)
>>> pct(series, 2)
timeseries(data  = [-- -- 0.0 200.0],
           dates = [2005 ... 2008],
           freq  = A-DEC)
"""
    newdata = masked_array(numeric.empty(series.shape, dtype=series.dtype),
                           mask=True)
    if nper < newdata.size:
        newdata[nper:] = 100*(series._series[nper:]/series._series[:-nper] - 1)
    newseries = newdata.view(type(series))
    newseries._dates = series._dates
    newseries._update_from(series)
    return newseries
TimeSeries.pct = pct
#...............................................................................
def fill_missing_dates(data, dates=None, freq=None, fill_value=None):
    """Finds and fills the missing dates in a time series. The data
corresponding to the initially missing dates are masked, or filled to
`fill_value`.

*Parameters*:
    data : {TimeSeries, ndarray}
        Initial array of data.
    dates : {DateArray} (optional)
        Initial array of dates. Specify this if you are passing a plain ndarray
        for the data instead of a TimeSeries.
    freq : {freq_spec} (optional)
        Frequency of result. If not specified, the initial frequency is used.
    fill_value : {scalar of type data.dtype} (optional)
        Default value for missing data. If Not specified, the data are just
        masked.
"""
    # Check the frequency ........
    orig_freq = freq
    freq = check_freq(freq)
    if orig_freq is not None and freq == _c.FR_UND:
        freqstr = check_freq_str(freq)
        raise ValueError,\
              "Unable to define a proper date resolution (found %s)." % freqstr
    # Check the dates .............
    if dates is None:
        if not isTimeSeries(data):
            raise InsufficientDateError
        dates = data._dates
    else:
        if not isinstance(dates, DateArray):
            dates = DateArray(dates, freq)
    dflat = dates.asfreq(freq).ravel()
    if not dflat.has_missing_dates():
        if isinstance(data, TimeSeries):
            return data
        data = data.view(TimeSeries)
        data._dates = dflat
        return data
    # Check the data ..............
    if isinstance(data, MaskedArray):
        datad = data._data
        datam = data._mask
        if isinstance(data, TimeSeries):
            datat = type(data)
        else:
            datat = TimeSeries
    else:
        datad = numpy.asarray(data)
        datam = nomask
        datat = TimeSeries
    # Check whether we need to flatten the data
    if dates.ndim > 1 and dates.ndim == datad.ndim:
        datad.shape = -1
    # ...and now, fill it ! ......
    (tstart, tend) = dflat[[0,-1]]
    newdates = date_array(start_date=tstart, end_date=tend)
    (osize, nsize) = (dflat.size, newdates.size)
    #.............................
    # Get the steps between consecutive data.
    delta = dflat.get_steps()-1
    gap = delta.nonzero()
    slcid = numpy.r_[[0,], numpy.arange(1,osize)[gap], [osize,]]
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
    newshape = list(datad.shape)
    newshape[0] = nsize
    newdatad = numeric.empty(newshape, data.dtype)
    newdatam = numeric.ones(newshape, bool_)
    #....
    if datam is nomask:
        for (new,old) in zip(newslc,oldslc):
            newdatad[new] = datad[old]
            newdatam[new] = False
    else:
        for (new,old) in zip(newslc,oldslc):
            newdatad[new] = datad[old]
            newdatam[new] = datam[old]
    if fill_value is None:
        fill_value = getattr(data, 'fill_value', None)
    newdata = maskedarray.masked_array(newdatad, mask=newdatam, fill_value=fill_value)
    _data = newdata.view(datat)
    _data._dates = newdates
    return _data
#..............................................................................
def stack(*series):
    """Performs a column_stack on the data from each series, and the
resulting series has the same dates as each individual series. All series
must be date compatible.

*Parameters*:
    series : the series to be stacked
"""
    _timeseriescompat_multiple(*series)
    return time_series(maskedarray.column_stack(series), series[0]._dates,
                       **_attrib_dict(series[0]))
#...............................................................................
def concatenate_series(*series, **kwargs):
    msg = """The use of this function is deprecated.
Please use concatenate instead.
Note: Please pay attention to the order of the series!"""
    raise NameError(msg)


def concatenate(series, axis=0, remove_duplicates=True, fill_missing=False):
    """Joins series together.

The series are joined in chronological order. Duplicated dates are handled with
the `remove_duplicates` parameter. If remove_duplicate=False, duplicated dates are
saved. Otherwise, only the first occurence of the date is conserved.

Example
>>> a = time_series([1,2,3], start_date=now('D'))
>>> b = time_series([10,20,30], start_date=now('D')+1)
>>> c = concatenate((a,b))
>>> c._series
masked_array(data = [ 1  2  3 30],
      mask = False,
      fill_value=999999)


*Parameters*:
    series : {sequence}
        Sequence of time series to join
    axis : {integer}
        Axis along which to join
    remove_duplicates : boolean
        Whether to remove duplicated dates.
    fill_missing : {boolean}
        Whether to fill the missing dates with missing values.
    """
    # Get the common frequency, raise an error if incompatibility
    common_f = _compare_frequencies(*series)
    # Concatenate the order of series
    sidx = numpy.concatenate([numpy.repeat(i,len(s))
                              for (i,s) in enumerate(series)], axis=axis)
    # Concatenate the dates and data
    ndates = numpy.concatenate([s._dates for s in series], axis=axis)
    ndata = maskedarray.concatenate([s._series for s in series], axis=axis)
    # Resort the data chronologically
    norder = ndates.argsort(kind='mergesort')
    ndates = ndates[norder]
    ndata = ndata[norder]
    sidx = sidx[norder]
    #
    if not remove_duplicates:
        ndates = date_array_fromlist(ndates, freq=common_f)
        result = time_series(ndata, dates=ndates)
    else:
        # Find the original dates
        orig = numpy.concatenate([[True],(numpy.diff(ndates) != 0)])
        result = time_series(ndata.compress(orig),
                             dates=ndates.compress(orig),freq=common_f)
    if fill_missing:
        result = fill_missing_dates(result)
    return result


#...............................................................................
def empty_like(series):
    """Returns an empty series with the same dtype, mask and dates as series."""
    result = numpy.empty_like(series).view(type(series))
    result._dates = series._dates
    result._mask = series._mask
    return result

################################################################################
if __name__ == '__main__':
    from maskedarray.testutils import assert_equal, assert_array_equal
    if 0:
        dlist = ['2007-01-%02i' % i for i in range(1,16)]
        dates = date_array_fromlist(dlist)
        data = masked_array(numeric.arange(15), mask=[1,0,0,0,0]*3)
        series = time_series(data, dlist)
        #
        aseries = time_series(data, dates+10)
        bseries = time_series(data, dates-10)
        (a, b) = align_with(series, aseries, bseries)
        assert_equal(a._dates, series._dates)
        assert_equal(b._dates, series._dates)
        assert_equal(a[-5:], series[:5])
        assert_equal(b[:5], series[-5:])
    #
    if 0:
        data = numpy.arange(5*24).reshape(5,24)
        datelist = ['2007-07-01','2007-07-02','2007-07-03','2007-07-05','2007-07-06']
        dates = date_array_fromlist(datelist, 'D')
        dseries = time_series(data, dates)
        ndates = date_array_fromrange(start_date=dates[0],end_date=dates[-2])
        #
        fseries = fill_missing_dates(dseries)
        assert_equal(fseries.shape, (6,24))
        assert_equal(fseries._mask[:,0], [0,0,0,1,0,0])
        #
        fseries = fill_missing_dates(dseries[:,0])
        assert_equal(fseries.shape, (6,))
        assert_equal(fseries._mask, [0,0,0,1,0,0])
        #
        series = time_series(data.ravel()[:4].reshape(2,2),dates=dates[:-1])
        fseries = fill_missing_dates(series)
        assert_equal(fseries.shape, (5,))
        assert_equal(fseries._mask, [0,0,0,1,0,])
        #
        fseries = fill_missing_dates(data, date_array_fromlist(datelist,'D'))
    #
    if 0:
        "Make sure we're not losing the fill_value"
        dlist = ['2007-01-%02i' % i for i in range(1,16)]
        dates = date_array_fromlist(dlist)
        series = time_series(maskedarray.zeros(dates.shape), dates=dates, fill_value=-9999)
        assert_equal(series.fill_value, -9999)
    if 0:
        "Check time_series w/ an existing time series"
        dlist = ['2007-01-%02i' % i for i in range(1,16)]
        dates = date_array_fromlist(dlist)
        series = time_series(maskedarray.zeros(dates.shape), dates=dates, fill_value=-9999)
        newseries = time_series(series, fill_value=+9999)
        assert_equal(newseries._data, series._data)
        assert_equal(newseries._mask, series._mask)
        assert_equal(newseries.fill_value, +9999)
    if 1:
        "Check that the fill_value is kept"
        data = [0,1,2,3,4,]
        datelist = ['2007-07-01','2007-07-02','2007-07-03','2007-07-05','2007-07-06']
        dates = date_array_fromlist(datelist, 'D')
        dseries = time_series(data, dates, fill_value=-999)
        ndates = date_array_fromrange(start_date=dates[0],end_date=dates[-2])
        fseries = fill_missing_dates(dseries)
        assert_equal(dseries.fill_value, fseries.fill_value)

    #
    if 0:
        dlist = ['2007-01-%02i' % i for i in (3,2,1)]
        data = [10,20,30]
#        series = time_series(data, dlist, mask=[1,0,0])
#        data = masked_array([10,20,30],mask=[1,0,0])
#        series = time_series(data, dlist)
        series = time_series(data, dlist, mask=[1,0,0])
        assert_equal(series._mask,[0,0,1])
    if 0:
        dlist = ['2007-01-%02i' % i for i in range(1,16)]
        dates = date_array_fromlist(dlist)
        data = masked_array(numeric.arange(15), mask=[1,0,0,0,0]*3)
        series = time_series(data, dlist)

        empty_series = time_series([], freq='d')
        a, b = align_series(series, empty_series)

    if 1:
        "Check concatenate..."
        import dates
        tt = time_series([.2,.2,.3],start_date=dates.Date('T',string='2007-10-10 01:10'))
        tt._dates += [0, 9, 18]
