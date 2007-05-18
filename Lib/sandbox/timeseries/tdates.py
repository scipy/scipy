"""
Classes definition for the support of individual dates and array of dates.

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknox_ca_at_hotmail_dot_com
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import datetime as dt

import itertools
import warnings
import types


import numpy
from numpy import bool_, float_, int_, object_
from numpy import ndarray
import numpy.core.numeric as numeric
import numpy.core.fromnumeric as fromnumeric
import numpy.core.numerictypes as ntypes
from numpy.core.numerictypes import generic

import maskedarray as MA

from parser import DateFromString, DateTimeFromString

import tcore as corelib
import const as _c
import cseries

cseries.set_callback_DateFromString(DateFromString)
cseries.set_callback_DateTimeFromString(DateTimeFromString)

from cseries import Date, thisday, check_freq, check_freq_str, get_freq_group,\
                    DateCalc_Error, DateCalc_RangeError
today = thisday

__all__ = [
'Date', 'DateArray','isDate','isDateArray',
'DateError', 'ArithmeticDateError', 'FrequencyDateError','InsufficientDateError',
'datearray','date_array', 'date_array_fromlist', 'date_array_fromrange',
'day_of_week','day_of_year','day','month','quarter','year','hour','minute',
'second','thisday','today','prevbusday','period_break', 'check_freq',
'check_freq_str','get_freq_group', 'DateCalc_Error', 'DateCalc_RangeError'
           ]


#####---------------------------------------------------------------------------
#---- --- Date Exceptions ---
#####---------------------------------------------------------------------------
class DateError(Exception):
    "Defines a generic DateArrayError."
    def __init__ (self, value=None):
        "Creates an exception."
        self.value = value
    def __str__(self):
        "Calculates the string representation."
        return str(self.value)
    __repr__ = __str__

class InsufficientDateError(DateError):
    """Defines the exception raised when there is not enough information
    to create a Date object."""
    def __init__(self, msg=None):
        if msg is None:
            msg = "Insufficient parameters given to create a date at the given frequency"
        DateError.__init__(self, msg)

class FrequencyDateError(DateError):
    """Defines the exception raised when the frequencies are incompatible."""
    def __init__(self, msg, freql=None, freqr=None):
        msg += " : Incompatible frequencies!"
        if not (freql is None or freqr is None):
            msg += " (%s<>%s)" % (freql, freqr)
        DateError.__init__(self, msg)

class ArithmeticDateError(DateError):
    """Defines the exception raised when dates are used in arithmetic expressions."""
    def __init__(self, msg=''):
        msg += " Cannot use dates for arithmetics!"
        DateError.__init__(self, msg)


#####---------------------------------------------------------------------------
#---- --- Functions ---
#####---------------------------------------------------------------------------

def prevbusday(day_end_hour=18, day_end_min=0):
    """Returns the previous business day (Monday-Friday) at business frequency.

:Parameters:
    - day_end_hour : (int, *[18]* )
    - day_end_min : (int, *[0]*)

:Return values:
    If it is currently Saturday or Sunday, then the preceding Friday will be
    returned. If it is later than the specified day_end_hour and day_end_min,
    thisday('b') will be returned. Otherwise, thisday('b')-1 will be returned.
"""
    tempDate = dt.datetime.now()
    dateNum = tempDate.hour + float(tempDate.minute)/60
    checkNum = day_end_hour + float(day_end_min)/60
    if dateNum < checkNum:
        return thisday(_c.FR_BUS) - 1
    else:
        return thisday(_c.FR_BUS)


def isDate(data):
    "Returns whether `data` is an instance of Date."
    return isinstance(data, Date) or \
           (hasattr(data,'freq') and hasattr(data,'value'))


#####---------------------------------------------------------------------------
#---- --- DateArray ---
#####---------------------------------------------------------------------------
ufunc_dateOK = ['add','subtract',
                'equal','not_equal','less','less_equal', 'greater','greater_equal',
                'isnan']

class _datearithmetics(object):
    """Defines a wrapper for arithmetic methods.
Instead of directly calling a ufunc, the corresponding method of  the `array._data`
object is called instead.
If `asdates` is True, a DateArray object is returned , else a regular ndarray
is returned.
    """
    def __init__ (self, methodname, asdates=True):
        """
:Parameters:
    - `methodname` (String) : Method name.
        """
        self.methodname = methodname
        self._asdates = asdates
        self.__doc__ = getattr(methodname, '__doc__')
        self.obj = None
    #
    def __get__(self, obj, objtype=None):
        self.obj = obj
        return self
    #
    def __call__ (self, other, *args, **kwargs):
        "Execute the call behavior."
        instance = self.obj
        freq = instance.freq
        if 'context' not in kwargs:
            kwargs['context'] = 'DateOK'
        method = getattr(super(DateArray,instance), self.methodname)
        if isinstance(other, DateArray):
            if other.freq != freq:
                raise FrequencyDateError("Cannot operate on dates", \
                                         freq, other.freq)
        elif isinstance(other, Date):
            if other.freq != freq:
                raise FrequencyDateError("Cannot operate on dates", \
                                         freq, other.freq)
            other = other.value
        elif isinstance(other, ndarray):
            if other.dtype.kind not in ['i','f']:
                raise ArithmeticDateError
        if self._asdates:
            return instance.__class__(method(other, *args),
                                      freq=freq)
        else:
            return method(other, *args)

class DateArray(ndarray):
    """Defines a ndarray of dates, as ordinals.

When viewed globally (array-wise), DateArray is an array of integers.
When viewed element-wise, DateArray is a sequence of dates.
For example, a test such as :
>>> DateArray(...) = value
will be valid only if value is an integer, not a Date
However, a loop such as :
>>> for d in DateArray(...):
accesses the array element by element. Therefore, `d` is a Date object.
    """
    def __new__(cls, dates=None, freq=None, copy=False):
        # Get the frequency ......
        if freq is None:
            _freq = getattr(dates, 'freq', _c.FR_UND)
        else:
            _freq = check_freq(freq)
        # Get the dates ..........
        _dates = numeric.array(dates, copy=copy, dtype=int_, subok=1)
        if _dates.ndim == 0:
            _dates.shape = (1,)
        _dates = _dates.view(cls)
        _dates.freq = _freq
        _dates._unsorted = None
        return _dates

    def __array_wrap__(self, obj, context=None):
        if context is None:
            return self
        elif context[0].__name__ not in ufunc_dateOK:
            raise ArithmeticDateError, "(function %s)" % context[0].__name__

    def __array_finalize__(self, obj):
        self.freq = getattr(obj, 'freq', _c.FR_UND)
        self._unsorted = getattr(obj,'_unsorted',None)
        self._cachedinfo = dict(toobj=None, tostr=None, toord=None,
                                steps=None, full=None, hasdups=None)
        if hasattr(obj,'_cachedinfo'):
            self._cachedinfo.update(obj._cachedinfo)
        return

    def __getitem__(self, indx):
        reset_full = True
        if isinstance(indx, Date):
            indx = self.find_dates(indx)
            reset_full = False
        elif numeric.asarray(indx).dtype.kind == 'O':
            try:
                indx = self.find_dates(indx)
            except AttributeError:
                pass
        r = ndarray.__getitem__(self, indx)
        if isinstance(r, (generic, int)):
            return Date(self.freq, value=r)
        elif hasattr(r, 'size') and r.size == 1:
            # need to check if it has a size attribute for situations
            # like when the datearray is the data for a maskedarray
            # or some other subclass of ndarray with wierd getitem
            # behaviour
            return Date(self.freq, value=r.item())
        else:
            if hasattr(r, '_cachedinfo'):
                _cache = r._cachedinfo
                _cache.update(dict([(k,_cache[k][indx])
                                    for k in ('toobj', 'tostr', 'toord')
                                    if _cache[k] is not None]))
                _cache['steps'] = None
            if reset_full:
                _cache['full'] = None
                _cache['hasdups'] = None
                
            return r

    def __getslice__(self, i, j):
        r = ndarray.__getslice__(self, i, j)
        if hasattr(r, '_cachedinfo'):
            _cache = r._cachedinfo
            _cache.update(dict([(k,_cache[k][i:j])
                                for k in ('toobj', 'tostr', 'toord')
                                if _cache[k] is not None]))
            _cache['steps'] = None
        return r

    def __repr__(self):
        return ndarray.__repr__(self)[:-1] + \
               ",\n          freq='%s')" % self.freqstr
    #......................................................
    __add__ = _datearithmetics('__add__', asdates=True)
    __radd__ = _datearithmetics('__add__', asdates=True)
    __sub__ = _datearithmetics('__sub__', asdates=True)
    __rsub__ = _datearithmetics('__rsub__', asdates=True)
    __le__ = _datearithmetics('__le__', asdates=False)
    __lt__ = _datearithmetics('__lt__', asdates=False)
    __ge__ = _datearithmetics('__ge__', asdates=False)
    __gt__ = _datearithmetics('__gt__', asdates=False)
    __eq__ = _datearithmetics('__eq__', asdates=False)
    __ne__ = _datearithmetics('__ne__', asdates=False)
    #......................................................
    @property
    def freqstr(self):
        "Returns the frequency string code."
        return check_freq_str(self.freq)
    @property
    def day(self):
        "Returns the day of month."
        return self.__getdateinfo__('D')
    @property
    def day_of_week(self):
        "Returns the day of week."
        return self.__getdateinfo__('W')
    @property
    def day_of_year(self):
        "Returns the day of year."
        return self.__getdateinfo__('R')
    @property
    def month(self):
        "Returns the month."
        return self.__getdateinfo__('M')
    @property
    def quarter(self):
        "Returns the quarter."
        return self.__getdateinfo__('Q')
    @property
    def year(self):
        "Returns the year."
        return self.__getdateinfo__('Y')
    @property
    def qyear(self):
        """For quarterly frequency dates, returns the year corresponding to the
year end (start) month. When using QTR or QTR-E based quarterly
frequencies, this is the fiscal year in a financial context.

For non-quarterly dates, this simply returns the year of the date."""

        return self.__getdateinfo__('F')
    @property
    def second(self):
        "Returns the seconds."
        return self.__getdateinfo__('S')
    @property
    def minute(self):
        "Returns the minutes."
        return self.__getdateinfo__('T')
    @property
    def hour(self):
        "Returns the hour."
        return self.__getdateinfo__('H')
    @property
    def week(self):
        "Returns the week."
        return self.__getdateinfo__('I')

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

    def __getdateinfo__(self, info):
        return numeric.asarray(cseries.DA_getDateInfo(numeric.asarray(self),
                                                      self.freq, info),
                               dtype=int_)
    __getDateInfo = __getdateinfo__
    #.... Conversion methods ....................
    #
    def tovalue(self):
        "Converts the dates to integer values."
        return numeric.asarray(self)
    #
    def toordinal(self):
        "Converts the dates from values to ordinals."
        # Note: we better try to cache the result
        if self._cachedinfo['toord'] is None:
#            diter = (Date(self.freq, value=d).toordinal() for d in self)
            if self.freq == _c.FR_UND:
                diter = (d.value for d in self)
            else:
                diter = (d.toordinal() for d in self)
            toord = numeric.fromiter(diter, dtype=float_)
            self._cachedinfo['toord'] = toord
        return self._cachedinfo['toord']
    #
    def tostring(self):
        "Converts the dates to strings."
        # Note: we better cache the result
        if self._cachedinfo['tostr'] is None:
            firststr = str(self[0])
            if self.size > 0:
                ncharsize = len(firststr)
                tostr = numpy.fromiter((str(d) for d in self),
                                        dtype='|S%i' % ncharsize)
            else:
                tostr = firststr
            self._cachedinfo['tostr'] = tostr
        return self._cachedinfo['tostr']
    #
    def asfreq(self, freq=None, relation="AFTER"):
        "Converts the dates to another frequency."
        # Note: As we define a new object, we don't need caching
        if freq is None or freq == _c.FR_UND:
            return self
        tofreq = check_freq(freq)
        if tofreq == self.freq:
            return self
        _rel = relation.upper()[0]
        fromfreq = self.freq
        if fromfreq == _c.FR_UND:
             fromfreq = _c.FR_DAY
        new = cseries.DA_asfreq(numeric.asarray(self), fromfreq, tofreq, _rel)
        return DateArray(new, freq=freq)

    #......................................................
    def find_dates(self, *dates):
        "Returns the indices corresponding to given dates, as an array."
        ifreq = self.freq
        c = numpy.zeros(self.shape, bool_)
        for d in corelib.flatargs(*dates):
            if d.freq != ifreq:
                d = d.asfreq(ifreq)
            c += (self == d.value)
        c = c.nonzero()
        if fromnumeric.size(c) == 0:
            raise IndexError, "Date out of bounds!"
        return c

    def date_to_index(self, date):
        "Returns the index corresponding to one given date, as an integer."
        if self.isvalid():
            index = date.value - self[0].value
            if index < 0 or index > self.size:
                raise IndexError, "Date out of bounds!"
            return index
        else:
            index_asarray = (self == date.value).nonzero()
            if fromnumeric.size(index_asarray) == 0:
                raise IndexError, "Date out of bounds!"
            return index_asarray[0][0]
    #......................................................
    def get_steps(self):
        """Returns the time steps between consecutive dates.
    The timesteps have the same unit as the frequency of the series."""
        if self.freq == _c.FR_UND:
            warnings.warn("Undefined frequency: assuming integers!")
        if self._cachedinfo['steps'] is None:
            _cached = self._cachedinfo
            val = numeric.asarray(self).ravel()
            if val.size > 1:
                steps = val[1:] - val[:-1]
                if _cached['full'] is None:
                    _cached['full'] = (steps.max() == 1)
                if _cached['hasdups'] is None:
                    _cached['hasdups'] = (steps.min() == 0)
            else:
                _cached['full'] = True
                _cached['hasdups'] = False
                steps = numeric.array([], dtype=int_)
            self._cachedinfo['steps'] = steps
        return self._cachedinfo['steps']

    def has_missing_dates(self):
        "Returns whether the DateArray have missing dates."
        if self._cachedinfo['full'] is None:
            steps = self.get_steps()
        return not(self._cachedinfo['full'])

    def isfull(self):
        "Returns whether the DateArray has no missing dates."
        if self._cachedinfo['full'] is None:
            steps = self.get_steps()
        return self._cachedinfo['full']

    def has_duplicated_dates(self):
        "Returns whether the DateArray has duplicated dates."
        if self._cachedinfo['hasdups'] is None:
            steps = self.get_steps()
        return self._cachedinfo['hasdups']

    def isvalid(self):
        "Returns whether the DateArray is valid: no missing/duplicated dates."
        return  (self.isfull() and not self.has_duplicated_dates())
    #......................................................

#............................


#####---------------------------------------------------------------------------
#---- --- DateArray functions ---
#####---------------------------------------------------------------------------
def isDateArray(a):
    "Tests whether an array is a DateArray object."
    return isinstance(a,DateArray)

def guess_freq(dates):
    """Tries to estimate the frequency of a list of dates, by checking the steps
    between consecutive dates The steps should be in days.
    Returns a frequency code (alpha character)."""
    ddif = numeric.asarray(numpy.diff(dates))
    ddif.sort()
    if ddif.size == 0:
        fcode = _c.FR_UND
    elif ddif[0] == ddif[-1] == 1.:
        fcode = _c.FR_DAY
    elif (ddif[0] == 1.) and (ddif[-1] == 3.):
        fcode = _c.FR_BUS
    elif (ddif[0] > 3.) and  (ddif[-1] == 7.):
        fcode = _c.FR_WK
    elif (ddif[0] >= 28.) and (ddif[-1] <= 31.):
        fcode = _c.FR_MTH
    elif (ddif[0] >= 90.) and (ddif[-1] <= 92.):
        fcode = _c.FR_QTR
    elif (ddif[0] >= 365.) and (ddif[-1] <= 366.):
        fcode = _c.FR_ANN
    elif numpy.abs(24.*ddif[0] - 1) <= 1e-5 and \
         numpy.abs(24.*ddif[-1] - 1) <= 1e-5:
        fcode = _c.FR_HR
    elif numpy.abs(1440.*ddif[0] - 1) <= 1e-5 and \
         numpy.abs(1440.*ddif[-1] - 1) <= 1e-5:
        fcode = _c.FR_MIN
    elif numpy.abs(86400.*ddif[0] - 1) <= 1e-5 and \
         numpy.abs(86400.*ddif[-1] - 1) <= 1e-5:
        fcode = _c.FR_SEC
    else:
        warnings.warn("Unable to estimate the frequency! %.3f<>%.3f" %\
                      (ddif[0], ddif[-1]))
        fcode = _c.FR_UND
    return fcode


def _listparser(dlist, freq=None):
    "Constructs a DateArray from a list."
    dlist = numeric.asarray(dlist)
    idx = dlist.argsort()
    dlist = dlist[idx]
    if dlist.ndim == 0:
        dlist.shape = (1,)
    # Case #1: dates as strings .................
    if dlist.dtype.kind in 'SU':
        #...construct a list of ordinals
        ords = numpy.fromiter((DateTimeFromString(s).toordinal() for s in dlist),
                               float_)
        ords += 1
        #...try to guess the frequency
        if freq is None or freq == _c.FR_UND:
            freq = guess_freq(ords)
        #...construct a list of dates
        for s in dlist:
            x = Date(freq, string=s)
        dates = [Date(freq, string=s) for s in dlist]
    # Case #2: dates as numbers .................
    elif dlist.dtype.kind in 'if':
        #...hopefully, they are values
        if freq is None or freq == _c.FR_UND:
            freq = guess_freq(dlist)
        dates = dlist
    # Case #3: dates as objects .................
    elif dlist.dtype.kind == 'O':
        template = dlist[0]
        #...as Date objects
        if isinstance(template, Date):
            dates = numpy.fromiter((d.value for d in dlist), int_)
        #...as mx.DateTime objects
        elif hasattr(template,'absdays'):
            # no freq given: try to guess it from absdays
            if freq == _c.FR_UND:
                ords = numpy.fromiter((s.absdays for s in dlist), float_)
                ords += 1
                freq = guess_freq(ords)
            dates = [Date(freq, datetime=m) for m in dlist]
        #...as datetime objects
        elif hasattr(template, 'toordinal'):
            ords = numpy.fromiter((d.toordinal() for d in dlist), float_)
            if freq == _c.FR_UND:
                freq = guess_freq(ords)
            dates = [Date(freq, datetime=dt.datetime.fromordinal(a)) for a in ords]
    #
    result = DateArray(dates, freq)
    result._unsorted = idx
    return result


def date_array(dlist=None, start_date=None, end_date=None, length=None,
               freq=None):
    """Constructs a DateArray from:
    - a starting date and either an ending date or a given length.
    - a list of dates.
    """
    freq = check_freq(freq)
    # Case #1: we have a list ...................
    if dlist is not None:
        # Already a DateArray....................
        if isinstance(dlist, DateArray):
            if (freq != _c.FR_UND) and (dlist.freq != check_freq(freq)):
                return dlist.asfreq(freq)
            else:
                return dlist
        # Make sure it's a sequence, else that's a start_date
        if hasattr(dlist,'__len__'):
            return _listparser(dlist, freq)
        elif start_date is not None:
            if end_date is not None:
                dmsg = "What starting date should be used ? '%s' or '%s' ?"
                raise DateError, dmsg % (dlist, start_date)
            else:
                (start_date, end_date) = (dlist, start_date)
        else:
            start_date = dlist
    # Case #2: we have a starting date ..........
    if start_date is None:
        if length == 0:
            return DateArray([], freq=freq)
        raise InsufficientDateError
    if not isDate(start_date):
        dmsg = "Starting date should be a valid Date instance! "
        dmsg += "(got '%s' instead)" % type(start_date)
        raise DateError, dmsg
    # Check if we have an end_date
    if end_date is None:
        if length is None:
#            raise ValueError,"No length precised!"
            length = 1
    else:
        if not isDate(end_date):
            raise DateError, "Ending date should be a valid Date instance!"
        length = int(end_date - start_date) + 1
#    dlist = [(start_date+i).value for i in range(length)]
    dlist = numeric.arange(length, dtype=int_)
    dlist += start_date.value
    if freq == _c.FR_UND:
        freq = start_date.freq
    return DateArray(dlist, freq=freq)
datearray = date_array

def date_array_fromlist(dlist, freq=None):
    "Constructs a DateArray from a list of dates."
    return date_array(dlist=dlist, freq=freq)

def date_array_fromrange(start_date, end_date=None, length=None,
                         freq=None):
    """Constructs a DateArray from a starting date and either an ending date or
    a length."""
    return date_array(start_date=start_date, end_date=end_date,
                      length=length, freq=freq)

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
            return getattr(DateArray, self._methodname).__doc__
        except AttributeError:
            return "???"
    #
    def __call__(self, caller, *args, **params):
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


def period_break(dates, period):
    """Returns the indices where the given period changes.

:Parameters:
    dates : DateArray
        Array of dates to monitor.
    period : string
        Name of the period to monitor.
    """
    current = getattr(dates, period)
    previous = getattr(dates-1, period)
    return (current - previous).nonzero()[0]


################################################################################

if __name__ == '__main__':
    import maskedarray.testutils
    from maskedarray.testutils import assert_equal
    if 0:
        dlist = ['2007-%02i' % i for i in range(1,5)+range(7,13)]
        mdates = date_array_fromlist(dlist, 'M')
        # Using an integer
        assert_equal(mdates[0].value, 24073)
        assert_equal(mdates[-1].value, 24084)
        # Using a date
        lag = mdates.find_dates(mdates[0])
        print mdates[lag]
        assert_equal(mdates[lag], mdates[0])
    if 0:
        hodie = today('D')
        D = DateArray(today('D'))
        assert_equal(D.freq, 6000)
    if 0:
        freqs = [x[0] for x in corelib.freq_dict.values() if x[0] != 'U']
        print freqs
        for f in freqs:
            print f
            today = thisday(f)
            assert(Date(freq=f, value=today.value) == today)
    if 0:
        D = date_array(freq='U', start_date=Date('U',1), length=10)
    if 0:
        dlist = ['2007-01-%02i' % i for i in (1,2,4,5,7,8,10,11,13)]
        ords = numpy.fromiter((DateTimeFromString(s).toordinal() for s in dlist),
                               float_)
    if 0:
        "Tests the automatic sorting of dates."
        D = date_array_fromlist(dlist=['2006-01','2005-01','2004-01'],freq='M')
        assert_equal(D.view(ndarray), [24037, 24049, 24061])

    if 1:
        dlist = ['2007-%02i' % i for i in range(1,5)+range(7,13)]
        mdates = date_array_fromlist(dlist, 'M')
        
        print mdates.tostr()
        