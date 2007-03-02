"""
Classes definition for the support of individual dates and array of dates.

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknox_ca_at_hotmail_dot_com
:version: $Id: tdates.py 2805 2007-03-01 21:40:02Z mattknox_ca $
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author: mattknox_ca $)"
__version__ = '1.0'
__revision__ = "$Revision: 2805 $"
__date__     = '$Date: 2007-03-01 16:40:02 -0500 (Thu, 01 Mar 2007) $'

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

try:
    from mx.DateTime import DateTimeType
except ImportError:
    class DateTimeType: pass

from parser import DateFromString, DateTimeFromString    

import tcore as corelib
import tcore as _c
import cseries



__all__ = [
'Date', 'DateArray','isDate','isDateArray',
'DateError', 'ArithmeticDateError', 'FrequencyDateError','InsufficientDateError',
'datearray','date_array', 'date_array_fromlist', 'date_array_fromrange',
'day_of_week','day_of_year','day','month','quarter','year','hour','minute','second',
'truncateDate','monthToQuarter','thisday','today','prevbusday','asfreq',
'period_break'
           ]

#####---------------------------------------------------------------------------
#---- --- Date Info ---
#####---------------------------------------------------------------------------

OriginDate = dt.datetime(1970, 1, 1)
secondlyOriginDate = OriginDate - dt.timedelta(seconds=1)
minutelyOriginDate = OriginDate - dt.timedelta(minutes=1)
hourlyOriginDate = OriginDate - dt.timedelta(hours=1)


#####---------------------------------------------------------------------------
#---- --- Date Exceptions ---
#####---------------------------------------------------------------------------
class DateError(Exception):
    """Defines a generic DateArrayError."""
    def __init__ (self, args=None):
        "Create an exception"
        Exception.__init__(self)
        self.args = args
    def __str__(self):
        "Calculate the string representation"
        return str(self.args)
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
#---- --- Date Class ---
#####---------------------------------------------------------------------------

class Date:
    """Defines a Date object, as the combination of a date and a frequency.
    Several options are available to construct a Date object explicitly:

    - Give appropriate values to the `year`, `month`, `day`, `quarter`, `hours`, 
      `minutes`, `seconds` arguments.
      
      >>> td.Date(freq='Q',year=2004,quarter=3)
      >>> td.Date(freq='D',year=2001,month=1,day=1)
      
    - Use the `string` keyword. This method calls the `mx.DateTime.Parser`
      submodule, more information is available in its documentation.
      
      >>> ts.Date('D', '2007-01-01')
      
    - Use the `datetime` keyword with an existing datetime.datetime object.

      >>> td.Date('D', datetime=datetime.datetime.now())
      """
    default_fmtstr = {'A': "%Y",
                      'Q': "%YQ%q",
                      'M': "%b-%Y",
                      'W': "%d-%b-%Y",
                      'B': "%d-%b-%Y",
                      'D': "%d-%b-%Y",
                      'U': "%d-%b-%Y",
                      'H': "%d-%b-%Y %H:00",
                      'T': "%d-%b-%Y %H:%M",
                      'S': "%d-%b-%Y %H:%M:%S"
                      }
      
    def __init__(self, freq, value=None, string=None,
                 year=None, month=None, day=None, quarter=None, 
                 hour=None, minute=None, second=None, 
                 datetime=None):
        
        if hasattr(freq, 'freq'):
            self.freq = corelib.check_freq(freq.freq)
        else:
            self.freq = corelib.check_freq(freq)
        self.freqstr = corelib.freq_tostr(self.freq)
        
        
        if value is not None:
            if isinstance(value, ndarray):
                value = int(value)
        
            if isinstance(value, str):
                if self.freq in (_c.FR_HR, _c.FR_MIN, _c.FR_SEC):
                    self.datetime = DateTimeFromString(value)
                else:
                    self.datetime = DateFromString(value)

            elif self.freq == _c.FR_SEC:
                self.datetime = secondlyOriginDate + dt.timedelta(seconds=value)
            elif self.freq == _c.FR_MIN:
                self.datetime = minutelyOriginDate + dt.timedelta(minutes=value)
            elif self.freq == _c.FR_HR:
                self.datetime = hourlyOriginDate + dt.timedelta(hours=value)
            elif self.freq in (_c.FR_DAY, _c.FR_UND):
                self.datetime = dt.datetime.fromordinal(value)
            elif self.freq == _c.FR_BUS:
                valtmp = (value - 1)//5
                self.datetime = dt.datetime.fromordinal(value + valtmp*2)
            elif self.freq == _c.FR_WK:
                self.datetime = dt.datetime(1,1,7) + \
                                dt.timedelta(days=(value-1)*7)
            elif self.freq == _c.FR_MTH:
                year = (value - 1)//12 + 1
                month = value - (year - 1)*12
                self.datetime = dt.datetime(year, month, 1)
            elif self.freq == _c.FR_QTR:
                year = (value - 1)//4 + 1
                month = (value - (year - 1)*4)*3
                self.datetime = dt.datetime(year, month, 1)
            elif self.freq == _c.FR_ANN:
                self.datetime = dt.datetime(value, 1, 1)
            else:
                raise ValueError("unrecognized frequency: "+str(self.freq))
        
        elif string is not None:
            if self.freq in (_c.FR_HR, _c.FR_MIN, _c.FR_SEC):
                self.datetime = DateTimeFromString(string)
            else:
                self.datetime = DateFromString(string)
            
        elif datetime is not None:
            if isinstance(datetime, DateTimeType):
                datetime = mx_to_datetime(datetime)
            self.datetime = truncateDate(self.freq, datetime)
            
        else:
            # First, some basic checks.....
            if year is None:
                raise InsufficientDateError            
            if self.freq in (_c.FR_BUS, _c.FR_DAY, _c.FR_WK, _c.FR_UND):
                if month is None or day is None: 
                    raise InsufficientDateError
            elif self.freq == _c.FR_MTH:
                if month is None: 
                    raise InsufficientDateError
                day = 1
            elif self.freq == _c.FR_QTR:
                if quarter is None: 
                    raise InsufficientDateError
                month = quarter * 3
                day = 1
            elif self.freq == _c.FR_ANN:
                month = 1
                day = 1
            elif self.freq == _c.FR_SEC:
                if month is None or day is None or second is None: 
                    raise InsufficientDateError
                
            if self.freq in (_c.FR_BUS, _c.FR_DAY, _c.FR_WK,
                             _c.FR_MTH, _c.FR_QTR, _c.FR_ANN):
                self.datetime = truncateDate(self.freq, dt.datetime(year, month, day))
                if self.freq == _c.FR_BUS:
                    if self.datetime.isoweekday() in [6,7]:
                        raise ValueError("Weekend passed as business day")
            elif self.freq in (_c.FR_HR, _c.FR_MIN, _c.FR_SEC):
                if hour is None:
                    if minute is None:
                        if second is None:
                            hour = 0
                        else:
                            hour = second//3600
                    else:
                        hour = minute // 60
                if minute is None:
                    if second is None:
                        minute = 0
                    else:
                        minute = (second-hour*3600)//60
                if second is None:
                    second = 0
                else:
                    second = second % 60
                self.datetime = truncateDate(self.freqstr,
                                             dt.datetime(year, month, day, 
                                                         hour, minute, second))
            else:
                raise ValueError("unrecognized frequency: "+str(self.freq))

        self.value = self.__value()

    def __getitem__(self, indx):
        return self

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
        
    def __getdateinfo__(self, info):
        return int(cseries.getDateInfo(numpy.asarray(self.value), 
                                       self.freq, info))
    __getDateInfo = __getdateinfo__
 
    def __add__(self, other):
        if isinstance(other, Date):
            raise FrequencyDateError("Cannot add dates", 
                                     self.freqstr, other.freqstr)
        return Date(freq=self.freq, value=int(self) + other)
    
    def __radd__(self, other): 
        return self+other
    
    def __sub__(self, other):
        if isinstance(other, Date):
            if self.freq != other.freq:
                raise FrequencyDateError("Cannot subtract dates", \
                                         self.freqstr, other.freqstr)
            else:
                return int(self) - int(other) 
        else:
            return self + (-1) * int(other)
    
    def __eq__(self, other):
        if not hasattr(other, 'freq'):
            return False
        elif self.freq != other.freq:
            raise FrequencyDateError("Cannot compare dates", \
                                     self.freqstr, other.freqstr)
        return int(self) == int(other) 
    
    def __cmp__(self, other): 
        if not hasattr(other, 'freq'):
            return False
        elif self.freq != other.freq:
            raise FrequencyDateError("Cannot compare dates", \
                                     self.freqstr, other.freqstr)
        return int(self)-int(other)    
        
    def __hash__(self): 
        return hash(int(self)) ^ hash(self.freq)
    
    def __int__(self):
        return self.value
    
    def __float__(self):
        return float(self.value)
    
    def __value(self):   
        "Converts the date to an integer, depending on the current frequency."
        # Secondly......
        if self.freq == _c.FR_SEC:
            delta = (self.datetime - secondlyOriginDate)
            val = delta.days*86400 + delta.seconds
        # Minutely......
        elif self.freq == _c.FR_MIN:
            delta = (self.datetime - minutelyOriginDate)
            val = delta.days*1440 + delta.seconds/(60)
        # Hourly........
        elif self.freq == _c.FR_HR:
            delta = (self.datetime - hourlyOriginDate)
            val = delta.days*24 + delta.seconds/(3600)
        # Daily/undefined
        elif self.freq in (_c.FR_DAY, _c.FR_UND):
            val = self.datetime.toordinal()
        # Business days.
        elif self.freq == _c.FR_BUS:
            days = self.datetime.toordinal()
            weeks = days // 7
            val = days - weeks*2  
        # Weekly........
        elif self.freq == _c.FR_WK:
            val = self.datetime.toordinal()//7
        # Monthly.......
        elif self.freq == _c.FR_MTH:
            val = (self.datetime.year-1)*12 + self.datetime.month
        # Quarterly.....
        elif self.freq == _c.FR_QTR:
            val = (self.datetime.year-1)*4 + self.datetime.month//3
        # Annual .......
        elif self.freq == _c.FR_ANN:
            val = self.datetime.year

        return int(val)
    #......................................................        
    def strfmt(self, fmt):
        "Formats the date"
        if fmt is None:
            fmt = self.default_fmtstr[self.freqstr]
        return cseries.strfmt(self.datetime, fmt)
            
    def __str__(self):
        return self.strfmt(self.default_fmtstr[self.freqstr])

    def __repr__(self): 
        return "<%s : %s>" % (str(self.freqstr), str(self))
    #......................................................
    def toordinal(self):
        "Returns the date as an ordinal."
        return self.datetime.toordinal()

    def fromordinal(self, ordinal):
        "Returns the date as an ordinal."
        return Date(self.freq, datetime=dt.datetime.fromordinal(ordinal))
    
    def tostring(self):
        "Returns the date as a string."
        return str(self)
    
    def toobject(self):
        "Returns the date as itself."
        return self
    
    def isvalid(self):
        "Returns whether the DateArray is valid: no missing/duplicated dates."
        # A date is always valid by itself, but we need the object to support the function
        # when we're working with singletons.
        return True
    #......................................................
        
    
#####---------------------------------------------------------------------------
#---- --- Functions ---
#####---------------------------------------------------------------------------

def mx_to_datetime(mxDate):
    microsecond = 1000000*(mxDate.second % 1)
    return dt.datetime(mxDate.year, mxDate.month,
                       mxDate.day, mxDate.hour,
                       mxDate.minute,
                       int(mxDate.second), microsecond)


def truncateDate(freq, datetime):
    "Chops off the irrelevant information from the datetime object passed in."
    freq = corelib.check_freq(freq)
    if freq == _c.FR_MIN:
        return dt.datetime(datetime.year, datetime.month, datetime.day, \
                           datetime.hour, datetime.minute)
    elif freq == _c.FR_HR:
        return dt.datetime(datetime.year, datetime.month, datetime.day, \
                           datetime.hour)
    elif freq in (_c.FR_BUS, _c.FR_DAY):
        if freq == _c.FR_BUS and datetime.isoweekday() in (6,7):
            raise ValueError("Weekend passed as business day")
        return dt.datetime(datetime.year, datetime.month, datetime.day)
    elif freq == _c.FR_WK:
        d = datetime.toordinal()
        return dt.datetime.fromordinal(d + (7 - d % 7) % 7)
    elif freq == _c.FR_MTH:
        return dt.datetime(datetime.year, datetime.month, 1)
    elif freq == _c.FR_QTR:
        return dt.datetime(datetime.year, monthToQuarter(datetime.month)*3, 1)
    elif freq == _c.FR_ANN:
        return dt.datetime(datetime.year, 1, 1)
    else:
        return datetime
    
def monthToQuarter(monthNum):
    """Returns the quarter corresponding to the month `monthnum`.
    For example, December is the 4th quarter, Januray the first."""
    return int((monthNum-1)/3)+1

def thisday(freq):
    "Returns today's date, at the given frequency `freq`."
    freq = corelib.check_freq(freq)
    tempDate = dt.datetime.now()
    # if it is Saturday or Sunday currently, freq==B, then we want to use Friday
    if freq == _c.FR_BUS and tempDate.isoweekday() >= 6:
        tempDate = tempDate - dt.timedelta(days=(tempDate.isoweekday() - 5))
    if freq in (_c.FR_BUS,_c.FR_DAY,
                _c.FR_HR,_c.FR_SEC,_c.FR_MIN,
                _c.FR_WK,_c.FR_UND):
        return Date(freq, datetime=tempDate)
    elif freq == _c.FR_MTH:
        return Date(freq, year=tempDate.year, month=tempDate.month)
    elif freq == _c.FR_QTR:
        return Date(freq, year=tempDate.year, quarter=monthToQuarter(tempDate.month))
    elif freq == _c.FR_ANN:
        return Date(freq, year=tempDate.year)
today = thisday

def prevbusday(day_end_hour=18, day_end_min=0):
    "Returns the previous business day."
    tempDate = dt.datetime.now()
    dateNum = tempDate.hour + float(tempDate.minute)/60
    checkNum = day_end_hour + float(day_end_min)/60
    if dateNum < checkNum: 
        return thisday(_c.FR_BUS) - 1
    else: 
        return thisday(_c.FR_BUS)
                
def asfreq(date, toFreq, relation="BEFORE"):
    """Returns a date converted to another frequency `toFreq`, according to the
    relation `relation` ."""
    tofreq = corelib.check_freq(toFreq)
    _rel = relation.upper()[0]
    if _rel not in ['B', 'A']:
        msg = "Invalid relation '%s': Should be in ['before', 'after']"
        raise ValueError, msg % relation

    if not isinstance(date, Date):
        raise DateError, "Date should be a valid Date instance!"

    if date.freq == _c.FR_UND:
        warnings.warn("Undefined frequency: assuming daily!")
        fromfreq = _c.FR_DAY
    else:
        fromfreq = date.freq
    
    if fromfreq == tofreq:
        return date
    else:
        value = cseries.asfreq(numeric.asarray(date.value), fromfreq, tofreq, _rel)
        if value > 0:
            return Date(freq=tofreq, value=value)
        else:
            return None
Date.asfreq = asfreq
            
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
    _defcachedinfo = dict(toobj=None, tostr=None, toord=None, 
                          steps=None, full=None, hasdups=None)
    def __new__(cls, dates=None, freq=None, copy=False):
        # Get the frequency ......
        if freq is None:
            _freq = getattr(dates, 'freq', -9999)
        else:
            _freq = corelib.check_freq(freq)
        cls._defaultfreq = corelib.check_freq(_freq)
        # Get the dates ..........
        _dates = numeric.array(dates, copy=copy, dtype=int_, subok=1)
        if _dates.ndim == 0:
            _dates.shape = (1,)
        _dates = _dates.view(cls)
        _dates.freq = _freq
        return _dates
    
    def __array_wrap__(self, obj, context=None):
        if context is None:
            return self
        elif context[0].__name__ not in ufunc_dateOK:
            raise ArithmeticDateError, "(function %s)" % context[0].__name__
    
    def __array_finalize__(self, obj):
        self.freq = getattr(obj, 'freq', -9999)
        self._cachedinfo = dict(toobj=None, tostr=None, toord=None, 
                                steps=None, full=None, hasdups=None)
        if hasattr(obj,'_cachedinfo'):
            self._cachedinfo.update(obj._cachedinfo)
        return
    
    def __getitem__(self, indx):
        if isinstance(indx, Date):
            indx = self.find_dates(indx)
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
                r._cachedinfo.update(dict(steps=None, full=None, hasdups=None))
                for attr in ('tostr','toobj','toord'):
                    if r._cachedinfo[attr] is not None:
                        r._cachedinfo[attr] = r._cachedinfo[attr][indx]
            return r
        
    def __repr__(self):
        return ndarray.__repr__(self)
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
        return corelib.freq_tostr(self.freq)
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
        return numeric.asarray(cseries.getDateInfo(numeric.asarray(self), 
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
    def asfreq(self, freq=None, relation="BEFORE"):
        "Converts the dates to another frequency."
        # Note: As we define a new object, we don't need caching
        if freq is None or freq == -9999:
            return self
        tofreq = corelib.check_freq(freq)
        if tofreq == self.freq:
            return self        
        fromfreq = self.freq
        _rel = relation.upper()[0]
        new = cseries.asfreq(numeric.asarray(self), fromfreq, tofreq, _rel)
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
        if self.freq == 'U':
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
    if ddif[0] == ddif[-1] == 1.:
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
    dlist.sort()
    if dlist.ndim == 0:
        dlist.shape = (1,)
    # Case #1: dates as strings .................
    if dlist.dtype.kind == 'S':
        #...construct a list of ordinals
        ords = numpy.fromiter((DateTimeFromString(s).toordinal() for s in dlist),
                               float_)
        ords += 1
        #...try to guess the frequency
        if freq is None:
            freq = guess_freq(ords)
        #...construct a list of dates
        dates = [Date(freq, string=s) for s in dlist]
    # Case #2: dates as numbers .................
    elif dlist.dtype.kind in 'if':
        #...hopefully, they are values
        if freq is None:
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
            if freq is None:
                ords = numpy.fromiter((s.absdays for s in dlist), float_)
                ords += 1
                freq = guess_freq(ords)
            dates = [Date(freq, datetime=m) for m in dlist]
        #...as datetime objects
        elif hasattr(template, 'toordinal'):
            ords = numpy.fromiter((d.toordinal() for d in dlist), float_)
            if freq is None:
                freq = guess_freq(ords)
            dates = [Date(freq, datetime=dt.datetime.fromordinal(a)) for a in ords]
    #
    result = DateArray(dates, freq)
    return result


def date_array(dlist=None, start_date=None, end_date=None, length=None, 
               include_last=True, freq=None):
    """Constructs a DateArray from:
    - a starting date and either an ending date or a given length.
    - a list of dates.
    """
    freq = corelib.check_freq(freq)
    # Case #1: we have a list ...................
    if dlist is not None:
        # Already a DateArray....................
        if isinstance(dlist, DateArray):
            if (freq is not None) and (dlist.freq != corelib.check_freq(freq)):
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
        length = int(end_date - start_date)
        if include_last:
            length += 1
#    dlist = [(start_date+i).value for i in range(length)]
    dlist = numeric.arange(length, dtype=int_)
    dlist += start_date.value
    if freq is None:
        freq = start_date.freq
    return DateArray(dlist, freq=freq)
datearray = date_array

def date_array_fromlist(dlist, freq=None):
    "Constructs a DateArray from a list of dates."
    return date_array(dlist=dlist, freq=freq)

def date_array_fromrange(start_date, end_date=None, length=None, 
                         include_last=True, freq=None):
    """Constructs a DateArray from a starting date and either an ending date or 
    a length."""
    return date_array(start_date=start_date, end_date=end_date, 
                      length=length, include_last=include_last, freq=freq)    

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
    if 1:
        dlist = ['2007-%02i' % i for i in range(1,5)+range(7,13)]
        mdates = date_array_fromlist(dlist, 'M')
        # Using an integer
        assert_equal(mdates[0].value, 24073)
        assert_equal(mdates[-1].value, 24084)
        # Using a date
        lag = mdates.find_dates(mdates[0])
        print mdates[lag]
        assert_equal(mdates[lag], mdates[0])
    if 1:
        hodie = today('D')
        D = DateArray(today('D'))
        assert_equal(D.freq, 6000)
    
    if 1:
        freqs = [x[0] for x in corelib.freq_dict.values() if x[0] != 'U']
        print freqs
        for f in freqs:
            print f
            today = thisday(f)
            assert(Date(freq=f, value=today.value) == today)
    
    if 1:
        D = date_array(start_date=thisday('D'), length=5)
        Dstr = D.tostring()
        assert_equal(D.tostring(), Dstr)
        DL = D[[0,-1]]
        assert_equal(DL.tostring(), Dstr[[0,-1]])
    