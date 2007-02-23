"""
Classes definition for the support of individual dates and array of dates.

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknow_ca_at_hotmail_dot_com
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import datetime
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

import mx.DateTime as mxD
from mx.DateTime.Parser import DateFromString as mxDFromString

import tcore as corelib
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
OriginDate = mxD.Date(1970)
secondlyOriginDate = OriginDate - mxD.DateTimeDeltaFrom(seconds=1)
minutelyOriginDate = OriginDate - mxD.DateTimeDeltaFrom(minutes=1)
hourlyOriginDate = OriginDate - mxD.DateTimeDeltaFrom(hours=1)

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
      
    - Use the `mxDate` keyword with an existing mx.DateTime.DateTime object, or 
      even a datetime.datetime object.
      
      >>> td.Date('D', mxDate=mx.DateTime.now())
      >>> td.Date('D', mxDate=datetime.datetime.now())
      """
    default_fmtstr = {'A': "%Y",
                      'Q': "%YQ%q",
                      'M': "%b-%Y",
                      'W': "%d-%b-%y",
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
                 mxDate=None):
        
        if hasattr(freq, 'freq'):
            self.freq = corelib.check_freq(freq.freq)
        else:
            self.freq = corelib.check_freq(freq)
        self.freqstr = corelib.freq_tostr(self.freq)
        
        
        if value is not None:
            if isinstance(value, str):
                self.mxDate = mxDFromString(value)
            elif self.freqstr == 'A':
                self.mxDate = mxD.Date(value, -1, -1)
            elif self.freqstr == 'B':
                valtmp = (value - 1)//5
                #... (value + valtmp*7 - valtmp*5)
                self.mxDate = mxD.DateTimeFromAbsDateTime(value + valtmp*2)
            elif self.freqstr in ['D','U']:
                self.mxDate = mxD.DateTimeFromAbsDateTime(value)
            elif self.freqstr == 'H':
                self.mxDate = hourlyOriginDate + mxD.DateTimeDeltaFrom(hours=value)
            elif self.freqstr == 'M':
                self.mxDate = mxD.DateTimeFromAbsDateTime(1) + \
                              mxD.RelativeDateTime(months=value-1, day=-1)
            elif self.freqstr == 'Q':
                self.mxDate = mxD.DateTimeFromAbsDateTime(1) + \
                              mxD.RelativeDateTime(years=(value // 4), 
                                                   month=((value * 3) % 12), day=-1)
            elif self.freqstr == 'S':
                self.mxDate = secondlyOriginDate + mxD.DateTimeDeltaFromSeconds(value)
            elif self.freqstr == 'T':
                self.mxDate = minutelyOriginDate + mxD.DateTimeDeltaFrom(minutes=value)
            elif self.freqstr == 'W':
                self.mxDate = mxD.Date(1,1,7) + \
                              mxD.RelativeDateTime(weeks=value-1)
        
        elif string is not None:
            self.mxDate = mxDFromString(string)  
            
        elif mxDate is not None:
            if isinstance(mxDate, datetime.datetime):
                mxDate = mxD.strptime(mxDate.isoformat()[:19], "%Y-%m-%dT%H:%M:%S")
            self.mxDate = truncateDate(self.freq, mxDate)
            
        else:
            # First, some basic checks.....
            if year is None:
                raise InsufficientDateError            
            if self.freqstr in 'BDWU':
                if month is None or day is None: 
                    raise InsufficientDateError
            elif self.freqstr == 'M':
                if month is None: 
                    raise InsufficientDateError
                day = -1
            elif self.freqstr == 'Q':
                if quarter is None: 
                    raise InsufficientDateError
                month = quarter * 3
                day = -1
            elif self.freqstr == 'A':
                month = -1
                day = -1
            elif self.freqstr == 'S':
                if month is None or day is None or second is None: 
                    raise InsufficientDateError
                
            if self.freqstr in ['A','B','D','M','Q','W']:
                self.mxDate = truncateDate(self.freq, mxD.Date(year, month, day))
                if self.freqstr == 'B':
                    if self.mxDate.day_of_week in [5,6]:
                        raise ValueError("Weekend passed as business day")
            elif self.freqstr in 'HTS':
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
                self.mxDate = truncateDate(self.freqstr,
                                           mxD.Date(year, month, day, 
                                                    hour, minute, second))
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
        # Annual .......
        if self.freqstr == 'A':
            val = self.mxDate.year
        # Business days.
        elif self.freqstr == 'B':
            days = self.mxDate.absdate
            weeks = days // 7
            val = days - weeks*2  
            # (weeks*5) + (days - weeks*7)
        # Daily/undefined
        elif self.freqstr in 'DU':
            val = self.mxDate.absdate
        # Hourly........
        elif self.freqstr == 'H':
            val = (self.mxDate - hourlyOriginDate).hours
        # Monthly.......
        elif self.freqstr == 'M':
            val = (self.mxDate.year-1)*12 + self.mxDate.month
        # Quarterly.....
        elif self.freqstr == 'Q':
            val = (self.mxDate.year-1)*4 + self.mxDate.month//3
        # Secondly......
        elif self.freqstr == 'S':
            val = (self.mxDate - secondlyOriginDate).seconds
        # Minutely......
        elif self.freqstr == 'T':
            val = (self.mxDate - minutelyOriginDate).minutes
        # Weekly........
        elif self.freqstr == 'W':
            val = self.mxDate.absdate//7
        return int(val)
    #......................................................        
    def strfmt(self, fmt):
        "Formats the date"
        if fmt is None:
            fmt = self.default_fmtstr[self.freqstr]
        qFmt = fmt.replace("%q", "XXXX")
        tmpStr = self.mxDate.strftime(qFmt)
        if "XXXX" in tmpStr: 
            tmpStr = tmpStr.replace("XXXX", str(self.quarter))
        return tmpStr
            
    def __str__(self):
        return self.strfmt(self.default_fmtstr[self.freqstr])

    def __repr__(self): 
        return "<%s : %s>" % (str(self.freqstr), str(self))
    #......................................................
    def toordinal(self):
        "Returns the date as an ordinal."
        return self.mxDate.absdays

    def fromordinal(self, ordinal):
        "Returns the date as an ordinal."
        return Date(self.freq, mxDate=mxD.DateTimeFromAbsDays(ordinal))
    
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
def truncateDate(freq, mxDate):
    """Chops off the irrelevant information from the mxDate passed in."""
    freqstr = corelib.check_freqstr(freq)
    if freqstr == 'A':
        return mxD.Date(mxDate.year)
    elif freqstr == 'Q':
        return mxD.Date(mxDate.year, monthToQuarter(mxDate.month)*3)
    elif freqstr == 'M':
        return mxD.Date(mxDate.year, mxDate.month)
    elif freqstr == 'W':
        d = mxDate.absdate
        return mxD.DateTimeFromAbsDateTime(d + (7 - d % 7) % 7)
    elif freqstr in 'BD':
        if freq == 'B' and mxDate.day_of_week in [5,6]:
            raise ValueError("Weekend passed as business day")
        return mxD.Date(mxDate.year, mxDate.month, mxDate.day)
    elif freqstr == 'H':
        return mxD.Date(mxDate.year, mxDate.month, mxDate.day, \
                        mxDate.hour)
    elif freqstr == 'T':
        return mxD.Date(mxDate.year, mxDate.month, mxDate.day, \
                        mxDate.hour, mxDate.minute)
    else:
        return mxDate
    
def monthToQuarter(monthNum):
    """Returns the quarter corresponding to the month `monthnum`.
    For example, December is the 4th quarter, Januray the first."""
    return int((monthNum-1)/3)+1

def thisday(freq):
    "Returns today's date, at the given frequency `freq`."
    freqstr = corelib.check_freqstr(freq)
    tempDate = mxD.now()
    # if it is Saturday or Sunday currently, freq==B, then we want to use Friday
    if freqstr == 'B' and tempDate.day_of_week >= 5:
        tempDate -= (tempDate.day_of_week - 4)
    if freqstr in ('B','D','H','S','T','W','U'):
        return Date(freq, mxDate=tempDate)
    elif freqstr == 'M':
        return Date(freq, year=tempDate.year, month=tempDate.month)
    elif freqstr == 'Q':
        return Date(freq, year=tempDate.year, quarter=monthToQuarter(tempDate.month))
    elif freqstr == 'A':
        return Date(freq, year=tempDate.year)
today = thisday

def prevbusday(day_end_hour=18, day_end_min=0):
    "Returns the previous business day."
    tempDate = mxD.localtime()
    dateNum = tempDate.hour + float(tempDate.minute)/60
    checkNum = day_end_hour + float(day_end_min)/60
    if dateNum < checkNum: 
        return thisday('B') - 1
    else: 
        return thisday('B')
                
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

    if date.freqstr == 'U':
        warnings.warn("Undefined frequency: assuming daily!")
        fromfreq = corelib.freq_revdict['D']
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
    (_tostr, _toord, _steps) = (None, None, None)
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
        self.freq = getattr(obj, 'freq', self._defaultfreq)
        self.freqstr = getattr(obj, 'freqstr', corelib.freq_tostr(self.freq))
        for attr in ('_toobj','_toord','_tostr',
                     '_steps','_full','_hasdups'):
            setattr(self, attr, getattr(obj, attr, None))
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
        return numeric.asarray(cseries.getDateInfo(numeric.asarray(self), self.freq, info), dtype=int_)
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
        if self._toord is None:
#            diter = (Date(self.freq, value=d).toordinal() for d in self)
            diter = (d.toordinal() for d in self)
            toord = numeric.fromiter(diter, dtype=float_)
            self._toord = toord
        return self._toord
    #
    def tostring(self):
        "Converts the dates to strings."
        # Note: we better cache the result
        if self._tostr is None:
            firststr = str(self[0])
            if self.size > 0:
                ncharsize = len(firststr)
                tostr = numpy.fromiter((str(d) for d in self),
                                        dtype='|S%i' % ncharsize)
            else:
                tostr = firststr
            self._tostr = tostr
        return self._tostr
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
            raise ValueError, "Date out of bounds!"
        return c  

    def date_to_index(self, date):
        "Returns the index corresponding to one given date, as an integer."
        if self.isvalid():
            index = date.value - self[0].value
            if index < 0 or index > self.size:
                raise ValueError, "Date out of bounds!"
            return index
        else:
            index_asarray = (self == date.value).nonzero()
            if fromnumeric.size(index_asarray) == 0:
                raise ValueError, "Date out of bounds!" 
            return index_asarray[0][0]
    #......................................................        
    def get_steps(self):
        """Returns the time steps between consecutive dates.
    The timesteps have the same unit as the frequency of the series."""
        if self.freq == 'U':
            warnings.warn("Undefined frequency: assuming integers!")
        if self._steps is None:
            val = numeric.asarray(self).ravel()
            if val.size > 1:
                steps = val[1:] - val[:-1]
                if self._full is None:
                    self._full = (steps.max() == 1)
                if self._hasdups is None:
                    self._hasdups = (steps.min() == 0)
            else:
                self._full = True
                self._hasdups = False
                steps = numeric.array([], dtype=int_)
            self._steps = steps
        return self._steps
    
    def has_missing_dates(self):
        "Returns whether the DateArray have missing dates."
        if self._full is None:
            steps = self.get_steps()
        return not(self._full)
    
    def isfull(self):
        "Returns whether the DateArray has no missing dates."
        if self._full is None:
            steps = self.get_steps()
        return self._full
    
    def has_duplicated_dates(self):
        "Returns whether the DateArray has duplicated dates."
        if self._hasdups is None:
            steps = self.get_steps()
        return self._hasdups
    
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
        fcode = 'D'
    elif (ddif[0] == 1.) and (ddif[-1] == 3.):
        fcode = 'B'
    elif (ddif[0] > 3.) and  (ddif[-1] == 7.):
        fcode = 'W'
    elif (ddif[0] >= 28.) and (ddif[-1] <= 31.):
        fcode = 'M'
    elif (ddif[0] >= 90.) and (ddif[-1] <= 92.):
        fcode = 'Q'
    elif (ddif[0] >= 365.) and (ddif[-1] <= 366.):
        fcode = 'A'
    elif numpy.abs(24.*ddif[0] - 1) <= 1e-5 and \
         numpy.abs(24.*ddif[-1] - 1) <= 1e-5:
        fcode = 'H'
    elif numpy.abs(1440.*ddif[0] - 1) <= 1e-5 and \
         numpy.abs(1440.*ddif[-1] - 1) <= 1e-5:
        fcode = 'T'
    elif numpy.abs(86400.*ddif[0] - 1) <= 1e-5 and \
         numpy.abs(86400.*ddif[-1] - 1) <= 1e-5:
        fcode = 'S'
    else:
        warnings.warn("Unable to estimate the frequency! %.3f<>%.3f" %\
                      (ddif[0], ddif[-1]))
        fcode = 'U'
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
        ords = numpy.fromiter((mxDFromString(s).absdays for s in dlist),
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
            dates = [Date(freq, mxDate=m) for m in dlist]
        #...as datetime objects
        elif hasattr(template, 'toordinal'):
            ords = numpy.fromiter((d.toordinal() for d in dlist), float_)
            if freq is None:
                freq = guess_freq(ords)
            dates = [Date(freq, mxDate=mxD.DateTimeFromAbsDays(a)) for a in ords]
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
    