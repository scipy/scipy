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


:author: Pierre GF Gerard-Marchant
:contact: pierregm_at_uga_edu
:date: $Date: 2006-12-05 20:40:46 -0500 (Tue, 05 Dec 2006) $
:version: $Id: timeseries.py 25 2006-12-06 01:40:46Z backtopop $
"""
__author__ = "Pierre GF Gerard-Marchant ($Author: backtopop $)"
__version__ = '1.0'
__revision__ = "$Revision: 25 $"
__date__     = '$Date: 2006-12-05 20:40:46 -0500 (Tue, 05 Dec 2006) $'


#import numpy as N

#from numpy.core import bool_, float_, int_
import numpy.core.umath as umath
import numpy.core.fromnumeric as fromnumeric
import numpy.core.numeric as numeric
from numpy import ndarray
from numpy.core.records import recarray
from numpy.core.records import fromarrays as recfromarrays

#from cPickle import dump, dumps

import core
reload(core)
from core import *


import maskedarray as MA
#reload(MA)
from maskedarray.core import domain_check_interval, \
    domain_greater_equal, domain_greater, domain_safe_divide, domain_tan  
from maskedarray.core import MaskedArray, MAError, masked_array, isMaskedArray,\
    getmask, getmaskarray, filled, mask_or, make_mask
from maskedarray.core import convert_typecode, masked_print_option, \
    masked_singleton
from maskedarray.core import load, loads, dump, dumps

import addons.numpyaddons
reload(addons.numpyaddons)
from addons.numpyaddons import ascondition

#...............................................................................
talog = logging.getLogger('TimeArray')
tslog = logging.getLogger('TimeSeries')
mtslog = logging.getLogger('MaskedTimeSeries')

nomask = MA.nomask
masked = MA.masked
masked_array = MA.masked_array
ufunc_domain = {}
ufunc_fills = {}


#### --------------------------------------------------------------------------
#--- ... Date Tools ...
#### --------------------------------------------------------------------------
dtmdefault = dtm.datetime(2001,1,1)

def isDate(obj):
    """Returns *True* if the argument is a `datetime.date` or `datetime.datetime`
instance."""
    return isinstance(obj,dtm.date) or isinstance(obj,dtm.datetime)
    
def fake_dates(count, start=dtmdefault):
    """fake_dates(count, start=dtmdefault)
Returns a *count x 1* array of daily `datetime` objects.
:Parameters:
    - `count` (Integer) : Number of dates to output
    - `start` (datetime object) : Starting date *[NOW]*."""
    dobj = list(dtmrule.rrule(DAILY,count=count,dtstart=start))
    return N.array(dobj)

#### --------------------------------------------------------------------------
#--- ... TimeArray class ...
#### --------------------------------------------------------------------------
class TimeArrayError(Exception):
    """Defines a generic DateArrayError."""
    def __init__ (self, args=None):
        "Create an exception"
        Exception.__init__(self)
        self.args = args
    def __str__(self):
        "Calculate the string representation"
        return str(self.args)
    __repr__ = __str__

class TimeArray(ndarray):
    """Stores datetime.dateime objects in an array."""
    def __new__(subtype, series, copy=False):
        vlev = 3
        if isinstance(series, TimeArray):
#            verbose.report("DA __new__: data isDA types %s, %s, %i" % \
#                           (type(series), series, len(series)), vlev)
            if not copy:            
                return series
            else: 
                return series.copy()
        else:
            data = N.array(series)
            if data.dtype.str == "|O8":
                new = data #.view(subtype)
                talog.debug("__new__: data isnotDA types %s, %s, %s, %i" % \
                            (type(series), data.dtype, series, data.size),)
            else:
                new = N.empty(data.shape, dtype="|O8")
                newflat = new.flat
                talog.debug("DA __new__: data isnotDA types %s, %s, %s, %i" % \
                            (type(series), data.dtype, series, data.size),)
                # Forces years (as integers) to be read as characters 
                if N.all(data < 9999):
                    data = data.astype("|S4")
                # Parse dates as characters
                if data.dtype.char == "S":
                    for (k,s) in enumerate(data.flat):
                        newflat[k] = dtmparser.parse(s, default=dtmdefault)
                else:
                    # Parse dates as ordinals
                    try:
                        for k,s in enumerate(data.flat):
                            newflat[k] = dtm.datetime.fromordinal(s)
                    except TypeError:
                        raise TypeError, "unable to process series !"
                subtype._asobject = new
        return new.view(subtype)
#    #............................................
#    def __array_wrap__(self, obj):
#        return TimeArray(obj)
    #............................................
    def __array_finalize__(self, obj):
        self._resolution = None
        (self._years, self._months, self._days, self._yeardays) = [None]*4
        self._asobjects = None
        self._asstrings = None
        self._asordinals = None
        return
    #............................................
    def __len__(self):
        "Returns the length of the object. Or zero if it fails."
        if self.ndim == 0:
            return 0
        return ndarray.__len__(self)
        
        
    def __getitem__(self,i):
        # Force a singleton to DateArray
        try:
            obj = ndarray.__getitem__(self,i)
        except IndexError:
            obj = self
        return TimeArray(obj)
#        if isDate(obj):
#            return DateArray(obj)
#        else:
#            return obj
    #............................................
    def __eq__(self, other):
        return N.asarray(self).__eq__(N.asarray(other))
    
    def __add__(self, other):
        if not isinstance(other, dtm.timedelta):
            raise TypeError, "Unsupported type for add operation"
        return TimeArray(N.asarray(self).__add__(other))
    
    def __sub__(self, other):
        if not isinstance(other, dtm.timedelta) and \
           not isinstance(other, TimeArray):
            raise TypeError, "Unsupported type for sub operation"
        return N.asarray(self).__sub__(N.asarray(other))       
    
    def __mul__(self, other):
        raise TimeArrayError, "TimeArray objects cannot be multiplied!"
    __rmul__ = __mul__
    __imul__ = __mul__
#    def shape(self):
#        if self.size == 1:
#            return (1,)
#        else:
#            return self._data.shape
    #............................................
    def __str__(self):
        """x.__str__() <=> str(x)"""
        return str(self.asstrings())
    def __repr__(self):
        """Calculate the repr representation, using masked for fill if
           it is enabled. Otherwise fill with fill value.
        """
        template = """\
timearray( 
 %(data)s,)"""
        template_short = """\
timearray(%(data)s)"""
        if self.ndim <= 1:
            return template_short % {'data': str(self)}
        else:
            return template % {'data': str(self)} 
    #............................................
    @property
    def resolution(self):
        """Calculates the initial resolution, as the smallest difference
between consecutive values (in days)."""
        if self._resolution is None:
            self._resolution = N.diff(self.asordinals().ravel()).min()
            if self._resolution <= 3:
                self._resolution = 1
            elif self._resolution <= 10:
                self._resolution = 7
            elif self._resolution <= 270:
                self._resolution = 30
            else:
                self._resolution = 365
        return self._resolution
    timestep = resolution
    #............................................
    def has_missing_dates(self, resolution=None):
        """Returns *True* there's a gap in the series, assuming a regular series
with a constant timestep of `resolution`.
If `resolution` is None, uses the default resolution of `self`.
"""
        if self.size < 2:
            return False
        dt = N.diff(self.asordinals().ravel())
        if resolution is None:
            resolution = self.resolution
        if resolution == 1:
            return (dt.max() > 1)
        elif resolution == 30:
            return (dt.max() > 31)
        elif resolution == 365:
            return (dt.max() > 730)
        else:
            #FIXME: OK, there's probly something wrong here...
            return True
    #............................................
    def asobjects(self):
        """Transforms an array of ordinal dates to an array of date objects."""
        if self._asobjects is None:
            if self.size == 1:
                self._asobjects = self.item()
            else:
                self._asobjects = self
        return self._asobjects
    #............................................
    def asordinals(self):
        """Transforms an array of dates to an array of date objects."""
        if self._asordinals is None:
            # Build list of datetime objects
            self._asordinals = N.empty(self.shape,dtype=int_)
            _asordinalsflat = self._asordinals.flat
            if self.size == 0:
                self._asordinals = dtm.datetime.toordinal(dtmdefault)
            elif self.size == 1:
                self._asordinals = dtm.datetime.toordinal(self.item())
            else:
                itr = (dtm.datetime.toordinal(val) for val in self.flat)
                self._asordinals = N.fromiter(itr, float_)
        return self._asordinals
    #............................................
    def asstrings(self, stringsize=10):
        """Transforms a *N*-array of ordinal dates to a *N*-list of  
datestrings `YYYY-MM-DD`.

:Parameters:
    `stringsize` : Integer *[10]*
        String size. 
        
        - `< 4' outputs 4
        - `< 8' outputs 8
        - anything else outputs 10.
        """ 
        if not stringsize:
            strsize = 10
        elif stringsize <= 4:
            strsize = 4
        elif stringsize <= 8:
            strsize = 7
        else:
            strsize = 10      
        if self._asstrings is None:
            deftype = "|S10"
            if self.size == 0:
                self._asstrings = N.array(dtmdefault.isoformat(), dtype=deftype)
            elif self.size == 1:
                self._asstrings = N.array(self.item().isoformat(), dtype=deftype)
            else:
                itr = (val.isoformat() for val in self.flat)
                self._asstrings = N.fromiter(itr,dtype=deftype).reshape(self.shape) 
        return self._asstrings.getfield('S%i' % strsize)
    #............................................  
    def astype(self, newtype):
        """Outputs the array in the type provided in argument."""
        newtype = N.dtype(newtype)
        if newtype.type == int_:
            return self.asordinals()
        elif newtype.type == str_:
            n = newtype.itemsize
            if n > 0 :
                return self.asstrings(stringsize=n)
            else:
                return self.asstrings(stringsize=n)
        elif newtype.type == object_:
            return self
        else:
            raise ValueError, "Invalid type: %s" % newtype
    #------------------------------------------------------
    @property
    def years(self):
        """Returns the corresponding year (as integer)."""
        if self._years is None:
            self._years = self.asstrings(1).astype(int_)
        return self._years
    @property
    def months(self):
        """Returns the corresponding month (as integer)."""
        if self._months is None:
            self._months = N.empty(self.shape, dtype=int_)
            _monthsflat = self._months.flat
            if self.size == 1:
                try:
                    _monthsflat[:] = self.asobjects().month
                except AttributeError:
                    _monthsflat[:] = 1
            else:
                for (k,val) in enumerate(self.asobjects().flat):
                    _monthsflat[k] = val.month
        return self._months    
    @property
    def days(self):
        """Returns the corresponding day of month (as integer)."""
        if self._days is None:
            self._days = N.empty(self.shape, dtype=int_)
            _daysflat = self._days.flat
            if self.size == 1:
                try:
                    _daysflat[:] = self.asobjects().day
                except AttributeError:
                    _daysflat[:] = 1
            else:
                for (k,val) in enumerate(self.asobjects().flat):
                    _daysflat[k] = val.day
        return self._days
    @property
    def yeardays(self):
        """Returns the corresponding day of year (as integer)."""
        if self._yeardays is None:
            self._yeardays = N.empty(self.size, dtype=int_)
            _doyflat = self._yeardays.flat
            for (k,val,yyy) in N.broadcast(N.arange(len(self)),
                                           self.asordinals(),
                                           self.years.flat):
                _doyflat[k] = val - dtm.datetime(yyy-1,12,31).toordinal()
            self._yeardays.reshape(self.shape).astype(int_)
        return self._yeardays
#................................................
def isTimeArray(x):
    """Checks whether `x` is n array of times/dates (an instance of `TimeArray` )."""
    return isinstance(x,TimeArray)

def time_array(t, copy=False):
    """Creates a date array from the series `t`.
    """
    if isinstance(t,TimeArray):
        dobj = t
    else:
        t = N.asarray(t)
        if t.dtype.str == "|O8":
            dobj = t
        else:
            dobj = N.empty(t.shape, dtype="|O8")
            dobjflat = dobj.flat
            if t.dtype.char == 'S':
                for k,s in enumerate(t.flat):
                    dobjflat[k] = dtmparser.parse(s)
            else:
                try:
                    for k,s in enumerate(t.flat):
                        dobjflat[k] = dtm.datetime.fromordinal(s)
                except TypeError:
                    raise TypeError, "unable to process series !"
    return TimeArray(dobj, copy=copy)
            
#### --------------------------------------------------------------------------
#--- ... Time Series ...
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
#......................................
#TODO: __eq__, __cmp__ classes
class TimeSeries(ndarray, object): 
    """Base class for the definition of time series.
A time series is here defined as the combination of two arrays:
    
    - an array storing the time information (as a `TimeArray` instance);
    - an array storing the data (as a `MaskedArray` instance.
    """
    def __new__(cls, data, dates=None, dtype=None, copy=True):
        tslog.debug("__new__: data types %s, %s" % (type(data), dtype))
#        if isinstance(data, TimeSeries):
        if hasattr(data,"_series") and hasattr(data,"_dates"):
            tslog.debug("__new__: data has _series and _dates")
            tslog.debug("__new__: setting basedates",)
            cls._basedates = data._dates
            tslog.debug("__new__: setting basedates done",)
            if (not copy) and (dtype == data.dtype):
                return N.asarray(data).view(cls)
            else:
                return N.array(data, dtype=dtype).view(cls)
        #..........
        dsize = N.size(data) 
        tslog.debug("__new__: args dates type %s, %s" % (type(dates), dates),)
        if dates is None:
            _dates = TimeArray(fake_dates(N.size(data)))
            tslog.debug("__new__: args dates FAKED type %s, %s, %i" % \
                       (type(dates), _dates, _dates.size),)
        else:
            _dates = TimeArray(dates)  
        tslog.debug("__new__: dates set to %s" % _dates,)
        if _dates.size != dsize:
            msg = "Incompatible sizes between dates (%i) and data (%i) !"
            raise ValueError, msg % (_dates.size, dsize)
        _dates.shape = N.shape(data)
        cls._basedates = _dates
        return N.array(data, copy=copy, dtype=dtype).view(cls)
    #..................................
    def __array_wrap__(self, obj):
        return TimeSeries(obj, dates=self._dates)
    #..................................
    def __array_finalize__(self,obj):
        if not hasattr(self,"_series"):
            try:
                self._series = obj._series
                tslog.debug("__array_finalize__: obj._data => : %s - %s" % \
                           (id(self._series), self._series.ravel() ),)
            except AttributeError:
                self._series = obj
                tslog.debug("__array_finalize__: obj       => : %s - %s" % \
                           (id(self._series), self._series.ravel() ),)
        if not hasattr(self,"_dates"):
            self._dates = self._basedates
            tslog.debug("__array_finalize__: dates set!  => : %s - %s" % \
                       (id(self._dates), self._dates.ravel() ))
        return
    #........................
    def _get_flat(self):
        """Calculate the flat value.
        """
        return self.__class__(self._series.ravel(), 
                              dates = self._dates.ravel(), copy=False)
    #....
    def _set_flat (self, value):
        "x.flat = value"
        y = self.ravel()
        y[:] = value
    flat = property(fget=_get_flat,fset=_set_flat,doc='Access array in flat form.')   
    #........................
    def _get_shape(self):
        "Return the current shape."
        return self._series.shape
    #....
    def _set_shape (self, newshape):
        "Set the array's shape."
        self._series.shape = self._dates.shape = newshape
    #....
    shape = property(fget=_get_shape, fset=_set_shape, doc="Array shape")
    #....
    @property
    def ndim(self):
        """Returns the number of dimensions."""
        return self._series.ndim
    @property
    def size(self):
        """Returns the total number of elements."""
        return self._series.size
    #........................
    def __getitem__(self, i):
        "Get item described by i. Not a copy as in previous versions."
        (tout, dout) = (self._dates[i], self._series[i])
        return self.__class__(dout, dates=tout)
    #....
    def __setitem__(self, index, value):
        "Gets item described by i. Not a copy as in previous versions."
        self._series[index] = value
        return self
    #........................
    def __getslice__(self, i, j):
        "Gets slice described by i, j"
        (tout, dout) = (self._dates[i:j], self._series[i:j])
        return self.__class__(dout, dates=tout,)
    #....
    def __setslice__(self, i, j, value):
        "Gets item described by i. Not a copy as in previous versions."
        #TODO: Here, we should have a better test for self._dates==value._dates
        if hasattr(value,"_dates"):
            tslog.debug("__setslice__: value: %s" % str(value))
            tslog.debug("__setslice__: dates: %s" % str(value._dates),)
            if not (self._dates == value._dates).all():
                raise TimeSeriesError,"Can't force a change of dates !"
        self._series[i:j] = value
        return self       
    #........................
    def ravel (self):
        """Returns a 1-D view of self."""
        return self.__class__(self._series.ravel(), 
                              dates=self._dates.ravel())
    #........................
    def __len__(self):
        if self.ndim == 0:
            return 0
        return ndarray.__len__(self)
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
 %(time)s, )
"""
        desc_short = """\
timeseries(data  = %(data)s,
           dates = %(time)s,)
"""
        if self.ndim <= 1:
            return desc_short % {
                'data': str(self._series),
                'time': str(self.dates),
                }
        return desc % {
            'data': str(self._series),
            'time': str(self.dates),
            }
 
    def ids (self):
        """Return the ids of the data, dates and mask areas"""
        return (id(self._series), id(self.dates),)
    
    @property
    def series(self):
        "Returnhs the series."
        return self._series
    
    #------------------------------------------------------
    @property
    def dates(self):
        """Returns the dates"""
        return self._dates
###    def _set_dates(self, object):
###        """Returns the dates"""
###        self._dates = object
###    dates = property(fget=_get_dates, fset=_set_dates, doc="Dates")
    @property
    def years(self):
        """Returns the corresponding years."""
        return self.dates.years
    @property
    def months(self):
        """Returns the corresponding months."""
        return self._dates.months
    @property
    def yeardays(self):
        """Returns the corresponding days of yuear."""
        return self._dates.yeardays

    def has_missing_dates(self):
        """Returns whether there's a date gap in the series."""
        return self._dates.has_missing_dates()

#### --------------------------------------------------------------------------
#--- ... Additional methods ...
#### --------------------------------------------------------------------------
class _tsmathmethod(object):
    """Defines a wrapper for arithmetic array methods (add, mul...).
When called, returns a new TimeSeries object, with the new series the result of
the method applied on the original series.
The `_dates` part remains unchanged.
    """
    def __init__ (self, binop):
        """abfunc(fillx, filly) must be defined.
           abinop(x, filly) = x for all x to enable reduce.
        """
        self.f = binop
    #
    def __get__(self, obj, objtype=None):
        self.obj = obj
        return self
    #
    def __call__ (self, other, *args):
        "Execute the call behavior."
        instance = self.obj
        _dates = instance._dates
        tslog.debug("_tsmathmethod: series: %s" % instance,)
        tslog.debug("_tsmathmethod: other  : %s" % other,)
        func = getattr(instance._series, self.f)
        if hasattr(instance,"_dates") and hasattr(other,"_dates"):
            tslog.debug("_tsmathmethod: instance  : %s" % instance,)
            tslog.debug("_tsmathmethod: other  : %s" % other,)
            if N.any(instance._dates != other._dates):
                tslog.debug("_tsmathmethod %s on %s and %s" % \
                               (self.f, instance, other),)
                raise TimeSeriesError, "Method %s not yet implemented !" % self.f
        return instance.__class__(func(other, *args), dates=_dates)  
#......................................
TimeSeries.__add__ = _tsmathmethod('__add__')
TimeSeries.__radd__ = _tsmathmethod('__add__')
TimeSeries.__sub__ = _tsmathmethod('__sub__')
TimeSeries.__rsub__ = _tsmathmethod('__rsub__')
TimeSeries.__pow__ = _tsmathmethod('__pow__')
TimeSeries.__mul__ = _tsmathmethod('__mul__')
TimeSeries.__rmul__ = _tsmathmethod('__mul__')
TimeSeries.__div__ = _tsmathmethod('__div__')
TimeSeries.__rdiv__ = _tsmathmethod('__rdiv__')
TimeSeries.__truediv__ = _tsmathmethod('__truediv__')
TimeSeries.__rtruediv__ = _tsmathmethod('__rtruediv__')
TimeSeries.__floordiv__ = _tsmathmethod('__floordiv__')
TimeSeries.__rfloordiv__ = _tsmathmethod('__rfloordiv__')
#................................................
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
        func_series = getattr(instance._series, _name)
        if self._ondates:
            func_dates = getattr(instance._dates, _name)
            return instance.__class__(func_series(*args), 
                                      dates=func_dates(*args))
        else:
            return instance.__class__(func_series(*args), 
                                      dates=instance._dates)  
#......................................
class _tsaxismethod(object):
    """Defines a wrapper for array methods working on an axis (mean...).
When called, returns a ndarray, as the result of the method applied on the original series.
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
        func = getattr(self.obj._series, self._name)
        return func(*args, **params)
#.......................................
TimeSeries.astype = _tsarraymethod('astype')
TimeSeries.reshape = _tsarraymethod('reshape', ondates=True)
TimeSeries.transpose = _tsarraymethod('transpose', ondates=True)
TimeSeries.swapaxes = _tsarraymethod('swapaxes', ondates=True)
TimeSeries.copy = _tsarraymethod('copy', ondates=True)
TimeSeries.compress = _tsarraymethod('compress', ondates=True)
#
TimeSeries.sum = _tsaxismethod('sum')
TimeSeries.cumsum = _tsaxismethod('cumsum')
TimeSeries.prod = _tsaxismethod('prod')
TimeSeries.cumprod = _tsaxismethod('cumprod')
TimeSeries.mean = _tsaxismethod('mean')
TimeSeries.var = _tsaxismethod('var')
TimeSeries.varu = _tsaxismethod('varu')
TimeSeries.std = _tsaxismethod('std')
TimeSeries.stdu = _tsaxismethod('stdu')


#..............................................................................
def concatenate(arrays, axis=0):
    """Concatenates a sequence of time series."""
    concatenate.__doc__ = N.concatenate.__doc__
    for a in arrays:
        if hasattr(a,"_dates"):
            raise TimeSeriesError, "Not yet implemented for TimeSeries !"

#### ---------------------------------------------------------------------------
#--- ... Pickling ...
#### ---------------------------------------------------------------------------
#FIXME: We're kinda stuck with forcing the mask to have the same shape as the data
def _tsreconstruct(baseclass, datesclass, baseshape, basetype):
    """Internal function that builds a new MaskedArray from the information stored
in a pickle."""
    _series = ndarray.__new__(ndarray, baseshape, basetype)
    _dates = ndarray.__new__(datesclass, baseshape, basetype)
    return TimeSeries.__new__(baseclass, _series, dates=_dates, dtype=basetype)

def _tsgetstate(a):
    "Returns the internal state of the TimeSeries, for pickling purposes."
    #TODO: We should prolly go through a recarray here as well.
    state = (1,
             a.shape, 
             a.dtype,
             a.flags.fnc,
             (a._series).__reduce__()[-1][-1],
             (a._dates).__reduce__()[-1][-1])
    return state
    
def _tssetstate(a, state):
    """Restores the internal state of the TimeSeries, for pickling purposes.
`state` is typically the output of the ``__getstate__`` output, and is a 5-tuple:

    - class name
    - a tuple giving the shape of the data
    - a typecode for the data
    - a binary string for the data
    - a binary string for the mask.
        """
    (ver, shp, typ, isf, raw, dti) = state
    (a._series).__setstate__((shp, typ, isf, raw))
    (a._dates).__setstate__((shp, N.dtype('|O8'), isf, dti))
    (a._dates)._asstrings = None
        
def _tsreduce(a):
    """Returns a 3-tuple for pickling a MaskedArray."""
    return (_tsreconstruct,
            (a.__class__, a._dates.__class__, (0,), 'b', ),
            a.__getstate__())

TimeSeries.__getstate__ = _tsgetstate
TimeSeries.__setstate__ = _tssetstate
TimeSeries.__reduce__ = _tsreduce
TimeSeries.__dump__ = dump
TimeSeries.__dumps__ = dumps

#................................................
def tofile(self, output, sep='\t', format='%s', format_dates=None):
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
    oformat = "%%s%s%s" % (sep,format)
    for (_dates,_data) in N.broadcast(self._dates.ravel().asstrings(), 
                                      filled(self)):
        ofile.write('%s\n' % sep.join([oformat % (_dates, _data) ]))
    ofile.close()
TimeSeries.tofile = tofile

#### --------------------------------------------------------------------------
#--- ... Shortcuts ...
#### --------------------------------------------------------------------------  
def isTimeSeries(x):
    """Checks whether `x` is a time series (an instance of `TimeSeries` )."""
    return isinstance(x, TimeSeries)        

#### --------------------------------------------------------------------------
#--- ... MaskedTimeSeries class ...
#### --------------------------------------------------------------------------
class MaskedTimeSeries(MaskedArray, TimeSeries): 
    """Base class for the definition of time series.
A time series is here defined as the combination of two arrays:
    
    - an array storing the time information (as a `TimeArray` instance);
    - an array storing the data (as a `MaskedArray` instance.
    """
    def __new__(cls, data, dates=None, mask=nomask, 
                dtype=None, copy=True, fill_value=-9999):
        mtslog.log(5, "__new__: data types %s, %s" % (type(data), dtype))
#        if isinstance(data, TimeSeries):
        #....................
        if isinstance(data, TimeSeries):  
            if isinstance(data, MaskedTimeSeries):
                _data = data._data
            else:
                _data = data
            _dates = data._dates
            _series = data._series
            mtslog.log(5, "__new__ from TS: data %i - %s - %s" % \
                          (id(_data._series), type(_data._series), _data.ravel()))
            mtslog.log(5,"__new__ from TS: dates %i - %s - %s" % \
                          (id(_dates), type(_dates), _dates.ravel()))
        elif isinstance(data, recarray):
            assert(data.dtype.names == ('_dates', '_series', '_mask'),
                   "Invalid fields names (got %s)" % (data.dtype.names,))
            _dates = data['_dates']
            _series = data['_series']
            _mask = data['_mask']
        else:
            if hasattr(data, "_data"):
                _data = TimeSeries(data._data, dates=dates, 
                                       dtype=dtype, copy=copy)
            else:
                _data = TimeSeries(data, dates=dates, 
                                       dtype=dtype, copy=copy)
            _dates = _data._dates    
            _series = _data._series
            mtslog.log(5,"__new__ from scratch: data %i - %s - %s" % \
                         (id(_data._series), type(_data._series), _data.ravel()))
            mtslog.log(5,"__new__ from TS: dates %i - %s - %s" % \
                         (id(_dates), type(_dates), _dates.ravel()))
        #.....................
        if mask is nomask:
            if hasattr(data, "_mask"):
                _mask = data._mask
            else:
                _mask = nomask
        else:
            _mask = make_mask(mask, copy=copy, flag=True)
        #....Check shapes compatibility
        if _mask is not nomask:
            (nd, nm) = (_data.size, _mask.size)
            if (nm != nd):
                if nm == 1:
                    _mask = N.resize(_mask, _data.shape)
                elif nd == 1:
                    _data = N.resize(_data, _mask.shape)
                else:
                    msg = "Mask and data not compatible (size issues: %i & %i)."
                    raise MAError, msg % (nm, nd)
            elif (_mask.shape != _data.shape):
                mtslog.log(5,"__new__ from scratch: force _mask shape %s > %s" % \
                             (_mask.shape, _data.shape))
                _mask.shape = _data.shape
        #....
        cls._fill_value = fill_value
        cls._basemask = _mask
        cls._basedates = _dates
        cls._baseseries = _series
        return _data.view(cls)
#        return _series.view(cls)
    #..............
    def __array_wrap__(self, obj, context=None):
        """Special hook for ufuncs.
Wraps the numpy array and sets the mask according to context.
        """
#        return MaskedArray.__array_wrap__(obj, context=None)
        return MaskedTimeSeries(obj, dates=self._dates, mask=self._mask,
                                fill_value=self._fill_value)
        
    #..............
    def __array_finalize__(self,obj):
        mtslog.log(5, "__array_finalize__: obj is %s" % (type(obj), ))
        if not hasattr(self, "_data"):
            self._data = obj
        if not hasattr(self, "_dates"):
            self._dates = self._basedates
            mtslog.log(5, "__array_finalize__: set dates to: %s - %s" % \
                          (id(self._dates), self._dates.ravel() ))
        if not hasattr(self, "_mask"):
            self._mask = self._basemask
            mtslog.log(5, "__array_finalize__: set mask to: %s - %s" % \
                          (id(self._mask), self._mask.ravel() ))
        if not hasattr(self, "_series"):
            if hasattr(obj, "_series"):
                self._series = obj._series
            else:
                self._series = obj
        self.fill_value = self._fill_value
        return

    #------------------------------------------------------
#    def __mul__(self):
    #------------------------------------------------------
    def __str__(self):
        """Calculate the str representation, using masked for fill if
           it is enabled. Otherwise fill with fill value.
        """
        if masked_print_option.enabled():
            f = masked_print_option
            # XXX: Without the following special case masked
            # XXX: would print as "[--]", not "--". Can we avoid
            # XXX: checks for masked by choosing a different value
            # XXX: for the masked singleton? 2005-01-05 -- sasha
            if self is masked:
                return str(f)
            m = self._mask
            if m is nomask:
                res = self._data
            else:
                if m.shape == () and m:
                    return str(f)
                # convert to object array to make filled work
                res = (self._series).astype("|O8")
                res[self._mask] = f
        else:
            res = self.filled(self.fill_value)
        return str(res)
    
    def __repr__(self):
        """Calculate the repr representation, using masked for fill if
           it is enabled. Otherwise fill with fill value.
        """
        desc = """\
timeseries(data =
 %(data)s,
           mask =
 %(mask)s, 
           date = 
 %(time)s, )
"""
        desc_short = """\
timeseries(data = %(data)s,
           mask = %(mask)s,
           date = %(time)s,)
"""
#        if (self._mask is nomask) and (not self._mask.any()):
#            if self.ndim <= 1:
#                return without_mask1 % {'data':str(self.filled()),
#                                        'time':str(self._dates.asstrings())}
#            return without_mask % {'data':str(self.filled()),
#                                   'time':str(self._dates.asstrings())}
#        else:
        if self.ndim <= 1:
            return desc_short % {
                'data': str(self),
                'mask': str(self._mask),
                'time': str(self.dates),
                }
        return desc % {
            'data': str(self),
            'mask': str(self._mask),
            'time': str(self.dates),
            }
    #............................................
    def ids (self):
        """Return the ids of the data, dates and mask areas"""
        return (id(self._series), id(self.dates), id(self._mask))
    #............................................
    @property
    def maskedseries(self):
        """Returns a masked array of the series (dates are omitteed)."""
        return masked_array(self._series, mask=self._mask)
    _mseries = maskedseries
    #............................................
    def filled(self, fill_value=None):
        """A numeric array with masked values filled. If fill_value is None,
           use self.fill_value().

           If mask is nomask, copy data only if not contiguous.
           Result is always a contiguous, numeric array.
# Is contiguous really necessary now?
        """
        (d, m) = (self._data, self._mask)
        if m is nomask:
            return d
        #
        if fill_value is None:
            value = self._fill_value
        else:
            value = fill_value
        #
        if self is masked_singleton:
            return numeric.array(value)
        #
        result = d.copy()
        try:
            result.__setitem__(m, value)
        except (TypeError, AttributeError):
            #ok, can't put that value in here
            value = numeric.array(value, dtype=object)
            d = d.astype(object)
            result = fromnumeric.choose(m, (d, value))
        except IndexError:
            #ok, if scalar
            if d.shape:
                raise
            elif m:
                result = numeric.array(value, dtype=d.dtype)
            else:
                result = d
        return result
    #............................................
    def sum(self, axis=None, dtype=None):
        """a.sum(axis=None, dtype=None) 
Sums the array `a` over the given axis `axis`.
Masked values are set to 0.
If `axis` is None, applies to a flattened version of the array.
    """
        if self._mask is nomask:
            return self._data.sum(axis, dtype=dtype)
        else:
            if axis is None:
                return self.filled(0).sum(None, dtype=dtype)
            return MaskedArray(self.filled(0).sum(axis, dtype=dtype),
                               mask=self._mask.all(axis))
            
    def cumsum(self, axis=None, dtype=None):
        """a.cumprod(axis=None, dtype=None)
Returns the cumulative sum of the elements of array `a` along the given axis `axis`. 
Masked values are set to 0.
If `axis` is None, applies to a flattened version of the array.
        """
        if self._mask is nomask:
            return self._data.cumsum(axis=axis, dtype=dtype)
        else:
            if axis is None:
                return self.filled(0).cumsum(None, dtype=dtype)
            return MaskedArray(self.filled(0).cumsum(axis=axis, dtype=dtype),
                               mask=self._mask)
        
    def prod(self, axis=None, dtype=None):
        """a.prod(axis=None, dtype=None)
Returns the product of the elements of array `a` along the given axis `axis`. 
Masked elements are set to 1.
If `axis` is None, applies to a flattened version of the array.
        """
        if self._mask is nomask:
            return self._data.prod(axis=axis, dtype=dtype)
        else:
            if axis is None:
                return self.filled(1).prod(None, dtype=dtype)
            return MaskedArray(self.filled(1).prod(axis=axis, dtype=dtype),
                               mask=self._mask.all(axis))
    product = prod
            
    def cumprod(self, axis=None, dtype=None):
        """a.cumprod(axis=None, dtype=None)
Returns the cumulative product of ethe lements of array `a` along the given axis `axis`. 
Masked values are set to 1.
If `axis` is None, applies to a flattened version of the array.
        """
        if self._mask is nomask:
            return self._data.cumprod(axis=axis, dtype=dtype)
        else:
            if axis is None:
                return self.filled(1).cumprod(None, dtype=dtype)
            return MaskedArray(self.filled(1).cumprod(axis=axis, dtype=dtype),
                               mask=self._mask)        
            
    def mean(self, axis=None, dtype=None):
        """mean(a, axis=None, dtype=None)
Returns the arithmetic mean.

The mean is the sum of the elements divided by the number of elements.
        """
        if self._mask is nomask:
            return self._data.mean(axis=axis, dtype=dtype)
        else:
            sum = N.sum(self.filled(0), axis=axis, dtype=dtype)
            cnt = self.count(axis=axis)
            if axis is None:
                if self._mask.all(None):
                    return masked
                else:
                    return sum*1./cnt
            return MaskedArray(sum*1./cnt, mask=self._mask.all(axis))

    def anom(self, axis=None, dtype=None):
        """a.anom(axis=None, dtype=None)
Returns the anomalies, or deviation from the average.
        """       
        m = self.mean(axis, dtype)
        if not axis:
            return (self - m)
        else:
            return (self - N.expand_dims(m,axis))
 
    def var(self, axis=None, dtype=None):
        """a.var(axis=None, dtype=None)
Returns the variance, a measure of the spread of a distribution.

The variance is the average of the squared deviations from the mean,
i.e. var = mean((x - x.mean())**2).
        """
        if self._mask is nomask:
            return MaskedArray(self._data.var(axis=axis, dtype=dtype),
                               mask=nomask)
        else:
            cnt = self.count(axis=axis)
            anom = self.anom(axis=axis, dtype=dtype)
            anom *= anom
            dvar = anom.sum(axis)
            dvar /= cnt
            if axis is None:
                return dvar
            return MaskedArray(dvar, mask=mask_or(self._mask.all(axis), (cnt==1)))
            
    def std(self, axis=None, dtype=None):
        """a.std(axis=None, dtype=None)
Returns the standard deviation, a measure of the spread of a distribution.

The standard deviation is the square root of the average of the squared
deviations from the mean, i.e. std = sqrt(mean((x - x.mean())**2)).
        """
        var = self.var(axis,dtype)
        if axis is None:
            if var is masked:
                return masked
            else:
                return N.sqrt(var)
        return MaskedArray(N.sqrt(var._data), mask=var._mask)
    
    def varu(self, axis=None, dtype=None):
        """a.var(axis=None, dtype=None)
Returns an unbiased estimate of the variance.

Instead of dividing the sum of squared anomalies by n, the number of elements,
this sum is divided by n-1.
        """
        cnt = self.count(axis=axis)
        anom = self.anom(axis=axis, dtype=dtype)
        anom *= anom
        var = anom.sum(axis)
        var /= (cnt-1)
        if axis is None:
            return var
        return MaskedArray(var, mask=mask_or(self._mask.all(axis), (cnt==1)))
            
    def stdu(self, axis=None, dtype=None):
        """a.var(axis=None, dtype=None)
Returns an unbiased estimate of the standard deviation.
        """
        var = self.varu(axis,dtype)
        if axis is None:
            if var is masked:
                return masked
            else:
                return N.sqrt(var)
        return MaskedArray(N.sqrt(var._data), mask=var._mask)
    #............................................
    def asrecords(self):
        """Returns the masked time series as a recarray.
Fields are `_dates`, `_data` and _`mask`.
        """
        desctype = [('_dates','|O8'), ('_series',self.dtype), ('_mask',N.bool_)]
        flat = self.ravel()
        if flat.size > 0:
            return recfromarrays([flat._dates, flat._series, getmaskarray(flat)],
                                 dtype=desctype,
                                 shape = (flat.size,),  
                                 )
        else:
            return recfromarrays([[], [], []], dtype=desctype, 
                                 shape = (flat.size,),  
                                 )
            
            
#    def reshape (self, *s):
#        """This array reshaped to shape s"""
#        self._data = self._data.reshape(*s)
#        self._dates = self._dates.reshape(*s)
#        if self._mask is not nomask:
#            self._mask = self._mask.reshape(*s)
#        return self.view()
#### --------------------------------------------------------------------------
#--- ... Pickling ...
#### --------------------------------------------------------------------------
def _mtsreconstruct(baseclass, datesclass, baseshape, basetype, fill_value):
    """Internal function that builds a new MaskedArray from the information stored
in a pickle."""
#    raise NotImplementedError,"Please use timeseries.archive/unarchive instead."""
    _series = ndarray.__new__(ndarray, baseshape, basetype)
    _dates = ndarray.__new__(datesclass, baseshape, '|O8')
    _mask = ndarray.__new__(ndarray, baseshape, '|O8')
    return baseclass.__new__(baseclass, _series, dates=_dates, mask=_mask, 
                             dtype=basetype, fill_value=fill_value)
#    
def _mtsgetstate(a):
    "Returns the internal state of the TimeSeries, for pickling purposes."
#    raise NotImplementedError,"Please use timeseries.archive/unarchive instead."""
    records = a.asrecords()
    state = (1,
             a.shape, 
             a.dtype,
             records.flags.fnc,
             a.fill_value,
             records
             )
    return state
#    
def _mtssetstate(a, state):
    """Restores the internal state of the TimeSeries, for pickling purposes.
`state` is typically the output of the ``__getstate__`` output, and is a 5-tuple:

    - class name
    - a tuple giving the shape of the data
    - a typecode for the data
    - a binary string for the data
    - a binary string for the mask.
        """
    (ver, shp, typ, isf, flv, rec) = state
    a.fill_value = flv
    a._data._series = a._series = N.asarray(rec['_series'])
    a._data._series.shape = a._series.shape = shp
    a._data._dates = a._dates = a._dates.__class__(rec['_dates'])
    a._data._dates.shape = a._dates.shape = shp
    (a._dates)._asstrings = None
    a._mask = N.array(rec['_mask'], dtype=MA.MaskType)
    a._mask.shape = shp
#        
def _mtsreduce(a):
    """Returns a 3-tuple for pickling a MaskedArray."""
    return (_mtsreconstruct,
            (a.__class__, a.dates.__class__, (0,), 'b', -9999),
            a.__getstate__())
#    
MaskedTimeSeries.__getstate__ = _mtsgetstate
MaskedTimeSeries.__setstate__ = _mtssetstate
MaskedTimeSeries.__reduce__ = _mtsreduce
MaskedTimeSeries.__dump__ = dump
MaskedTimeSeries.__dumps__ = dumps

##### -------------------------------------------------------------------------
#---- --- TimeSeries creator ---
##### -------------------------------------------------------------------------
def time_series(data, dates=None, mask=nomask, copy=False, fill_value=None):
    """Creates a TimeSeries object
    
:Parameters:
    `dates` : ndarray
        Array of dates.
    `data` : 
        Array of data.
    """
    if isinstance(data, MaskedTimeSeries):
        if not copy:
            data._mask = mask_or(data._mask, mask)
            return data
        _data = data._data
        _mask = mask_or(data._mask, mask)
        _dates = _data.dates
    elif isinstance(data, TimeSeries):
        _data = data._series
        _mask = make_mask(mask)
        _dates = data.dates    
    else:
        data = masked_array(data, copy=False)
        _data = data._data
        _mask = mask_or(data._mask, mask)
        if dates is None:
            _dates = fake_dates(data.size)
        else:
            _dates = time_array(dates)
    return MaskedTimeSeries(_data, dates=_dates, mask=_mask, copy=copy, 
                            fill_value=fill_value)


#

#### --------------------------------------------------------------------------
#--- ... Additional functions ...
#### --------------------------------------------------------------------------
def check_dates(a,b):
    """Returns the array of dates from the two objects `a` or `b` (or None)."""
    if isTimeSeries(a):
        if isTimeSeries(b) and (a._dates == b._dates).all() is False:
            raise ValueError, "Incompatible dates !"
        return a._dates
    elif isTimeSeries(b):
        return b._dates
    else:
        return
         
def parse_period(period):
    """Returns a TimeArray couple (starting date; ending date) from the arguments."""
####    print "........DEBUG PARSE DATES: period %s is %s" % (period, type(period))
#    if isinstance(period,TimeArray) or isinstance(period,Dates):
####        print "........DEBUG PARSE_PERIOD: OK"
    if isinstance(period,TimeArray):
        return (period[0],period[-1])
    elif hasattr(period,"__len__"):
        if not isinstance(period[0], TimeArray):
            tstart = TimeArray(period[0])
        else:
            tstart = period[0] 
        if not isinstance(period[-1], TimeArray):
            tend = TimeArray(period[-1])
        else:
            tend = period[-1] 
        return (tstart, tend)
    else:
        p = N.asarray(period)
        if N.all(p < 9999):
            p = N.array(period,dtype="|S4")
        p = time_array(p)
        return (p[0], p[-1])

def where_period(period, dates, *choices):
    """Returns choices fro True/False, whether dates fall during a given period.
If no choices are given, outputs the array indices  for the dates falling in the
period.

:Parameters:
    `period` : Sequence
        Selection period, as a sequence (starting date, ending date).
    `dates` : TimeArray
        Array of dates.
    `choices` : *(optional)*
        Arrays to select from when the condition is True/False.
    """
    (tstart, tend) = parse_period(period)
    condition = ascondition((dates>=tstart)&(dates<=tend))
    condition = (dates>=tstart)&(dates<=tend)
    return N.where(condition, *choices)

def masked_inside_period(data, period, dates=None):
    """Returns x as an array masked where dates fall inside the selection period,
as well as where data are initially missing (masked)."""
    (tstart, tend) = parse_period(period)
    # Get dates ..................
    if hasattr(data, "_dates"):
        dates = data._dates
    elif dates is None:
        raise ValueError,"Undefined dates !"
    else:
        assert(N.size(dates)==N.size(data), 
               "Inconsistent data and dates sizes!")
    # where_period yields True inside the period, when mask should yield False
    condition = ascondition(N.logical_and((dates>=tstart), (dates<=tend)))
    cm = filled(condition,True).reshape(data.shape)
    mask = mask_or(MA.getmaskarray(data), cm, copy=True)
    if isinstance(data, MaskedTimeSeries):
        return data.__class__(data._data, dates=dates, mask=mask, copy=True)
    if isinstance(data, TimeSeries):
        return MaskedTimeSeries(data, dates=dates, mask=mask, copy=True)
    else:
        return masked_array(data, mask=mask, copy=True) 

def masked_outside_period(data, period, dates=None):
    """Returns x as an array masked where dates fall outside the selection period,
as well as where data are initially missing (masked)."""
    (tstart, tend) = parse_period(period)
    if hasattr(data, "_dates"):
        dates = data._dates
    elif dates is None:
        raise ValueError,"Undefined dates !"
    else:
        assert(N.size(dates)==N.size(data), 
               "Inconsistent data and dates sizes!")
    #................
    condition = ascondition(N.logical_or((dates<tstart),(dates>tend)))
    cm = filled(condition,True).reshape(data.shape)
    mask = mask_or(MA.getmaskarray(data), cm, copy=True)
    if isinstance(data, MaskedTimeSeries):
        return data.__class__(data._data, dates=dates, mask=mask, copy=True)
    if isinstance(data, TimeSeries):
        return MaskedTimeSeries(data, dates=dates, mask=mask, copy=True)
    else:
        return masked_array(data, mask=mask, copy=True) 

#..............................................................................
def fill_missing_dates(dates,data,resolution=None,fill_value=None):
    """Finds and fills the missing dates in a time series, and allocates a 
default filling value to the data corresponding to the newly added dates.
    
:Parameters:
    `dates` 
        Initial array of dates.
    `data`
        Initial array of data.
    `resolution` : float *[None]*
        New date resolutions, in years. For example, a value of 1/365.25 indicates
        a daily resolution. If *None*, the initial resolution is used instead.
    `fill_value` : float
        Default value for missing data.
    """
    if not isinstance(dates, TimeArray):
        print "DEBUG FILL_MISSING_DATES: dates was %s" % type(dates)
        dates = TimeArray(dates)
    print "DEBUG FILL_MISSING_DATES: dates is %s" % type(dates)
    dflat = dates.ravel()
    n = len(dflat)
    # Get data ressolution .......
    if resolution is None:
        resolution = dflat.resolution
        if resolution >= 28 and resolution <= 31:
            resolution = 30
    else:
        resolution = int(1./float(resolution))
    # Find on what to fill .......
    if resolution == 1:
        (resol, freq, refdelta) = (DAILY, 1, dtmdelta.relativedelta(days=+1))    
    elif resolution == 7:
        (resol, freq, refdelta) = (WEEKLY, 1, dtmdelta.relativedelta(days=+7))
    elif resolution == 30:
        (resol, freq, refdelta) = (MONTHLY, 1, dtmdelta.relativedelta(months=+1))
    elif resolution == 365:
        (resol, freq, refdelta) = (YEARLY, 1, dtmdelta.relativedelta(years=+1))
    else:
        raise ValueError,\
              "Unable to define a proper date resolution (found %s)." % resolution
    # ...and now, fill it ! ......
    (tstart, tend) = dflat.asobjects()[[0,-1]].tolist()
    gaprule = dtmrule.rrule(resol, interval=freq, dtstart=tstart, until=tend)
    newdates = dates.__class__(list(gaprule))
    #.............................
    # Get the steps between consecutive data. We need relativedelta to deal w/ months
    delta = N.array([dtmdelta.relativedelta(b,a) 
                         for (b,a) in N.broadcast(dflat[1:],dflat[:-1])])
    dOK = N.equal(delta,refdelta)
    slcid = N.r_[[0,], N.arange(1,n).compress(-dOK), [n,]]
    oldslc = N.array([slice(i,e) for (i,e) in N.broadcast(slcid[:-1],slcid[1:])])
    if resolution == 1:
        addidx = N.cumsum([d.days for d in N.diff(dflat).compress(-dOK)])
    elif resolution == 30:
        addidx = N.cumsum([d.years*12+d.months for d in delta.compress(-dOK)])
    elif resolution == 365:
        addidx = N.cumsum([d.years for d in delta.compress(-dOK)])
    addidx -= N.arange(len(addidx))
    newslc = N.r_[[oldslc[0]], 
                  [slice(i+d-1,e+d-1) for (i,e,d) in \
                       N.broadcast(slcid[1:-1],slcid[2:],addidx)] 
                 ]
#    misslc = [slice(i,i+d-1) 
#                       for (i,d) in N.broadcast(slcid[1:-1],addidx)]
    #.............................
    # Just a quick check
    for (osl,nsl) in zip(oldslc,newslc):
        assert N.equal(dflat[osl],newdates[nsl]).all(),\
            "Slicing mishap ! Please check %s (old) and %s (new)" % (osl,nsl)
    #.............................
    data = MA.asarray(data)
    oldmask = MA.getmaskarray(data)
    newdata = N.empty(newdates.size,data.dtype)
    newmask = N.ones(newdates.size, bool_)
    if fill_value is None:
        if hasattr(data,'fill_value'):
            fill_value = data.fill_value
        else:
            fill_value = MA.default_fill_value(data)
    data = data.filled(fill_value)
    newdata.fill(fill_value)
    #....
    for (new,old) in zip(newslc,oldslc):
        newdata[new] = data[old]
        newmask[new] = oldmask[old]
#    for mis in misslc:
#        newdata[mis].fill(fill_value)
    # Get new shape ..............
    if data.ndim == 1:
        nshp = (newdates.size,)
    else:
        nshp = tuple([-1,] + list(data.shape[1:]))
    return MaskedTimeSeries(newdata.reshape(nshp),
                            dates=newdates.reshape(nshp),
                            mask=newmask.reshape(nshp),
                            fill_value=fill_value)

######--------------------------------------------------------------------------
##---- --- Archiving ---
######--------------------------------------------------------------------------
#import iodata.iotools as iotools
#def archive(timeseries,filename,compression=None):
#    data = timeseries.asrecords()
#    iotools.archive(data, filename, compression)
#
#def unarchive(filename):
#    raise NotImplementedError
    

###############################################################################