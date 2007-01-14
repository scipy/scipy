import numpy
from numpy import ma

import corelib
import cseries
import tsdate
import copy as copytools

from types import MethodType

masked = ma.masked
nomask = ma.nomask

def ts_compatible(a, b):
    if a.freq != b.freq:
        raise ValueError("Both TimeSeries must have same freq!")
    elif a.start_date() != b.start_date():
        raise ValueError("Both TimeSeries must have same start_date!")
    elif a.shape != b.shape:
        raise ValueError("Both TimeSeries must be of the same size!")


class ts_unary_operation:
    def __init__ (self, abfunc):
        self.f = abfunc
        self.__doc__ = getattr(abfunc, "__doc__", str(abfunc))

    def __call__ (self, a, *args, **kwargs):
        "Execute the call behavior."
        if isinstance(a, TimeSeries):
            return TimeSeries(self.f(a, *args, **kwargs), freq=a.freq, observed=a.observed, start_date=a.start_date())
        else:
            return self.f(a, *args, **kwargs)
        
        
class ts_binary_operation:
    def __init__ (self, abfunc):
        self.f = abfunc
        self.__doc__ = getattr(abfunc, "__doc__", str(abfunc))

    def __call__ (self, a, b, *args, **kwargs):
        "Execute the call behavior."

        if isinstance(a, TimeSeries) and isinstance(b, TimeSeries):
            ts_compatible(a, b)
            return TimeSeries(self.f(a, b, *args, **kwargs), freq=a.freq, observed=a.observed, start_date=a.start_date())
        elif isinstance(a, TimeSeries):
            if corelib.isDateType(a.tstype):
                return TimeSeries(self.f(a, b, *args, **kwargs), dtype=a.tstype, freq=a.freq, observed=a.observed, start_date=a.start_date())
            else:
                return TimeSeries(self.f(a, b, *args, **kwargs), freq=a.freq, observed=a.observed, start_date=a.start_date())
        elif isinstance(b, TimeSeries):
            if corelib.isDateType(b.tstype):
                return TimeSeries(self.f(a, b, *args, **kwargs), dtype=b.tstype, freq=b.freq, observed=b.observed, start_date=b.start_date())
            else:
                return TimeSeries(self.f(a, b, *args, **kwargs), freq=b.freq, observed=b.observed, start_date=b.start_date())
        else:
            return self.f(a, b, *args, **kwargs)
            
    def reduce (self, target, axis=0, dtype=None):
        return self.f.reduce(target, axis, dtype)

    def outer (self, a, b):
        return self.f.outer(a, b)

    def accumulate (self, target, axis=0):
        return datawrap(self.f.accumulate(target, axis), target)

    def __str__ (self):
        return "Masked version of " + str(self.f)            


class TimeSeries(ma.MaskedArray):

    __array_priority__ = 10.2

    def __init__(self, data, dtype=None, freq=None, start_date=None, observed=None, copy=True, order=False, mask=ma.nomask, fill_value=None):
    
        if isinstance(data, TimeSeries):
            if freq is None: freq = data.freq
            if start_date is None: start_date = data.start_date()
            if observed is None: observed = data.observed
        else:
            if observed is None: observed = 'END'
        
        self.freq = corelib.fmtFreq(freq)

        if isinstance(start_date, tsdate.Date):
            if start_date.freq != self.freq: raise ValueError("frequency of start_date must match frequency of series")
            else: self.__start_date = start_date
        else:
            self.__start_date = tsdate.Date(freq=self.freq, value=start_date)

        self.observed = corelib.fmtObserv(observed)

        self.tstype = None

        if corelib.isDateType(dtype) or (isinstance(data, TimeSeries) and corelib.isDateType(data.tstype)):
            self.tstype = dtype
            dtype = numpy.int_

        super(TimeSeries, self).__init__(data=data, dtype=dtype, copy=copy, order=order, mask=mask, fill_value=fill_value)
        
        if self.tstype is None: self.tstype = self.dtype


    def __getitem__(self, key):
        return super(TimeSeries, self).__getitem__(self.__prepKey(key))
        
    def __setitem__(self, key, value):
        super(TimeSeries, self).__setitem__(self.__prepKey(key), value)

    def __prepKey(self, key):
    
        if isinstance(key, tsdate.Date):
            key = int(key - self.start_date())
            if key < 0: raise ValueError("Date out of bounds")
            else: return key

        elif isinstance(key, TimeSeries):
            if corelib.isDateType(key.tstype):
                if key.tstype.freq != self.freq:
                    raise ValueError("series of frequency "+str(self.freq)+" given date expression of type "+str(key.tstype.freq))

                if key.mask is ma.nomask: key = numpy.asarray(key) - int(self.start_date())
                else: key = numpy.asarray(key[key.mask == False]) - int(self.start_date())
                
                if len(numpy.where(key < 0)[0]) > 0: raise ValueError("Indices out of bounds")
                
                return key
                
            else:

                # frequency, size, and start_date of key must all match self
                # when the data type is note a date
                ts_compatible(key, self)

                if key.tstype is numpy.bool_:
                    key = key.filled(False)
                elif numpy.ravel(key.mask).any():
                    raise ValueError("masked values cannot be used as indices!")

                return numpy.asarray(key)
        
        elif isinstance(key, ma.MaskedArray):

            if key.dtype is numpy.bool_:
                key = key.filled(False)
            elif numpy.ravel(key.mask).any():
                raise ValueError("masked values cannot be used as indices!")

            return numpy.asarray(key)
        
        else: return key

    
    def convert(self, freq, func='auto', position='END', interp=None):
        """
        return self converted to freq.
        
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
        placed at the end of the month). Interp is the method that will be used
        to fill in the gaps. Valid values are "CUBIC", "LINEAR", "CONSTANT", "DIVIDED",
        and None.
        
        Note: interp currently not implemented
        
        """
        
        if position.upper() not in ('END','START'): raise ValueError("invalid value for position argument: (%s)",str(position))
        
        toFreq = corelib.fmtFreq(freq)
        fromFreq = self.freq
        
        if fromFreq != toFreq:
        
            if func == 'auto':
                func = corelib.obsDict[self.observed]

            if self.size == 0:
                return TimeSeries(self, freq=toFreq, start_date=tsdate.dateOf(self.start_date(), toFreq))

            tempData = self.filled()

            if self.mask is ma.nomask:
                tempMask = numpy.empty(tempData.shape, dtype=numpy.bool_)
                tempMask[:] = False
            else: tempMask = self.mask

            cRetVal = cseries.convert(tempData, fromFreq, toFreq, position, int(self.start_date()), tempMask)

            _values = cRetVal['values']
            _mask = cRetVal['mask']
            startIndex = cRetVal['startindex']
            
            tempData = ma.array(_values)
            tempMask = ma.make_mask(_mask)
            tempData[tempMask] = ma.masked

            if func is not None and tempData.ndim == 2:
                tempData = corelib.apply_along_axis(func, 1, tempData)

            return TimeSeries(tempData, freq=toFreq, observed=self.observed, start_date=startIndex)
            
        else:
            return copytools.deepcopy(self)


    def adjust_endpoints(self, start_date=None, end_date=None):
        self.__init__(adjust_endpoints(self, start_date=start_date, end_date=end_date))

        
    def __str__(self):
        retVal = ""

        if self.shape[0] > 0:
            for i in range(self.shape[0]):
                index = str(self.start_date() + i)
                index = index + (" " * (6-len(index)))
                retVal += index + " --> " + str(self[i])+"\n"
            return retVal
        else:
            return "<no data>"
            
            
    def first_value(self, asDate=False):
        firstIndex = corelib.first_unmasked(self)
        if asDate:
            return self.start_date() + firstIndex
        else:
            return firstIndex
        
    def last_value(self, asDate=False):
        lastIndex = corelib.last_unmasked(self)
        if asDate:
            return self.start_date() + lastIndex
        else:
            return lastIndex
            
    def start_date(self): return self.__start_date
    def end_date(self): return self.__start_date + (self.shape[0] - 1)
            
    def date_to_index(self, date):
        if date.freq != self.freq: raise ValueError("date.freq != self.freq")
        return date - self.start_date()
            
            
    # built-in methods

    def __and__(self, other): return bitwise_and(self, other)
    def __or__(self, other): return bitwise_or(self, other)
    def __xor__(self, other): return bitwise_xor(self, other)
    __rand__ = __and__
    __ror__ = __or__
    __rxor__ = __xor__
    def __abs__(self): return absolute(self)
    def __neg__(self): return negative(self)
    def __pos__(self): return TimeSeries(self)
    def __add__(self, other): return add(self, other)
    __radd__ = __add__
    def __mod__ (self, other): return remainder(self, other)
    def __rmod__ (self, other): return remainder(other, self)
    def __lshift__ (self, n): return left_shift(self, n)
    def __rshift__ (self, n): return right_shift(self, n)
    def __sub__(self, other): return subtract(self, other)
    def __rsub__(self, other): return subtract(other, self)
    def __mul__(self, other): return multiply(self, other)
    __rmul__ = __mul__
    def __div__(self, other): return divide(self, other)
    def __rdiv__(self, other): return divide(other, self)
    def __truediv__(self, other): return true_divide(self, other)
    def __rtruediv__(self, other): return true_divide(other, self)
    def __floordiv__(self, other): return floor_divide(self, other)
    def __rfloordiv__(self, other): return floor_divide(other, self)
    def __pow__(self, other, third=None): return power(self, other, third)
    def __sqrt__(self): return sqrt(self)

    def __iadd__(self, other):
        return self + other

    def __imul__(self, other):
        return self * other

    def __isub__(self, other):
        return self - other

    def __idiv__(self, other):
        return self / other

    def __eq__(self, other): return equal(self,other)
    def __ne__(self, other): return not_equal(self,other)
    def __lt__(self, other): return less(self,other)
    def __le__(self, other): return less_equal(self,other)
    def __gt__(self, other): return greater(self,other)
    def __ge__(self, other): return greater_equal(self,other)

    def astype (self, tc):
        "return self as array of given type."
        d = self._data.astype(tc)
        return datawrap(ma.array(d, mask=self._mask), self)

    def filled (self, fill_value=None, ts=False):
        d = super(TimeSeries, self).filled(fill_value)
        if ts: return datawrap(d, self)
        else: return d


def datawrap(data, ts): return TimeSeries(data, freq=ts.freq, observed=ts.observed, start_date=ts.start_date())

## wrappers for numpy.ma funcs

sqrt = ts_unary_operation(ma.sqrt)
log = ts_unary_operation(ma.log)
log10 = ts_unary_operation(ma.log10)
exp = ts_unary_operation(ma.exp)
sin = ts_unary_operation(ma.sin)
cos = ts_unary_operation(ma.cos)
tan = ts_unary_operation(ma.tan)
arcsin = ts_unary_operation(ma.arcsin)
arccos = ts_unary_operation(ma.arccos)
arctan = ts_unary_operation(ma.arctan)

def cumprod(self, axis=0, dtype=None, out=None): return datawrap(ma._cumprod(self, axis, dtype, out), self)
def cumsum(self, axis=0, dtype=None, out=None): return datawrap(ma._cumsum(self, axis, dtype, out), self)

arcsinh = ts_unary_operation(ma.arcsinh)
arccosh = ts_unary_operation(ma.arccosh)
arctanh = ts_unary_operation(ma.arctanh)
sinh = ts_unary_operation(ma.sinh)
cosh = ts_unary_operation(ma.cosh)
tanh = ts_unary_operation(ma.tanh)
absolute = ts_unary_operation(ma.absolute)
fabs = ts_unary_operation(ma.fabs)
negative = ts_unary_operation(ma.negative)
nonzero = ts_unary_operation(ma.nonzero)

around = ts_unary_operation(ma.around)
floor = ts_unary_operation(ma.floor)
ceil = ts_unary_operation(ma.ceil)
logical_not = ts_unary_operation(ma.logical_not)

def zeros(shape, dtype=float, freq=None, start_date=None, observed=None):
    return TimeSeries(ma.zeros(shape, dtype), freq=freq, start_date=start_date, observed=observed)
def ones(shape, dtype=float, freq=None, start_date=None, observed=None):
    return TimeSeries(ma.ones(shape, dtype), freq=freq, start_date=start_date, observed=observed)


# functions from ma that we want to return scalars or masked arrays
count = ma.count
sum = ma.sum
product = ma.product
average = ma.average
compress = ma.compress
minimum = ma.minimum
maximum = ma.maximum
alltrue = ma.alltrue
allclose = ma.allclose
allequal = ma.allequal
sometrue = ma.sometrue
std = ma._std
var = ma._var

def argmin (x, axis = -1, out=None, fill_value=None):
    # same as argmin for ma, but returns a date instead of integer
    return x.start_date() + ma.argmin(x, axis, out, fill_value)    

def argmax (x, axis = -1, out=None, fill_value=None):
    # same as argmax for ma, but returns a date instead of integer
    return x.start_date() + ma.argmax(x, axis, out, fill_value)
    

# binary operators
add = ts_binary_operation(ma.add)
subtract = ts_binary_operation(ma.subtract)
multiply = ts_binary_operation(ma.multiply)
divide = ts_binary_operation(ma.divide)
power = ts_binary_operation(ma.power)
true_divide = ts_binary_operation(ma.true_divide)
floor_divide = ts_binary_operation(ma.floor_divide)
remainder = ts_binary_operation(ma.remainder)
fmod = ts_binary_operation(ma.fmod)
hypot = ts_binary_operation(ma.hypot)
arctan2 = ts_binary_operation(ma.arctan2)
equal = ts_binary_operation(ma.equal)
not_equal = ts_binary_operation(ma.not_equal)
less_equal = ts_binary_operation(ma.less_equal)
greater_equal = ts_binary_operation(ma.greater_equal)
less = ts_binary_operation(ma.less)
greater = ts_binary_operation(ma.greater)
logical_and = ts_binary_operation(ma.logical_and)
logical_or = ts_binary_operation(ma.logical_or)
logical_xor = ts_binary_operation(ma.logical_xor)
bitwise_and = ts_binary_operation(ma.bitwise_and)
bitwise_or = ts_binary_operation(ma.bitwise_or)
bitwise_xor = ts_binary_operation(ma.bitwise_xor)

def left_shift (a, n): return datawrap(ma.left_shift(a, n), a)
def right_shift (a, n): return datawrap(ma.right_shift(a, n), a)

def masked_where(condition, x, copy=1): return datawrap(ma.masked_where(condition, x, copy), x)
def masked_greater(x, value, copy=1): return datawrap(ma.masked_greater(x, value, copy), x)
def masked_greater_equal(x, value, copy=1): return datawrap(ma.masked_greater_equal(x, value, copy), x)
def masked_less(x, value, copy=1): return datawrap(ma.masked_less(x, value, copy), x)
def masked_less_equal(x, value, copy=1): return datawrap(ma.masked_less_equal(x, value, copy), x)
def masked_not_equal(x, value, copy=1): return datawrap(ma.masked_not_equal(x, value, copy), x)
def masked_equal(x, value, copy=1): return datawrap(ma.masked_equal(x, value, copy), x)
def masked_inside(x, v1, v2, copy=1): return datawrap(ma.masked_inside(x, v1, v2, copy), x)
def masked_outside(x, v1, v2, copy=1): return datawrap(ma.masked_outside(x, v1, v2, copy), x)

def clip(self,a_min,a_max,out=None): return datawrap(ma._clip(self, a_min, a_max, out=None), self)


array = TimeSeries

def _m(f):
    return MethodType(f, None, array)
    
array.clip = _m(clip)
array.argmax = _m(argmax)
array.argmin = _m(argmin)
array.cumprod = _m(cumprod)
array.cumsum = _m(cumsum)


# time series specific functions

def tser(start, end):
    if start.freq != end.freq:
        raise ValueError("start and end dates must have same frequency!")
    return TimeSeries(numpy.arange(int(start), int(end)+1), dtype=corelib.freqTypeMapping[start.freq], freq=start.freq, start_date=start)

def year(dateSer):
    return __getDateInfo(dateSer,'Y')

def quarter(dateSer):
    return __getDateInfo(dateSer,'Q')
    
def month(dateSer):
    return __getDateInfo(dateSer,'M')
    
def day(dateSer):
    return __getDateInfo(dateSer,'D')
    
def day_of_week(dateSer):
    return __getDateInfo(dateSer,'W')

def __getDateInfo(dateSer,infoCode):
    newData = ma.array(cseries.getDateInfo(dateSer.filled(), dateSer.tstype.freq, infoCode))
    if dateSer.mask is not ma.nomask:
        newData[dateSer.mask] = ma.masked
    return datawrap(newData, dateSer)


def adjust_endpoints(a, start_date=None, end_date=None):
    """adjust_endpoints(a, start_date=None, end_date=None) returns a new
    TimeSeries going from start_date to end_date"""
    
    if start_date is None: start_date = a.start_date()
    if end_date is None: end_date = a.end_date()

    tmpShape = list(a.shape)
    tmpShape[0] = max(end_date - start_date + 1, 0)
    tmpShape = tuple(tmpShape)
    
    tmpSer = TimeSeries(ma.resize(a, tmpShape), freq=a.freq, observed=a.observed, start_date=start_date)
    
    setStart, setEnd = max(start_date, a.start_date()), min(end_date, a.end_date())
    setLen = setEnd - setStart
    
    tmpSer[:] = ma.masked
    
    if setLen >= 0:
        tmpSer[tmpSer.date_to_index(setStart):tmpSer.date_to_index(setEnd)+1] = a[a.date_to_index(setStart):a.date_to_index(setEnd)+1]
            
    return tmpSer


def aligned(*series, **kwargs):
    
    if len(series) < 2:
        return series
        
    freq = series[0].freq
    
    if len(set([x.freq for x in series])) > 1: raise ValueError("All series must have same frequency!")
    
    if 'start_date' in kwargs: start_date = kwargs['start_date']
    else: start_date = min([x.start_date() for x in series])
    
    if 'end_date' in kwargs: end_date = kwargs['end_date']
    else: end_date = max([x.end_date() for x in series])
    
    return [adjust_endpoints(x, start_date=start_date, end_date=end_date) for x in series]
    