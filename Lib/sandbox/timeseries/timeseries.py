import numpy
from numpy import ma

import corelib
import shiftingarray as sa
from shiftingarray import doFunc, doFunc_oneseries
import cseries
import tsdate
import copy

class TimeSeries(sa.ShiftingArray):
    def __init__(self, values=[], dtype=None, freq=None, observed='END', startIndex=None, mask=ma.nomask):
    
        if freq is None: raise ValueError("freq not specified")
        
        if dtype is None: dtype = values.dtype

        super(TimeSeries, self).__init__(values, dtype, startIndex, mask)
        self.freq = corelib.fmtFreq(freq)
        self.observed = corelib.fmtObserv(observed)
        self.dtype = dtype
        
    def __getitem__(self, key):
        if isinstance(key,tsdate.Date):
            if self.freq != key.freq:
                raise "series of frequency "+str(self.freq)+" given date expression of type "+str(key.freq)
            else:
                key = int(key)
        return super(TimeSeries, self).__getitem__(key)
        
    def __setitem__(self, key, value):
        if isinstance(key, tsdate.Date):
            key = int(key)
        super(TimeSeries, self).__setitem__(key, value)

    
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

            firstIndex = corelib.first_unmasked(self.data)
            if firstIndex is None:
                return TimeSeries([], dtype=self.dtype, freq=toFreq, observed=self.observed)

            startIndexAdj = self.firstValue()

            lastIndex = corelib.last_unmasked(self.data)

            tempData = copy.deepcopy(self.data[firstIndex:lastIndex+1])
            tempMask = tempData.mask
            tempData = tempData.filled()

            cRetVal = cseries.reindex(tempData, fromFreq, toFreq, position, startIndexAdj, tempMask)

            _values = cRetVal['values']
            _mask = cRetVal['mask']
            
            tempData = ma.array(_values)
            tempMask = ma.make_mask(_mask)
            tempData[tempMask] = ma.masked

            if func is not None and tempData.ndim == 2:
                tempData = corelib.apply_along_axis(func, 1, tempData)
                
            startIndex = cseries.convert(startIndexAdj, fromFreq, toFreq)

            return TimeSeries(tempData, dtype=self.data.dtype, freq=toFreq, observed=self.observed, startIndex=startIndex)
            
        else:
            return copy.deepcopy(self)


        
    def __str__(self):
        retVal = ""
        if self.firstValue() is not None:
            for i in range(self.firstValue(),self.lastValue()+1):
                index = str(tsdate.Date(freq=self.freq,value=i))
                index = index + (" " * (6-len(index)))
                retVal += index + "---> " + str(super(TimeSeries, self).__getitem__(i)) + "\n"
            return retVal
        else:
            return "<no data>"
            
            
    def firstValue(self, asDate=False):
        value = super(TimeSeries, self).firstValue()
        if asDate:
            return tsdate.Date(freq=self.freq, value=value)
        else:
            return value
        
    def lastValue(self, asDate=False):
        value = super(TimeSeries, self).lastValue()
        if asDate:
            return tsdate.Date(freq=self.freq, value=value)
        else:
            return value

    ### DATA 
    
    def __add__(self, other):
        validOpInputs(self, other)
        return SAtoTS(super(TimeSeries, self).__add__(other), self.freq, self.observed)
        
    def __radd__(self, other):
        validOpInputs(self, other)
        return SAtoTS(super(TimeSeries, self).__add__(other), self.freq, self.observed)
        
    def __sub__(self, other):
        validOpInputs(self, other)
        return SAtoTS(super(TimeSeries, self).__sub__(other), self.freq, self.observed)
        
    def __rsub__(self, other):
        validOpInputs(self, other)
        return SAtoTS(super(TimeSeries, self).__rsub__(other), self.freq, self.observed)
        
    def __mul__(self, other):
        validOpInputs(self, other)
        return SAtoTS(super(TimeSeries, self).__mul__(other), self.freq, self.observed)
        
    def __rmul__(self, other):
        validOpInputs(self, other)
        return SAtoTS(super(TimeSeries, self).__rmul__(other), self.freq, self.observed)
        
    def __div__(self, other):
        validOpInputs(self, other)
        return SAtoTS(super(TimeSeries, self).__div__(other), self.freq, self.observed)
        
    def __rdiv__(self, other):
        validOpInputs(self, other)
        return SAtoTS(super(TimeSeries, self).__rdiv__(other), self.freq, self.observed)
        
    def __pow__(self, other):
        validOpInputs(self, other)
        return SAtoTS(super(TimeSeries, self).__pow__(other), self.freq, self.observed)
        
    ### IN PLACE
    
    def __iadd__(self, other):
        validOpInputs(self, other)
        self = SAtoTS(super(TimeSeries, self).__add__(other), self.freq, self.observed)
        return self
    
    def __isub__(self, other):
        validOpInputs(self, other)
        self = SAtoTS(super(TimeSeries, self).__sub__(other), self.freq, self.observed)
        return self
    
    def __imul__(self, other):
        validOpInputs(self, other)
        self = SAtoTS(super(TimeSeries, self).__mul__(other), self.freq, self.observed)
        return self
    
    def __idiv__(self, other):
        validOpInputs(self, other)
        self = SAtoTS(super(TimeSeries, self).__div__(other), self.freq, self.observed)
        return self
        
    # this overrides & and should only be used by boolean series
    def __and__(self, other):
        validOpInputs(self, other)
        return self * other

    # this overrides | and should only be used by boolean series
    def __or__(self, other):
        validOpInputs(self, other)
        return ~(~self & ~other)
            
    # this overrides ~ and should only be used by boolean series
    # it is our "not" operator
    def __invert__(self):
        return self == False
    
    ### COMPARISON
    
    def __eq__(self, other):
        validOpInputs(self, other)
        return SAtoTS(super(TimeSeries, self).__eq__(other), self.freq, self.observed)
        
    def __le__(self, other):
        validOpInputs(self, other)
        return SAtoTS(super(TimeSeries, self).__le__(other), self.freq, self.observed)
        
    def __lt__(self, other):
        validOpInputs(self, other)
        return SAtoTS(super(TimeSeries, self).__lt__(other), self.freq, self.observed)
        
    def __ge__(self, other):
        validOpInputs(self, other)
        return SAtoTS(super(TimeSeries, self).__ge__(other), self.freq, self.observed)
        
    def __gt__(self, other):
        validOpInputs(self, other)
        return SAtoTS(super(TimeSeries, self).__gt__(other), self.freq, self.observed)

def tser(start, end):
    if start.freq != end.freq:
        raise ValueError("start and end dates must have same frequency!")
    return TimeSeries(numpy.arange(int(start), int(end)+1), dtype=corelib.freqTypeMapping[start.freq], freq=start.freq, observed='END', startIndex=int(start))

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
    newData = ma.array(cseries.getDateInfo(dateSer.data.filled(), dateSer.dtype.freq, infoCode))
    newData[dateSer.data.mask] = ma.masked
    newSer = copy.deepcopy(dateSer)
    newSer.data = newData
    newSer.dtype = numpy.int_
    return newSer

        
def validOpInputs(ser1, ser2):
    if isinstance(ser1, TimeSeries) and isinstance(ser2, TimeSeries) and ser1.freq != ser2.freq:
        raise "operation cannot be performed on series with different frequencies ("+str(ser1.freq) + " and " + str(ser2.freq)+")"

    
def SAtoTS(values, freq, observed, dtype=None):
    if dtype is None: _dtype = values.dtype
    else: _dtype = dtype
    return TimeSeries(values.data, dtype=_dtype, freq=freq, observed=observed, startIndex=values.indexZeroRepresents)


# math functions (two series)
def add(ser1, ser2, fill_value=ma.masked):
    return apply_func_twoseries(ma.add, ser1, ser2, fill_value)

def multiply(ser1, ser2, fill_value=ma.masked):
    return apply_func_twoseries(ma.multiply, ser1, ser2, fill_value)

def divide(ser1, ser2, fill_value=ma.masked):
    return apply_func_twoseries(ma.divide, ser1, ser2, fill_value)
    
def subtract(ser1, ser2, fill_value=ma.masked):
    return apply_func_twoseries(ma.subtract, ser1, ser2, fill_value)
    
# math functions (one series, return series)
def sqrt(ser):
    return apply_func_oneseries(ma.sqrt, ser)
    
# math functions (one series, return scalar)
def sum(ser):
    return ma.sum(ser.data)

def product(ser):
    return ma.product(ser.data)
    
def average(ser):
    return ma.average(ser.data)
    
def where(condition, x, y):
    tempResult = ma.where(condition.data, x, y)
    return TimeSeries(tempResult, dtype=numpy.bool_, freq=condition.freq, observed=condition.observed, startIndex=condition.indexZeroRepresents)

# generic functions
def apply_func_twoseries(func, ser1, ser2, fill_value=ma.masked):
    validOpInputs(ser1, ser2)
    return SAtoTS(doFunc(ser1, ser2, func, fill_value=fill_value), ser1.freq, ser1.observed)
    
def apply_func_oneseries(func, ser):
    return SAtoTS(doFunc_oneseries(ser, func),ser.freq, ser.observed)
    
