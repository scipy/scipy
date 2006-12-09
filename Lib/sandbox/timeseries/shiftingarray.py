import numpy, types , corelib
from numpy import ma

class ShiftingArray(object):
    def __init__(self, values, dtype=None, startIndex=None, mask=ma.nomask):

        # hack to convert our fake date data types to real data types
        if corelib.isDateType(dtype):
            self.dtype = numpy.int_
        else:
            self.dtype = dtype
            
        if self.dtype is None:
            self.dtype = values.dtype

        # need to use the empty function instead of passing an empty list
        # because that won't work when type=numpy.object_
        if len(values) == 0 and dtype is numpy.object_: 
            tempData = ma.array(numpy.empty((0,), self.dtype))
        else:
            tempData = ma.array(values, self.dtype)
        
        newSize = tempData.size*2       

        firstIndex = newSize//4
        lastIndex = firstIndex + tempData.size - 1
        if startIndex is None:
            self.indexZeroRepresents = None
        else:    
            self.indexZeroRepresents = int(startIndex)-firstIndex
        
        if mask is not ma.nomask:
            tempMask = ma.make_mask(mask)
            tempData[tempMask] = ma.masked

        self.data = ma.array(numpy.empty(newSize,self.dtype))

        if firstIndex > 0:
            self.data[0:firstIndex] = ma.masked
            if self.data.size > lastIndex+1: self.data[lastIndex+1:self.data.size] = ma.masked

        self.data[firstIndex:lastIndex+1] = tempData[:]


    def shift(self, n):
        self.indexZeroRepresents += n


    #DATA ACCESS        

    def __setitem__(self, index, value):
        self.__expandToFit(self.__minIndex(index),self.__maxIndex(index))
        convIndex = self.__convIndex(index)
        self.data[convIndex] = value


    def __getitem__(self, index):
        self.__expandToFit(self.__minIndex(index),self.__maxIndex(index))
        convIndex = self.__convIndex(index)
        return self.data[convIndex]

    def _shift(self, startIndex, endIndex):
        self.__expandToFit(startIndex, endIndex)
        return self.data[startIndex-self.indexZeroRepresents:endIndex-self.indexZeroRepresents+1]

        
    #PRIVATE FUNCTIONS

    def __convIndex(self,index):

        if self.indexZeroRepresents is not None:
            if isinstance(index,ShiftingArray):

                if index.indexZeroRepresents > self.indexZeroRepresents:
                    #expand index to the left
                    originalSize = index.data.size
                    shiftAmt = index.indexZeroRepresents - self.indexZeroRepresents
                    newSize = originalSize + shiftAmt
                    temp = ma.array(numpy.empty(newSize, index.data.dtype))
                    temp[newSize-originalSize:] = index.data
                    temp[0:shiftAmt] = False
                    temp = temp.filled(False)
                else:
                    #chop off first portion of data
                    temp = index.data[self.indexZeroRepresents - index.indexZeroRepresents:].filled(False)

                # chop off extra values on right hand side
                if temp.size > self.data.size: return temp[:self.data.size]
                else: return temp

            elif type(index) == types.SliceType:
                if index.start is None: tempStart = None
                else: tempStart = index.start - self.indexZeroRepresents
                if index.stop is None: tempStop = None
                else: tempStop = index.stop - self.indexZeroRepresents
                tempStep = index.step
                
                return slice(tempStart,tempStop,tempStep)
            else:
                return index - self.indexZeroRepresents
                
        else:
            return index

    def __maxIndex(self,index):
        if type(index) == types.IntType: return index
        if type(index) == types.SliceType: return index.stop
        elif isinstance(index,ShiftingArray): return index.lastValue()
        elif hasattr(index,'__len__'): return max(index)
        else: return int(index)
    
    def __minIndex(self,index):
        if type(index) == types.IntType: return index
        if type(index) == types.SliceType: return index.start
        elif isinstance(index,ShiftingArray): return index.firstValue()
        elif hasattr(index,'__len__'): return min(index)
        else: return int(index)

    def __expandToFit(self, minRange, maxRange=None):

        if self.indexZeroRepresents is None:
            self.indexZeroRepresents = minRange

        if maxRange is None:
            maxRange = minRange
        if maxRange < minRange:
            raise ValueError("invalid range: " + str(minRange) + " to " + str(maxRange))

        minRange -= self.indexZeroRepresents
        maxRange -= self.indexZeroRepresents

        if maxRange > self.data.size-1:   #expand to the right
            originalSize = self.data.size
            newSize = originalSize
            while maxRange > newSize-1:
                newSize = expandAmt(newSize)
            
            self.data = self.data.resize(numpy.shape(numpy.empty(newSize)))
            self.data[originalSize:] = ma.masked


        if minRange < 0:                  #expand to the left
            originalSize = self.data.size
            newSize = originalSize
            shiftAmt = int(0)
            while minRange + shiftAmt < 0:
                newSize = expandAmt(newSize)
                shiftAmt = int(newSize - originalSize)

            temp = ma.array(numpy.empty(newSize, self.data.dtype))
            temp[newSize-originalSize:] = self.data
            self.data = temp
            self.data[0:shiftAmt] = ma.masked

            self.indexZeroRepresents -= shiftAmt



    #MATH FUNCTIONS

    def __add__(self, other,fill_value=ma.masked):  return doFunc(self,other, ma.add,fill_value=fill_value)
    def __radd__(self, other):                      return self+other
    def __sub__(self, other,fill_value=ma.masked):  return doFunc(self,other, ma.subtract,fill_value=fill_value)
    def __rsub__(self, other):                      return doFunc((self*-1),other, ma.add)
    def __mul__(self, other,fill_value=ma.masked):  return doFunc(self,other, ma.multiply,fill_value=fill_value)
    def __rmul__(self, other):                      return self*other
    def __div__(self, other,fill_value=ma.masked):  return doFunc(self,other, ma.divide,fill_value=fill_value)
    def __rdiv__(self, other):                      return doFunc(pow(self,-1),other, ma.multiply)
    def __pow__(self, other,fill_value=ma.masked):  return doFunc(self,other, ma.power,fill_value=fill_value)

    def __eq__(self, other):                        return doFunc(self,other, ma.equal)
    def __le__(self, other):                        return doFunc(self,other, ma.less_equal)
    def __lt__(self, other):                        return doFunc(self,other, ma.less)
    def __ge__(self, other):                        return doFunc(self,other, ma.greater_equal)
    def __gt__(self, other):                        return doFunc(self,other, ma.greater)

    def max(self,other):                            return doFunc(self,other, ma.maximum)
    def min(self,other):                            return doFunc(self,other, ma.minimum)

    #INFO FUNCTIONS

    def __len__(self):
        fv = self.firstValue()
        if fv is not None:
            return self.lastValue()-fv+1
        else:
            return 0

    def firstValue(self):
        firstIndex = first_unmasked(self.data)
        if self.indexZeroRepresents is None or firstIndex is None:
            return None
        else:
            return firstIndex+self.indexZeroRepresents


    def lastValue(self):
        lastIndex = last_unmasked(self.data)
        if self.indexZeroRepresents is None or lastIndex is None:
            return None
        else:
            return lastIndex+self.indexZeroRepresents


    #DISPLAY FUNCTIONS

    def __str__(self):
        retVal = ""
        if self.firstValue() is not None:
            for i in range(first_unmasked(self.data), last_unmasked(self.data)+1):
                index = str(i+self.indexZeroRepresents)
                index = index + (" " * (6-len(index)))
                retVal += index + "---> " + str(self.data[i]) + "\n"
            return retVal
        else:
            return "<no data>"



#apply func to ser1 and ser2, replacing masked values with fill_value
def doFunc(ser1, ser2, func,fill_value=ma.masked):
    if not isinstance(ser2, ShiftingArray):
        if ser1.indexZeroRepresents is None:
            return ShiftingArray([],ser1.data.dtype)        
        else:
            ser2 = ShiftingArray([ser2]*len(ser1),ser1.data.dtype, ser1.firstValue())
    
    sFV, sLV = ser1.firstValue(), ser1.lastValue()
    oFV, oLV = ser2.firstValue(), ser2.lastValue()

    if ser1.indexZeroRepresents is not None and ser2.indexZeroRepresents is not None:
        if fill_value is not ma.masked:
            minVal = min(sFV, oFV)
            maxVal = max(sLV, oLV)
        else:
            minVal = max(sFV, oFV)
            maxVal = min(sLV, oLV)
    elif ser1.indexZeroRepresents is None and ser2.indexZeroRepresents is None:
        return ShiftingArray([],ser1.data.dtype)
    elif ser1.indexZeroRepresents is None:
        minVal = oFV
        maxVal = oLV
    else: #ser2.indexZeroRepresents is None:
        minVal = sFV
        maxVal = sLV

    if maxVal < minVal:
        return ShiftingArray([],ser1.data.dtype, startIndex=minVal)

    data1 = ser1._shift(minVal, maxVal)
    data2 = ser2._shift(minVal, maxVal)

    if fill_value is ma.masked:
        mask = data1.mask | data2.mask
    else:
        mask = data1.mask & data2.mask

        data1 = data1.filled(fill_value)
        data2 = data2.filled(fill_value)

    data = func(data1,data2)

    return ShiftingArray(data,data.dtype,minVal, mask)
    
    
def doFunc_oneseries(ser,func):
    
    sFV = ser.firstValue()
    
    if sFV is None:
        return ser
    else:
        result = func(ser.data)
        return ShiftingArray(result,result.dtype,sFV,result.mask)


#MISC GLOBAL FUNCTIONS

def expandAmt(size):
    EXPAND_MULT = 1.2
    EXPAND_ADD  = 5
    return round(size*EXPAND_MULT) + EXPAND_ADD


def first_unmasked(m):
    idx = numpy.where(m.mask == False)
    if len(idx) != 0 and len(idx[0]) != 0:
        return idx[0][0]
    else:
        return None
    
def last_unmasked(m):
    idx = numpy.where(m.mask == False)
    if len(idx) != 0 and len(idx[0]) != 0:
        return idx[0][-1]
    else:
        return None