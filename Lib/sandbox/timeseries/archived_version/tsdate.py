import corelib
import cseries
import numpy as np
import mx.DateTime

class Date:
    def __init__(self, freq, year=None, month=None, day=None, seconds=None, quarter=None, mxDate=None, value=None):
        
        if hasattr(freq, 'freq'):
            self.freq = corelib.fmtFreq(freq.freq)
        else:
            self.freq = corelib.fmtFreq(freq)
        self.type = corelib.freqToType(self.freq)
        
        if value is not None:
            if self.freq == 'D':
                self.mxDate = mx.DateTime.DateTimeFromAbsDays(value-1)
            elif self.freq == 'B':
                value = value - 1
                self.mxDate = mx.DateTime.DateTimeFromAbsDays(value + (value//5)*7 - (value//5)*5)
            elif self.freq == 'S':
                self.mxDate = secondlyOriginDate + mx.DateTime.DateTimeDeltaFromSeconds(value)
            elif self.freq == 'M':
                self.mxDate = (mx.DateTime.Date(0)) + mx.DateTime.RelativeDateTime(months=value-1, day=-1)
            elif self.freq == 'A':
                self.mxDate = mx.DateTime.Date(value, -1, -1)
            elif self.freq == 'Q':
                self.mxDate = (mx.DateTime.Date(0)) + mx.DateTime.RelativeDateTime(years=(value // 4), month=((value * 3) % 12), day=-1)
        elif mxDate is not None:
            self.mxDate = mxDate
        else:
            error = ValueError("Insufficient parameters given to create a date at the given frequency")

            if year is None:
                raise error            
            
            if self.freq in ('B', 'D'):
                if month is None or day is None: raise error
            elif self.freq == 'M':
                if month is None: raise error
                day = -1
            elif self.freq == 'Q':
                if quarter is None: raise error
                month = quarter * 3
                day = -1
            elif self.freq == 'A':
                month = -1
                day = -1
            elif self.freq == 'S':
                if month is None or day is None or seconds is None: raise error
                
            if self.freq != 'S':
                self.mxDate = mx.DateTime.Date(year, month, day)
                if self.freq == 'B':
                    if self.mxDate.day_of_week == 5 or self.mxDate.day_of_week == 6:
                        raise ValueError("Weekend passed as business day")
            else:
                _hours = int(seconds/3600)
                _minutes = int((seconds - _hours*3600)/60)
                _seconds = seconds % 60
                
                self.mxDate = mx.DateTime.Date(year, month, day, _hours, _minutes, _seconds)
                
        self.value = self.__value()
                
    def day(self):          return self.mxDate.day
    def day_of_week(self):  return self.mxDate.day_of_week
    def month(self):        return self.mxDate.month
    def quarter(self):      return monthToQuarter(self.mxDate.month)
    def year(self):         return self.mxDate.year
    def seconds(self):      return int(self.mxDate.second)
    def minute(self):       return int(self.mxDate.minute)
    def hour(self):         return int(self.mxDate.hour)
    
    def strfmt(self, fmt):
        qFmt = fmt.replace("%q", "XXXX")
        tmpStr = self.mxDate.strftime(qFmt)
        return tmpStr.replace("XXXX", str(self.quarter()))
            
    def __str__(self):
        return self.strfmt(self.default_fmtstr())
    
    def default_fmtstr(self):
        if self.freq in ("B", "D"):
            return "%d-%b-%y"
        elif self.freq == "S":
            return "%d-%b-%Y %H:%M:%S"
        elif self.freq == "M":
            return "%b-%Y"
        elif self.freq == "Q":
            return "%Yq%q"
        elif self.freq == "A":
            return "%Y"
        else:
            return "%d-%b-%y"
        
    def __add__(self, other):
        if isinstance(other, Date):
            raise TypeError("Cannot add dates")
        return Date(freq=self.freq, value=int(self) + other)
    
    def __radd__(self, other): return self+other
    
    def __sub__(self, other):
        if isinstance(other, Date):
            if self.freq != other.freq:
                raise ValueError("Cannont subtract dates of different frequency (" + str(self.freq) + " != " + str(other.freq) + ")")
            else:
                return int(self) - int(other) 
        else:
            return self + (-1) * int(other)


    def __repr__(self): return "<" + str(self.freq) + ":" + str(self) + ">"
    
    def __eq__(self, other):
        if self.freq != other.freq:
            raise TypeError("frequencies are not equal!")
        return int(self) == int(other) 
    
    def __cmp__(self, other): 
        if self.freq != other.freq:
            raise TypeError("frequencies are not equal!")
        return int(self)-int(other)    
        
    def __hash__(self): return hash(int(self)) ^ hash(self.freq)
    
    def __int__(self):
        return self.value
            
    def __value(self):
        
        if self.freq == 'D':
            return self.mxDate.absdate
        elif self.freq == 'B':
            days = self.mxDate.absdate
            weeks = days // 7
            return int((weeks*5) + (days - weeks*7))
        elif self.freq == 'M':
            return self.mxDate.year*12 + self.mxDate.month
        elif self.freq == 'S':
            return int((self.mxDate - secondlyOriginDate).seconds)
        elif self.freq == 'A':
            return int(self.mxDate.year)
        elif self.freq == 'Q':
            return int(self.mxDate.year*4 + self.mxDate.month/3)

    
secondlyOriginDate = mx.DateTime.Date(1980) - mx.DateTime.DateTimeDeltaFromSeconds(1)

    
#######################
# FUNCTIONS
#######################
def monthToQuarter(monthNum):
    return int((monthNum-1)/3)+1

def thisday(freq):

    freq = corelib.fmtFreq(freq)

    tempDate = mx.DateTime.now()
    
    # if it is Saturday or Sunday currently, freq==B, then we want to use Friday
    if freq == 'B' and tempDate.day_of_week >= 5:
        tempDate -= (tempDate.day_of_week - 4)
    if freq == 'B' or freq == 'D' or freq == 'S':
        return Date(freq, mxDate=tempDate)
    elif freq == 'M':
        return Date(freq, year=tempDate.year, month=tempDate.month)
    elif freq == 'Q':
        return Date(freq, year=tempDate.year, quarter=monthToQuarter(tempDate.month))
    elif freq == 'A':
        return Date(freq, year=tempDate.year)


def prevbusday(day_end_hour=18, day_end_min=0):
    tempDate = mx.DateTime.localtime()

    dateNum = tempDate.hour + float(tempDate.minute)/60
    checkNum = day_end_hour + float(day_end_min)/60

    if dateNum < checkNum: return thisday('B') - 1
    else: return thisday('B')

                
#    returns date converted to a date of toFreq according to relation
#    relation = "BEFORE" or "AFTER" (not case sensitive)
def dateOf(date, toFreq, relation="BEFORE"):

    toFreq = corelib.fmtFreq(toFreq)
    _rel = relation.upper()[0]

    if date.freq == toFreq:
        return date
    else:
        return Date(freq=toFreq, value=cseries.asfreq(np.asarray(date.value), date.freq, toFreq, _rel))
