import corelib
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
                #originDate + val + (val//5)*7 - (val//5)*5
                value -= 1
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
        if self.freq in ("B", "D"):
            return self.strfmt("%d-%b-%y")
        elif self.freq == "S":
            return self.strfmt("%d-%b-%Y %H:%M:%S")
        elif self.freq == "M":
            return self.strfmt("%b-%Y")
        elif self.freq == "Q":
            return self.strfmt("%Yq%q")
        elif self.freq == "A":
            return self.strfmt("%Y")
        else:
            return self.strfmt("%d-%b-%y")

        
    def __add__(self, other):
        if isinstance(other, Date):
            raise TypeError("Cannot add dates")
        return Date(freq=self.freq, value=int(self) + other)
    
    def __radd__(self, other): return self+other
    
    def __sub__(self, other):
        try: return self + (-1) * other
        except: pass
        try:
            if self.freq != other.freq:
                raise ValueError("Cannont subtract dates of different frequency (" + str(self.freq) + " != " + str(other.freq) + ")")
            return int(self) - int(other)
        except TypeError: 
            raise TypeError("Could not subtract types " + str(type(self)) + " and " + str(type(other)))

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
    elif date.freq == 'D':

        if toFreq == 'B':
            # BEFORE result: preceeding Friday if date is a weekend, same day otherwise
            # AFTER result: following Monday if date is a weekend, same day otherwise
            tempDate = date.mxDate
            if _rel == 'B':
                if tempDate.day_of_week >= 5: tempDate -= (tempDate.day_of_week - 4)
            elif _rel == 'A':
                if tempDate.day_of_week >= 5: tempDate += 7 - tempDate.day_of_week
            return Date(freq='B', mxDate=tempDate)
            
        elif toFreq == 'M': return Date(freq='M', year=date.year(), month=date.month())

        elif toFreq == 'S':
            if _rel == 'B': return Date(freq='S', year=date.year(), month=date.month(), day=date.day(), seconds=0)
            elif _rel == "A": return Date(freq='S', year=date.year(), month=date.month(), day=date.day(), seconds=24*60*60-1)
            
        elif toFreq == 'Q': return Date(freq='Q', year=date.year(), quarter=date.quarter())
        
        elif toFreq == 'A': return Date(freq='A', year=date.year())
        
    elif date.freq == 'B':

        if toFreq == 'D': return Date(freq='D', year=date.year(), month=date.month(), day=date.day())

        elif toFreq == 'M': return Date(freq='M', year=date.year(), month=date.month())

        elif toFreq == 'S':
            if _rel == 'B': return Date(freq='S', year=date.year(), month=date.month(), day=date.day(), seconds=0)
            elif _rel == 'A': return Date(freq='S', year=date.year(), month=date.month(), dday=ate.day(), seconds=24*60*60-1)
            
        elif toFreq == 'Q': return Date(freq='Q', year=date.year(), quarter=date.quarter())
                
        elif toFreq == 'A': return Date(freq='A', year=date.year())

    elif date.freq == 'M':

        if toFreq == 'D':
            tempDate = date.mxDate
            if _rel == 'B':
                return Date(freq='D', year=date.year(), month=date.month(), day=1)
            elif _rel == 'A':
                if date.month() == 12:
                    tempMonth = 1
                    tempYear = date.year() + 1
                else:
                    tempMonth = date.month() + 1
                    tempYear = date.year()
                return Date('D', year=tempYear, month=tempMonth, day=1)-1

        elif toFreq == 'B':
            if _rel == 'B': return dateOf(dateOf(date, 'D', "BEFORE"), 'B', "AFTER")
            elif _rel == 'A': return dateOf(dateOf(date, 'D', "AFTER"), 'B', "BEFORE")

        elif toFreq == 'S':
            if _rel == 'B': return dateOf(dateOf(date, 'D', "BEFORE"), 'S', "BEFORE")
            elif _rel == 'A': return dateOf(dateOf(date, 'D', "AFTER"), 'S', "AFTER")
            
        elif toFreq == 'Q': return Date(freq='Q', year=date.year(), quarter=date.quarter())
                        
        elif toFreq == 'A': return Date(freq='A', year=date.year())
    
    elif date.freq == 'S':

        if toFreq == 'D':
            return Date('D', year=date.year(), month=date.month(), day=date.day())
        elif toFreq == 'B':
            if _rel == 'B': return dateOf(dateOf(date, 'D'), 'B', "BEFORE")
            elif _rel == 'A': return dateOf(dateOf(date, 'D'), 'B', "AFTER")
        elif toFreq == 'M':
            return Date(freq='M', year=date.year(), month=date.month())
            
    elif date.freq == 'Q':
    
        if toFreq == 'D':
            if _rel == 'B': return dateOf(date-1, 'D', "AFTER")+1
            elif _rel == 'A': return Date(freq='D', year=date.year(), month=date.month(), day=date.day())
        elif toFreq == 'B':
            if _rel == 'B': return dateOf(dateOf(date, 'D'), 'B', "AFTER")
            if _rel == 'A': return dateOf(dateOf(date, 'D', "AFTER"), 'B', "BEFORE")
        elif toFreq == 'M':
            if _rel == 'B': return dateOf(date-1, 'M', "AFTER")+1
            elif _rel == 'A': return Date(freq='M', year=date.year(), month=date.month())
        elif toFreq == 'A': return Date(freq='A', year=date.year())
    elif date.freq == 'A':
        
        if toFreq == 'D':
            if _rel == 'B': return Date(freq='D', year=date.year(), month=1, day=1)
            elif _rel == 'A': return Date(freq='D', year=date.year(), month=12, day=31)            
        elif toFreq == 'B':
            if _rel == 'B': return dateOf(dateOf(date, 'D'), 'B', "AFTER")
            if _rel == 'A': return dateOf(dateOf(date, 'D', "AFTER"), 'B', "BEFORE")
        elif toFreq == 'M':
            if _rel == 'B': return Date(freq='M', year=date.year(), month=1)
            elif _rel == 'A': return Date(freq='M', year=date.year(), month=12)
        elif toFreq == 'Q':
            if _rel == 'B': return Date(freq='Q', year=date.year(), quarter=1)
            elif _rel == 'A': return Date(freq='Q', year=date.year(), quarter=4)