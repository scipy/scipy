import corelib
import mx.DateTime
import numpy

class Date:
    def __init__(self,freq,year=None, month=None, day=None, seconds=None,quarter=None, date=None, val=None):
        
        if hasattr(freq,'freq'):
            self.freq = corelib.fmtFreq(freq.freq)
        else:
            self.freq = corelib.fmtFreq(freq)
        self.type = corelib.freqToType(self.freq)
        
        if val is not None:
            if self.freq == 'D':
                self.__date = val+originDate
            elif self.freq == 'B':
                self.__date = originDate + val + (val//5)*7 - (val//5)*5
            elif self.freq == 'S':
                self.__date = secondlyOriginDate + mx.DateTime.DateTimeDeltaFromSeconds(val)
            elif self.freq == 'M':
                self.__date = originDate + mx.DateTime.RelativeDateTime(months=val, day=-1)
            elif self.freq == 'A':
                self.__date = originDate + mx.DateTime.RelativeDateTime(years=val, month=-1, day=-1)
            elif self.freq == 'Q':
                self.__date = originDate + 1 + mx.DateTime.RelativeDateTime(years=int(val/4), month=int(12 * (float(val)/4 - val/4)), day=-1)
        elif date is not None:
            self.__date = date
        else:
            error = ValueError("Insufficient parameters given to create a date at the given frequency")

            if year is None:
                raise error            
            
            if self.freq in ('B','D'):
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
                self.__date = mx.DateTime.Date(year, month, day)
                if self.freq == 'B':
                    if self.__date.day_of_week == 5 or self.__date.day_of_week == 6:
                        raise ValueError("Weekend passed as business day")
            else:
                _hours = int(seconds/3600)
                _minutes = int((seconds - _hours*3600)/60)
                _seconds = seconds % 60
                
                self.__date = mx.DateTime.Date(year, month, day, _hours, _minutes, _seconds)
                        
                
    def day(self):          return self.getDate().day
    def day_of_week(self):  return self.getDate().day_of_week
    def month(self):        return self.getDate().month
    def quarter(self):      return monthToQuarter(self.getDate().month)
    def year(self):         return self.getDate().year
    def seconds(self):      return int(self.getDate().second)
    def minute(self):       return int(self.getDate().minute)
    def hour(self):         return int(self.getDate().hour)
    
    def strfmt(self,fmt):
        qFmt = fmt.replace("%q","XXXX")
        tmpStr = self.__date.strftime(qFmt)
        return tmpStr.replace("XXXX",str(self.quarter()))
            
    def __str__(self):
        if self.freq in ("B","D"):
            return self.__date.strftime("%d-%b-%y")
        elif self.freq == "S":
            return self.__date.strftime("%d-%b-%Y %H:%M:%S")
        elif self.freq == "M":
            return self.__date.strftime("%b-%Y")
        elif self.freq == "Q":
            return str(self.year())+"q"+str(self.quarter())
        elif self.freq == "A":
            return str(self.year())
        else:
            return self.__date.strftime("%d-%b-%y")

        
    def __add__(self, other):
        if isinstance(other, Date):
            raise TypeError("Cannot add dates")
        return Date(freq=self.freq, val=int(self) + other)
    
    def __radd__(self, other): return self+other
    
    def __sub__(self, other):
        try: return self + (-1) * other
        except: pass
        try:
            if self.freq <> other.freq:
                raise ValueError("Cannont subtract dates of different frequency (" + str(self.freq) + " <> " + str(other.freq) + ")")
            return int(self) - int(other)
        except TypeError: 
            raise TypeError("Could not subtract types " + str(type(self)) + " and " + str(type(other)))

    def __repr__(self): return "<" + str(self.freq) + ":" + str(self) + ">"
    
    def __eq__(self, other):
        if self.freq <> other.freq:
            raise TypeError("frequencies are not equal!")
        return int(self) == int(other) 
    
    def __cmp__(self, other): 
        if self.freq <> other.freq:
            raise TypeError("frequencies are not equal!")
        return int(self)-int(other)    
        
    def __hash__(self): return hash(int(self)) ^ hash(self.freq)
    
    def __int__(self):
        return self.value()
            
    def value(self):
        if self.freq == 'D':
            return int((self.__date-originDate).days)
        elif self.freq == 'B':
            days = (self.__date-originDate).days
            weeks = days // 7
            return int((weeks*5) + (days - weeks*7))
        elif self.freq == 'M':
            return (self.__date.year - originDate.year)*12 + (self.__date.month - originDate.month)
        elif self.freq == 'S':
            return int((self.__date - secondlyOriginDate).seconds)
        elif self.freq == 'A':
            return int(self.__date.year - originDate.year + 1)
        elif self.freq == 'Q':
            return int ((self.__date.year - originDate.year)*4 + (self.__date.month - originDate.month)/3)
            
    def mxDate(self):
        return self.__date
    
originDate = mx.DateTime.Date(1850)-1
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
        return Date(freq, date=tempDate)
    elif freq == 'M':
        return Date(freq,tempDate.year,tempDate.month)
    elif freq == 'Q':
        return Date(freq,tempDate.year,quarter=monthToQuarter(tempDate.month))
    elif freq == 'A':
        return Date(freq,tempDate.year)

def prevbusday(day_end_hour=18,day_end_min=0):
    tempDate = mx.DateTime.localtime()

    dateNum = tempDate.hour + float(tempDate.minute)/60
    checkNum = day_end_hour + float(day_end_min)/60

    if dateNum < checkNum: return thisday('B') - 1
    else: return thisday('B')

                
#    returns _date converted to a date of _destFreq according to _relation
#    _relation = "BEFORE" or "AFTER" (not case sensitive)
def dateOf(_date,_destFreq,_relation="BEFORE"):

    _destFreq = corelib.fmtFreq(_destFreq)
    _rel = _relation.upper()[0]

    if _date.freq == _destFreq:
        return _date
    elif _date.freq == 'D':

        if _destFreq == 'B':
            # BEFORE result: preceeding Friday if _date is a weekend, same day otherwise
            # AFTER result: following Monday if _date is a weekend, same day otherwise
            tempDate = _date.mxDate()
            if _rel == "B":
                if tempDate.day_of_week >= 5: tempDate -= (tempDate.day_of_week - 4)
            elif _rel == "A":
                if tempDate.day_of_week >= 5: tempDate += 7 - tempDate.day_of_week
            return Date(freq='B',date=tempDate)
            
        elif _destFreq == 'M': return Date('M',_date.mxDate().year,_date.mxDate().month)

        elif _destFreq == 'S':
            if _rel == "B": return Date('S',_date.mxDate().year,_date.mxDate().month,_date.mxDate().day,0)
            elif _rel == "A": return Date('S',_date.mxDate().year,_date.mxDate().month,_date.mxDate().day,24*60*60-1)
            
        elif _destFreq == 'Q': return Date('Q',_date.mxDate().year,quarter=monthToQuarter(_date.mxDate().month))
        
        elif _destFreq == 'A': return Date('A',_date.mxDate().year)
        
    elif _date.freq == 'B':

        if _destFreq == 'D': return Date('D',_date.mxDate().year,_date.mxDate().month,_date.mxDate().day)

        elif _destFreq == 'M': return Date('M',_date.mxDate().year,_date.mxDate().month)

        elif _destFreq == 'S':
            if _rel == "B": return Date('S',_date.mxDate().year,_date.mxDate().month,_date.mxDate().day,0)
            elif _rel == "A": return Date('S',_date.mxDate().year,_date.mxDate().month,_date.mxDate().day,24*60*60-1)
            
        elif _destFreq == 'Q': return Date('Q',_date.mxDate().year,quarter=monthToQuarter(_date.mxDate().month))
                
        elif _destFreq == 'A': return Date('A',_date.mxDate().year)

    elif _date.freq == 'M':

        if _destFreq == 'D':
            tempDate = _date.mxDate()
            if _rel == "B":
                return Date('D',_date.mxDate().year,_date.mxDate().month,1)
            elif _rel == "A":
                if _date.mxDate().month == 12:
                    tempMonth = 1
                    tempYear = _date.mxDate().year + 1
                else:
                    tempMonth = _date.mxDate().month + 1
                    tempYear = _date.mxDate().year
                return Date('D',tempYear,tempMonth,1)-1

        elif _destFreq == 'B':
            if _rel == "B": return dateOf(dateOf(_date,'D',"BEFORE"),'B',"AFTER")
            elif _rel == "A": return dateOf(dateOf(_date,'D',"AFTER"),'B',"BEFORE")

        elif _destFreq == 'S':
            if _rel == "B": return dateOf(dateOf(_date,'D',"BEFORE"),'S',"BEFORE")
            elif _rel == "A": return dateOf(dateOf(_date,'D',"AFTER"),'S',"AFTER")
            
        elif _destFreq == 'Q': return Date('Q',_date.mxDate().year,quarter=monthToQuarter(_date.mxDate().month))
                        
        elif _destFreq == 'A': return Date('A',_date.mxDate().year)
    
    elif _date.freq == 'S':

        if _destFreq == 'D':
            return Date('D',_date.mxDate().year,_date.mxDate().month,_date.mxDate().day)
        elif _destFreq == 'B':
            if _rel == "B": return dateOf(dateOf(_date,'D'),'B',"BEFORE")
            elif _rel == "A": return dateOf(dateOf(_date,'D'),'B',"AFTER")
        elif _destFreq == 'M':
            return Date('M',_date.mxDate().year,_date.mxDate().month)
            
    elif _date.freq == 'Q':
    
        if _destFreq == 'D':
            if _rel == "B": return dateOf(_date-1,'D',"AFTER")+1
            elif _rel == "A": return Date('D',_date.mxDate().year,_date.mxDate().month,_date.mxDate().day)
        elif _destFreq == 'B':
            if _rel == "B": return dateOf(dateOf(_date,'D'),'B',"AFTER")
            if _rel == "A": return dateOf(dateOf(_date,'D',"AFTER"),'B',"BEFORE")
        elif _destFreq == 'M':
            if _rel == "B": return dateOf(_date-1,'M',"AFTER")+1
            elif _rel == "A": return Date('M',_date.mxDate().year,_date.mxDate().month)
        elif _destFreq == 'A': return Date('A',_date.mxDate().year)
    elif _date.freq == 'A':
        
        if _destFreq == 'D':
            if _rel == "B": return Date('D',_date.mxDate().year, 1, 1)
            elif _rel == "A": return Date('D',_date.mxDate().year,12,31)            
        elif _destFreq == 'B':
            if _rel == "B": return dateOf(dateOf(_date,'D'),'B',"AFTER")
            if _rel == "A": return dateOf(dateOf(_date,'D',"AFTER"),'B',"BEFORE")
        elif _destFreq == 'M':
            if _rel == "B": return Date('M',_date.mxDate().year,1)
            elif _rel == "A": return Date('M',_date.mxDate().year,12)
        elif _destFreq == 'Q':
            if _rel == "B": return Date('Q',_date.mxDate().year,quarter=1)
            elif _rel == "A": return Date('Q',_date.mxDate().year,quarter=4)