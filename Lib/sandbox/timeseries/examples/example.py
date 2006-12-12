import numpy as np
import timeseries as ts

 
# create a time series at business frequency and fill it with random data
bSer = ts.TimeSeries(np.random.uniform(-100,100,600),dtype=np.float64,freq='B',observed='SUMMED',startIndex=ts.thisday('B')-600)


# Set negative values to zero.
bSer[bSer < 0] = 0


# Set values occurring on Fridays to 100.
weekdays = ts.day_of_week(ts.tser(bSer.firstValue(asDate=True),bSer.lastValue(asDate=True)))
bSer[weekdays == 4] = 100


"""
Convert bSer to a monthly frequency series.

The optional observed argument to the convert method specifies what
method will be used to perform the frequency conversion. If it is
not specified, the observed attribute of the series will be used to
determine the method.
"""
mSer1 = bSer.convert('M',observed='AVERAGED')


# create another monthly frequency series
mSer2 = ts.TimeSeries(np.random.uniform(-100,100,100),dtype=np.float64,freq='m',observed='END',startIndex=ts.thisday('M')-110)


"""
Slicing also supported. The intention is to have indexing behave
largely in the same manner as regular numpy arrays. It sure would be
nice if we could slice with the dates directly, but as it stands we
shall have to cast the dates to integers
"""
mSer2[int(ts.thisday('m')-60):int(ts.thisday('m')-45)] = 12


# Mask a value. series.lastValue() returns the index of the last
# unmasked value in the series (as an integer, not a Date object)
mSer2[mSer2.lastValue()-40] = ts.masked #ts.masked is the same thing as numpy.ma.masked


"""
Only series of the same frequency can be used in the basic operations.
The results are the same as you would expect for masked arrays with the
basic operations.

Notice that the start and end indices of mSer1 and mSer2 do not need to
line up. This conversion is done implicitly.    
"""
mAdd1 = mSer1 + mSer2


"""
if you want more control over behaviour of masked values, use ts.add
(or multiply, etc) instead.

if a fill_value is specified, both TimeSeries objects are filled from
min(mSer1.firstValue(),mSer2.firstValue()) to max(mSer1.lastValue(),mSer2.lastValue())
wherever the series are masked before performing the operation
"""
mAdd2 = ts.add(mSer1,mSer2,fill_value=0)


# calculate the average value in the series. Behaves the same as in ma
bAverage = ts.average(bSer)


# Take the sqrt root of each element in the series (returns a TimeSeries object).
# Not all functions from ma supported yet, but they are easy to implement
# for the most part.
bSqrt = ts.sqrt(bSer)


# get the last day of this year, at daily frequency
dLastDayOfYear = ts.dateOf(ts.thisday('A'),'D','AFTER')

# get the first day of this year, at business frequency
bFirstDayOfYear = ts.dateOf(ts.thisday('A'),'B','BEFORE')

# get the last day of the previous quarter, business frequency
bFirstDayOfLastQuarter = ts.dateOf(ts.thisday('Q')-1,'B','AFTER')

# dateOf can also go from high frequency to low frequency. In this case, the third parameter has no impact
aTrueValue = (ts.thisday('Q') == ts.dateOf(ts.thisday('b'),'Q'))

# dates of the same frequency can be subtracted (but not added obviously)
numberOfBusinessDaysPassedThisYear = ts.thisday('b') - bFirstDayOfYear

# integers can be added/substracted to/from dates
fiveDaysFromNow = ts.thisday('d') + 5


# get the previous business day, where business day is considered to
# end at day_end_hour and day_end_min
pbd = ts.prevbusday(day_end_hour=18,day_end_min=0)


# construct a date object explicitly
myDateQ = ts.Date(freq='Q',year=2004,quarter=3)
myDateD = ts.Date(freq='D',year=1985,month=10,day=4)