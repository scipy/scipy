import numpy as np
import timeseries as ts
from numpy import ma

 
# create a time series at business frequency and fill it with random data
bSer = ts.TimeSeries(np.random.uniform(-100,100,600),dtype=np.float64,freq='B',observed='SUMMED',start_date=ts.thisday('B')-600)


# Set negative values to zero.
bSer[bSer < 0] = 0


# Set values occurring on Fridays to 100.
weekdays = ts.day_of_week(ts.tser(bSer.start_date(),bSer.end_date()))
bSer[weekdays == 4] = 100


"""
Convert bSer to a monthly frequency series.

The optional func argument to the convert method specifies is a 
function that acts on a 1-dimension masked array and returns a single
value.
"""
mSer1 = bSer.convert('M',func=ts.average)


"""
If func is None, a "2 dimensional" time series will be returned. In this
example, the value for each month will be a length 23 masked array (23
being the max number of business days in a month)
"""
mSer1_2d = bSer.convert('M',func=None)


"""
If func is not specified, the observed attribute of the series
will be used to determine the method. (SUMMED for this example)
"""
mSer1_default = bSer.convert('M')


"""
Convert mSer to a business frequency series.

when converting from a lower frequency to a higher frequency, position is one
of 'START' or 'END', and determines where the data point will be placed in the
period. In the future, interpolation methods will be supported to fill in the
resulting masked values.
"""
mToB = bSer.convert('M',position='START')


# create another monthly frequency series
mSer2 = ts.TimeSeries(np.random.uniform(-100,100,100),dtype=np.float64,freq='m',observed='END',start_date=ts.thisday('M')-110)


"""
Slicing also supported. The intention is to have indexing behave
largely in the same manner as regular numpy arrays.

series.adjust_date  convert a date object into the corresponding
integer for indexing the series
"""
sixtyMonthsAgoIdx = mSer2.date_to_index(ts.thisday('m')-60)
mSer2[sixtyMonthsAgoIdx:sixtyMonthsAgoIdx+10] = 12


# Mask the last value in the series
mSer2[-1] = ts.masked #ts.masked is the same thing as numpy.ma.masked


# dates can be used as indices as well
mSer2[ts.thisday('m')-55] = 400


"""
the tser function makes it easy to index a series over a range of dates
without worrying about converting the dates to appropriate integers first
"""
mSer2[ts.tser(ts.thisday('m')-59, ts.thisday('m')-45)] = 25


"""
Only series of the same frequency and size and same start date
can be used in the basic operations.

The results are the same as you would expect for masked arrays with the
basic operations.

start_date and end_date are optional parameters to the aligned function.
If omitted, the min start_date() and end_date() of all series is used as
the new boundaries for each series.
"""
mSer1, mSer2 = ts.aligned(mSer1, mSer2, start_date=ts.thisday('m')-100, end_date=ts.thisday('m'))
mAdd1 = mSer1 + mSer2


# add the two series together, first filling in masked values with zeros
mAdd1_filled = mSer1.filled(fill_value=0, ts=True) + mSer2.filled(fill_value=0, ts=True)

# adjust the start and end dates of a series
newSer = ts.adjust_endpoints(mSer1, start_date=ts.Date(freq='M', year=1954, month=5),  end_date=ts.Date(freq='M', year=2000, month=6))

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
bLastDayOfLastQuarter = ts.dateOf(ts.thisday('Q')-1,'B','AFTER')


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
