
"""
=== TimeSeries ===

A TimeSeries object is the combination of three ndarrays:
    
    - `dates`: DateArray object.
    - `data` : ndarray.
    - `mask` : Boolean ndarray, indicating missing or invalid data.
    

==== Construction ====

To construct a TimeSeries, you can use the class constructor:

>>> TimeSeries(data, dates=None, mask=nomask, 
               freq=None, observed=None, start_date=None, 
               dtype=None, copy=False, fill_value=None,
               keep_mask=True, small_mask=True, hard_mask=False)
               
where `data` is a sequence.
If `dates` is None, a DateArray of the same length as the data is constructed at 
a `freq` frequency, starting at `start_date`.

Alternatively, you can use the `time_series` function:


time_series(data, dates=None, freq=None, 
            start_date=None, end_date=None, length=None, include_last=True,
            mask=nomask, dtype=None, copy=False, fill_value=None,
            keep_mask=True, small_mask=True, hard_mask=False)    


Let us construct a series of 600 random elements, starting 600 business days ago, 
at  a business daily frequency

>>> import numpy as np
>>> import tseries as ts
>>> import tsdate as td
>>> data = np.random.uniform(-100,100,600)
>>> today = td.thisday('B')
>>> series = ts.time_series(data, dtype=np.float_, freq='B', observed='SUMMED',
                            start_date=today-600)
                            
Let us set negative values to zero...

>>> series[series<0] = 0

... and the values falling on Fridays to 100
>>> series[series.day_of_week == 4] = 100

Note that we could also create a temporary array of 'day_of weeks' for the 
corresponding period, and use it as condition.

>>> weekdays = td.day_of_week(series)
>>> series[weekdays == 4] = 100

==== Slicing ====

Accessing elements of a TimeSeries works just like slicing
>>> series[-30:]

But you can also use a date:
>>> thirtydaysago = today-30
>>> series[thirtydaysago:]

Or even a string
>>> series[thirtydaysago.tostring():]


==== Conversions ====

To convert a TimeSeries to another frequency, use the `convert` method or function.
The optional argument `func` must be a function that acts on a 1D masked array 
and returns a scalar. 

>>> mSer1 = series.convert('M',func=ma.average)

If `func` is not specified, the default value `'auto'` is used instead. In that case,
the conversion function is estimated from the `observed` attribute of the series.
For example, if `observed='SUMMED'`, then `func='auto'` is in fact `func=sum`.

>>> mSer1_default = series.convert('M')

If `func` is None, the convert method/function returns a 2D array, where each row 
corresponds to the new frequency, and the columns to the original data. In our 
example, convert will return a 2D array with 23 columns, as there are at most
23 business days per month.

>>> mSer1_2d = series.convert('M',func=None)

When converting from a lower frequency to a higher frequency, an extra argument
`position` is required. The value of the argument is either 'START' or 'END', 
and determines where the data point will be placed in the
period. In the future, interpolation methods will be supported to fill in the
resulting masked values.

Let us create a second series, this time with a monthly frequency, starting 110
months ago.
>>> data = np.random.uniform(-100,100,100).astype(np.float_)
>>> today = today.asfreq('M') - 110
>>> mSer2 = ts.TimeSeries(data, freq='m', observed='END',start_date=today)
>>> sixtymonthsago = today-60
>>> mSer2[sixtymonthsago:sixtymonthsago+10] = 12

"""
import numpy as np
import tseries as ts
reload(ts)
import tsdate as td
reload(td)
#from numpy import ma
import maskedarray as ma

data = np.random.uniform(-100,100,600)
today = td.thisday('B')
series = ts.time_series(data, dtype=np.float_, freq='B', observed='SUMMED',
                        start_date=today-600)
series[series < 0] = 0

#WARNING: The ORIGINAL day_of_week version seemed bugged. 
#It told me that 2006/12/28 was a Friday.
weekdays = td.day_of_week(series)
series[weekdays == 4] = 100

mSer1 = series.convert('M',func=ma.average)
mSer1_2d = series.convert('M',func=None)
mSer1_default = series.convert('M')
mToB = series.convert('M',position='START')


# Create another monthly frequency series
data = np.random.uniform(-100,100,100).astype(np.float_)
today = today.asfreq('M') - 110
mSer2 = ts.TimeSeries(data, freq='m', observed='END',start_date=today)
sixtymonthsago = today-60
mSer2[sixtymonthsago:sixtymonthsago+10] = 12


"""
==== Operations on TimeSeries ====

If you work with only one TimeSeries, you can use regular commands to process
the data. For example:

>>> mSer2_log = np.log(mSer2)

Note that invalid values (negative, in that case), are automatically masked.
Note as well that trying to use the maskedarray command directly is unlikely to 
work: you would end up with a regular MaskedArray.

When working with multiple series, only series of the same frequency, size and
starting date can be used in basic operations. The function `align_series` forces
series to have matching starting and ending dates. By default, the starting date
will be set to the smallest starting date of the series, and the ending date to
the largest. [An alias to `align_series` is aligned]

Let's construct a list of months, starting on Jan 2005 and ending on Dec 2006,
with a gap from Oct 2005 to Dec 2006. 

>>> mlist = ['2005-%02i' % i for i in range(1,10)]
>>> mlist += ['2006-%02i' % i for i in range(2,13)]
>>> mdata = numpy.arange(len(mlist))
>>> mser1 = time_series(mdata, mlist, observed='SUMMED')

Note that the frequency is 'U', for undefined. In fact, you have to specify what
kind of data is actually missing by forcing a given frequency.

>>> mser1 = mser1.asfreq('M')

Let us check whether there are some duplicated dates (no):

>>> mser1.has_duplicated_dates()

...or missing dates (yes):

>>> mser1.has_missing_dates()

Let us construct a second monthly series, this time without missing dates

>>> mlist2 = ['2004-%02i' % i for i in range(1,13)]
>>> mlist2 += ['2005-%02i' % i for i in range(1,13)]
>>> mser2 = time_series(mdata, mlist2, observed='SUMMED')

Let's try to add the two series:

>>> mser3 = mser1 + mser2

That doesn't work, as the series have different starting dates. We need to align 
them first.

>>> (malg1,malg2) = aligned(mser1, mser2) 

That still doesnt' work, as `malg1` has missing dates. Let us fill it, then: 
data falling on a date that was missing will be masked.

>>> mser1_filled = fill_missing_dates(mser1)
>>> (malg1,malg2) = align_series(mser1_filled, mser2) 

Now we can add the two series. Only the data that fall on dates common to the
original, non-aligned series will be actually added, the others will be masked.
After all, we are adding masked arrays.

>>> mser3 = malg1 + malg2

We could have filled the initial series first:
>>> mser3 = filled(malg1,0) + filled(malg2,0)

Alternatively, we can force the series to start/end at some given dates

>>> (malg1,malg2) = aligned(mser1_filled, mser2, 
>>>                         start_date='2004-06', end_date='2006-06')

"""
mlist = ['2005-%02i' % i for i in range(1,10)]
mlist += ['2006-%02i' % i for i in range(2,13)]
mdata = np.arange(len(mlist))
mser1 = ts.time_series(mdata, mlist, observed='SUMMED')
mser1 = mser1.asfreq('M')
#
mlist2 = ['2004-%02i' % i for i in range(1,13)]
mlist2 += ['2005-%02i' % i for i in range(1,13)]
mser2 = ts.time_series(np.arange(len(mlist2)), mlist2, observed='SUMMED')
#
mser1_filled = ts.fill_missing_dates(mser1)
(malg1,malg2) = ts.aligned(mser1_filled, mser2) 
mser3 = malg1 + malg2 

"""
# add the two series together, first filling in masked values with zeros
mAdd1_filled = mSer1.filled(fill_value=0, ts=True) + mSer2.filled(fill_value=0, ts=True)

# adjust the start and end dates of a series
newSer = ts.adjust_endpoints(mSer1, start_date=td.Date(freq='M', year=1954, month=5),  end_date=td.Date(freq='M', year=2000, month=6))

# calculate the average value in the series. Behaves the same as in ma
bAverage = ma.average(series)





# Get the last day of this year, at daily frequency
dLastDayOfYear = td.dateOf(td.thisday('A'),'D','AFTER')


# Get the first day of this year, at business frequency
bFirstDayOfYear = td.dateOf(td.thisday('A'),'B','BEFORE')

# Get the last day of the previous quarter, business frequency
bLastDayOfLastQuarter = td.dateOf(td.thisday('Q')-1,'B','AFTER')

# dateOf can also go from high frequency to low frequency. In this case, the third parameter has no impact
aTrueValue = (td.thisday('Q') == td.dateOf(td.thisday('b'),'Q'))

# Dates of the same frequency can be subtracted (but not added obviously)
numberOfBusinessDaysPassedThisYear = td.thisday('b') - bFirstDayOfYear

# Integers can be added/substracted to/from dates
fiveDaysFromNow = td.thisday('d') + 5


# Get the previous business day, where business day is considered to
# end at day_end_hour and day_end_min
pbd = td.prevbusday(day_end_hour=18,day_end_min=0)
"""
# ------------------------------------------------------------------------------
"""
=== Date construction ===
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

# Construct a sequence of dates at a given frequency
