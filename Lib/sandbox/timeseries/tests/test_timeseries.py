# pylint: disable-msg=W0611, W0612, W0511,R0201
"""Tests suite for MaskedArray.
Adapted from the original test_ma by Pierre Gerard-Marchant

:author: Pierre Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu & mattknox_ca_at_hotmail_dot_com
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import numpy as N
from numpy import bool_, complex_, float_, int_, object_
import numpy.core.fromnumeric  as fromnumeric
import numpy.core.numeric as numeric
from numpy.testing import NumpyTest, NumpyTestCase
from numpy.testing.utils import build_err_msg

import maskedarray
import maskedarray as MA
from maskedarray import masked_array, masked, nomask

import maskedarray.testutils
from maskedarray.testutils import assert_equal, assert_array_equal

from timeseries import tseries
from timeseries import Date, date_array_fromlist, date_array, thisday
from timeseries import time_series, TimeSeries, adjust_endpoints, \
    mask_period, align_series, fill_missing_dates, tsmasked, concatenate_series,\
    stack, split

class test_creation(NumpyTestCase):
    "Base test class for MaskedArrays."
    def __init__(self, *args, **kwds):
        NumpyTestCase.__init__(self, *args, **kwds)
        dlist = ['2007-01-%02i' % i for i in range(1,16)]
        dates = date_array_fromlist(dlist)
        data = masked_array(numeric.arange(15), mask=[1,0,0,0,0]*3)
        self.d = (dlist, dates, data)

    def test_fromlist (self):
        "Base data definition."
        (dlist, dates, data) = self.d
        series = time_series(data, dlist)
        assert(isinstance(series, TimeSeries))
        assert_equal(series._mask, [1,0,0,0,0]*3)
        assert_equal(series._series, data)
        assert_equal(series._dates, date_array_fromlist(dlist))
        assert_equal(series.freqstr, 'D')

    def test_fromrange (self):
        "Base data definition."
        (dlist, dates, data) = self.d
        series = time_series(data, start_date=dates[0], length=15)
        assert(isinstance(series, TimeSeries))
        assert_equal(series._mask, [1,0,0,0,0]*3)
        assert_equal(series._series, data)
        assert_equal(series._dates, dates)
        assert_equal(series.freqstr, 'D')

    def test_fromseries (self):
        "Base data definition."
        (dlist, dates, data) = self.d
        series = time_series(data, dlist)
        dates = dates+15
        series = time_series(series, dates)
        assert(isinstance(series, TimeSeries))
        assert_equal(series._mask, [1,0,0,0,0]*3)
        assert_equal(series._series, data)
        assert_equal(series._dates, dates)
        assert_equal(series.freqstr, 'D')


    def test_fromdatearray(self):
        "Tests the creation of a series from a datearray"
        _, dates, _ = self.d
        data = dates
        #
        series = time_series(data, dates)
        assert(isinstance(series, TimeSeries))
        assert_equal(series._dates, dates)
        assert_equal(series._data, data)
        assert_equal(series.freqstr, 'D')
        #
        series[5] = MA.masked
        # ensure that series can be represented by a string after masking a value
        # (there was a bug before that prevented this from working when using a
        # DateArray for the data)
        strrep = str(series)


    def test_datafromlist(self):
        "Check the creation of a time series from a list of data."
        (_, dates, _) = self.d
        data = list(range(15))
        series = time_series(data, dates)
        assert_equal(series._data.size, 15)
        
    def test_unsorted(self):
        "Tests that the data are properly sorted along the dates."
        dlist = ['2007-01-%02i' % i for i in (3,2,1)]
        data = [10,20,30]
        series = time_series(data,dlist)
        assert_equal(series._data,[30,20,10])
        #
        series = TimeSeries(data, dlist)
        assert_equal(series._data,[30,20,10])
        #
        series = TimeSeries(data, dlist, mask=[1,0,0])
        assert_equal(series._mask,[0,0,1])
        #
        data = masked_array([10,20,30],mask=[1,0,0])
        series = TimeSeries(data, dlist)
        assert_equal(series._mask,[0,0,1])
#...............................................................................

class test_arithmetics(NumpyTestCase):
    "Some basic arithmetic tests"
    def __init__(self, *args, **kwds):
        NumpyTestCase.__init__(self, *args, **kwds)
        dlist = ['2007-01-%02i' % i for i in range(1,16)]
        dates = date_array_fromlist(dlist)
        data = masked_array(numeric.arange(15), mask=[1,0,0,0,0]*3)
        self.d = (time_series(data, dlist), data)

    def test_intfloat(self):
        "Test arithmetic timeseries/integers"
        (series, data) =self.d
        #
        nseries = series+1
        assert(isinstance(nseries, TimeSeries))
        assert_equal(nseries._mask, [1,0,0,0,0]*3)
        assert_equal(nseries._series, data+1)
        assert_equal(nseries._dates, series._dates)
        #
        nseries = series-1
        assert(isinstance(nseries, TimeSeries))
        assert_equal(nseries._mask, [1,0,0,0,0]*3)
        assert_equal(nseries._series, data-1)
        assert_equal(nseries._dates, series._dates)
        #
        nseries = series*1
        assert(isinstance(nseries, TimeSeries))
        assert_equal(nseries._mask, [1,0,0,0,0]*3)
        assert_equal(nseries._series, data*1)
        assert_equal(nseries._dates, series._dates)
        #
        nseries = series/1.
        assert(isinstance(nseries, TimeSeries))
        assert_equal(nseries._mask, [1,0,0,0,0]*3)
        assert_equal(nseries._series, data/1.)
        assert_equal(nseries._dates, series._dates)

    def test_intfloat_inplace(self):
        "Test int/float arithmetics in place."
        (series, data) =self.d
        nseries = series.astype(float_)
        idini = id(nseries)
        data = data.astype(float_)
        #
        nseries += 1.
        assert(isinstance(nseries, TimeSeries))
        assert_equal(nseries._mask, [1,0,0,0,0]*3)
        assert_equal(nseries._series, data+1.)
        assert_equal(nseries._dates, series._dates)
        assert_equal(id(nseries),idini)
        #
        nseries -= 1.
        assert(isinstance(nseries, TimeSeries))
        assert_equal(nseries._mask, [1,0,0,0,0]*3)
        assert_equal(nseries._series, data)
        assert_equal(nseries._dates, series._dates)
        assert_equal(id(nseries),idini)
        #
        nseries *= 2.
        assert(isinstance(nseries, TimeSeries))
        assert_equal(nseries._mask, [1,0,0,0,0]*3)
        assert_equal(nseries._series, data*2.)
        assert_equal(nseries._dates, series._dates)
        assert_equal(id(nseries),idini)
        #
        nseries /= 2.
        assert(isinstance(nseries, TimeSeries))
        assert_equal(nseries._mask, [1,0,0,0,0]*3)
        assert_equal(nseries._series, data)
        assert_equal(nseries._dates, series._dates)
        assert_equal(id(nseries),idini)
    #
    def test_updatemask(self):
        "Checks modification of mask."
        (series, data) =self.d
        assert_equal(series._mask, [1,0,0,0,0]*3)
        series.mask = nomask
        assert(series._mask is nomask)
        assert(series._series._mask is nomask)
        #series._series.mask = [1,0,0]*5
        series.mask = [1,0,0]*5
        assert_equal(series._mask, [1,0,0]*5)
        assert_equal(series._series._mask, [1,0,0]*5)
        series[2] = masked
        assert_equal(series._mask, [1,0,1]+[1,0,0]*4)
        assert_equal(series._series._mask, [1,0,1]+[1,0,0]*4)
    #
    def test_ismasked(self):
        "Checks checks on masked"
        (series, data) =self.d
        assert(series[0] is tsmasked)
        assert(tsmasked._series is masked)
        assert(series._series[0] is masked)
        assert(series[0]._series is masked)


#...............................................................................

class test_getitem(NumpyTestCase):
    "Some getitem tests"
    def __init__(self, *args, **kwds):
        NumpyTestCase.__init__(self, *args, **kwds)
        dlist = ['2007-01-%02i' % i for i in range(1,16)]
        dates = date_array_fromlist(dlist)
        data = masked_array(numeric.arange(15), mask=[1,0,0,0,0]*3, dtype=float_)
        self.d = (time_series(data, dlist), data, dates)

    def test_wdate(self):
        "Tests  getitem with date as index"
        (series, data, dates) = self.d
        last_date = series._dates[-1]
        assert_equal(series[-1], series[last_date])
        assert_equal(series._dates[-1], dates[-1])
        assert_equal(series[-1]._dates[0], dates[-1])
        assert_equal(series[last_date]._dates[0], dates[-1])
        assert_equal(series._series[-1], data._data[-1])
        assert_equal(series[-1]._series, data._data[-1])
        assert_equal(series._mask[-1], data._mask[-1])
        #
        series['2007-01-06'] = 999
        assert_equal(series[5], 999)
        #
    def test_wtimeseries(self):
        "Tests getitem w/ TimeSeries as index"
        (series, data, dates) = self.d
        # Testing a basic condition on data
        cond = (series<8).filled(False)
        dseries = series[cond]
        assert_equal(dseries._data, [1,2,3,4,6,7])
        assert_equal(dseries._dates, series._dates[[1,2,3,4,6,7]])
        assert_equal(dseries._mask, nomask)
        # Testing a basic condition on dates
        series[series._dates < Date('D',string='2007-01-06')] = masked
        assert_equal(series[:5]._series._mask, [1,1,1,1,1])

    def test_wslices(self):
        "Test get/set items."
        (series, data, dates) = self.d
        # Basic slices
        assert_equal(series[3:7]._series._data, data[3:7]._data)
        assert_equal(series[3:7]._series._mask, data[3:7]._mask)
        assert_equal(series[3:7]._dates, dates[3:7])
        # Ditto
        assert_equal(series[:5]._series._data, data[:5]._data)
        assert_equal(series[:5]._series._mask, data[:5]._mask)
        assert_equal(series[:5]._dates, dates[:5])
        # With set
        series[:5] = 0
        assert_equal(series[:5]._series, [0,0,0,0,0])
        dseries = N.log(series)
        series[-5:] = dseries[-5:]
        assert_equal(series[-5:], dseries[-5:])
        # Now, using dates !
        dseries = series[series.dates[3]:series.dates[7]]
        assert_equal(dseries, series[3:7])

    def test_on2d(self):
        "Tests getitem on a 2D series"
        (a,b,d) = ([1,2,3],[3,2,1], date_array(thisday('M'),length=3))
        ser_x = time_series(N.column_stack((a,b)), dates=d)
        assert_equal(ser_x[0,0], time_series(a[0],d[0]))
        assert_equal(ser_x[0,:], time_series([(a[0],b[0])], d[0]))
        assert_equal(ser_x[:,0], time_series(a, d))
        assert_equal(ser_x[:,:], ser_x)

    def test_onnd(self):
        "Tests getitem on a nD series"
        hodie = thisday('D')
        # Case 1D
        series = time_series(N.arange(5), mask=[1,0,0,0,0], start_date=hodie)
        assert_equal(series[0], 0)
        # Case 1D + mask
        series = time_series(N.arange(5), mask=[1,0,0,0,0], start_date=hodie)
        assert series[0] is tsmasked
        # Case 2D
        series = time_series(N.arange(10).reshape(5,2), start_date=hodie)
        assert_equal(len(series), 5)
        assert_equal(series[0], [[0,1]])
        assert_equal(series[0]._dates[0], (hodie))
        assert_equal(series[:,0], [0,2,4,6,8])
        assert_equal(series[:,0]._dates, series._dates)
        # Case 2D + mask
        series = time_series(N.arange(10).reshape(5,2), start_date=hodie,
                             mask=[[1,1],[0,0],[0,0],[0,0],[0,0]])
        assert_equal(len(series), 5)
        assert_equal(series[0], [[0,1]])
        assert_equal(series[0]._mask, [[1,1]])
        assert_equal(series[0]._dates[0], (hodie))
        assert_equal(series[:,0]._data, [0,2,4,6,8])
        assert_equal(series[:,0]._mask, [1,0,0,0,0])
        assert_equal(series[:,0]._dates, series._dates)
        # Case 3D
        series = time_series(N.arange(30).reshape(5,3,2), start_date=hodie)
        x = series[0]
        assert_equal(len(series), 5)
        assert_equal(series[0], [[[0,1],[2,3],[4,5]]])
        assert_equal(series[0]._dates[0], (hodie))
        assert_equal(series[:,0], series._data[:,0])
        assert_equal(series[:,0]._dates, series._dates)
        x = series[:,:,0]
        assert_equal(series[:,:,0], series._data[:,:,0])
        assert_equal(series[:,:,0]._dates, series._dates)

class test_functions(NumpyTestCase):
    "Some getitem tests"
    def __init__(self, *args, **kwds):
        NumpyTestCase.__init__(self, *args, **kwds)
        dlist = ['2007-01-%02i' % i for i in range(1,16)]
        dates = date_array_fromlist(dlist)
        data = masked_array(numeric.arange(15), mask=[1,0,0,0,0]*3)
        self.d = (time_series(data, dlist), data, dates)
    #
    def test_adjustendpoints(self):
        "Tests adjust_endpoints"
        (series, data, dates) = self.d
        dseries = adjust_endpoints(series, series.dates[0], series.dates[-1])
        assert_equal(dseries, series)
        dseries = adjust_endpoints(series, series.dates[3], series.dates[-3])
        assert_equal(dseries, series[3:-2])
        dseries = adjust_endpoints(series, end_date=Date('D', string='2007-01-31'))
        assert_equal(dseries.size, 31)
        assert_equal(dseries._mask, N.r_[series._mask, [1]*16])
        dseries = adjust_endpoints(series, end_date=Date('D', string='2007-01-06'))
        assert_equal(dseries.size, 6)
        assert_equal(dseries, series[:6])
        dseries = adjust_endpoints(series,
                                   start_date=Date('D', string='2007-01-06'),
                                   end_date=Date('D', string='2007-01-31'))
        assert_equal(dseries.size, 26)
        assert_equal(dseries._mask, N.r_[series._mask[5:], [1]*16])
    #
    def test_alignseries(self):
        "Tests align_series & align_with"
        (series, data, dates) = self.d
        #
        empty_series = time_series([], freq='d')
        a, b = align_series(series, empty_series)
        assert_equal(a.start_date, b.start_date)
        assert_equal(a.end_date, b.end_date)
        #
        aseries = time_series(data, dates+10)
        bseries = time_series(data, dates-10)
        (a, b) = align_with(series, aseries, bseries)
        assert_equal(a._dates, series._dates)
        assert_equal(b._dates, series._dates)
        assert_equal(a[-5:], series[:5])
        assert_equal(b[:5], series[-5:])
    #
    def test_tshift(self):
        "Test tshift function"
        series = self.d[0]
        shift_negative = series.tshift(-1)
        result_data = [999] + [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
        result_mask = [1  ] + [1,0,0,0,0,1,0,0,0,0,1, 0, 0, 0 ]
        shift_negative_result = time_series(result_data, series._dates, mask=result_mask)

        shift_positive = series.tshift(1)
        result_data = [1,2,3,4,5,6,7,8,9,10,11,12,13,14] + [999]
        result_mask = [0,0,0,0,1,0,0,0,0,1, 0, 0, 0, 0 ] + [1  ]
        shift_positive_result = time_series(result_data, series._dates, mask=result_mask)

        assert_array_equal(shift_negative, shift_negative_result)
        assert_array_equal(shift_positive, shift_positive_result)
    #
    def test_split(self):
        """Test the split function."""
        ms = time_series(N.arange(62).reshape(31,2),
                         start_date=Date(freq='d', year=2005, month=7, day=1))
        d1,d2 = split(ms)
        assert_array_equal(d1.data, ms.data[:,0])
        assert_array_equal(d1.dates, ms.dates)
        assert_array_equal(d2.data, ms.data[:,1])

        series = self.d[0]
        ss = split(series)[0]
        assert_array_equal(series, ss)
    #
    def test_convert(self):
        """Test convert function

Just check basic functionality. The details of the actual
date conversion algorithms already tested by asfreq in the
test_dates test suite.
        """
        lowFreqSeries = time_series(N.arange(10),
                                    start_date=Date(freq='m', year=2005, month=6))
        highFreqSeries = time_series(N.arange(100),
                                    start_date=Date(freq='b', year=2005, month=6, day=1))
        ndseries = time_series(N.arange(124).reshape(62,2), 
                             start_date=Date(freq='d', year=2005, month=7, day=1))

        lowToHigh_start = lowFreqSeries.convert('B', position='START')

        assert_equal(lowToHigh_start.start_date,
                     Date(freq='m', year=2005, month=6).asfreq("B", relation="BEFORE"))
        assert_equal(lowToHigh_start.end_date,
                     (Date(freq='m', year=2005, month=6) + 9).asfreq("B", relation="AFTER"))

        assert_equal(lowToHigh_start._mask[0], False)
        assert_equal(lowToHigh_start._mask[-1], True)

        lowToHigh_end = lowFreqSeries.convert('B', position='END')

        assert_equal(lowToHigh_end.start_date,
                     Date(freq='m', year=2005, month=6).asfreq("B", relation="BEFORE"))
        assert_equal(lowToHigh_end.end_date,
                     (Date(freq='m', year=2005, month=6) + 9).asfreq("B", relation="AFTER"))

        assert_equal(lowToHigh_end._mask[0], True)
        assert_equal(lowToHigh_end._mask[-1], False)


        highToLow = highFreqSeries.convert('M', func=None)

        assert_equal(highToLow.ndim, 2)
        assert_equal(highToLow.shape[1], 23)
        assert_equal(highToLow.start_date,
                     Date(freq='b', year=2005, month=6, day=1).asfreq('M'))
        assert_equal(highToLow.end_date,
                     (Date(freq='b', year=2005, month=6, day=1) + 99).asfreq('M'))

        assert_array_equal(lowFreqSeries, lowFreqSeries.convert("M"))
                
        assert_equal(ndseries.convert('M',sum), [[930,961],[2852,2883]])
    #
    def test_fill_missing_dates(self):
        """Test fill_missing_dates function"""
        _start = Date(freq='m', year=2005, month=1)
        _end = Date(freq='m', year=2005, month=4)

        dates = date_array([_start, _end], freq='M')
        series = time_series([1, 2], dates)
        filled_ser = fill_missing_dates(series)

        assert_equal(filled_ser.start_date, _start)
        assert_equal(filled_ser.end_date, _end)
        assert(filled_ser.isfull())
        assert(not filled_ser.has_duplicated_dates())
        assert_equal(filled_ser.size, _end - _start + 1)
    #
    def test_maskperiod(self):
        "Test mask_period"
        (series, data, dates) = self.d
        series.mask = nomask
        (start, end) = ('2007-01-06', '2007-01-12')
        mask = mask_period(series, start, end, inside=True, include_edges=True,
                           inplace=False)
        assert_equal(mask._mask, N.array([0,0,0,0,0,1,1,1,1,1,1,1,0,0,0]))
        mask = mask_period(series, start, end, inside=True, include_edges=False,
                           inplace=False)
        assert_equal(mask._mask, [0,0,0,0,0,0,1,1,1,1,1,0,0,0,0])
        mask = mask_period(series, start, end, inside=False, include_edges=True,
                           inplace=False)
        assert_equal(mask._mask, [1,1,1,1,1,1,0,0,0,0,0,1,1,1,1])
        mask = mask_period(series, start, end, inside=False, include_edges=False,
                           inplace=False)
        assert_equal(mask._mask, [1,1,1,1,1,0,0,0,0,0,0,0,1,1,1])
        # Now w/ multivariables
        data = masked_array(numeric.arange(30).reshape(-1,2), dtype=float_)
        series = time_series(data, dates=dates)
        mask = mask_period(series, start, end, inside=True, include_edges=True,
                           inplace=False)
        result = N.array([0,0,0,0,0,1,1,1,1,1,1,1,0,0,0])
        assert_equal(mask._mask, result.repeat(2).reshape(-1,2))
    #
    def test_pickling(self):
        "Tests pickling/unpickling"
        (series, data, dates) = self.d
        import cPickle
        series_pickled = cPickle.loads(series.dumps())
        assert_equal(series_pickled._dates, series._dates)
        assert_equal(series_pickled._data, series._data)
        assert_equal(series_pickled._mask, series._mask)
        #
        data = masked_array(N.matrix(range(10)).T, mask=[1,0,0,0,0]*2)
        dates = date_array(start_date=thisday('D'), length=10)
        series = time_series(data,dates=dates)
        series_pickled = cPickle.loads(series.dumps())
        assert_equal(series_pickled._dates, series._dates)
        assert_equal(series_pickled._data, series._data)
        assert_equal(series_pickled._mask, series._mask)
        assert(isinstance(series_pickled._data, N.matrix))


    def test_empty_timeseries(self):
        "Tests that empty TimeSeries are  handled properly"
        empty_ts = time_series([], freq='b')
        assert_array_equal(empty_ts, empty_ts + 1)
        assert_array_equal(empty_ts, empty_ts + empty_ts)
        assert_equal(empty_ts.start_date, None)
        assert_equal(empty_ts.end_date, None)

    def test__timeseriescompat_multiple(self):
        "Tests the compatibility of multiple time series."
        seriesM_10 = time_series(N.arange(10),
                                    date_array(
                                      start_date=Date(freq='m', year=2005, month=1),
                                      length=10)
                                )

        seriesD_10 = time_series(N.arange(10),
                                    date_array(
                                      start_date=Date(freq='d', year=2005, month=1, day=1),
                                      length=10)
                                )

        seriesD_5 = time_series(N.arange(5),
                                    date_array(
                                      start_date=Date(freq='d', year=2005, month=1, day=1),
                                      length=5)
                                )

        seriesD_5_apr = time_series(N.arange(5),
                                    date_array(
                                      start_date=Date(freq='d', year=2005, month=4, day=1),
                                      length=5)
                                )

        assert(tseries._timeseriescompat_multiple(seriesM_10, seriesM_10, seriesM_10))

        try:
            tseries._timeseriescompat_multiple(seriesM_10, seriesD_10)
            exception = False
        except:
            exception = True
        assert(exception)

        try:
            tseries._timeseriescompat_multiple(seriesD_5, seriesD_10)
            exception = False
        except:
            exception = True
        assert(exception)

        try:
            tseries._timeseriescompat_multiple(seriesD_5, seriesD_5_apr)
            exception = False
        except:
            exception = True
        assert(exception)

    def test_compressed(self):
        "Tests compress"
        dlist = ['2007-01-%02i' % i for i in range(1,16)]
        dates = date_array_fromlist(dlist)
        data = masked_array(numeric.arange(15), mask=[1,0,0,0,0]*3, dtype=float_)
        series = time_series(data, dlist)
        #
        keeper = N.array([0,1,1,1,1]*3, dtype=bool_)
        c_series = series.compressed()
        assert_equal(c_series._data, [1,2,3,4,6,7,8,9,11,12,13,14])
        assert_equal(c_series._mask, nomask)
        assert_equal(c_series._dates, dates[keeper])
        #
        series_st = time_series(MA.column_stack((data,data[::-1])),
                                dates=dates)
        c_series = series_st.compressed()
        d = [1,2,3,6,7,8,11,12,13]
        assert_equal(c_series._data, N.c_[(d,list(reversed(d)))])
        assert_equal(c_series._mask, nomask)
        assert_equal(c_series._dates, dates[d])

    def test_concatenate(self):
        "Tests concatenate"
        dlist = ['2007-%02i' % i for i in range(1,6)]
        dates = date_array_fromlist(dlist)
        data = masked_array(numeric.arange(5), mask=[1,0,0,0,0], dtype=float_)
        #
        ser_1 = time_series(data, dates)
        ser_2 = time_series(data, dates=dates+10)
        newseries = concatenate_series(ser_1, ser_2)
        assert_equal(newseries._data,[0,1,2,3,4,0,0,0,0,0,0,1,2,3,4])
        assert_equal(newseries._mask,[1,0,0,0,0]+[1]*5+[1,0,0,0,0])
         #
        ser_1 = time_series(data, dates)
        ser_2 = time_series(data, dates=dates+10)
        newseries = concatenate_series(ser_1, ser_2, keep_gap=False)
        assert_equal(newseries._data,[0,1,2,3,4,0,1,2,3,4])
        assert_equal(newseries._mask,[1,0,0,0,0]+[1,0,0,0,0])
        assert newseries.has_missing_dates()
        #
        ser_2 = time_series(data, dates=dates+3)
        newseries = concatenate_series(ser_1, ser_2)
        assert_equal(newseries._data,[0,1,2,0,1,2,3,4])
        assert_equal(newseries._mask,[1,0,0,1,0,0,0,0])
        #



###############################################################################
#------------------------------------------------------------------------------
if __name__ == "__main__":
    NumpyTest().run()