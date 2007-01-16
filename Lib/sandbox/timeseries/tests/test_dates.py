# pylint: disable-msg=W0611, W0612, W0511,R0201
"""Tests suite for timeseries.tdates.

:author: Pierre Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknow_ca_at_hotmail_dot_com
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import types
import datetime

import numpy
import numpy.core.fromnumeric  as fromnumeric
import numpy.core.numeric as numeric
from numpy.testing import NumpyTest, NumpyTestCase
from numpy.testing.utils import build_err_msg

import maskedarray
from maskedarray import masked_array

import maskedarray.testutils
from maskedarray.testutils import assert_equal, assert_array_equal

import tdates
import tcore
reload(tdates)
from tdates import date_array_fromlist, Date, DateArray, mxDFromString

class test_creation(NumpyTestCase):
    "Base test class for MaskedArrays."
    
    def __init__(self, *args, **kwds):
        NumpyTestCase.__init__(self, *args, **kwds)
    
    def test_fromstrings(self):
        "Tests creation from list of strings"
        dlist = ['2007-01-%02i' % i for i in range(1,15)]
        # A simple case: daily data
        dates = date_array_fromlist(dlist, 'D')
        assert_equal(dates.freq,'D')
        assert(dates.isfull())
        assert(not dates.has_duplicated_dates())
        assert_equal(dates, 732677+numpy.arange(len(dlist)))
        # as simple, but we need to guess the frequency this time
        dates = date_array_fromlist(dlist, 'D')
        assert_equal(dates.freq,'D')
        assert(dates.isfull())
        assert(not dates.has_duplicated_dates())
        assert_equal(dates, 732677+numpy.arange(len(dlist)))
        # Still daily data, that we force to month
        dates = date_array_fromlist(dlist, 'M')
        assert_equal(dates.freq,'M')
        assert(not dates.isfull())
        assert(dates.has_duplicated_dates())
        assert_equal(dates, [24073]*len(dlist))
        # Now, for monthly data
        dlist = ['2007-%02i' % i for i in range(1,13)]
        dates = date_array_fromlist(dlist, 'M')
        assert_equal(dates.freq,'M')
        assert(dates.isfull())
        assert(not dates.has_duplicated_dates())
        assert_equal(dates, 24073 + numpy.arange(12))
        # Monthly data  w/ guessing
        dlist = ['2007-%02i' % i for i in range(1,13)]
        dates = date_array_fromlist(dlist, )
        assert_equal(dates.freq,'M')
        assert(dates.isfull())
        assert(not dates.has_duplicated_dates())
        assert_equal(dates, 24073 + numpy.arange(12))
        
    def test_fromstrings_wmissing(self):
        "Tests creation from list of strings w/ missing dates"
        dlist = ['2007-01-%02i' % i for i in (1,2,4,5,7,8,10,11,13)]
        dates = date_array_fromlist(dlist)
        assert_equal(dates.freq,'U')
        assert(not dates.isfull())
        assert(not dates.has_duplicated_dates())
        assert_equal(dates.tovalue(),732676+numpy.array([1,2,4,5,7,8,10,11,13]))
        #
        ddates = date_array_fromlist(dlist, 'D')
        assert_equal(ddates.freq,'D')
        assert(not ddates.isfull())
        assert(not ddates.has_duplicated_dates())
        #
        mdates = date_array_fromlist(dlist, 'M')
        assert_equal(mdates.freq,'M')
        assert(not dates.isfull())
        assert(mdates.has_duplicated_dates())
        #
    
    def test_fromsobjects(self):
        "Tests creation from list of objects."
        dlist = ['2007-01-%02i' % i for i in (1,2,4,5,7,8,10,11,13)]
        dates = date_array_fromlist(dlist)
        dobj = [datetime.datetime.fromordinal(d) for d in dates.toordinal()]
        odates = date_array_fromlist(dobj)
        assert_equal(dates,odates)
        dobj = [mxDFromString(d) for d in dlist]
        odates = date_array_fromlist(dobj)
        assert_equal(dates,odates)

    def test_consistent_value(self):
        "Tests that values don't get mutated when constructing dates from a value"
        freqs = [x for x in list(tcore.fmtfreq_dict) if x != 'U']
        for f in freqs:
            today = tdates.thisday(f)
            assert(tdates.Date(freq=f, value=today.value) == today)


class test_freq_conversion(NumpyTestCase):
    "Test frequency conversion of date objects"
    
    def __init__(self, *args, **kwds):
        NumpyTestCase.__init__(self, *args, **kwds)       
        
    def test_conv_annual(self):
        "frequency conversion tests: from Annual Frequency"

        date_A = Date(freq='A', year=2007)
        date_A_to_Q_before = Date(freq='Q', year=2007, quarter=1)
        date_A_to_Q_after = Date(freq='Q', year=2007, quarter=4)
        date_A_to_M_before = Date(freq='M', year=2007, month=1)
        date_A_to_M_after = Date(freq='M', year=2007, month=12)
        date_A_to_W_before = Date(freq='W', year=2007, month=1, day=1)
        date_A_to_W_after = Date(freq='W', year=2007, month=12, day=31)
        date_A_to_B_before = Date(freq='B', year=2007, month=1, day=1)
        date_A_to_B_after = Date(freq='B', year=2007, month=12, day=31)
        date_A_to_D_before = Date(freq='D', year=2007, month=1, day=1)
        date_A_to_D_after = Date(freq='D', year=2007, month=12, day=31)
        date_A_to_H_before = Date(freq='H', year=2007, month=1, day=1, 
                                  hours=0)
        date_A_to_H_after = Date(freq='H', year=2007, month=12, day=31, 
                                 hours=23)
        date_A_to_T_before = Date(freq='T', year=2007, month=1, day=1, 
                                  hours=0, minutes=0)
        date_A_to_T_after = Date(freq='T', year=2007, month=12, day=31, 
                                 hours=23, minutes=59)
        date_A_to_S_before = Date(freq='S', year=2007, month=1, day=1, 
                                  hours=0, minutes=0, seconds=0)
        date_A_to_S_after = Date(freq='S', year=2007, month=12, day=31, 
                                 hours=23, minutes=59, seconds=59)
        
        assert_equal(date_A.asfreq('Q', "BEFORE"), date_A_to_Q_before)
        assert_equal(date_A.asfreq('Q', "AFTER"), date_A_to_Q_after)
        assert_equal(date_A.asfreq('M', "BEFORE"), date_A_to_M_before)
        assert_equal(date_A.asfreq('M', "AFTER"), date_A_to_M_after)
        assert_equal(date_A.asfreq('W', "BEFORE"), date_A_to_W_before)
        assert_equal(date_A.asfreq('W', "AFTER"), date_A_to_W_after)
        assert_equal(date_A.asfreq('B', "BEFORE"), date_A_to_B_before)
        assert_equal(date_A.asfreq('B', "AFTER"), date_A_to_B_after)
        assert_equal(date_A.asfreq('D', "BEFORE"), date_A_to_D_before)
        assert_equal(date_A.asfreq('D', "AFTER"), date_A_to_D_after)
        assert_equal(date_A.asfreq('H', "BEFORE"), date_A_to_H_before)
        assert_equal(date_A.asfreq('H', "AFTER"), date_A_to_H_after)
        assert_equal(date_A.asfreq('T', "BEFORE"), date_A_to_T_before)
        assert_equal(date_A.asfreq('T', "AFTER"), date_A_to_T_after)
        assert_equal(date_A.asfreq('S', "BEFORE"), date_A_to_S_before)
        assert_equal(date_A.asfreq('S', "AFTER"), date_A_to_S_after)

        
    def test_conv_quarterly(self):
        "frequency conversion tests: from Quarterly Frequency"

        date_Q = Date(freq='Q', year=2007, quarter=1)
        date_Q_end_of_year = Date(freq='Q', year=2007, quarter=4)
        date_Q_to_A = Date(freq='A', year=2007)
        date_Q_to_M_before = Date(freq='M', year=2007, month=1)
        date_Q_to_M_after = Date(freq='M', year=2007, month=3)
        date_Q_to_W_before = Date(freq='W', year=2007, month=1, day=1)
        date_Q_to_W_after = Date(freq='W', year=2007, month=3, day=31)
        date_Q_to_B_before = Date(freq='B', year=2007, month=1, day=1)
        date_Q_to_B_after = Date(freq='B', year=2007, month=3, day=30)
        date_Q_to_D_before = Date(freq='D', year=2007, month=1, day=1)
        date_Q_to_D_after = Date(freq='D', year=2007, month=3, day=31)
        date_Q_to_H_before = Date(freq='H', year=2007, month=1, day=1, 
                                  hours=0)
        date_Q_to_H_after = Date(freq='H', year=2007, month=3, day=31, 
                                 hours=23)
        date_Q_to_T_before = Date(freq='T', year=2007, month=1, day=1, 
                                  hours=0, minutes=0)
        date_Q_to_T_after = Date(freq='T', year=2007, month=3, day=31, 
                                 hours=23, minutes=59)
        date_Q_to_S_before = Date(freq='S', year=2007, month=1, day=1, 
                                  hours=0, minutes=0, seconds=0)
        date_Q_to_S_after = Date(freq='S', year=2007, month=3, day=31, 
                                 hours=23, minutes=59, seconds=59)
        
        assert_equal(date_Q.asfreq('A'), date_Q_to_A)
        assert_equal(date_Q_end_of_year.asfreq('A'), date_Q_to_A)
        
        assert_equal(date_Q.asfreq('M', "BEFORE"), date_Q_to_M_before)
        assert_equal(date_Q.asfreq('M', "AFTER"), date_Q_to_M_after)
        assert_equal(date_Q.asfreq('W', "BEFORE"), date_Q_to_W_before)
        assert_equal(date_Q.asfreq('W', "AFTER"), date_Q_to_W_after)
        assert_equal(date_Q.asfreq('B', "BEFORE"), date_Q_to_B_before)
        assert_equal(date_Q.asfreq('B', "AFTER"), date_Q_to_B_after)
        assert_equal(date_Q.asfreq('D', "BEFORE"), date_Q_to_D_before)
        assert_equal(date_Q.asfreq('D', "AFTER"), date_Q_to_D_after)
        assert_equal(date_Q.asfreq('H', "BEFORE"), date_Q_to_H_before)
        assert_equal(date_Q.asfreq('H', "AFTER"), date_Q_to_H_after)
        assert_equal(date_Q.asfreq('T', "BEFORE"), date_Q_to_T_before)
        assert_equal(date_Q.asfreq('T', "AFTER"), date_Q_to_T_after)
        assert_equal(date_Q.asfreq('S', "BEFORE"), date_Q_to_S_before)
        assert_equal(date_Q.asfreq('S', "AFTER"), date_Q_to_S_after)
        

    def test_conv_monthly(self):
        "frequency conversion tests: from Monthly Frequency"
        
        date_M = Date(freq='M', year=2007, month=1)
        date_M_end_of_year = Date(freq='M', year=2007, month=12)
        date_M_end_of_quarter = Date(freq='M', year=2007, month=3)
        date_M_to_A = Date(freq='A', year=2007)
        date_M_to_Q = Date(freq='Q', year=2007, quarter=1)
        date_M_to_W_before = Date(freq='W', year=2007, month=1, day=1)
        date_M_to_W_after = Date(freq='W', year=2007, month=1, day=31)
        date_M_to_B_before = Date(freq='B', year=2007, month=1, day=1)
        date_M_to_B_after = Date(freq='B', year=2007, month=1, day=31)
        date_M_to_D_before = Date(freq='D', year=2007, month=1, day=1)
        date_M_to_D_after = Date(freq='D', year=2007, month=1, day=31)
        date_M_to_H_before = Date(freq='H', year=2007, month=1, day=1, 
                                  hours=0)
        date_M_to_H_after = Date(freq='H', year=2007, month=1, day=31, 
                                 hours=23)
        date_M_to_T_before = Date(freq='T', year=2007, month=1, day=1, 
                                  hours=0, minutes=0)
        date_M_to_T_after = Date(freq='T', year=2007, month=1, day=31, 
                                 hours=23, minutes=59)
        date_M_to_S_before = Date(freq='S', year=2007, month=1, day=1, 
                                  hours=0, minutes=0, seconds=0)
        date_M_to_S_after = Date(freq='S', year=2007, month=1, day=31, 
                                 hours=23, minutes=59, seconds=59)
        
        assert_equal(date_M.asfreq('A'), date_M_to_A)
        assert_equal(date_M_end_of_year.asfreq('A'), date_M_to_A)
        assert_equal(date_M.asfreq('Q'), date_M_to_Q)
        assert_equal(date_M_end_of_quarter.asfreq('Q'), date_M_to_Q)

        assert_equal(date_M.asfreq('W', "BEFORE"), date_M_to_W_before)
        assert_equal(date_M.asfreq('W', "AFTER"), date_M_to_W_after)
        assert_equal(date_M.asfreq('B', "BEFORE"), date_M_to_B_before)
        assert_equal(date_M.asfreq('B', "AFTER"), date_M_to_B_after)
        assert_equal(date_M.asfreq('D', "BEFORE"), date_M_to_D_before)
        assert_equal(date_M.asfreq('D', "AFTER"), date_M_to_D_after)
        assert_equal(date_M.asfreq('H', "BEFORE"), date_M_to_H_before)
        assert_equal(date_M.asfreq('H', "AFTER"), date_M_to_H_after)
        assert_equal(date_M.asfreq('T', "BEFORE"), date_M_to_T_before)
        assert_equal(date_M.asfreq('T', "AFTER"), date_M_to_T_after)
        assert_equal(date_M.asfreq('S', "BEFORE"), date_M_to_S_before)
        assert_equal(date_M.asfreq('S', "AFTER"), date_M_to_S_after)

        
    def test_conv_weekly(self):
        "frequency conversion tests: from Weekly Frequency"
        
        date_W = Date(freq='W', year=2007, month=1, day=1)
        date_W_end_of_year = Date(freq='W', year=2007, month=12, day=31)
        date_W_end_of_quarter = Date(freq='W', year=2007, month=3, day=31)
        date_W_end_of_month = Date(freq='W', year=2007, month=1, day=31)
        date_W_to_A = Date(freq='A', year=2007)
        date_W_to_Q = Date(freq='Q', year=2007, quarter=1)
        date_W_to_M = Date(freq='M', year=2007, month=1)

        if Date(freq='D', year=2007, month=12, day=31).day_of_week() == 6:
            date_W_to_A_end_of_year = Date(freq='A', year=2007)
        else:
            date_W_to_A_end_of_year = Date(freq='A', year=2008)

        if Date(freq='D', year=2007, month=3, day=31).day_of_week() == 6:
            date_W_to_Q_end_of_quarter = Date(freq='Q', year=2007, quarter=1)
        else:
            date_W_to_Q_end_of_quarter = Date(freq='Q', year=2007, quarter=2)

        if Date(freq='D', year=2007, month=1, day=31).day_of_week() == 6:
            date_W_to_M_end_of_month = Date(freq='M', year=2007, month=1)
        else:
            date_W_to_M_end_of_month = Date(freq='M', year=2007, month=2)

        date_W_to_B_before = Date(freq='B', year=2007, month=1, day=1)
        date_W_to_B_after = Date(freq='B', year=2007, month=1, day=5)
        date_W_to_D_before = Date(freq='D', year=2007, month=1, day=1)
        date_W_to_D_after = Date(freq='D', year=2007, month=1, day=7)
        date_W_to_H_before = Date(freq='H', year=2007, month=1, day=1, 
                                  hours=0)
        date_W_to_H_after = Date(freq='H', year=2007, month=1, day=7, 
                                 hours=23)
        date_W_to_T_before = Date(freq='T', year=2007, month=1, day=1, 
                                  hours=0, minutes=0)
        date_W_to_T_after = Date(freq='T', year=2007, month=1, day=7, 
                                 hours=23, minutes=59)
        date_W_to_S_before = Date(freq='S', year=2007, month=1, day=1, 
                                  hours=0, minutes=0, seconds=0)
        date_W_to_S_after = Date(freq='S', year=2007, month=1, day=7, 
                                 hours=23, minutes=59, seconds=59)
        
        assert_equal(date_W.asfreq('A'), date_W_to_A)
        assert_equal(date_W_end_of_year.asfreq('A'), date_W_to_A_end_of_year)
        assert_equal(date_W.asfreq('Q'), date_W_to_Q)
        assert_equal(date_W_end_of_quarter.asfreq('Q'), date_W_to_Q_end_of_quarter)
        assert_equal(date_W.asfreq('M'), date_W_to_M)
        assert_equal(date_W_end_of_month.asfreq('M'), date_W_to_M_end_of_month)

        assert_equal(date_W.asfreq('B', "BEFORE"), date_W_to_B_before)
        assert_equal(date_W.asfreq('B', "AFTER"), date_W_to_B_after)
        assert_equal(date_W.asfreq('D', "BEFORE"), date_W_to_D_before)
        assert_equal(date_W.asfreq('D', "AFTER"), date_W_to_D_after)
        assert_equal(date_W.asfreq('H', "BEFORE"), date_W_to_H_before)
        assert_equal(date_W.asfreq('H', "AFTER"), date_W_to_H_after)
        assert_equal(date_W.asfreq('T', "BEFORE"), date_W_to_T_before)
        assert_equal(date_W.asfreq('T', "AFTER"), date_W_to_T_after)
        assert_equal(date_W.asfreq('S', "BEFORE"), date_W_to_S_before)
        assert_equal(date_W.asfreq('S', "AFTER"), date_W_to_S_after)
        
    def test_conv_business(self):
        "frequency conversion tests: from Business Frequency"
        
        date_B = Date(freq='B', year=2007, month=1, day=1)
        date_B_end_of_year = Date(freq='B', year=2007, month=12, day=31)
        date_B_end_of_quarter = Date(freq='B', year=2007, month=3, day=30)
        date_B_end_of_month = Date(freq='B', year=2007, month=1, day=31)
        date_B_end_of_week = Date(freq='B', year=2007, month=1, day=5)
        
        date_B_to_A = Date(freq='A', year=2007)
        date_B_to_Q = Date(freq='Q', year=2007, quarter=1)
        date_B_to_M = Date(freq='M', year=2007, month=1)
        date_B_to_W = Date(freq='W', year=2007, month=1, day=7)
        date_B_to_D = Date(freq='D', year=2007, month=1, day=1)
        date_B_to_H_before = Date(freq='H', year=2007, month=1, day=1, 
                                  hours=0)
        date_B_to_H_after = Date(freq='H', year=2007, month=1, day=1, 
                                 hours=23)
        date_B_to_T_before = Date(freq='T', year=2007, month=1, day=1, 
                                  hours=0, minutes=0)
        date_B_to_T_after = Date(freq='T', year=2007, month=1, day=1, 
                                 hours=23, minutes=59)
        date_B_to_S_before = Date(freq='S', year=2007, month=1, day=1, 
                                  hours=0, minutes=0, seconds=0)
        date_B_to_S_after = Date(freq='S', year=2007, month=1, day=1, 
                                 hours=23, minutes=59, seconds=59)
        
        assert_equal(date_B.asfreq('A'), date_B_to_A)
        assert_equal(date_B_end_of_year.asfreq('A'), date_B_to_A)
        assert_equal(date_B.asfreq('Q'), date_B_to_Q)
        assert_equal(date_B_end_of_quarter.asfreq('Q'), date_B_to_Q)
        assert_equal(date_B.asfreq('M'), date_B_to_M)
        assert_equal(date_B_end_of_month.asfreq('M'), date_B_to_M)
        assert_equal(date_B.asfreq('W'), date_B_to_W)
        assert_equal(date_B_end_of_week.asfreq('W'), date_B_to_W)

        assert_equal(date_B.asfreq('D'), date_B_to_D)

        assert_equal(date_B.asfreq('H', "BEFORE"), date_B_to_H_before)
        assert_equal(date_B.asfreq('H', "AFTER"), date_B_to_H_after)
        assert_equal(date_B.asfreq('T', "BEFORE"), date_B_to_T_before)
        assert_equal(date_B.asfreq('T', "AFTER"), date_B_to_T_after)
        assert_equal(date_B.asfreq('S', "BEFORE"), date_B_to_S_before)
        assert_equal(date_B.asfreq('S', "AFTER"), date_B_to_S_after)

    def test_conv_daily(self):
        "frequency conversion tests: from Business Frequency"
        
        date_D = Date(freq='D', year=2007, month=1, day=1)
        date_D_end_of_year = Date(freq='D', year=2007, month=12, day=31)
        date_D_end_of_quarter = Date(freq='D', year=2007, month=3, day=31)
        date_D_end_of_month = Date(freq='D', year=2007, month=1, day=31)
        date_D_end_of_week = Date(freq='D', year=2007, month=1, day=7)
        
        date_D_friday = Date(freq='D', year=2007, month=1, day=5)
        date_D_saturday = Date(freq='D', year=2007, month=1, day=6)
        date_D_sunday = Date(freq='D', year=2007, month=1, day=7)
        date_D_monday = Date(freq='D', year=2007, month=1, day=8)
        
        date_B_friday = Date(freq='B', year=2007, month=1, day=5)
        date_B_monday = Date(freq='B', year=2007, month=1, day=8)
        
        date_D_to_A = Date(freq='A', year=2007)
        date_D_to_Q = Date(freq='Q', year=2007, quarter=1)
        date_D_to_M = Date(freq='M', year=2007, month=1)
        date_D_to_W = Date(freq='W', year=2007, month=1, day=7)

        date_D_to_H_before = Date(freq='H', year=2007, month=1, day=1, 
                                  hours=0)
        date_D_to_H_after = Date(freq='H', year=2007, month=1, day=1, 
                                 hours=23)
        date_D_to_T_before = Date(freq='T', year=2007, month=1, day=1, 
                                  hours=0, minutes=0)
        date_D_to_T_after = Date(freq='T', year=2007, month=1, day=1, 
                                 hours=23, minutes=59)
        date_D_to_S_before = Date(freq='S', year=2007, month=1, day=1, 
                                  hours=0, minutes=0, seconds=0)
        date_D_to_S_after = Date(freq='S', year=2007, month=1, day=1, 
                                 hours=23, minutes=59, seconds=59)
        
        assert_equal(date_D.asfreq('A'), date_D_to_A)
        assert_equal(date_D_end_of_year.asfreq('A'), date_D_to_A)
        assert_equal(date_D.asfreq('Q'), date_D_to_Q)
        assert_equal(date_D_end_of_quarter.asfreq('Q'), date_D_to_Q)
        assert_equal(date_D.asfreq('M'), date_D_to_M)
        assert_equal(date_D_end_of_month.asfreq('M'), date_D_to_M)
        assert_equal(date_D.asfreq('W'), date_D_to_W)
        assert_equal(date_D_end_of_week.asfreq('W'), date_D_to_W)

        assert_equal(date_D_friday.asfreq('B'), date_B_friday)
        assert_equal(date_D_saturday.asfreq('B', "BEFORE"), date_B_friday)
        assert_equal(date_D_saturday.asfreq('B', "AFTER"), date_B_monday)
        assert_equal(date_D_sunday.asfreq('B', "BEFORE"), date_B_friday)
        assert_equal(date_D_sunday.asfreq('B', "AFTER"), date_B_monday)

        assert_equal(date_D.asfreq('H', "BEFORE"), date_D_to_H_before)
        assert_equal(date_D.asfreq('H', "AFTER"), date_D_to_H_after)
        assert_equal(date_D.asfreq('T', "BEFORE"), date_D_to_T_before)
        assert_equal(date_D.asfreq('T', "AFTER"), date_D_to_T_after)
        assert_equal(date_D.asfreq('S', "BEFORE"), date_D_to_S_before)
        assert_equal(date_D.asfreq('S', "AFTER"), date_D_to_S_after)

    def test_conv_hourly(self):
        "frequency conversion tests: from Hourly Frequency"
        
        date_H = Date(freq='H', year=2007, month=1, day=1, hours=0)
        date_H_end_of_year = Date(freq='H', year=2007, month=12, day=31, 
                                  hours=23)
        date_H_end_of_quarter = Date(freq='H', year=2007, month=3, day=31, 
                                     hours=23)
        date_H_end_of_month = Date(freq='H', year=2007, month=1, day=31, 
                                   hours=23)
        date_H_end_of_week = Date(freq='H', year=2007, month=1, day=7, 
                                  hours=23)
        date_H_end_of_day = Date(freq='H', year=2007, month=1, day=1, 
                                 hours=23)
        date_H_end_of_bus = Date(freq='H', year=2007, month=1, day=1, 
                                 hours=23)
        
        date_H_to_A = Date(freq='A', year=2007)
        date_H_to_Q = Date(freq='Q', year=2007, quarter=1)
        date_H_to_M = Date(freq='M', year=2007, month=1)
        date_H_to_W = Date(freq='W', year=2007, month=1, day=7)
        date_H_to_D = Date(freq='D', year=2007, month=1, day=1)
        date_H_to_B = Date(freq='B', year=2007, month=1, day=1)
        
        date_H_to_T_before = Date(freq='T', year=2007, month=1, day=1, 
                                  hours=0, minutes=0)
        date_H_to_T_after = Date(freq='T', year=2007, month=1, day=1, 
                                 hours=0, minutes=59)
        date_H_to_S_before = Date(freq='S', year=2007, month=1, day=1, 
                                  hours=0, minutes=0, seconds=0)
        date_H_to_S_after = Date(freq='S', year=2007, month=1, day=1, 
                                 hours=0, minutes=59, seconds=59)
        
        assert_equal(date_H.asfreq('A'), date_H_to_A)
        assert_equal(date_H_end_of_year.asfreq('A'), date_H_to_A)
        assert_equal(date_H.asfreq('Q'), date_H_to_Q)
        assert_equal(date_H_end_of_quarter.asfreq('Q'), date_H_to_Q)
        assert_equal(date_H.asfreq('M'), date_H_to_M)
        assert_equal(date_H_end_of_month.asfreq('M'), date_H_to_M)
        assert_equal(date_H.asfreq('W'), date_H_to_W)
        assert_equal(date_H_end_of_week.asfreq('W'), date_H_to_W)
        assert_equal(date_H.asfreq('D'), date_H_to_D)
        assert_equal(date_H_end_of_day.asfreq('D'), date_H_to_D)
        assert_equal(date_H.asfreq('B'), date_H_to_B)
        assert_equal(date_H_end_of_bus.asfreq('B'), date_H_to_B)

        assert_equal(date_H.asfreq('T', "BEFORE"), date_H_to_T_before)
        assert_equal(date_H.asfreq('T', "AFTER"), date_H_to_T_after)
        assert_equal(date_H.asfreq('S', "BEFORE"), date_H_to_S_before)
        assert_equal(date_H.asfreq('S', "AFTER"), date_H_to_S_after)

    def test_conv_minutely(self):
        "frequency conversion tests: from Minutely Frequency"
        
        date_T = Date(freq='T', year=2007, month=1, day=1, 
                      hours=0, minutes=0)
        date_T_end_of_year = Date(freq='T', year=2007, month=12, day=31, 
                                  hours=23, minutes=59)
        date_T_end_of_quarter = Date(freq='T', year=2007, month=3, day=31, 
                                     hours=23, minutes=59)
        date_T_end_of_month = Date(freq='T', year=2007, month=1, day=31, 
                                   hours=23, minutes=59)
        date_T_end_of_week = Date(freq='T', year=2007, month=1, day=7, 
                                  hours=23, minutes=59)
        date_T_end_of_day = Date(freq='T', year=2007, month=1, day=1, 
                                 hours=23, minutes=59)
        date_T_end_of_bus = Date(freq='T', year=2007, month=1, day=1, 
                                 hours=23, minutes=59)
        date_T_end_of_hour = Date(freq='T', year=2007, month=1, day=1, 
                                  hours=0, minutes=59)
        
        date_T_to_A = Date(freq='A', year=2007)
        date_T_to_Q = Date(freq='Q', year=2007, quarter=1)
        date_T_to_M = Date(freq='M', year=2007, month=1)
        date_T_to_W = Date(freq='W', year=2007, month=1, day=7)
        date_T_to_D = Date(freq='D', year=2007, month=1, day=1)
        date_T_to_B = Date(freq='B', year=2007, month=1, day=1)
        date_T_to_H = Date(freq='H', year=2007, month=1, day=1, hours=0)
        
        date_T_to_S_before = Date(freq='S', year=2007, month=1, day=1, 
                                  hours=0, minutes=0, seconds=0)
        date_T_to_S_after = Date(freq='S', year=2007, month=1, day=1, 
                                 hours=0, minutes=0, seconds=59)
        
        assert_equal(date_T.asfreq('A'), date_T_to_A)
        assert_equal(date_T_end_of_year.asfreq('A'), date_T_to_A)
        assert_equal(date_T.asfreq('Q'), date_T_to_Q)
        assert_equal(date_T_end_of_quarter.asfreq('Q'), date_T_to_Q)
        assert_equal(date_T.asfreq('M'), date_T_to_M)
        assert_equal(date_T_end_of_month.asfreq('M'), date_T_to_M)
        assert_equal(date_T.asfreq('W'), date_T_to_W)
        assert_equal(date_T_end_of_week.asfreq('W'), date_T_to_W)
        assert_equal(date_T.asfreq('D'), date_T_to_D)
        assert_equal(date_T_end_of_day.asfreq('D'), date_T_to_D)
        assert_equal(date_T.asfreq('B'), date_T_to_B)
        assert_equal(date_T_end_of_bus.asfreq('B'), date_T_to_B)
        assert_equal(date_T.asfreq('H'), date_T_to_H)
        assert_equal(date_T_end_of_hour.asfreq('H'), date_T_to_H)

        assert_equal(date_T.asfreq('S', "BEFORE"), date_T_to_S_before)
        assert_equal(date_T.asfreq('S', "AFTER"), date_T_to_S_after)


    def test_conv_secondly(self):
        "frequency conversion tests: from Secondly Frequency"
        
        date_S = Date(freq='S', year=2007, month=1, day=1, 
                      hours=0, minutes=0, seconds=0)
        date_S_end_of_year = Date(freq='S', year=2007, month=12, day=31, 
                                  hours=23, minutes=59, seconds=59)
        date_S_end_of_quarter = Date(freq='S', year=2007, month=3, day=31, 
                                     hours=23, minutes=59, seconds=59)
        date_S_end_of_month = Date(freq='S', year=2007, month=1, day=31, 
                                   hours=23, minutes=59, seconds=59)
        date_S_end_of_week = Date(freq='S', year=2007, month=1, day=7, 
                                  hours=23, minutes=59, seconds=59)
        date_S_end_of_day = Date(freq='S', year=2007, month=1, day=1, 
                                 hours=23, minutes=59, seconds=59)
        date_S_end_of_bus = Date(freq='S', year=2007, month=1, day=1, 
                                 hours=23, minutes=59, seconds=59)
        date_S_end_of_hour = Date(freq='S', year=2007, month=1, day=1, 
                                  hours=0, minutes=59, seconds=59)
        date_S_end_of_minute = Date(freq='S', year=2007, month=1, day=1, 
                                    hours=0, minutes=0, seconds=59)
        
        date_S_to_A = Date(freq='A', year=2007)
        date_S_to_Q = Date(freq='Q', year=2007, quarter=1)
        date_S_to_M = Date(freq='M', year=2007, month=1)
        date_S_to_W = Date(freq='W', year=2007, month=1, day=7)
        date_S_to_D = Date(freq='D', year=2007, month=1, day=1)
        date_S_to_B = Date(freq='B', year=2007, month=1, day=1)
        date_S_to_H = Date(freq='H', year=2007, month=1, day=1, 
                           hours=0)
        date_S_to_T = Date(freq='T', year=2007, month=1, day=1, 
                           hours=0, minutes=0)
        
        assert_equal(date_S.asfreq('A'), date_S_to_A)
        assert_equal(date_S_end_of_year.asfreq('A'), date_S_to_A)
        assert_equal(date_S.asfreq('Q'), date_S_to_Q)
        assert_equal(date_S_end_of_quarter.asfreq('Q'), date_S_to_Q)
        assert_equal(date_S.asfreq('M'), date_S_to_M)
        assert_equal(date_S_end_of_month.asfreq('M'), date_S_to_M)
        assert_equal(date_S.asfreq('W'), date_S_to_W)
        assert_equal(date_S_end_of_week.asfreq('W'), date_S_to_W)
        assert_equal(date_S.asfreq('D'), date_S_to_D)
        assert_equal(date_S_end_of_day.asfreq('D'), date_S_to_D)
        assert_equal(date_S.asfreq('B'), date_S_to_B)
        assert_equal(date_S_end_of_bus.asfreq('B'), date_S_to_B)
        assert_equal(date_S.asfreq('H'), date_S_to_H)
        assert_equal(date_S_end_of_hour.asfreq('H'), date_S_to_H)
        assert_equal(date_S.asfreq('T'), date_S_to_T)
        assert_equal(date_S_end_of_minute.asfreq('T'), date_S_to_T)
        

class test_methods(NumpyTestCase):
    "Base test class for MaskedArrays."
    
    def __init__(self, *args, **kwds):
        NumpyTestCase.__init__(self, *args, **kwds)       
        
    def test_getitem(self):
        "Tests getitem"
        dlist = ['2007-%02i' % i for i in range(1,5)+range(7,13)]
        mdates = date_array_fromlist(dlist, 'M')
        # Using an integer
        assert_equal(mdates[0].value, 24073)
        assert_equal(mdates[-1].value, 24084)
        # Using a date
        lag = mdates.find_dates(mdates[0])
        assert_equal(mdates[lag], mdates[0])
        lag = mdates.find_dates(Date('M',value=24080))
        assert_equal(mdates[lag], mdates[5])
        # Using several dates
        lag = mdates.find_dates(Date('M',value=24073), Date('M',value=24084))
        assert_equal(mdates[lag], 
                     DateArray([mdates[0], mdates[-1]], freq='M'))
        assert_equal(mdates[[mdates[0],mdates[-1]]], mdates[lag])
        #
        assert_equal(mdates>=mdates[-4], [0,0,0,0,0,0,1,1,1,1])
        dlist = ['2006-%02i' % i for i in range(1,5)+range(7,13)]
        mdates = date_array_fromlist(dlist).asfreq('M')

        
    def test_getsteps(self):
        dlist = ['2007-01-%02i' %i for i in (1,2,3,4,8,9,10,11,12,15)]
        ddates = date_array_fromlist(dlist)
        assert_equal(ddates.get_steps(), [1,1,1,4,1,1,1,1,3])
        
        


###############################################################################
#------------------------------------------------------------------------------
if __name__ == "__main__":
    NumpyTest().run()        