# pylint: disable-msg=W0611, W0612, W0511,R0201
"""Tests suite for Date handling.

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

import timeseries as ts
from timeseries import const as C
from timeseries.parser import DateFromString, DateTimeFromString
from timeseries import Date, DateArray,\
    thisday, today, date_array, date_array_fromlist
from timeseries.cseries import freq_dict


class TestCreation(NumpyTestCase):
    "Base test class for MaskedArrays."

    def __init__(self, *args, **kwds):
        NumpyTestCase.__init__(self, *args, **kwds)

    def test_fromstrings(self):
        "Tests creation from list of strings"
        print "starting test_fromstrings..."
        dlist = ['2007-01-%02i' % i for i in range(1,15)]
        # A simple case: daily data
        dates = date_array_fromlist(dlist, 'D')
        assert_equal(dates.freqstr,'D')
        assert(dates.isfull())
        assert(not dates.has_duplicated_dates())
        assert_equal(dates, 732677+numpy.arange(len(dlist)))
        # as simple, but we need to guess the frequency this time
        dates = date_array_fromlist(dlist, 'D')
        assert_equal(dates.freqstr,'D')
        assert(dates.isfull())
        assert(not dates.has_duplicated_dates())
        assert_equal(dates, 732677+numpy.arange(len(dlist)))
        # Still daily data, that we force to month
        dates = date_array_fromlist(dlist, 'M')
        assert_equal(dates.freqstr,'M')
        assert(not dates.isfull())
        assert(dates.has_duplicated_dates())
        assert_equal(dates, [24073]*len(dlist))
        # Now, for monthly data
        dlist = ['2007-%02i' % i for i in range(1,13)]
        dates = date_array_fromlist(dlist, 'M')
        assert_equal(dates.freqstr,'M')
        assert(dates.isfull())
        assert(not dates.has_duplicated_dates())
        assert_equal(dates, 24073 + numpy.arange(12))
        # Monthly data  w/ guessing
        dlist = ['2007-%02i' % i for i in range(1,13)]
        dates = date_array_fromlist(dlist, )
        assert_equal(dates.freqstr,'M')
        assert(dates.isfull())
        assert(not dates.has_duplicated_dates())
        assert_equal(dates, 24073 + numpy.arange(12))
        print "finished test_fromstrings"

    def test_fromstrings_wmissing(self):
        "Tests creation from list of strings w/ missing dates"
        print "starting test_fromstrings_wmissing..."
        dlist = ['2007-01-%02i' % i for i in (1,2,4,5,7,8,10,11,13)]
        dates = date_array_fromlist(dlist)
        assert_equal(dates.freqstr,'U')
        assert(not dates.isfull())
        assert(not dates.has_duplicated_dates())
        assert_equal(dates.tovalue(),732676+numpy.array([1,2,4,5,7,8,10,11,13]))
        #
        ddates = date_array_fromlist(dlist, 'D')
        assert_equal(ddates.freqstr,'D')
        assert(not ddates.isfull())
        assert(not ddates.has_duplicated_dates())
        #
        mdates = date_array_fromlist(dlist, 'M')
        assert_equal(mdates.freqstr,'M')
        assert(not dates.isfull())
        assert(mdates.has_duplicated_dates())
        print "finished test_fromstrings_wmissing"
        #

    def test_fromsobjects(self):
        "Tests creation from list of objects."
        print "starting test_fromsobjects..."
        dlist = ['2007-01-%02i' % i for i in (1,2,4,5,7,8,10,11,13)]
        dates = date_array_fromlist(dlist)
        dobj = [datetime.datetime.fromordinal(d) for d in dates.toordinal()]
        odates = date_array_fromlist(dobj)
        assert_equal(dates,odates)
        dobj = [DateFromString(d) for d in dlist]
        odates = date_array_fromlist(dobj)
        assert_equal(dates,odates)
        #
        D = date_array_fromlist(dlist=['2006-01'])
        assert_equal(D.tovalue(), [732312, ])
        assert_equal(D.freq, C.FR_UND)
        print "finished test_fromsobjects"

    def test_consistent_value(self):
        "Tests that values don't get mutated when constructing dates from a value"
        print "starting test_consistent_value..."
        freqs = [x[0] for x in freq_dict.values() if x[0] != 'U']

        for f in freqs:
            today = thisday(f)
            assert_equal(Date(freq=f, value=today.value), today)
        print "finished test_consistent_value"

    def test_shortcuts(self):
        "Tests some creation shortcuts. Because I'm lazy like that."
        print "starting test_shortcuts..."
        # Dates shortcuts
        assert_equal(Date('D','2007-01'), Date('D',string='2007-01'))
        assert_equal(Date('D','2007-01'), Date('D', value=732677))
        assert_equal(Date('D',732677), Date('D', value=732677))
        # DateArray shortcuts
        n = today('M')
        d = date_array(start_date=n, length=3)
        assert_equal(date_array(n,length=3), d)
        assert_equal(date_array(n, n+2), d)
        print "finished test_shortcuts"

class TestDateProperties(NumpyTestCase):
    "Test properties such as year, month, day_of_week, etc...."

    def __init__(self, *args, **kwds):
        NumpyTestCase.__init__(self, *args, **kwds)

    def test_properties(self):

        a_date = Date(freq='A', year=2007)

        q_date = Date(freq=C.FR_QTREDEC, year=2007, quarter=1)

        qedec_date = Date(freq=C.FR_QTREDEC, year=2007, quarter=1)
        qejan_date = Date(freq=C.FR_QTREJAN, year=2007, quarter=1)
        qejun_date = Date(freq=C.FR_QTREJUN, year=2007, quarter=1)

        qsdec_date = Date(freq=C.FR_QTREDEC, year=2007, quarter=1)
        qsjan_date = Date(freq=C.FR_QTREJAN, year=2007, quarter=1)
        qsjun_date = Date(freq=C.FR_QTREJUN, year=2007, quarter=1)

        m_date = Date(freq='M', year=2007, month=1)
        w_date = Date(freq='W', year=2007, month=1, day=7)
        b_date = Date(freq='B', year=2007, month=1, day=1)
        d_date = Date(freq='D', year=2007, month=1, day=1)
        h_date = Date(freq='H', year=2007, month=1, day=1,
                                       hour=0)
        t_date = Date(freq='T', year=2007, month=1, day=1,
                                       hour=0, minute=0)
        s_date = Date(freq='T', year=2007, month=1, day=1,
                                       hour=0, minute=0, second=0)

        assert_equal(a_date.year, 2007)

        for x in range(3):
            for qd in (qedec_date, qejan_date, qejun_date,
                       qsdec_date, qsjan_date, qsjun_date):
                assert_equal((qd+x).qyear, 2007)
                assert_equal((qd+x).quarter, x+1)

        for x in range(11):

            m_date_x = m_date+x
            assert_equal(m_date_x.year, 2007)

            if   1  <= x + 1 <= 3:  assert_equal(m_date_x.quarter, 1)
            elif 4  <= x + 1 <= 6:  assert_equal(m_date_x.quarter, 2)
            elif 7  <= x + 1 <= 9:  assert_equal(m_date_x.quarter, 3)
            elif 10 <= x + 1 <= 12: assert_equal(m_date_x.quarter, 4)

            assert_equal(m_date_x.month, x+1)

        assert_equal(w_date.year, 2007)
        assert_equal(w_date.quarter, 1)
        assert_equal(w_date.month, 1)
        assert_equal(w_date.week, 1)
        assert_equal((w_date-1).week, 52)

        assert_equal(b_date.year, 2007)
        assert_equal(b_date.quarter, 1)
        assert_equal(b_date.month, 1)
        assert_equal(b_date.day, 1)
        assert_equal(b_date.day_of_week, 0)
        assert_equal(b_date.day_of_year, 1)

        assert_equal(d_date.year, 2007)
        assert_equal(d_date.quarter, 1)
        assert_equal(d_date.month, 1)
        assert_equal(d_date.day, 1)
        assert_equal(d_date.day_of_week, 0)
        assert_equal(d_date.day_of_year, 1)

        assert_equal(h_date.year, 2007)
        assert_equal(h_date.quarter, 1)
        assert_equal(h_date.month, 1)
        assert_equal(h_date.day, 1)
        assert_equal(h_date.day_of_week, 0)
        assert_equal(h_date.day_of_year, 1)
        assert_equal(h_date.hour, 0)

        assert_equal(t_date.year, 2007)
        assert_equal(t_date.quarter, 1)
        assert_equal(t_date.month, 1)
        assert_equal(t_date.day, 1)
        assert_equal(t_date.day_of_week, 0)
        assert_equal(t_date.day_of_year, 1)
        assert_equal(t_date.hour, 0)
        assert_equal(t_date.minute, 0)

        assert_equal(s_date.year, 2007)
        assert_equal(s_date.quarter, 1)
        assert_equal(s_date.month, 1)
        assert_equal(s_date.day, 1)
        assert_equal(s_date.day_of_week, 0)
        assert_equal(s_date.day_of_year, 1)
        assert_equal(s_date.hour, 0)
        assert_equal(s_date.minute, 0)
        assert_equal(s_date.second, 0)


def dArrayWrap(date):
    "wrap a date into a DateArray of length 1"
    return date_array(start_date=date,length=1)

def noWrap(item): return item

class TestFreqConversion(NumpyTestCase):
    "Test frequency conversion of date objects"

    def __init__(self, *args, **kwds):
        NumpyTestCase.__init__(self, *args, **kwds)
        self.dateWrap = [(dArrayWrap, assert_array_equal),
                         (noWrap, assert_equal)]

    def test_conv_annual(self):
        "frequency conversion tests: from Annual Frequency"

        for dWrap, assert_func in self.dateWrap:
            date_A = dWrap(Date(freq='A', year=2007))

            date_AJAN = dWrap(Date(freq=C.FR_ANNJAN, year=2007))
            date_AJUN = dWrap(Date(freq=C.FR_ANNJUN, year=2007))
            date_ANOV = dWrap(Date(freq=C.FR_ANNNOV, year=2007))

            date_A_to_Q_before = dWrap(Date(freq='Q', year=2007, quarter=1))
            date_A_to_Q_after = dWrap(Date(freq='Q', year=2007, quarter=4))
            date_A_to_M_before = dWrap(Date(freq='M', year=2007, month=1))
            date_A_to_M_after = dWrap(Date(freq='M', year=2007, month=12))
            date_A_to_W_before = dWrap(Date(freq='W', year=2007, month=1, day=1))
            date_A_to_W_after = dWrap(Date(freq='W', year=2007, month=12, day=31))
            date_A_to_B_before = dWrap(Date(freq='B', year=2007, month=1, day=1))
            date_A_to_B_after = dWrap(Date(freq='B', year=2007, month=12, day=31))
            date_A_to_D_before = dWrap(Date(freq='D', year=2007, month=1, day=1))
            date_A_to_D_after = dWrap(Date(freq='D', year=2007, month=12, day=31))
            date_A_to_H_before = dWrap(Date(freq='H', year=2007, month=1, day=1,
                                      hour=0))
            date_A_to_H_after = dWrap(Date(freq='H', year=2007, month=12, day=31,
                                     hour=23))
            date_A_to_T_before = dWrap(Date(freq='T', year=2007, month=1, day=1,
                                      hour=0, minute=0))
            date_A_to_T_after = dWrap(Date(freq='T', year=2007, month=12, day=31,
                                     hour=23, minute=59))
            date_A_to_S_before = dWrap(Date(freq='S', year=2007, month=1, day=1,
                                      hour=0, minute=0, second=0))
            date_A_to_S_after = dWrap(Date(freq='S', year=2007, month=12, day=31,
                                     hour=23, minute=59, second=59))

            date_AJAN_to_D_after = dWrap(Date(freq='D', year=2007, month=1, day=31))
            date_AJAN_to_D_before = dWrap(Date(freq='D', year=2006, month=2, day=1))
            date_AJUN_to_D_after = dWrap(Date(freq='D', year=2007, month=6, day=30))
            date_AJUN_to_D_before = dWrap(Date(freq='D', year=2006, month=7, day=1))
            date_ANOV_to_D_after = dWrap(Date(freq='D', year=2007, month=11, day=30))
            date_ANOV_to_D_before = dWrap(Date(freq='D', year=2006, month=12, day=1))

            assert_func(date_A.asfreq('Q', "BEFORE"), date_A_to_Q_before)
            assert_func(date_A.asfreq('Q', "AFTER"), date_A_to_Q_after)
            assert_func(date_A.asfreq('M', "BEFORE"), date_A_to_M_before)
            assert_func(date_A.asfreq('M', "AFTER"), date_A_to_M_after)
            assert_func(date_A.asfreq('W', "BEFORE"), date_A_to_W_before)
            assert_func(date_A.asfreq('W', "AFTER"), date_A_to_W_after)
            assert_func(date_A.asfreq('B', "BEFORE"), date_A_to_B_before)
            assert_func(date_A.asfreq('B', "AFTER"), date_A_to_B_after)
            assert_func(date_A.asfreq('D', "BEFORE"), date_A_to_D_before)
            assert_func(date_A.asfreq('D', "AFTER"), date_A_to_D_after)
            assert_func(date_A.asfreq('H', "BEFORE"), date_A_to_H_before)
            assert_func(date_A.asfreq('H', "AFTER"), date_A_to_H_after)
            assert_func(date_A.asfreq('T', "BEFORE"), date_A_to_T_before)
            assert_func(date_A.asfreq('T', "AFTER"), date_A_to_T_after)
            assert_func(date_A.asfreq('S', "BEFORE"), date_A_to_S_before)
            assert_func(date_A.asfreq('S', "AFTER"), date_A_to_S_after)

            assert_func(date_AJAN.asfreq('D', "BEFORE"), date_AJAN_to_D_before)
            assert_func(date_AJAN.asfreq('D', "AFTER"), date_AJAN_to_D_after)

            assert_func(date_AJUN.asfreq('D', "BEFORE"), date_AJUN_to_D_before)
            assert_func(date_AJUN.asfreq('D', "AFTER"), date_AJUN_to_D_after)

            assert_func(date_ANOV.asfreq('D', "BEFORE"), date_ANOV_to_D_before)
            assert_func(date_ANOV.asfreq('D', "AFTER"), date_ANOV_to_D_after)

    def test_conv_quarterly(self):
        "frequency conversion tests: from Quarterly Frequency"

        for dWrap, assert_func in self.dateWrap:
            date_Q = dWrap(Date(freq='Q', year=2007, quarter=1))
            date_Q_end_of_year = dWrap(Date(freq='Q', year=2007, quarter=4))

            date_QEJAN = dWrap(Date(freq=C.FR_QTREJAN, year=2007, quarter=1))
            date_QEJUN = dWrap(Date(freq=C.FR_QTREJUN, year=2007, quarter=1))

            date_QSJAN = dWrap(Date(freq=C.FR_QTRSJAN, year=2007, quarter=1))
            date_QSJUN = dWrap(Date(freq=C.FR_QTRSJUN, year=2007, quarter=1))
            date_QSDEC = dWrap(Date(freq=C.FR_QTRSDEC, year=2007, quarter=1))

            date_Q_to_A = dWrap(Date(freq='A', year=2007))
            date_Q_to_M_before = dWrap(Date(freq='M', year=2007, month=1))
            date_Q_to_M_after = dWrap(Date(freq='M', year=2007, month=3))
            date_Q_to_W_before = dWrap(Date(freq='W', year=2007, month=1, day=1))
            date_Q_to_W_after = dWrap(Date(freq='W', year=2007, month=3, day=31))
            date_Q_to_B_before = dWrap(Date(freq='B', year=2007, month=1, day=1))
            date_Q_to_B_after = dWrap(Date(freq='B', year=2007, month=3, day=30))
            date_Q_to_D_before = dWrap(Date(freq='D', year=2007, month=1, day=1))
            date_Q_to_D_after = dWrap(Date(freq='D', year=2007, month=3, day=31))
            date_Q_to_H_before = dWrap(Date(freq='H', year=2007, month=1, day=1,
                                      hour=0))
            date_Q_to_H_after = dWrap(Date(freq='H', year=2007, month=3, day=31,
                                     hour=23))
            date_Q_to_T_before = dWrap(Date(freq='T', year=2007, month=1, day=1,
                                      hour=0, minute=0))
            date_Q_to_T_after = dWrap(Date(freq='T', year=2007, month=3, day=31,
                                     hour=23, minute=59))
            date_Q_to_S_before = dWrap(Date(freq='S', year=2007, month=1, day=1,
                                      hour=0, minute=0, second=0))
            date_Q_to_S_after = dWrap(Date(freq='S', year=2007, month=3, day=31,
                                     hour=23, minute=59, second=59))

            date_QEJAN_to_D_before = dWrap(Date(freq='D', year=2006, month=2, day=1))
            date_QEJAN_to_D_after = dWrap(Date(freq='D', year=2006, month=4, day=30))

            date_QEJUN_to_D_before = dWrap(Date(freq='D', year=2006, month=7, day=1))
            date_QEJUN_to_D_after = dWrap(Date(freq='D', year=2006, month=9, day=30))

            date_QSJAN_to_D_before = dWrap(Date(freq='D', year=2007, month=2, day=1))
            date_QSJAN_to_D_after = dWrap(Date(freq='D', year=2007, month=4, day=30))

            date_QSJUN_to_D_before = dWrap(Date(freq='D', year=2007, month=7, day=1))
            date_QSJUN_to_D_after = dWrap(Date(freq='D', year=2007, month=9, day=30))

            date_QSDEC_to_D_before = dWrap(Date(freq='D', year=2007, month=1, day=1))
            date_QSDEC_to_D_after = dWrap(Date(freq='D', year=2007, month=3, day=31))

            assert_func(date_Q.asfreq('A'), date_Q_to_A)
            assert_func(date_Q_end_of_year.asfreq('A'), date_Q_to_A)

            assert_func(date_Q.asfreq('M', "BEFORE"), date_Q_to_M_before)
            assert_func(date_Q.asfreq('M', "AFTER"), date_Q_to_M_after)
            assert_func(date_Q.asfreq('W', "BEFORE"), date_Q_to_W_before)
            assert_func(date_Q.asfreq('W', "AFTER"), date_Q_to_W_after)
            assert_func(date_Q.asfreq('B', "BEFORE"), date_Q_to_B_before)
            assert_func(date_Q.asfreq('B', "AFTER"), date_Q_to_B_after)
            assert_func(date_Q.asfreq('D', "BEFORE"), date_Q_to_D_before)
            assert_func(date_Q.asfreq('D', "AFTER"), date_Q_to_D_after)
            assert_func(date_Q.asfreq('H', "BEFORE"), date_Q_to_H_before)
            assert_func(date_Q.asfreq('H', "AFTER"), date_Q_to_H_after)
            assert_func(date_Q.asfreq('T', "BEFORE"), date_Q_to_T_before)
            assert_func(date_Q.asfreq('T', "AFTER"), date_Q_to_T_after)
            assert_func(date_Q.asfreq('S', "BEFORE"), date_Q_to_S_before)
            assert_func(date_Q.asfreq('S', "AFTER"), date_Q_to_S_after)

            assert_func(date_QEJAN.asfreq('D', "BEFORE"), date_QEJAN_to_D_before)
            assert_func(date_QEJAN.asfreq('D', "AFTER"), date_QEJAN_to_D_after)
            assert_func(date_QEJUN.asfreq('D', "BEFORE"), date_QEJUN_to_D_before)
            assert_func(date_QEJUN.asfreq('D', "AFTER"), date_QEJUN_to_D_after)

            assert_func(date_QSJAN.asfreq('D', "BEFORE"), date_QSJAN_to_D_before)
            assert_func(date_QSJAN.asfreq('D', "AFTER"), date_QSJAN_to_D_after)
            assert_func(date_QSJUN.asfreq('D', "BEFORE"), date_QSJUN_to_D_before)
            assert_func(date_QSJUN.asfreq('D', "AFTER"), date_QSJUN_to_D_after)
            assert_func(date_QSDEC.asfreq('D', "BEFORE"), date_QSDEC_to_D_before)
            assert_func(date_QSDEC.asfreq('D', "AFTER"), date_QSDEC_to_D_after)

    def test_conv_monthly(self):
        "frequency conversion tests: from Monthly Frequency"

        for dWrap, assert_func in self.dateWrap:
            date_M = dWrap(Date(freq='M', year=2007, month=1))
            date_M_end_of_year = dWrap(Date(freq='M', year=2007, month=12))
            date_M_end_of_quarter = dWrap(Date(freq='M', year=2007, month=3))
            date_M_to_A = dWrap(Date(freq='A', year=2007))
            date_M_to_Q = dWrap(Date(freq='Q', year=2007, quarter=1))
            date_M_to_W_before = dWrap(Date(freq='W', year=2007, month=1, day=1))
            date_M_to_W_after = dWrap(Date(freq='W', year=2007, month=1, day=31))
            date_M_to_B_before = dWrap(Date(freq='B', year=2007, month=1, day=1))
            date_M_to_B_after = dWrap(Date(freq='B', year=2007, month=1, day=31))
            date_M_to_D_before = dWrap(Date(freq='D', year=2007, month=1, day=1))
            date_M_to_D_after = dWrap(Date(freq='D', year=2007, month=1, day=31))
            date_M_to_H_before = dWrap(Date(freq='H', year=2007, month=1, day=1,
                                      hour=0))
            date_M_to_H_after = dWrap(Date(freq='H', year=2007, month=1, day=31,
                                     hour=23))
            date_M_to_T_before = dWrap(Date(freq='T', year=2007, month=1, day=1,
                                      hour=0, minute=0))
            date_M_to_T_after = dWrap(Date(freq='T', year=2007, month=1, day=31,
                                     hour=23, minute=59))
            date_M_to_S_before = dWrap(Date(freq='S', year=2007, month=1, day=1,
                                      hour=0, minute=0, second=0))
            date_M_to_S_after = dWrap(Date(freq='S', year=2007, month=1, day=31,
                                     hour=23, minute=59, second=59))

            assert_func(date_M.asfreq('A'), date_M_to_A)
            assert_func(date_M_end_of_year.asfreq('A'), date_M_to_A)
            assert_func(date_M.asfreq('Q'), date_M_to_Q)
            assert_func(date_M_end_of_quarter.asfreq('Q'), date_M_to_Q)

            assert_func(date_M.asfreq('W', "BEFORE"), date_M_to_W_before)
            assert_func(date_M.asfreq('W', "AFTER"), date_M_to_W_after)
            assert_func(date_M.asfreq('B', "BEFORE"), date_M_to_B_before)
            assert_func(date_M.asfreq('B', "AFTER"), date_M_to_B_after)
            assert_func(date_M.asfreq('D', "BEFORE"), date_M_to_D_before)
            assert_func(date_M.asfreq('D', "AFTER"), date_M_to_D_after)
            assert_func(date_M.asfreq('H', "BEFORE"), date_M_to_H_before)
            assert_func(date_M.asfreq('H', "AFTER"), date_M_to_H_after)
            assert_func(date_M.asfreq('T', "BEFORE"), date_M_to_T_before)
            assert_func(date_M.asfreq('T', "AFTER"), date_M_to_T_after)
            assert_func(date_M.asfreq('S', "BEFORE"), date_M_to_S_before)
            assert_func(date_M.asfreq('S', "AFTER"), date_M_to_S_after)


    def test_conv_weekly(self):
        "frequency conversion tests: from Weekly Frequency"

        for dWrap, assert_func in self.dateWrap:
            date_W = dWrap(Date(freq='W', year=2007, month=1, day=1))

            date_WSUN = dWrap(Date(freq='W-SUN', year=2007, month=1, day=7))
            date_WSAT = dWrap(Date(freq='W-SAT', year=2007, month=1, day=6))
            date_WFRI = dWrap(Date(freq='W-FRI', year=2007, month=1, day=5))
            date_WTHU = dWrap(Date(freq='W-THU', year=2007, month=1, day=4))
            date_WWED = dWrap(Date(freq='W-WED', year=2007, month=1, day=3))
            date_WTUE = dWrap(Date(freq='W-TUE', year=2007, month=1, day=2))
            date_WMON = dWrap(Date(freq='W-MON', year=2007, month=1, day=1))

            date_WSUN_to_D_before = dWrap(Date(freq='D', year=2007, month=1, day=1))
            date_WSUN_to_D_after = dWrap(Date(freq='D', year=2007, month=1, day=7))
            date_WSAT_to_D_before = dWrap(Date(freq='D', year=2006, month=12, day=31))
            date_WSAT_to_D_after = dWrap(Date(freq='D', year=2007, month=1, day=6))
            date_WFRI_to_D_before = dWrap(Date(freq='D', year=2006, month=12, day=30))
            date_WFRI_to_D_after = dWrap(Date(freq='D', year=2007, month=1, day=5))
            date_WTHU_to_D_before = dWrap(Date(freq='D', year=2006, month=12, day=29))
            date_WTHU_to_D_after = dWrap(Date(freq='D', year=2007, month=1, day=4))
            date_WWED_to_D_before = dWrap(Date(freq='D', year=2006, month=12, day=28))
            date_WWED_to_D_after = dWrap(Date(freq='D', year=2007, month=1, day=3))
            date_WTUE_to_D_before = dWrap(Date(freq='D', year=2006, month=12, day=27))
            date_WTUE_to_D_after = dWrap(Date(freq='D', year=2007, month=1, day=2))
            date_WMON_to_D_before = dWrap(Date(freq='D', year=2006, month=12, day=26))
            date_WMON_to_D_after = dWrap(Date(freq='D', year=2007, month=1, day=1))

            date_W_end_of_year = dWrap(Date(freq='W', year=2007, month=12, day=31))
            date_W_end_of_quarter = dWrap(Date(freq='W', year=2007, month=3, day=31))
            date_W_end_of_month = dWrap(Date(freq='W', year=2007, month=1, day=31))
            date_W_to_A = dWrap(Date(freq='A', year=2007))
            date_W_to_Q = dWrap(Date(freq='Q', year=2007, quarter=1))
            date_W_to_M = dWrap(Date(freq='M', year=2007, month=1))

            if Date(freq='D', year=2007, month=12, day=31).day_of_week == 6:
                date_W_to_A_end_of_year = dWrap(Date(freq='A', year=2007))
            else:
                date_W_to_A_end_of_year = dWrap(Date(freq='A', year=2008))

            if Date(freq='D', year=2007, month=3, day=31).day_of_week == 6:
                date_W_to_Q_end_of_quarter = dWrap(Date(freq='Q', year=2007, quarter=1))
            else:
                date_W_to_Q_end_of_quarter = dWrap(Date(freq='Q', year=2007, quarter=2))

            if Date(freq='D', year=2007, month=1, day=31).day_of_week == 6:
                date_W_to_M_end_of_month = dWrap(Date(freq='M', year=2007, month=1))
            else:
                date_W_to_M_end_of_month = dWrap(Date(freq='M', year=2007, month=2))

            date_W_to_B_before = dWrap(Date(freq='B', year=2007, month=1, day=1))
            date_W_to_B_after = dWrap(Date(freq='B', year=2007, month=1, day=5))
            date_W_to_D_before = dWrap(Date(freq='D', year=2007, month=1, day=1))
            date_W_to_D_after = dWrap(Date(freq='D', year=2007, month=1, day=7))
            date_W_to_H_before = dWrap(Date(freq='H', year=2007, month=1, day=1,
                                      hour=0))
            date_W_to_H_after = dWrap(Date(freq='H', year=2007, month=1, day=7,
                                     hour=23))
            date_W_to_T_before = dWrap(Date(freq='T', year=2007, month=1, day=1,
                                      hour=0, minute=0))
            date_W_to_T_after = dWrap(Date(freq='T', year=2007, month=1, day=7,
                                     hour=23, minute=59))
            date_W_to_S_before = dWrap(Date(freq='S', year=2007, month=1, day=1,
                                      hour=0, minute=0, second=0))
            date_W_to_S_after = dWrap(Date(freq='S', year=2007, month=1, day=7,
                                     hour=23, minute=59, second=59))

            assert_func(date_W.asfreq('A'), date_W_to_A)
            assert_func(date_W_end_of_year.asfreq('A'), date_W_to_A_end_of_year)
            assert_func(date_W.asfreq('Q'), date_W_to_Q)
            assert_func(date_W_end_of_quarter.asfreq('Q'), date_W_to_Q_end_of_quarter)
            assert_func(date_W.asfreq('M'), date_W_to_M)
            assert_func(date_W_end_of_month.asfreq('M'), date_W_to_M_end_of_month)

            assert_func(date_W.asfreq('B', "BEFORE"), date_W_to_B_before)
            assert_func(date_W.asfreq('B', "AFTER"), date_W_to_B_after)

            assert_func(date_W.asfreq('D', "BEFORE"), date_W_to_D_before)
            assert_func(date_W.asfreq('D', "AFTER"), date_W_to_D_after)

            assert_func(date_WSUN.asfreq('D', "BEFORE"), date_WSUN_to_D_before)
            assert_func(date_WSUN.asfreq('D', "AFTER"), date_WSUN_to_D_after)
            assert_func(date_WSAT.asfreq('D', "BEFORE"), date_WSAT_to_D_before)
            assert_func(date_WSAT.asfreq('D', "AFTER"), date_WSAT_to_D_after)
            assert_func(date_WFRI.asfreq('D', "BEFORE"), date_WFRI_to_D_before)
            assert_func(date_WFRI.asfreq('D', "AFTER"), date_WFRI_to_D_after)
            assert_func(date_WTHU.asfreq('D', "BEFORE"), date_WTHU_to_D_before)
            assert_func(date_WTHU.asfreq('D', "AFTER"), date_WTHU_to_D_after)
            assert_func(date_WWED.asfreq('D', "BEFORE"), date_WWED_to_D_before)
            assert_func(date_WWED.asfreq('D', "AFTER"), date_WWED_to_D_after)
            assert_func(date_WTUE.asfreq('D', "BEFORE"), date_WTUE_to_D_before)
            assert_func(date_WTUE.asfreq('D', "AFTER"), date_WTUE_to_D_after)
            assert_func(date_WMON.asfreq('D', "BEFORE"), date_WMON_to_D_before)
            assert_func(date_WMON.asfreq('D', "AFTER"), date_WMON_to_D_after)

            assert_func(date_W.asfreq('H', "BEFORE"), date_W_to_H_before)
            assert_func(date_W.asfreq('H', "AFTER"), date_W_to_H_after)
            assert_func(date_W.asfreq('T', "BEFORE"), date_W_to_T_before)
            assert_func(date_W.asfreq('T', "AFTER"), date_W_to_T_after)
            assert_func(date_W.asfreq('S', "BEFORE"), date_W_to_S_before)
            assert_func(date_W.asfreq('S', "AFTER"), date_W_to_S_after)

    def test_conv_business(self):
        "frequency conversion tests: from Business Frequency"

        for dWrap, assert_func in self.dateWrap:
            date_B = dWrap(Date(freq='B', year=2007, month=1, day=1))
            date_B_end_of_year = dWrap(Date(freq='B', year=2007, month=12, day=31))
            date_B_end_of_quarter = dWrap(Date(freq='B', year=2007, month=3, day=30))
            date_B_end_of_month = dWrap(Date(freq='B', year=2007, month=1, day=31))
            date_B_end_of_week = dWrap(Date(freq='B', year=2007, month=1, day=5))

            date_B_to_A = dWrap(Date(freq='A', year=2007))
            date_B_to_Q = dWrap(Date(freq='Q', year=2007, quarter=1))
            date_B_to_M = dWrap(Date(freq='M', year=2007, month=1))
            date_B_to_W = dWrap(Date(freq='W', year=2007, month=1, day=7))
            date_B_to_D = dWrap(Date(freq='D', year=2007, month=1, day=1))
            date_B_to_H_before = dWrap(Date(freq='H', year=2007, month=1, day=1,
                                      hour=0))
            date_B_to_H_after = dWrap(Date(freq='H', year=2007, month=1, day=1,
                                     hour=23))
            date_B_to_T_before = dWrap(Date(freq='T', year=2007, month=1, day=1,
                                      hour=0, minute=0))
            date_B_to_T_after = dWrap(Date(freq='T', year=2007, month=1, day=1,
                                     hour=23, minute=59))
            date_B_to_S_before = dWrap(Date(freq='S', year=2007, month=1, day=1,
                                      hour=0, minute=0, second=0))
            date_B_to_S_after = dWrap(Date(freq='S', year=2007, month=1, day=1,
                                     hour=23, minute=59, second=59))

            assert_func(date_B.asfreq('A'), date_B_to_A)
            assert_func(date_B_end_of_year.asfreq('A'), date_B_to_A)
            assert_func(date_B.asfreq('Q'), date_B_to_Q)
            assert_func(date_B_end_of_quarter.asfreq('Q'), date_B_to_Q)
            assert_func(date_B.asfreq('M'), date_B_to_M)
            assert_func(date_B_end_of_month.asfreq('M'), date_B_to_M)
            assert_func(date_B.asfreq('W'), date_B_to_W)
            assert_func(date_B_end_of_week.asfreq('W'), date_B_to_W)

            assert_func(date_B.asfreq('D'), date_B_to_D)

            assert_func(date_B.asfreq('H', "BEFORE"), date_B_to_H_before)
            assert_func(date_B.asfreq('H', "AFTER"), date_B_to_H_after)
            assert_func(date_B.asfreq('T', "BEFORE"), date_B_to_T_before)
            assert_func(date_B.asfreq('T', "AFTER"), date_B_to_T_after)
            assert_func(date_B.asfreq('S', "BEFORE"), date_B_to_S_before)
            assert_func(date_B.asfreq('S', "AFTER"), date_B_to_S_after)

    def test_conv_daily(self):
        "frequency conversion tests: from Business Frequency"

        for dWrap, assert_func in self.dateWrap:
            date_D = dWrap(Date(freq='D', year=2007, month=1, day=1))
            date_D_end_of_year = dWrap(Date(freq='D', year=2007, month=12, day=31))
            date_D_end_of_quarter = dWrap(Date(freq='D', year=2007, month=3, day=31))
            date_D_end_of_month = dWrap(Date(freq='D', year=2007, month=1, day=31))
            date_D_end_of_week = dWrap(Date(freq='D', year=2007, month=1, day=7))

            date_D_friday = dWrap(Date(freq='D', year=2007, month=1, day=5))
            date_D_saturday = dWrap(Date(freq='D', year=2007, month=1, day=6))
            date_D_sunday = dWrap(Date(freq='D', year=2007, month=1, day=7))
            date_D_monday = dWrap(Date(freq='D', year=2007, month=1, day=8))

            date_B_friday = dWrap(Date(freq='B', year=2007, month=1, day=5))
            date_B_monday = dWrap(Date(freq='B', year=2007, month=1, day=8))

            date_D_to_A = dWrap(Date(freq='A', year=2007))

            date_Deoq_to_AJAN = dWrap(Date(freq='A-JAN', year=2008))
            date_Deoq_to_AJUN = dWrap(Date(freq='A-JUN', year=2007))
            date_Deoq_to_ADEC = dWrap(Date(freq='A-DEC', year=2007))

            date_D_to_QEJAN = dWrap(Date(freq=C.FR_QTREJAN, year=2007, quarter=4))
            date_D_to_QEJUN = dWrap(Date(freq=C.FR_QTREJUN, year=2007, quarter=3))
            date_D_to_QEDEC = dWrap(Date(freq=C.FR_QTREDEC, year=2007, quarter=1))

            date_D_to_QSJAN = dWrap(Date(freq=C.FR_QTRSJAN, year=2006, quarter=4))
            date_D_to_QSJUN = dWrap(Date(freq=C.FR_QTRSJUN, year=2006, quarter=3))
            date_D_to_QSDEC = dWrap(Date(freq=C.FR_QTRSDEC, year=2007, quarter=1))

            date_D_to_M = dWrap(Date(freq='M', year=2007, month=1))
            date_D_to_W = dWrap(Date(freq='W', year=2007, month=1, day=7))

            date_D_to_H_before = dWrap(Date(freq='H', year=2007, month=1, day=1,
                                      hour=0))
            date_D_to_H_after = dWrap(Date(freq='H', year=2007, month=1, day=1,
                                     hour=23))
            date_D_to_T_before = dWrap(Date(freq='T', year=2007, month=1, day=1,
                                      hour=0, minute=0))
            date_D_to_T_after = dWrap(Date(freq='T', year=2007, month=1, day=1,
                                     hour=23, minute=59))
            date_D_to_S_before = dWrap(Date(freq='S', year=2007, month=1, day=1,
                                      hour=0, minute=0, second=0))
            date_D_to_S_after = dWrap(Date(freq='S', year=2007, month=1, day=1,
                                     hour=23, minute=59, second=59))

            assert_func(date_D.asfreq('A'), date_D_to_A)

            assert_func(date_D_end_of_quarter.asfreq('A-JAN'), date_Deoq_to_AJAN)
            assert_func(date_D_end_of_quarter.asfreq('A-JUN'), date_Deoq_to_AJUN)
            assert_func(date_D_end_of_quarter.asfreq('A-DEC'), date_Deoq_to_ADEC)

            assert_func(date_D_end_of_year.asfreq('A'), date_D_to_A)
            assert_func(date_D_end_of_quarter.asfreq('Q'), date_D_to_QEDEC)
            assert_func(date_D.asfreq(C.FR_QTREJAN), date_D_to_QEJAN)
            assert_func(date_D.asfreq(C.FR_QTREJUN), date_D_to_QEJUN)
            assert_func(date_D.asfreq(C.FR_QTREDEC), date_D_to_QEDEC)
            assert_func(date_D.asfreq(C.FR_QTRSJAN), date_D_to_QSJAN)
            assert_func(date_D.asfreq(C.FR_QTRSJUN), date_D_to_QSJUN)
            assert_func(date_D.asfreq(C.FR_QTRSDEC), date_D_to_QSDEC)
            assert_func(date_D.asfreq('M'), date_D_to_M)
            assert_func(date_D_end_of_month.asfreq('M'), date_D_to_M)
            assert_func(date_D.asfreq('W'), date_D_to_W)
            assert_func(date_D_end_of_week.asfreq('W'), date_D_to_W)

            assert_func(date_D_friday.asfreq('B'), date_B_friday)
            assert_func(date_D_saturday.asfreq('B', "BEFORE"), date_B_friday)
            assert_func(date_D_saturday.asfreq('B', "AFTER"), date_B_monday)
            assert_func(date_D_sunday.asfreq('B', "BEFORE"), date_B_friday)
            assert_func(date_D_sunday.asfreq('B', "AFTER"), date_B_monday)

            assert_func(date_D.asfreq('H', "BEFORE"), date_D_to_H_before)
            assert_func(date_D.asfreq('H', "AFTER"), date_D_to_H_after)
            assert_func(date_D.asfreq('T', "BEFORE"), date_D_to_T_before)
            assert_func(date_D.asfreq('T', "AFTER"), date_D_to_T_after)
            assert_func(date_D.asfreq('S', "BEFORE"), date_D_to_S_before)
            assert_func(date_D.asfreq('S', "AFTER"), date_D_to_S_after)

    def test_conv_hourly(self):
        "frequency conversion tests: from Hourly Frequency"

        for dWrap, assert_func in self.dateWrap:
            date_H = dWrap(Date(freq='H', year=2007, month=1, day=1, hour=0))
            date_H_end_of_year = dWrap(Date(freq='H', year=2007, month=12, day=31,
                                      hour=23))
            date_H_end_of_quarter = dWrap(Date(freq='H', year=2007, month=3, day=31,
                                         hour=23))
            date_H_end_of_month = dWrap(Date(freq='H', year=2007, month=1, day=31,
                                       hour=23))
            date_H_end_of_week = dWrap(Date(freq='H', year=2007, month=1, day=7,
                                      hour=23))
            date_H_end_of_day = dWrap(Date(freq='H', year=2007, month=1, day=1,
                                     hour=23))
            date_H_end_of_bus = dWrap(Date(freq='H', year=2007, month=1, day=1,
                                     hour=23))

            date_H_to_A = dWrap(Date(freq='A', year=2007))
            date_H_to_Q = dWrap(Date(freq='Q', year=2007, quarter=1))
            date_H_to_M = dWrap(Date(freq='M', year=2007, month=1))
            date_H_to_W = dWrap(Date(freq='W', year=2007, month=1, day=7))
            date_H_to_D = dWrap(Date(freq='D', year=2007, month=1, day=1))
            date_H_to_B = dWrap(Date(freq='B', year=2007, month=1, day=1))

            date_H_to_T_before = dWrap(Date(freq='T', year=2007, month=1, day=1,
                                      hour=0, minute=0))
            date_H_to_T_after = dWrap(Date(freq='T', year=2007, month=1, day=1,
                                     hour=0, minute=59))
            date_H_to_S_before = dWrap(Date(freq='S', year=2007, month=1, day=1,
                                      hour=0, minute=0, second=0))
            date_H_to_S_after = dWrap(Date(freq='S', year=2007, month=1, day=1,
                                     hour=0, minute=59, second=59))

            assert_func(date_H.asfreq('A'), date_H_to_A)
            assert_func(date_H_end_of_year.asfreq('A'), date_H_to_A)
            assert_func(date_H.asfreq('Q'), date_H_to_Q)
            assert_func(date_H_end_of_quarter.asfreq('Q'), date_H_to_Q)
            assert_func(date_H.asfreq('M'), date_H_to_M)
            assert_func(date_H_end_of_month.asfreq('M'), date_H_to_M)
            assert_func(date_H.asfreq('W'), date_H_to_W)
            assert_func(date_H_end_of_week.asfreq('W'), date_H_to_W)
            assert_func(date_H.asfreq('D'), date_H_to_D)
            assert_func(date_H_end_of_day.asfreq('D'), date_H_to_D)
            assert_func(date_H.asfreq('B'), date_H_to_B)
            assert_func(date_H_end_of_bus.asfreq('B'), date_H_to_B)

            assert_func(date_H.asfreq('T', "BEFORE"), date_H_to_T_before)
            assert_func(date_H.asfreq('T', "AFTER"), date_H_to_T_after)
            assert_func(date_H.asfreq('S', "BEFORE"), date_H_to_S_before)
            assert_func(date_H.asfreq('S', "AFTER"), date_H_to_S_after)

    def test_conv_minutely(self):
        "frequency conversion tests: from Minutely Frequency"

        for dWrap, assert_func in self.dateWrap:
            date_T = dWrap(Date(freq='T', year=2007, month=1, day=1,
                          hour=0, minute=0))
            date_T_end_of_year = dWrap(Date(freq='T', year=2007, month=12, day=31,
                                      hour=23, minute=59))
            date_T_end_of_quarter = dWrap(Date(freq='T', year=2007, month=3, day=31,
                                         hour=23, minute=59))
            date_T_end_of_month = dWrap(Date(freq='T', year=2007, month=1, day=31,
                                       hour=23, minute=59))
            date_T_end_of_week = dWrap(Date(freq='T', year=2007, month=1, day=7,
                                      hour=23, minute=59))
            date_T_end_of_day = dWrap(Date(freq='T', year=2007, month=1, day=1,
                                     hour=23, minute=59))
            date_T_end_of_bus = dWrap(Date(freq='T', year=2007, month=1, day=1,
                                     hour=23, minute=59))
            date_T_end_of_hour = dWrap(Date(freq='T', year=2007, month=1, day=1,
                                      hour=0, minute=59))

            date_T_to_A = dWrap(Date(freq='A', year=2007))
            date_T_to_Q = dWrap(Date(freq='Q', year=2007, quarter=1))
            date_T_to_M = dWrap(Date(freq='M', year=2007, month=1))
            date_T_to_W = dWrap(Date(freq='W', year=2007, month=1, day=7))
            date_T_to_D = dWrap(Date(freq='D', year=2007, month=1, day=1))
            date_T_to_B = dWrap(Date(freq='B', year=2007, month=1, day=1))
            date_T_to_H = dWrap(Date(freq='H', year=2007, month=1, day=1, hour=0))

            date_T_to_S_before = dWrap(Date(freq='S', year=2007, month=1, day=1,
                                      hour=0, minute=0, second=0))
            date_T_to_S_after = dWrap(Date(freq='S', year=2007, month=1, day=1,
                                     hour=0, minute=0, second=59))

            assert_func(date_T.asfreq('A'), date_T_to_A)
            assert_func(date_T_end_of_year.asfreq('A'), date_T_to_A)
            assert_func(date_T.asfreq('Q'), date_T_to_Q)
            assert_func(date_T_end_of_quarter.asfreq('Q'), date_T_to_Q)
            assert_func(date_T.asfreq('M'), date_T_to_M)
            assert_func(date_T_end_of_month.asfreq('M'), date_T_to_M)
            assert_func(date_T.asfreq('W'), date_T_to_W)
            assert_func(date_T_end_of_week.asfreq('W'), date_T_to_W)
            assert_func(date_T.asfreq('D'), date_T_to_D)
            assert_func(date_T_end_of_day.asfreq('D'), date_T_to_D)
            assert_func(date_T.asfreq('B'), date_T_to_B)
            assert_func(date_T_end_of_bus.asfreq('B'), date_T_to_B)
            assert_func(date_T.asfreq('H'), date_T_to_H)
            assert_func(date_T_end_of_hour.asfreq('H'), date_T_to_H)

            assert_func(date_T.asfreq('S', "BEFORE"), date_T_to_S_before)
            assert_func(date_T.asfreq('S', "AFTER"), date_T_to_S_after)


    def test_conv_secondly(self):
        "frequency conversion tests: from Secondly Frequency"

        for dWrap, assert_func in self.dateWrap:
            date_S = dWrap(Date(freq='S', year=2007, month=1, day=1,
                          hour=0, minute=0, second=0))
            date_S_end_of_year = dWrap(Date(freq='S', year=2007, month=12, day=31,
                                      hour=23, minute=59, second=59))
            date_S_end_of_quarter = dWrap(Date(freq='S', year=2007, month=3, day=31,
                                         hour=23, minute=59, second=59))
            date_S_end_of_month = dWrap(Date(freq='S', year=2007, month=1, day=31,
                                       hour=23, minute=59, second=59))
            date_S_end_of_week = dWrap(Date(freq='S', year=2007, month=1, day=7,
                                      hour=23, minute=59, second=59))
            date_S_end_of_day = dWrap(Date(freq='S', year=2007, month=1, day=1,
                                     hour=23, minute=59, second=59))
            date_S_end_of_bus = dWrap(Date(freq='S', year=2007, month=1, day=1,
                                     hour=23, minute=59, second=59))
            date_S_end_of_hour = dWrap(Date(freq='S', year=2007, month=1, day=1,
                                      hour=0, minute=59, second=59))
            date_S_end_of_minute = dWrap(Date(freq='S', year=2007, month=1, day=1,
                                        hour=0, minute=0, second=59))

            date_S_to_A = dWrap(Date(freq='A', year=2007))
            date_S_to_Q = dWrap(Date(freq='Q', year=2007, quarter=1))
            date_S_to_M = dWrap(Date(freq='M', year=2007, month=1))
            date_S_to_W = dWrap(Date(freq='W', year=2007, month=1, day=7))
            date_S_to_D = dWrap(Date(freq='D', year=2007, month=1, day=1))
            date_S_to_B = dWrap(Date(freq='B', year=2007, month=1, day=1))
            date_S_to_H = dWrap(Date(freq='H', year=2007, month=1, day=1,
                               hour=0))
            date_S_to_T = dWrap(Date(freq='T', year=2007, month=1, day=1,
                               hour=0, minute=0))

            assert_func(date_S.asfreq('A'), date_S_to_A)
            assert_func(date_S_end_of_year.asfreq('A'), date_S_to_A)
            assert_func(date_S.asfreq('Q'), date_S_to_Q)
            assert_func(date_S_end_of_quarter.asfreq('Q'), date_S_to_Q)
            assert_func(date_S.asfreq('M'), date_S_to_M)
            assert_func(date_S_end_of_month.asfreq('M'), date_S_to_M)
            assert_func(date_S.asfreq('W'), date_S_to_W)
            assert_func(date_S_end_of_week.asfreq('W'), date_S_to_W)
            assert_func(date_S.asfreq('D'), date_S_to_D)
            assert_func(date_S_end_of_day.asfreq('D'), date_S_to_D)
            assert_func(date_S.asfreq('B'), date_S_to_B)
            assert_func(date_S_end_of_bus.asfreq('B'), date_S_to_B)
            assert_func(date_S.asfreq('H'), date_S_to_H)
            assert_func(date_S_end_of_hour.asfreq('H'), date_S_to_H)
            assert_func(date_S.asfreq('T'), date_S_to_T)
            assert_func(date_S_end_of_minute.asfreq('T'), date_S_to_T)


class TestMethods(NumpyTestCase):
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
        #CHECK : Oops, what were we supposed to do here ?


    def test_getsteps(self):
        "Tests the getsteps method"
        dlist = ['2007-01-%02i' %i for i in (1,2,3,4,8,9,10,11,12,15)]
        ddates = date_array_fromlist(dlist)
        assert_equal(ddates.get_steps(), [1,1,1,4,1,1,1,1,3])


    def test_empty_datearray(self):
        empty_darray = DateArray([], freq='b')
        assert_equal(empty_darray.isfull(), True)
        assert_equal(empty_darray.isvalid(), True)
        assert_equal(empty_darray.get_steps(), None)

    def test_cachedinfo(self):
        D = date_array(start_date=thisday('D'), length=5)
        Dstr = D.tostring()
        assert_equal(D.tostring(), Dstr)
        DL = D[[0,-1]]
        assert_equal(DL.tostring(), Dstr[[0,-1]])

###############################################################################
#------------------------------------------------------------------------------
if __name__ == "__main__":
    NumpyTest().run()
