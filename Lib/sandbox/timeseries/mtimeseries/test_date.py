# pylint: disable-msg=W0611, W0612, W0511,R0201
"""Tests suite for MaskedArray.
Adapted from the original test_ma by Pierre Gerard-Marchant

:author: Pierre Gerard-Marchant
:contact: pierregm_at_uga_dot_edu
:version: $Id: test_core.py 59 2006-12-22 23:58:11Z backtopop $
"""
__author__ = "Pierre GF Gerard-Marchant ($Author: backtopop $)"
__version__ = '1.0'
__revision__ = "$Revision: 59 $"
__date__     = '$Date: 2006-12-22 18:58:11 -0500 (Fri, 22 Dec 2006) $'

import types
import datetime

import numpy as N
import numpy.core.fromnumeric  as fromnumeric
import numpy.core.numeric as numeric
from numpy.testing import NumpyTest, NumpyTestCase
from numpy.testing.utils import build_err_msg

import maskedarray
from maskedarray import masked_array

import maskedarray.testutils
reload(maskedarray.testutils)
from maskedarray.testutils import assert_equal, assert_array_equal

import tsdate
reload(tsdate)
from tsdate import *

class test_creation(NumpyTestCase):
    "Base test class for MaskedArrays."
    
    def __init__(self, *args, **kwds):
        NumpyTestCase.__init__(self, *args, **kwds)
    
    def test_fromstrings(self):
        "Tests creation from list of strings"
        dlist = ['2007-01-%02i' % i for i in range(1,15)]
        # A simple case: daily data
        dates = datearray_fromlist(dlist, 'D')
        assert_equal(dates.freq,'D')
        assert(dates.isfull())
        assert(not dates.ispacked())
        assert_equal(dates, 732677+numpy.arange(len(dlist)))
        # as simple, but we need to guess the frequency this time
        dates = datearray_fromlist(dlist, 'D')
        assert_equal(dates.freq,'D')
        assert(dates.isfull())
        assert(not dates.ispacked())
        assert_equal(dates, 732677+numpy.arange(len(dlist)))
        # Still daily data, that we force to month
        dates = datearray_fromlist(dlist, 'M')
        assert_equal(dates.freq,'M')
        assert(not dates.isfull())
        assert(dates.ispacked())
        assert_equal(dates, [24085]*len(dlist))
        # Now, for monthly data
        dlist = ['2007-%02i' % i for i in range(1,13)]
        dates = datearray_fromlist(dlist, 'M')
        assert_equal(dates.freq,'M')
        assert(dates.isfull())
        assert(not dates.ispacked())
        assert_equal(dates, 24085 + numpy.arange(12))
        # Monthly data  w/ guessing
        dlist = ['2007-%02i' % i for i in range(1,13)]
        dates = datearray_fromlist(dlist, )
        assert_equal(dates.freq,'M')
        assert(dates.isfull())
        assert(not dates.ispacked())
        assert_equal(dates, 24085 + numpy.arange(12))
        
    def test_fromstrings_wmissing(self):
        "Tests creation from list of strings w/ missing dates"
        dlist = ['2007-01-%02i' % i for i in (1,2,4,5,7,8,10,11,13)]
        dates = datearray_fromlist(dlist)
        assert_equal(dates.freq,'U')
        assert(not dates.isfull())
        assert(not dates.ispacked())
        assert_equal(dates.tovalue(),732676+numpy.array([1,2,4,5,7,8,10,11,13]))
        #
        ddates = datearray_fromlist(dlist, 'D')
        assert_equal(ddates.freq,'D')
        assert(not ddates.isfull())
        assert(not ddates.ispacked())
        #
        mdates = datearray_fromlist(dlist, 'M')
        assert_equal(mdates.freq,'M')
        assert(not dates.isfull())
        assert(mdates.ispacked())
        #
    
    def test_fromsobjects(self):
        "Tests creation from list of objects."
        dlist = ['2007-01-%02i' % i for i in (1,2,4,5,7,8,10,11,13)]
        dates = datearray_fromlist(dlist)
        dobj = [datetime.datetime.fromordinal(d) for d in dates.toordinal()]
        odates = datearray_fromlist(dobj)
        assert_equal(dates,odates)
        dobj = [mxDFromString(d) for d in dlist]
        odates = datearray_fromlist(dobj)
        assert_equal(dates,odates)


class test_methods(NumpyTestCase):
    "Base test class for MaskedArrays."
    
    def __init__(self, *args, **kwds):
        NumpyTestCase.__init__(self, *args, **kwds)       
        
    def test_getitem(self):
        "Tests getitem"
        dlist = ['2007-%02i' % i for i in range(1,5)+range(7,13)]
        mdates = datearray_fromlist(dlist).asfreq('M')
        # Using an integer
        assert_equal(mdates[0].value, 24085)
        assert_equal(mdates[-1].value, 24096)
        # Using a date
        lag = mdates.find_dates(mdates[0])
        assert_equal(mdates[lag], mdates[0])
        lag = mdates.find_dates(Date('M',value=24092))
        assert_equal(mdates[lag], mdates[5])
        # Using several dates
        lag = mdates.find_dates(Date('M',value=24085), Date('M',value=24096))
        assert_equal(mdates[lag], 
                     DateArray([mdates[0], mdates[-1]], freq='M'))
        assert_equal(mdates[[mdates[0],mdates[-1]]], mdates[lag])
        #
        assert_equal(mdates>=mdates[-4], [0,0,0,0,0,0,1,1,1,1])
        dlist = ['2006-%02i' % i for i in range(1,5)+range(7,13)]
        mdates = datearray_fromlist(dlist).asfreq('M')

        
    def test_getsteps(self):
        dlist = ['2007-01-%02i' %i for i in (1,2,3,4,8,9,10,11,12,15)]
        ddates = datearray_fromlist(dlist)
        assert_equal(ddates.get_steps(), [1,1,1,4,1,1,1,1,3])

###############################################################################
#------------------------------------------------------------------------------
if __name__ == "__main__":
    NumpyTest().run()        