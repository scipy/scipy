# pylint: disable-msg=W0611, W0612, W0511,R0201
"""Tests suite for MaskedArray & subclassing.

:author: Pierre Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu & mattknox_ca_at_hotmail_dot_com
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import numpy as N
import numpy.core.numeric as numeric

from numpy.testing import NumpyTest, NumpyTestCase

import maskedarray.testutils
from maskedarray.testutils import *

import maskedarray.core as coremodule
from maskedarray.core import MaskedArray, masked

import tseries
from tseries import time_series, thisday

import addons.filters
from addons.filters import running_mean


class test_runningmean(NumpyTestCase):
    
    def __init__(self, *args, **kwds):
        NumpyTestCase.__init__(self, *args, **kwds)
        self.data = numeric.arange(25)
        self.maskeddata = MaskedArray(self.data)
        self.maskeddata[10] = masked
    #         
    def test_onregulararray(self):
        data = self.data
        for width in [3,5,7]:
            k = (width-1)/2
            ravg = running_mean(data,width)
            assert(isinstance(ravg, MaskedArray))
            assert_equal(ravg, data)
            assert_equal(ravg._mask, [1]*k+[0]*(len(data)-2*k)+[1]*k)
    #
    def test_onmaskedarray(self):
        data = self.maskeddata
        for width in [3,5,7]:
            k = (width-1)/2
            ravg = running_mean(data,width)
            assert(isinstance(ravg, MaskedArray))
            assert_equal(ravg, data)
            m = N.zeros(len(data), N.bool_)
            m[:k] = m[-k:] = m[10-k:10+k+1] = True
            assert_equal(ravg._mask, m)
    #
    def test_ontimeseries(self):
        data = time_series(self.maskeddata, start_date=thisday('D'))
        for width in [3,5,7]:
            k = (width-1)/2
            ravg = running_mean(data,width)
            assert(isinstance(ravg, MaskedArray))
            assert_equal(ravg, data)
            m = N.zeros(len(data), N.bool_)
            m[:k] = m[-k:] = m[10-k:10+k+1] = True
            assert_equal(ravg._mask, m)
            assert_equal(ravg._dates, data._dates)
    #
    def tests_onmultitimeseries(self):
        maskeddata = MaskedArray(N.random.random(75).reshape(25,3))
        maskeddata[10] = masked
        data = time_series(maskeddata, start_date=thisday('D'))
        for width in [3,5,7]:
            k = (width-1)/2
            ravg = running_mean(data,width)
            assert(isinstance(ravg, MaskedArray))
            assert_almost_equal(ravg[18].squeeze(), data[18-k:18+k+1].mean(0))
            m = N.zeros(data.shape, N.bool_)
            m[:k] = m[-k:] = m[10-k:10+k+1] = True
            assert_equal(ravg._mask, m)
            assert_equal(ravg._dates, data._dates)

#------------------------------------------------------------------------------
if __name__ == "__main__":
    NumpyTest().run()                