# pylint: disable-msg=W0611, W0612, W0511,R0201
"""Tests suite for MaskedArray & subclassing.

:author: Pierre Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu & mattknox_ca_at_hotmail_dot_com
:version: $Id: test_filters.py 2819 2007-03-03 23:00:20Z pierregm $
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author: pierregm $)"
__version__ = '1.0'
__revision__ = "$Revision: 2819 $"
__date__     = '$Date: 2007-03-03 18:00:20 -0500 (Sat, 03 Mar 2007) $'

import numpy as N
import numpy.core.numeric as numeric

from numpy.testing import NumpyTest, NumpyTestCase

import maskedarray.testutils
from maskedarray.testutils import *

import maskedarray as MA
import maskedarray.core as coremodule
from maskedarray.core import MaskedArray, masked
from maskedarray import mstats

import timeseries as TS
from timeseries import time_series, thisday

from timeseries.lib import moving_funcs as MF

class TestCMovAverage(NumpyTestCase):

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
            ravg = MF.cmov_average(data,width)
            assert(isinstance(ravg, MaskedArray))
            assert_equal(ravg, data)
            assert_equal(ravg._mask, [1]*k+[0]*(len(data)-2*k)+[1]*k)
    #
    def test_onmaskedarray(self):
        data = self.maskeddata
        for width in [3,5,7]:
            k = (width-1)/2
            ravg = MF.cmov_average(data,width)
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
            ravg = MF.cmov_average(data,width)
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
            ravg = MF.cmov_average(data,width)
            assert(isinstance(ravg, MaskedArray))
            assert_almost_equal(ravg[18].squeeze(), data[18-k:18+k+1].mean(0))
            m = N.zeros(data.shape, N.bool_)
            m[:k] = m[-k:] = m[10-k:10+k+1] = True
            assert_equal(ravg._mask, m)
            assert_equal(ravg._dates, data._dates)



class TestMovFuncs(NumpyTestCase):

    def __init__(self, *args, **kwds):
        NumpyTestCase.__init__(self, *args, **kwds)
        self.data = numeric.arange(25)
        self.maskeddata = MaskedArray(self.data)
        self.maskeddata[10] = masked
        self.func_pairs = [
            (MF.mov_average, MA.mean),
            (MF.mov_median, mstats.mmedian),
            ((lambda x, span : MF.mov_stddev(x, span, bias=True)), MA.std)]
    #
    def test_onregulararray(self):
        data = self.data
        for Mfunc, Nfunc in self.func_pairs:
            for k in [3,4,5]:
                result = Mfunc(data, k)
                assert(isinstance(result, MaskedArray))
                for x in range(len(data)-k+1):
                    assert_almost_equal(result[x+k-1], Nfunc(data[x:x+k]))
                assert_equal(result._mask, [1]*(k-1)+[0]*(len(data)-k+1))

    #
    def test_onmaskedarray(self):
        data = self.maskeddata

        for Mfunc, Nfunc in self.func_pairs:
            for k in [3,4,5]:
                result = Mfunc(data, k)
                assert(isinstance(result, MaskedArray))
                for x in range(len(data)-k+1):
                    if result[x+k-1] is not MA.masked:
                        assert_almost_equal(result[x+k-1], Nfunc(data[x:x+k]))
                result_mask = N.array([1]*(k-1)+[0]*(len(data)-k+1))
                result_mask[10:10+k] = 1
                assert_equal(result._mask, result_mask)

    #
    def test_ontimeseries(self):

        data = time_series(self.maskeddata, start_date=thisday('D'))

        for Mfunc, Nfunc in self.func_pairs:
            for k in [3,4,5]:
                result = Mfunc(data, k)
                assert(isinstance(result, MaskedArray))
                for x in range(len(data)-k+1):
                    if result[x+k-1] is not TS.tsmasked:
                        assert_almost_equal(
                                N.asarray(result[x+k-1]),
                                N.asarray(Nfunc(data[x:x+k])))
                result_mask = N.array([1]*(k-1)+[0]*(len(data)-k+1))
                result_mask[10:10+k] = 1
                assert_equal(result._mask, result_mask)
                assert_equal(result._dates, data._dates)

    def test_covar(self):
        # test that covariance of series with itself is equal to variance
        data = self.maskeddata
        for bias in [True, False]:
            covar = MF.mov_covar(data, data, 3, bias=bias)
            var = MF.mov_var(data, 3, bias=bias)
            assert_equal(covar, var)

#------------------------------------------------------------------------------
if __name__ == "__main__":
    NumpyTest().run()
