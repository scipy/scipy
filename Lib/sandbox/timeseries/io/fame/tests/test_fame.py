"""Tests suite for fame io submodule.

:author: Matt Knox
:contact: mattknox_ca_at_hotmail_dot_com
:version: $Id: test_fame.py 2578 2007-01-17 19:25:10Z mattknox_ca $
"""
__author__ = "Matt Knox ($Author: mattknox_ca $)"
__version__ = '1.0'
__revision__ = "$Revision: 2578 $"
__date__     = '$Date: 2007-01-17 14:25:10 -0500 (Wed, 17 Jan 2007) $'

import numpy as N
from numpy import bool_, complex_, float_, int_, object_
import numpy.core.fromnumeric  as fromnumeric
import numpy.core.numeric as numeric
from numpy.testing import NumpyTest, NumpyTestCase
from numpy.testing.utils import build_err_msg

from timeseries.io import fame
from timeseries import Report
import timeseries as ts
import maskedarray as ma
import numpy as np

import maskedarray
from maskedarray import masked_array, masked, nomask

import maskedarray.testutils
from maskedarray.testutils import assert_equal, assert_array_equal, approx, assert_mask_equal

# setup all the data to be used for reading and writing
data = {'dates':{}, 'darrays':{}, 'freqs':{}, 'cser':{}, 'tser':{}, 'scalars':{}}

data['dates']['a'] = ts.Date(freq='A', year=2004)
data['dates']['q'] = ts.Date(freq='Q', year=2004, quarter=1)
data['dates']['m'] = ts.Date(freq='M', year=2004, month=1)
data['dates']['w'] = ts.Date(freq='W', year=2004, month=1, day=1)
data['dates']['b'] = ts.Date(freq='b', year=2004, month=1, day=1)
data['dates']['d'] = ts.Date(freq='d', year=2004, month=1, day=1)
data['dates']['h'] = ts.Date(freq='h', year=2004, month=1, day=1, hour=0)
data['dates']['t'] = ts.Date(freq='t', year=2004, month=1, day=1, hour=0, minute=0)
data['dates']['s'] = ts.Date(freq='s', year=2004, month=1, day=1, hour=0, minute=0, second=0)

for freq in data['dates']:
    data['darrays'][freq] = ts.date_array(start_date=data['dates'][freq], length=10)
    data['cser']['date_'+freq] = data['darrays'][freq]

data['cser']['bool'] = [True, False, True, False, True, True]
data['cser']['int32'] = np.arange(6).astype(np.int32)
data['cser']['int64'] = np.arange(6).astype(np.int64)
data['cser']['float32'] = np.arange(6).astype(np.float32)
data['cser']['float64'] = np.arange(6).astype(np.float64)
data['cser']['str'] = ["asdf", "aasssssssss", "zzzzzzzzzzzz", "", "blah"]

for x in data['cser']:
    data['cser'][x] = ma.masked_array(data['cser'][x])
    data['tser'][x] = ts.time_series(data['cser'][x], start_date=data['dates']['a'])

for freq in data['dates']:
    data['freqs'][freq] = ts.time_series(np.arange(20).astype(np.float32), start_date=data['dates'][freq])

# test writing for all data types as time series and as case series
for x in data['tser']:
    data['tser'][x][1] = ma.masked
    data['cser'][x][1] = ma.masked

# series for testing appending data to an existing series
appendTSer = ts.time_series(np.arange(10, 15).astype(np.float32), freq='A', start_date=ts.Date(freq='A', year=2007))
appendCSer = np.arange(10, 15).astype(np.float32)

# series for testing writing over a specified range
rangeTSer = ts.time_series(np.arange(20).astype(np.float32), freq='A', start_date=ts.Date(freq='A', year=2004))
rangeCSer = np.arange(20).astype(np.float32)

data['scalars']['int32'] = np.int32(5)
data['scalars']['int64'] = np.int64(5)
data['scalars']['float32'] = np.float32(5)
data['scalars']['float64'] = np.float64(5)
data['scalars']['pyInt'] = 5
data['scalars']['pyFloat'] = 5234.6323
data['scalars']['string'] = "mystring"
data['scalars']['namelist'] = ["mystring", "$asdf","gggggggg"]
data['scalars']['boolean'] = True
for f in data['dates']:
    data['scalars']['date_'+f] = data['dates'][f]

class test_write(NumpyTestCase):
    
    def setUp(self):
        self.db = fame.FameDb("testdb.db",'o')
        
    def test_main(self):
        "execute all the tests. Order is important here"

        self._test_write_scalars()
        self._test_read_scalars()
        
        self._test_dict_scalars()

        self._test_write_freqs_tser()
        self._test_read_freqs_tser()

        self._test_write_dtypes_tser()
        self._test_read_dtypes_tser()
        
        self._test_read_range_tser()

        self._test_write_append_tser()
        self._test_read_append_tser()
        
        self._test_write_range_tser()
        self._test_verify_write_range_tser()
        
        self._test_write_empty_tser()
        self._test_read_empty_tser()
        
        self._test_overwrite_tser()
        
        self._test_assume_exists_tser()
        
        self._test_dict_tser()
        
        self._test_write_dtypes_cser()
        self._test_read_dtypes_cser()
        
        self._test_read_range_cser()

        self._test_write_append_cser()
        self._test_read_append_cser()
        
        self._test_write_range_cser()
        self._test_verify_write_range_cser()

        self._test_write_empty_cser()
        self._test_read_empty_cser()
        
        self._test_overwrite_cser()
        
        self._test_assume_exists_cser()
        
        self._test_dict_cser()
        
        self._test_whats()
        
        self._test_exists()
        
        self._test_remove()
        
        self._test_wildlist()
        
        

    def _test_write_scalars(self):
        "test writing all types of scalar values"
        for s in data['scalars']:
            self.db.write_scalar('$scalar_'+s, data['scalars'][s])
            
    def _test_dict_scalars(self):
        "test writing multiple scalars at once using write_scalar_dict"
        self.db.write_scalar_dict({'$scalar_1':data['scalars']['float32'],
                                   '$scalar_2':data['scalars']['float32']})
        result = self.db.read(['$scalar_1', '$scalar_2'])
        assert_equal(result['$scalar_1'], data['scalars']['float32'])
        assert_equal(result['$scalar_2'], data['scalars']['float32'])

    def _test_read_scalars(self):
        "read scalars of every data type"
        for s in data['scalars']:
            sclr = self.db.read('$scalar_'+s)
            orig = data['scalars'][s]

            if s == 'int32':
                assert_equal(sclr, orig.astype(np.float32))
            elif s in ('pyInt', 'pyFloat', 'int64'):
                assert_equal(sclr, np.float64(orig))
            elif s == 'namelist':
                assert_equal(sclr, [x.upper() for x in orig])
            else:
                assert_equal(sclr, orig)

    def _test_write_freqs_tser(self):
        "test writing time series for all frequencies"
        for x in data['freqs']:
            self.db.write_tser('$freq_'+x, data['freqs'][x])

    def _test_read_freqs_tser(self):
        """read series at every frequency and ensure they are the
        same as what was written"""
        for x in data['freqs']:
            ser = self.db.read('$freq_'+x)
            assert_mask_equal(ser.mask, data['freqs'][x].mask)
            assert((ser == data['freqs'][x]).all())

    def _test_write_dtypes_tser(self):
        "test writing for all dtypes for time series"
        for x in data['tser']:
            self.db.write_tser('$tser_'+x, data['tser'][x])

    def _test_read_dtypes_tser(self):
        "read time series of every data type"
        for x in data['tser']:
            ser = self.db.read('$tser_'+x)
            if str(ser.dtype)[:5] == 'float' and str(data['tser'][x].dtype)[:3] == 'int':
                ser = ser.astype(data['tser'][x].dtype)
                
            assert_mask_equal(ser.mask, data['tser'][x].mask)
            assert((ser == data['tser'][x]).all())
            
    def _test_read_range_tser(self):
        "test reading a time series over specified ranges"
        src = data['tser']['float32']
        s1 = src.start_date+2
        s2 = src.start_date-2
        e1 = src.end_date+2
        e2 = src.end_date-2
        
        dateList = [(s1, e1),
                    (s1, e2),
                    (s2, e1),
                    (s2, e2)]
                    
        for s, e in dateList:
            res = ts.adjust_endpoints(src, start_date=s, end_date=e)
            ser = self.db.read('$tser_float32', start_date=s, end_date=e)
            assert_array_equal(res, ser)
            

    def _test_write_append_tser(self):
        "test appending data to an existing time series"
        self.db.write_tser('$appendTSer', data['tser']['float32'])
        self.db.write_tser('$appendTSer', appendTSer)
        
    def _test_read_append_tser(self):
        "test reading of appended time series"
        result = ts.adjust_endpoints(data['tser']['float32'],
                                     start_date=data['tser']['float32'].start_date,
                                     end_date=appendTSer.end_date)
        result[appendTSer.start_date:appendTSer.end_date+1] = appendTSer
        
        ser = self.db.read('$appendTSer')
        
        assert_array_equal(result, ser)


    def _test_write_range_tser(self):
        "test writing a time series over a specified range"
        self.db.write_tser('$rangeTSer', rangeTSer,
                           start_date=ts.Date(freq='A', year=2008),
                           end_date=ts.Date(freq='A', year=2012))

    def _test_verify_write_range_tser(self):
        "verify that _test_write_range_write_tser worked as expected"
        
        ser = self.db.read('$rangeTSer')
        
        sDate = ts.Date(freq='A', year=2008)
        eDate = ts.Date(freq='A', year=2012)
        
        assert_array_equal(ser, rangeTSer[sDate:eDate+1])

    def _test_write_empty_tser(self):
        "test writing a time series with no data"
        emptySer = ts.time_series([], freq='A')
        self.db.write_tser('$emptyTSer', emptySer)

    def _test_read_empty_tser(self):
        "test reading a time series with no data"
        ser = self.db.read('$emptyTSer')
        assert(ser.start_date is None)
        
    def _test_overwrite_tser(self):
        "test overwriting a time series"
        self.db.write_tser('$tser_float32', data['tser']['bool'], overwrite=True)
        ser = self.db.read('$tser_float32')
        assert_array_equal(ser, data['tser']['bool'])

    def _test_assume_exists_tser(self):
        "check to see if the assume_exists flag works for write_tser"
        exception = False
        try:
            self.db.write_tser('$doesNotExist', appendTSer, assume_exists=True)
        except fame.DBError:
            exception = True
        assert(exception)
        
    def _test_dict_tser(self):
        "test writing multiple time series at once using write_tser_dict"
        self.db.write_tser_dict({'$tser_1':data['tser']['float32'],
                                   '$tser_2':data['tser']['float32']})
        result = self.db.read(['$tser_1', '$tser_2'])
        assert_array_equal(result['$tser_1'], data['tser']['float32'])
        assert_array_equal(result['$tser_2'], data['tser']['float32'])

    def _test_write_dtypes_cser(self):
        "test writing for all dtypes for case series"""
        for x in data['cser']:
            self.db.write_cser('$cser_'+x, data['cser'][x])

    def _test_read_dtypes_cser(self):
        "read case series of every data type"
        for x in data['cser']:
            ser = self.db.read('$cser_'+x)
            if str(ser.dtype)[:5] == 'float' and str(data['cser'][x].dtype)[:3] == 'int':
                ser = ser.astype(data['cser'][x].dtype)
                
            assert_mask_equal(ser.mask, data['cser'][x].mask)
            assert((ser == data['cser'][x]).all())

    def _test_read_range_cser(self):
        "test reading case series over specified ranges"
        src = data['cser']['float32']
        s1 = 3
        s2 = 1
        e1 = 8
        e2 = 4
        
        caseList = [(s1, e1),
                    (s1, e2),
                    (s2, e1),
                    (s2, e2)]
                    
        for s, e in caseList:
            size = (e - s + 1)
            res = ma.array([0]*size , np.float32, mask=[1]*size )

            if e < src.size: _e = size
            else: _e = size - max(e-size, 0, size - src.size)

            res[0:_e] = src[s-1:min(e, src.size)]
            ser = self.db.read('$cser_float32', start_case=s, end_case=e)

            assert_array_equal(res, ser)

    def _test_write_append_cser(self):
        "test appending to an existing case series"
        self.db.write_cser('$appendCSer', data['cser']['float32'])
        self.db.write_cser('$appendCSer', appendCSer, zero_represents=4)
        
    def _test_read_append_cser(self):
        "test reading of appended case series"
        
        result = ma.concatenate([data['cser']['float32'][:3], appendCSer])
        ser = self.db.read('$appendCSer')
        assert_array_equal(result, ser)
        
    def _test_write_range_cser(self):
        "test writing over a specified range"
        self.db.write_cser('$rangeCSer', rangeCSer,
                           start_case=5, end_case=9)

    def _test_verify_write_range_cser(self):
        "verify that _test_write_range_write_cser worked as expected"
        
        ser = self.db.read('$rangeCSer')
        result = ma.arange(9).astype(np.float32)
        result[:4] = ma.masked
        
        assert_array_equal(ser, result)

    def _test_write_empty_cser(self):
        "test writing a case series with no data"
        self.db.write_cser('$emptyCSer', ma.array([]))

    def _test_read_empty_cser(self):
        "test reading a case series with no data"
        ser = self.db.read('$emptyCSer')
        assert_equal(ser.size, 0)
    
    def _test_overwrite_cser(self):
        "test overwriting a case series"
        self.db.write_cser('$cser_float32', data['cser']['bool'], overwrite=True)
        ser = self.db.read('$cser_float32')
        assert_array_equal(ser, data['cser']['bool'])
        
    def _test_assume_exists_cser(self):
        "check to see if the assume_exists flag works for write_cser"
        exception = False
        try:
            self.db.write_cser('$doesNotExist', appendCSer, assume_exists=True)
        except fame.DBError:
            exception = True
        assert(exception)

    def _test_dict_cser(self):
        "test writing multiple case series at once using write_cser_dict"
        self.db.write_cser_dict({'$cser_1':data['cser']['float32'],
                                   '$cser_2':data['cser']['float32']})
        result = self.db.read(['$cser_1', '$cser_2'])
        assert_array_equal(result['$cser_1'], data['cser']['float32'])
        assert_array_equal(result['$cser_2'], data['cser']['float32'])
        
    def _test_whats(self):
        "test whats method"
        
        # just make sure it doesn't crash for now
        what_dict = self.db.whats('$tser_float32')
        
    def _test_exists(self):
        assert(self.db.exists('$cser_float32'))
        assert(not self.db.exists('$fake_series'))
        
    def _test_remove(self):
        assert(self.db.exists('$cser_1'))
        assert(self.db.exists('$cser_2'))
        self.db.remove(['$cser_1', '$cser_2'])
        assert(not self.db.exists('$cser_1'))
        assert(not self.db.exists('$cser_2'))
        
        self.db.remove('$cser_1', must_exist=False)

        
    def _test_wildlist(self):
        wl1 = self.db.wildlist("$cser_?")
        wl2 = self.db.wildlist("$cser_?", wildonly=True)
        
        res1 = sorted(["$CSER_"+x.upper() for x in list(data['cser'])])
        res2 = sorted([x.upper() for x in list(data['cser'])])
        
        assert_equal(wl1, res1)
        assert_equal(wl2, res2)

    
    def tearDown(self):
        self.db.close()

        
        
###############################################################################
#------------------------------------------------------------------------------
if __name__ == "__main__":
    NumpyTest().run()
