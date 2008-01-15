# pylint: disable-msg=W0611, W0612, W0511,R0201
"""Tests suite for trecords.

:author: Pierre Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu  & mattknox_ca_at_hotmail_dot_com
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import types

import numpy
import numpy.core.fromnumeric  as fromnumeric
from scipy.testing import *

from scipy.sandbox.maskedarray.testutils import assert_equal, assert_array_equal, assert_equal_records

import scipy.sandbox.maskedarray.core as MA
import scipy.sandbox.maskedarray.mrecords as MR
from scipy.sandbox.maskedarray.mrecords import addfield

from scipy.sandbox.maskedarray.core import getmaskarray, nomask, masked_array

from scipy.sandbox.timeseries.trecords import \
     TimeSeriesRecords, TimeSeries,\
    fromarrays, fromtextfile, fromrecords, \
    date_array, time_series


#..............................................................................
class TestMRecords(TestCase):
    "Base test class for MaskedArrays."
    def __init__(self, *args, **kwds):
        TestCase.__init__(self, *args, **kwds)
        self.setup()

    def setup(self):
        "Generic setup"
        d = numpy.arange(5)
        m = MA.make_mask([1,0,0,1,1])
        base_d = numpy.r_[d,d[::-1]].reshape(2,-1).T
        base_m = numpy.r_[[m, m[::-1]]].T
        base = MA.array(base_d, mask=base_m)
        mrec = MR.fromarrays(base.T,)
        dlist = ['2007-%02i' % (i+1) for i in d]
        dates = date_array(dlist)
        ts = time_series(mrec,dates)
        mts = TimeSeriesRecords(mrec,dates)
        self.data = [d, m, mrec, dlist, dates, ts, mts]

    def test_get(self):
        "Tests fields retrieval"
        [d, m, mrec, dlist, dates, ts, mts] = self.data
        assert(isinstance(mts['f0'], TimeSeries))
        assert_equal(mts['f0']._dates, dates)
        assert_equal(mts['f0']._data, d)
        assert_equal(mts['f0']._mask, m)
        #
        assert(isinstance(mts[0], TimeSeriesRecords))
        assert_equal(mts._data[0], mrec._data[0])
        # We can't use assert_equal here, as it tries to convert the tuple into a singleton
#        assert(mts[0]._data.view(numpyndarray) == mrec[0])
        assert_equal_records(mts._data[0], mrec[0])
        assert_equal(mts._dates[0], dates[0])
        assert_equal(mts[0]._dates, dates[0])
        #
        assert(isinstance(mts['2007-01'], TimeSeriesRecords))
        assert_equal_records(mts['2007-01']._data, mrec[0])
        assert_equal(mts['2007-01']._dates, dates[0])
        #
        assert(isinstance(mts.f0, TimeSeries))
        assert_equal(mts.f0, time_series(d, dates=dates, mask=m))
        assert_equal(mts.f1, time_series(d[::-1], dates=dates, mask=m[::-1]))
        assert((mts._fieldmask == numpy.core.records.fromarrays([m, m[::-1]])).all())
        assert_equal(mts._mask, numpy.r_[[m,m[::-1]]].all(0))
        assert_equal(mts.f0[1], mts[1].f0)
        #
        assert(isinstance(mts[:2], TimeSeriesRecords))
        assert_equal(mts[:2]._data.f0, mrec[:2].f0)
        assert_equal(mts[:2]._data.f1, mrec[:2].f1)
        assert_equal(mts[:2]._dates, dates[:2])

    def test_set(self):
        "Tests setting fields/attributes."
        [d, m, mrec, dlist, dates, ts, mts] = self.data
        mts.f0._data[:] = 5
        assert_equal(mts['f0']._data, [5,5,5,5,5])
        mts.f0 = 1
        assert_equal(mts['f0']._data, [1]*5)
        assert_equal(getmaskarray(mts['f0']), [0]*5)
        mts.f1 = MA.masked
        assert_equal(mts.f1.mask, [1]*5)
        assert_equal(getmaskarray(mts['f1']), [1]*5)
        mts._mask = MA.masked
        assert_equal(getmaskarray(mts['f1']), [1]*5)
        assert_equal(mts['f0']._mask, mts['f1']._mask)
        mts._mask = MA.nomask
        assert_equal(getmaskarray(mts['f1']), [0]*5)
        assert_equal(mts['f0']._mask, mts['f1']._mask)

    def test_setslices(self):
        "Tests setting slices."
        [d, m, mrec, dlist, dates, ts, mts] = self.data
        #
        mts[:2] = 5
        assert_equal(mts.f0._data, [5,5,2,3,4])
        assert_equal(mts.f1._data, [5,5,2,1,0])
        assert_equal(mts.f0._mask, [0,0,0,1,1])
        assert_equal(mts.f1._mask, [0,0,0,0,1])
        mts.harden_mask()
        mts[-2:] = 5
        assert_equal(mts.f0._data, [5,5,2,3,4])
        assert_equal(mts.f1._data, [5,5,2,5,0])
        assert_equal(mts.f0._mask, [0,0,0,1,1])
        assert_equal(mts.f1._mask, [0,0,0,0,1])

    def test_hardmask(self):
        "Test hardmask"
        [d, m, mrec, dlist, dates, ts, mts] = self.data
        mts.harden_mask()
        assert(mts._hardmask)
        mts._mask = nomask
        assert_equal(mts._mask, numpy.r_[[m,m[::-1]]].all(0))
        mts.soften_mask()
        assert(not mts._hardmask)
        mts._mask = nomask
        assert(mts['f1']._mask is nomask)
        assert_equal(mts['f0']._mask,mts['f1']._mask)

    def test_addfield(self):
        "Tests addfield"
        [d, m, mrec, dlist, dates, ts, mts] = self.data
        mts = addfield(mts, masked_array(d+10, mask=m[::-1]))
        assert_equal(mts.f2, d+10)
        assert_equal(mts.f2._mask, m[::-1])

    def test_fromrecords(self):
        "Test from recarray."
        [d, m, mrec, dlist, dates, ts, mts] = self.data
        nrec = numpy.core.records.fromarrays(numpy.r_[[d,d[::-1]]])
        mrecfr = fromrecords(nrec.tolist(), dates=dates)
        assert_equal(mrecfr.f0, mrec.f0)
        assert_equal(mrecfr.dtype, mrec.dtype)
        #....................
        altrec = [tuple([d,]+list(r)) for (d,r) in zip(dlist,nrec)]
        mrecfr = fromrecords(altrec, names='dates,f0,f1')
        assert_equal(mrecfr.f0, mrec.f0)
        assert_equal(mrecfr.dtype, mrec.dtype)
        #....................
        tmp = TimeSeriesRecords(mts._series[::-1], dates=mts.dates)
        mrecfr = fromrecords(tmp)
        assert_equal(mrecfr.f0, mrec.f0[::-1])

    def test_fromtextfile(self):
        "Tests reading from a text file."
        fcontent = """#
'Dates', 'One (S)','Two (I)','Three (F)','Four (M)','Five (-)','Six (C)'
'2007-01', 'strings',1,1.0,'mixed column',,1
'2007-02', 'with embedded "double quotes"',2,2.0,1.0,,1
'2007-03', 'strings',3,3.0E5,3,,1
'2007-05','strings',4,-1e-10,,,1
"""
        import os
        from datetime import datetime
        fname = 'tmp%s' % datetime.now().strftime("%y%m%d%H%M%S%s")
        f = open(fname, 'w')
        f.write(fcontent)
        f.close()
        mrectxt = fromtextfile(fname,delimitor=',',varnames='ABCDEFG',
                               dates_column=0)
        os.unlink(fname)
        #
        dlist = ['2007-%02i' % i for i in (1,2,3,5)]
        assert(isinstance(mrectxt, TimeSeriesRecords))
        assert_equal(mrectxt._dates, date_array(dlist,'M'))
        assert_equal(mrectxt.dtype.names, ['B','C','D','E','F','G'])
        assert_equal(mrectxt.G, [1,1,1,1])
        assert_equal(mrectxt.F._mask, [1,1,1,1])
        assert_equal(mrectxt.D, [1,2,3.e+5,-1e-10])

###############################################################################
#------------------------------------------------------------------------------
if __name__ == "__main__":
    nose.run(argv=['', __file__])
