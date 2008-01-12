# pylint: disable-msg=W0611, W0612, W0511,R0201
"""Tests suite for mrecarray.

:author: Pierre Gerard-Marchant
:contact: pierregm_at_uga_dot_edu
:version: $Id: test_mrecarray.py 68 2007-01-10 17:46:04Z backtopop $
"""
__author__ = "Pierre GF Gerard-Marchant ($Author: backtopop $)"
__version__ = '1.0'
__revision__ = "$Revision: 68 $"
__date__     = '$Date: 2007-01-10 12:46:04 -0500 (Wed, 10 Jan 2007) $'

import types

import numpy as N
import numpy.core.fromnumeric  as fromnumeric
from scipy.testing import *
from scipy.sandbox.maskedarray.testutils import *

import scipy.sandbox.maskedarray.core as MA
from scipy.sandbox.maskedarray.mrecords import mrecarray, fromarrays, fromtextfile, fromrecords


#..............................................................................
class TestMrecarray(TestCase):
    "Base test class for MaskedArrays."
    def __init__(self, *args, **kwds):
        TestCase.__init__(self, *args, **kwds)
        self.setup()

    def setup(self):
        "Generic setup"
        d = N.arange(5)
        m = MA.make_mask([1,0,0,1,1])
        base_d = N.r_[d,d[::-1]].reshape(2,-1).T
        base_m = N.r_[[m, m[::-1]]].T
        base = MA.array(base_d, mask=base_m)
        mrec = fromarrays(base.T,)
        self.data = [d, m, mrec]

    def test_get(self):
        "Tests fields retrieval"
        [d, m, mrec] = self.data
        assert_equal(mrec.f0, MA.array(d,mask=m))
        assert_equal(mrec.f1, MA.array(d[::-1],mask=m[::-1]))
        assert((mrec._fieldmask == N.core.records.fromarrays([m, m[::-1]])).all())
        assert_equal(mrec._mask, N.r_[[m,m[::-1]]].all(0))
        assert_equal(mrec.f0, mrec['f0'])

    def test_set(self):
        "Tests setting fields/attributes."
        [d, m, mrec] = self.data
        mrec.f0_data = 5
        assert_equal(mrec['f0_data'], [5,5,5,5,5])
        mrec.f0 = 1
        assert_equal(mrec['f0_data'], [1]*5)
        assert_equal(mrec['f0_mask'], [0]*5)
        mrec.f1 = MA.masked
        assert_equal(mrec.f1.mask, [1]*5)
        assert_equal(mrec['f1_mask'], [1]*5)
        mrec._mask = MA.masked
        assert_equal(mrec['f1_mask'], [1]*5)
        assert_equal(mrec['f0_mask'],mrec['f1_mask'])
        mrec._mask = MA.nomask
        assert_equal(mrec['f1_mask'], [0]*5)
        assert_equal(mrec['f0_mask'],mrec['f1_mask'])

    def test_hardmask(self):
        "Test hardmask"
        [d, m, mrec] = self.data
        print mrec._mask
        mrec.harden_mask()
        assert(mrec._hardmask)
        mrec._mask = nomask
        assert_equal(mrec._mask, N.r_[[m,m[::-1]]].all(0))
        mrec.soften_mask()
        assert(not mrec._hardmask)
        mrec._mask = nomask
        assert_equal(mrec['f1_mask'], [0]*5)
        assert_equal(mrec['f0_mask'],mrec['f1_mask'])

    def test_fromtextfile(self):
        "Tests reading from a text file."
        fcontent = """#
'One (S)','Two (I)','Three (F)','Four (M)','Five (-)','Six (C)'
'strings',1,1.0,'mixed column',,1
'with embedded "double quotes"',2,2.0,1.0,,1
'strings',3,3.0E5,3,,1
'strings',4,-1e-10,,,1
"""
        import os
        from datetime import datetime
        fname = 'tmp%s' % datetime.now().strftime("%y%m%d%H%M%S%s")
        f = open(fname, 'w')
        f.write(fcontent)
        f.close()
        mrectxt = fromtextfile(fname,delimitor=',',varnames='ABCDEFG')
        os.unlink(fname)
        #
        assert(isinstance(mrectxt, mrecarray))
        assert_equal(mrectxt.F, [1,1,1,1])
        assert_equal(mrectxt.E._mask, [1,1,1,1])
        assert_equal(mrectxt.C, [1,2,3.e+5,-1e-10])

    def test_fromrecords(self):
        "Test from recarray."
        [d, m, mrec] = self.data
        nrec = N.core.records.fromarrays(N.r_[[d,d[::-1]]])
        mrecfr = fromrecords(nrec.tolist())
        assert_equal(mrecfr.f0, mrec.f0)
        assert_equal(mrecfr.dtype, mrec.dtype)
        #....................
        mrecfr = fromrecords(nrec)
        assert_equal(mrecfr.f0, mrec.f0)
        assert_equal(mrecfr.dtype, mrec.dtype)
        #....................
        tmp = mrec[::-1] #.tolist()
        mrecfr = fromrecords(tmp)
        assert_equal(mrecfr.f0, mrec.f0[::-1])

###############################################################################
#------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
