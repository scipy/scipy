# pylint: disable-msg=W0611, W0612, W0511
"""Tests suite for MaskedArray.
Adapted from the original test_ma by Pierre Gerard-Marchant

:author: Pierre Gerard-Marchant
:contact: pierregm_at_uga_dot_edu
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import numpy as N
from numpy.testing import NumpyTest, NumpyTestCase
from numpy.testing.utils import build_err_msg

import maskedarray.testutils
reload(maskedarray.testutils)
from maskedarray.testutils import *

import maskedarray.core
reload(maskedarray.core)
from maskedarray.core import *
import maskedarray.extras
reload(maskedarray.extras)
from maskedarray.extras import *

       
class test_concatenator(NumpyTestCase):
    "Tests for mr_, the equivalent of r_ for masked arrays."
    def check_1d(self):
        "Tests mr_ on 1D arrays."
        assert_array_equal(mr_[1,2,3,4,5,6],array([1,2,3,4,5,6]))
        b = ones(5)
        m = [1,0,0,0,0]
        d = masked_array(b,mask=m)
        c = mr_[d,0,0,d]
        assert(isinstance(c,MaskedArray) or isinstance(c,core.MaskedArray))
        assert_array_equal(c,[1,1,1,1,1,0,0,1,1,1,1,1])
        assert_array_equal(c.mask, mr_[m,0,0,m])

    def check_2d(self):
        "Tests mr_ on 2D arrays."
        a_1 = rand(5,5)
        a_2 = rand(5,5)
        m_1 = N.round_(rand(5,5),0)
        m_2 = N.round_(rand(5,5),0)
        b_1 = masked_array(a_1,mask=m_1)
        b_2 = masked_array(a_2,mask=m_2)
        d = mr_['1',b_1,b_2]  # append columns
        assert(d.shape == (5,10))
        assert_array_equal(d[:,:5],b_1)
        assert_array_equal(d[:,5:],b_2)
        assert_array_equal(d.mask, N.r_['1',m_1,m_2])
        d = mr_[b_1,b_2]
        assert(d.shape == (10,5))
        assert_array_equal(d[:5,:],b_1)
        assert_array_equal(d[5:,:],b_2)
        assert_array_equal(d.mask, N.r_[m_1,m_2])

class test_notmasked(NumpyTestCase):
    "Tests notmasked_edges and notmasked_contiguous."
    def check_edges(self):
        "Tests unmasked_edges"
        a = masked_array(N.arange(24).reshape(3,8),
                         mask=[[0,0,0,0,1,1,1,0],
                               [1,1,1,1,1,1,1,1],
                               [0,0,0,0,0,0,1,0],])
        #
        assert_equal(notmasked_edges(a, None), [0,23])
        #
        tmp = notmasked_edges(a, 0)
        assert_equal(tmp[0], (array([0,0,0,0,2,2,0]), array([0,1,2,3,4,5,7])))
        assert_equal(tmp[1], (array([2,2,2,2,2,2,2]), array([0,1,2,3,4,5,7])))
        #
        tmp = notmasked_edges(a, 1)
        assert_equal(tmp[0], (array([0,2,]), array([0,0])))
        assert_equal(tmp[1], (array([0,2,]), array([7,7])))
    
    def check_contiguous(self):
        "Tests notmasked_contiguous"
        a = masked_array(N.arange(24).reshape(3,8),
                         mask=[[0,0,0,0,1,1,1,1],
                               [1,1,1,1,1,1,1,1],
                               [0,0,0,0,0,0,1,0],])
        tmp = notmasked_contiguous(a, None)
        assert_equal(tmp[-1], (6, (16,21)))
        assert_equal(tmp[-2], (4, (0,3)))
        #
        tmp = notmasked_contiguous(a, 0)
        assert(len(tmp[-1]) == 1)
        assert(tmp[-2] is None)
        assert_equal(tmp[-3],tmp[-1])
        assert(len(tmp[0]) == 2)
        #
        tmp = notmasked_contiguous(a, 1)
        assert_equal(tmp[0][-1], (4, (0,3)))
        assert(tmp[1] is None)
        assert_equal(tmp[2][-1], (6, (0,5)))

###############################################################################
#------------------------------------------------------------------------------
if __name__ == "__main__":
    NumpyTest().run()