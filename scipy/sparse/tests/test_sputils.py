"""unit tests for sparse utility functions"""

from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import TestCase, run_module_suite, assert_equal
from scipy.sparse import sputils


class TestSparseUtils(TestCase):

    def test_upcast(self):
        assert_equal(sputils.upcast('intc'),np.intc)
        assert_equal(sputils.upcast('int32','float32'),np.float64)
        assert_equal(sputils.upcast('bool',complex,float),np.complex128)
        assert_equal(sputils.upcast('i','d'),np.float64)

    def test_getdtype(self):
        A = np.array([1],dtype='int8')

        assert_equal(sputils.getdtype(None,default=float),np.float)
        assert_equal(sputils.getdtype(None,a=A),np.int8)

    def test_isscalarlike(self):
        assert_equal(sputils.isscalarlike(3.0),True)
        assert_equal(sputils.isscalarlike(-4),True)
        assert_equal(sputils.isscalarlike(2.5),True)
        assert_equal(sputils.isscalarlike(1 + 3j),True)
        assert_equal(sputils.isscalarlike(np.array(3)),True)
        assert_equal(sputils.isscalarlike("16"), True)

        assert_equal(sputils.isscalarlike(np.array([3])), False)
        assert_equal(sputils.isscalarlike([[3]]), False)
        assert_equal(sputils.isscalarlike((1,)), False)
        assert_equal(sputils.isscalarlike((1,2)), False)

    def test_isintlike(self):
        assert_equal(sputils.isintlike(3.0),True)
        assert_equal(sputils.isintlike(-4),True)
        assert_equal(sputils.isintlike(np.array(3)),True)
        assert_equal(sputils.isintlike(np.array([3])), False)

        assert_equal(sputils.isintlike(2.5),False)
        assert_equal(sputils.isintlike(1 + 3j),False)
        assert_equal(sputils.isintlike((1,)), False)
        assert_equal(sputils.isintlike((1,2)), False)

    def test_isshape(self):
        assert_equal(sputils.isshape((1,2)),True)
        assert_equal(sputils.isshape((5,2)),True)

        assert_equal(sputils.isshape((1.5,2)),False)
        assert_equal(sputils.isshape((2,2,2)),False)
        assert_equal(sputils.isshape(([2],2)),False)

    def test_issequence(self):
        assert_equal(sputils.issequence((1,)),True)
        assert_equal(sputils.issequence((1,2,3)),True)
        assert_equal(sputils.issequence([1]),True)
        assert_equal(sputils.issequence([1,2,3]),True)
        assert_equal(sputils.issequence(np.array([1,2,3])),True)

        assert_equal(sputils.issequence(np.array([[1],[2],[3]])),False)
        assert_equal(sputils.issequence(3),False)

    def test_isdense(self):
        assert_equal(sputils.isdense(np.array([1])),True)
        assert_equal(sputils.isdense(np.matrix([1])),True)

if __name__ == "__main__":
    run_module_suite()
