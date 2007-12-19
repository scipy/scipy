"""unit tests for sparse utility functions"""

import numpy

from numpy.testing import *
set_package_path()
from scipy.sparse.sputils import *
restore_path()



class TestSparseUtils(NumpyTestCase):

    def check_upcast(self):
        assert_equal(upcast('int32'),numpy.int32)
        assert_equal(upcast('int32','float32'),numpy.float64)
        assert_equal(upcast('bool',complex,float),numpy.complex128)
        assert_equal(upcast('i','d'),numpy.float64)

    def check_getdtype(self):
        A = numpy.array([1],dtype='int8')

        assert_equal(getdtype(None,default=float),numpy.float)
        assert_equal(getdtype(None,a=A),numpy.int8)

    def check_isscalarlike(self):
        assert_equal(isscalarlike(3.0),True)
        assert_equal(isscalarlike(-4),True)
        assert_equal(isscalarlike(2.5),True)
        assert_equal(isscalarlike(1 + 3j),True)
        assert_equal(isscalarlike(numpy.array(3)),True)
        assert_equal(isscalarlike( "16" ), True)

        assert_equal(isscalarlike( numpy.array([3])), False)
        assert_equal(isscalarlike( [[3]] ), False)
        assert_equal(isscalarlike( (1,) ), False)
        assert_equal(isscalarlike( (1,2) ), False)

    def check_isintlike(self):
        assert_equal(isintlike(3.0),True)
        assert_equal(isintlike(-4),True)
        assert_equal(isintlike(numpy.array(3)),True)
        assert_equal(isintlike( numpy.array([3])), True)

        assert_equal(isintlike(2.5),False)
        assert_equal(isintlike(1 + 3j),False)
        assert_equal(isintlike( (1,) ), False)
        assert_equal(isintlike( (1,2) ), False)

    def check_isshape(self):
        assert_equal(isshape( (1,2) ),True)
        assert_equal(isshape( (5,2) ),True)

        assert_equal(isshape( (-1,4) ),False)
        assert_equal(isshape( (1.5,2) ),False)
        assert_equal(isshape( (0,4) ),False)
        assert_equal(isshape( (2,2,2) ),False)

    def check_issequence(self):
        assert_equal(issequence( (1,) ),True)
        assert_equal(issequence( (1,2,3) ),True)
        assert_equal(issequence( [1] ),True)
        assert_equal(issequence( [1,2,3] ),True)
        assert_equal(issequence( numpy.array([1,2,3]) ),True)
        
        assert_equal(issequence( numpy.array([[1],[2],[3]]) ),False)
        assert_equal(issequence( 3 ),False)

    def check_isdense(self):
        assert_equal(isdense( numpy.array([1]) ),True)
        assert_equal(isdense( numpy.matrix([1]) ),True)
                
if __name__ == "__main__":
    NumpyTest().run()


