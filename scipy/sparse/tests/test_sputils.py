"""unit tests for sparse utility functions"""

import numpy as np

from np.testing import *
set_package_path()
from scipy.sparse.sputils import *
restore_path()



class TestSparseUtils(NumpyTestCase):

    def check_upcast(self):
        assert_equal(upcast('int32'),np.int32)
        assert_equal(upcast('int32','float32'),np.float64)
        assert_equal(upcast('bool',complex,float),np.complex128)
        assert_equal(upcast('i','d'),np.float64)

    def check_getdtype(self):
        A = np.array([1],dtype='int8')

        assert_equal(getdtype(None,default=float),np.float)
        assert_equal(getdtype(None,a=A),np.int8)

    def check_isscalarlike(self):
        assert_equal(isscalarlike(3.0),True)
        assert_equal(isscalarlike(-4),True)
        assert_equal(isscalarlike(2.5),True)
        assert_equal(isscalarlike(1 + 3j),True)
        assert_equal(isscalarlike(np.array(3)),True)
        assert_equal(isscalarlike( "16" ), True)

        assert_equal(isscalarlike( np.array([3])), False)
        assert_equal(isscalarlike( [[3]] ), False)
        assert_equal(isscalarlike( (1,) ), False)
        assert_equal(isscalarlike( (1,2) ), False)

    def check_isintlike(self):
        assert_equal(isintlike(3.0),True)
        assert_equal(isintlike(-4),True)
        assert_equal(isintlike(np.array(3)),True)
        assert_equal(isintlike( np.array([3])), True)

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
        assert_equal(issequence( np.array([1,2,3]) ),True)
        
        assert_equal(issequence( np.array([[1],[2],[3]]) ),False)
        assert_equal(issequence( 3 ),False)

    def check_isdense(self):
        assert_equal(isdense( np.array([1]) ),True)
        assert_equal(isdense( np.matrix([1]) ),True)
                
if __name__ == "__main__":
    NumpyTest().run()


