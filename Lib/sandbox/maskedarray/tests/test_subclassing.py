import numpy as N
import numpy.core.numeric as numeric

from numpy.testing import NumpyTest, NumpyTestCase

import maskedarray.testutils
reload(maskedarray.testutils)
from maskedarray.testutils import *

import maskedarray.core as coremodule
reload(coremodule)
from maskedarray.core import *



class SubArray(N.ndarray):
    """Defines a generic N.ndarray subclass, that stores some metadata
    in the  dictionary `info`."""
    def __new__(cls,arr,info={}):
        x = N.asanyarray(arr).view(cls)
        x.info = info
        return x
    def __array_finalize__(self, obj):
        self.info = getattr(obj,'info',{})
        return
    def __add__(self, other):
        result = N.ndarray.__add__(self, other)
        result.info.update({'added':result.info.pop('added',0)+1})
        return result
subarray = SubArray

class MSubArray(SubArray,MaskedArray):
    def __new__(cls, data, info=None, mask=nomask):
        subarr = SubArray(data, info)
        _data = MaskedArray.__new__(cls, data=subarr, mask=mask)
        _data.info = subarr.info
        return _data
    def __array_finalize__(self,obj):
        SubArray.__array_finalize__(self, obj)
        MaskedArray.__array_finalize__(self,obj)
        return
msubarray = MSubArray

class MMatrix(MaskedArray, N.matrix,):
    def __new__(cls, data, mask=nomask):
        mat = N.matrix(data)
        _data = MaskedArray.__new__(cls, data=mat, mask=mask)
        return _data
    def __array_finalize__(self,obj):
        N.matrix.__array_finalize__(self, obj)
        MaskedArray.__array_finalize__(self,obj)
        return     
mmatrix = MMatrix 
        


class test_subclassing(NumpyTestCase):
    """Test suite for masked subclasses of ndarray."""
                
    def check_data_subclassing(self):
        "Tests whether the subclass is kept."
        x = N.arange(5)
        m = [0,0,1,0,0]
        xsub = SubArray(x)
        xmsub = masked_array(xsub, mask=m)
        assert isinstance(xmsub, MaskedArray)
        assert_equal(xmsub._data, xsub)
        assert isinstance(xmsub._data, SubArray)
        
    def check_maskedarray_subclassing(self):
        "Tests subclassing MaskedArray"
        x = N.arange(5)
        mx = mmatrix(x,mask=[0,1,0,0,0])
        assert isinstance(mx._data, N.matrix)
        "Tests masked_unary_operation"
        assert isinstance(add(mx,mx), mmatrix)
        assert isinstance(add(mx,x), mmatrix)
        assert_equal(add(mx,x), mx+x)
        assert isinstance(add(mx,mx)._data, N.matrix)
        assert isinstance(add.outer(mx,mx), mmatrix)
        "Tests masked_binary_operation"
        assert isinstance(hypot(mx,mx), mmatrix)
        assert isinstance(hypot(mx,x), mmatrix)

    def check_attributepropagation(self):
        x = array(arange(5), mask=[0]+[1]*4)
        my = masked_array(subarray(x))
        ym = msubarray(x)
        #
        z = (my+1)
        assert isinstance(z,MaskedArray)
        assert not isinstance(z, MSubArray)
        assert isinstance(z._data, SubArray)
        assert_equal(z._data.info, {})
        #
        z = (ym+1)
        assert isinstance(z, MaskedArray)
        assert isinstance(z, MSubArray)
        assert isinstance(z._data, SubArray)
        assert z._data.info['added'] > 0
        
################################################################################
if __name__ == '__main__':
    NumpyTest().run()
    if 1:
        x = N.arange(5)
        mx = mmatrix(x,mask=[0,1,0,0,0])
        
        c = mx.ravel()
        
        a = add(mx,x)
        b = mx+x
        assert_equal(a,b)


