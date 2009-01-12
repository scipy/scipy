#!/usr/bin/env python
import numpy as np
import nose
from scipy.stsci.image import *
from numpy.testing import *

def test_average1():
    """
    average() nominally computes the average pixel value for a stack of
    identically shaped images.

    arrays     specifies a sequence of inputs arrays, which are nominally a
               stack of identically shaped images.

    output     may be used to specify the output array.  If none is specified,
               either arrays[0] is copied or a new array of type 'outtype'
               is created.

    outtype    specifies the type of the output array when no 'output' is
               specified.

    nlow       specifies the number of pixels to be excluded from average
               on the low end of the pixel stack.

    nhigh      specifies the number of pixels to be excluded from average
               on the high end of the pixel stack.

    badmasks   specifies boolean arrays corresponding to 'arrays', where true
               indicates that a particular pixel is not to be included in the
               average calculation.
    """
    a = np.arange(4)
    a = a.reshape((2,2))
    arrays = [a*16, a*4, a*2, a*8]
    result = average(arrays)
    test = np.array([[ 0,  7],
           [15, 22]])
    assert_equal(result,test)
    
def test_average2():
    a = np.arange(4)
    a = a.reshape((2,2))
    arrays = [a*16, a*4, a*2, a*8]
    result = average(arrays, nhigh=1)
    test = np.array([[ 0,  4],
           [ 9, 14]])
    assert_equal(result,test)

def test_average3():
    a = np.arange(4)
    a = a.reshape((2,2))
    arrays = [a*16, a*4, a*2, a*8]
    result = average(arrays, nlow=1)
    test = np.array([[ 0,  9],
           [18, 28]])
    assert_equal(result,test)

def test_average4():
    a = np.arange(4)
    a = a.reshape((2,2))
    arrays = [a*16, a*4, a*2, a*8]
    result = average(arrays, outtype=np.float32)
    test = np.array([[  0. ,   7.5],
           [ 15. ,  22.5]], dtype=np.float32)
    assert_equal(result,test)

def test_average5():
    a = np.arange(4)
    a = a.reshape((2,2))
    arrays = [a*16, a*4, a*2, a*8]
    bm = np.zeros((4,2,2), dtype=np.bool8)
    bm[2,...] = 1
    result = average(arrays, badmasks=bm)
    test = np.array([[ 0,  9],
           [18, 28]])
    assert_equal(result,test)

def test_average6():
    a = np.arange(4)
    a = a.reshape((2,2))
    arrays = [a*16, a*4, a*2, a*8]
    result = average(arrays, badmasks=threshhold(arrays, high=25))
    test = np.array([[ 0,  7],
           [ 9, 14]])
    assert_equal(result,test)

if __name__ == "__main__":
    run_module_suite()
