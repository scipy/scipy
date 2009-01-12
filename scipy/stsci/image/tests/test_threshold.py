#!/usr/bin/env python
import numpy as np
import nose
from scipy.stsci.image import *
from numpy.testing import *

def test_threshhold1():
    """
    threshhold() computes a boolean array 'outputs' with
    corresponding elements for each element of arrays.  The
    boolean value is true where each of the arrays values
    is < the low or >= the high threshholds.
    """

    a=np.arange(100)
    a=a.reshape((10,10))
    result = (threshhold(a, 1, 50)).astype(np.int8)
    test = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]], dtype=np.int8)
    assert_equal(result,test)

def test_threshold2():
    a=np.arange(100)
    a=a.reshape((10,10))
    result = (threshhold([ range(10)]*10, 3, 7)).astype(np.int8)
    test = np.array([[1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
           [1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
           [1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
           [1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
           [1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
           [1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
           [1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
           [1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
           [1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
           [1, 1, 1, 0, 0, 0, 0, 1, 1, 1]], dtype=np.int8)
    assert_equal(result,test)

def test_threshold3():
    a=np.arange(100)
    a=a.reshape((10,10))
    result = (threshhold(a, high=50)).astype(np.int8)
    test = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]], dtype=np.int8)
    assert_equal(result,test)
    
def test_threshold4():
    a=np.arange(100)
    a=a.reshape((10,10))
    result = (threshhold(a, low=50)).astype(np.int8)
    test = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype=np.int8)

if __name__ == "__main__":
    run_module_suite()
