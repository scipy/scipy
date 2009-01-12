#!/usr/bin/env python
import numpy as np
import nose
from scipy.stsci.convolve import *
from numpy.testing import *
import scipy.stsci.convolve._correlate as _correlate
import scipy.stsci.convolve.iraf_frame as iraf_frame
import numpy.fft as dft


def test_correlate1():
    """
    correlate(data, kernel, mode=FULL)
    """
    result = correlate(np.arange(8), [1, 2], mode=VALID)
    test = np.array([ 2,  5,  8, 11, 14, 17, 20])
    assert_equal(result,test)
    
def test_correlate2():
    result = correlate(np.arange(8), [1, 2], mode=SAME)
    test = np.array([ 0,  2,  5,  8, 11, 14, 17, 20])
    assert_equal(result,test)

def test_correlate3():
    result = correlate(np.arange(8), [1, 2], mode=FULL)
    test = np.array([ 0,  2,  5,  8, 11, 14, 17, 20,  7])
    assert_equal(result,test)

def test_correlate4():
    test = correlate(np.arange(8), [1, 2, 3], mode=VALID)
    result = np.array([ 8, 14, 20, 26, 32, 38])
    assert_equal(result,test)

def test_correlate5():
    test = correlate(np.arange(8), [1, 2, 3], mode=SAME)
    result = np.array([ 3,  8, 14, 20, 26, 32, 38, 20])
    assert_equal(result,test)

def test_correlate6():
    test = correlate(np.arange(8), [1, 2, 3], mode=FULL)
    result = np.array([ 0,  3,  8, 14, 20, 26, 32, 38, 20,  7])
    assert_equal(result,test)

def test_correlate7():
    test = correlate(np.arange(8), [1, 2, 3, 4, 5, 6], mode=VALID)
    result = np.array([ 70,  91, 112])
    assert_equal(result,test)
    
def test_correlate8():
    test = correlate(np.arange(8), [1, 2, 3, 4, 5, 6], mode=SAME)
    result = np.array([ 17,  32,  50,  70,  91, 112,  85,  60])
    assert_equal(result,test)
    
def test_correlate9():
    test = correlate(np.arange(8), [1, 2, 3, 4, 5, 6], mode=FULL)
    result = np.array([  0,   6,  17,  32,  50,  70,  91, 112,  85,  60,  38,  20,   7])
    assert_equal(result,test)
    
def test_correlate10():
    test = False
    try:
        result = correlate(np.arange(8), 1+1j)
    except TypeError:
        test=True
    assert_equal(test,True)

def test_convolve1():
    """
    convolve(data, kernel, mode=FULL)
    Returns the discrete, linear convolution of 1-D
    sequences a and v; mode can be 0 (VALID), 1 (SAME), or 2 (FULL)
    to specify size of the resulting sequence.
    """
    result = convolve(np.arange(8), [1, 2], mode=VALID)
    test = np.array([ 1,  4,  7, 10, 13, 16, 19])
    assert_equal(result,test)

def test_convolve2():
    result = convolve(np.arange(8), [1, 2], mode=SAME)
    test = np.array([ 0,  1,  4,  7, 10, 13, 16, 19])
    assert_equal(result,test)

def test_convolve3():
    result = convolve(np.arange(8), [1, 2], mode=FULL)
    test = np.array([ 0,  1,  4,  7, 10, 13, 16, 19, 14])
    assert_equal(result,test)

def test_convolve4():
    result = convolve(np.arange(8), [1, 2, 3], mode=VALID)
    test = np.array([ 4, 10, 16, 22, 28, 34])
    assert_equal(result,test)

def test_convolve5():
    result = convolve(np.arange(8), [1, 2, 3], mode=SAME)
    test = np.array([ 1,  4, 10, 16, 22, 28, 34, 32])
    assert_equal(result,test)

def test_convolve6():
    result = convolve(np.arange(8), [1, 2, 3], mode=FULL)
    test = np.array([ 0,  1,  4, 10, 16, 22, 28, 34, 32, 21])
    assert_equal(result,test)

def test_convolve7():
    result = convolve(np.arange(8), [1, 2, 3, 4, 5, 6], mode=VALID)
    test = np.array([35, 56, 77])
    assert_equal(result,test)

def test_convolve8():
    result = convolve(np.arange(8), [1, 2, 3, 4, 5, 6], mode=SAME)
    test = np.array([ 4, 10, 20, 35, 56, 77, 90, 94])
    assert_equal(result,test)

def test_convolve9():
    result = convolve(np.arange(8), [1, 2, 3, 4, 5, 6], mode=FULL)
    test = np.array([ 0,  1,  4, 10, 20, 35, 56, 77, 90, 94, 88, 71, 42])
    assert_equal(result,test)

def test_convolve10():
    result = convolve([1.,2.], np.arange(10.))
    test = np.array([  0.,   1.,   4.,   7.,  10.,  13.,  16.,  19.,  22.,  25.,  18.])
    assert_equal(result,test)

def test_correlate2d():
    """
    correlate2d does 2d correlation of 'data' with 'kernel', storing
    the result in 'output'.

    supported 'mode's include:
        'nearest'   elements beyond boundary come from nearest edge pixel.
        'wrap'      elements beyond boundary come from the opposite array edge.
        'reflect'   elements beyond boundary come from reflection on same array edge.
        'constant'  elements beyond boundary are set to 'cval'

    If fft is True,  the correlation is performed using the FFT, else the
    correlation is performed using the naive approach.
    """
    a = np.arange(20*20)
    a = a.reshape((20,20))
    b = np.ones((5,5), dtype=np.float64)
    rn = correlate2d(a, b, fft=0)
    rf = correlate2d(a, b, fft=1)
    result = np.alltrue(np.ravel(rn-rf<1e-10))
    test = True
    assert_equal(result,test)

def test_boxcar1():
    """
    boxcar computes a 1D or 2D boxcar filter on every 1D or 2D subarray of data.

    'boxshape' is a tuple of integers specifying the dimensions of the filter: e.g. (3,3)

    if 'output' is specified, it should be the same shape as 'data' and
    None will be returned.

    supported 'mode's include:
        'nearest'   elements beyond boundary come from nearest edge pixel.
        'wrap'      elements beyond boundary come from the opposite array edge.
        'reflect'   elements beyond boundary come from reflection on same array edge.
        'constant'  elements beyond boundary are set to 'cval'
    """
    result = boxcar(np.array([10, 0, 0, 0, 0, 0, 1000]), (3,), mode="nearest").astype(np.longlong)
    test = np.array([  6,   3,   0,   0,   0, 333, 666], dtype=np.int64)
    assert_equal(result,test)

def test_boxcar2():
    result = boxcar(np.array([10, 0, 0, 0, 0, 0, 1000]), (3,), mode="wrap").astype(np.longlong)
    test = np.array([336,   3,   0,   0,   0, 333, 336], dtype=np.int64)
    assert_equal(result,test)

def test_boxcar3():
    result = boxcar(np.array([10, 0, 0, 0, 0, 0, 1000]), (3,), mode="reflect").astype(np.longlong)
    test = np.array([  6,   3,   0,   0,   0, 333, 666], dtype=np.int64)
    assert_equal(result,test)

def test_boxcar4():
    result = boxcar(np.array([10, 0, 0, 0, 0, 0, 1000]), (3,), mode="constant").astype(np.longlong)
    test = np.array([  3,   3,   0,   0,   0, 333, 333], dtype=np.int64)
    assert_equal(result,test)

def test_boxcar5():
    a = np.zeros((10,10))
    a[0,0] = 100
    a[5,5] = 1000
    a[9,9] = 10000
    result = boxcar(a, (3,3)).astype(np.longlong)
    test = np.array([[  44,   22,    0,    0,    0,    0,    0,    0,    0,    0],
           [  22,   11,    0,    0,    0,    0,    0,    0,    0,    0],
           [   0,    0,    0,    0,    0,    0,    0,    0,    0,    0],
           [   0,    0,    0,    0,    0,    0,    0,    0,    0,    0],
           [   0,    0,    0,    0,  111,  111,  111,    0,    0,    0],
           [   0,    0,    0,    0,  111,  111,  111,    0,    0,    0],
           [   0,    0,    0,    0,  111,  111,  111,    0,    0,    0],
           [   0,    0,    0,    0,    0,    0,    0,    0,    0,    0],
           [   0,    0,    0,    0,    0,    0,    0,    0, 1111, 2222],
           [   0,    0,    0,    0,    0,    0,    0,    0, 2222, 4444]], dtype=np.int64)
    assert_equal(result,test)

def test_boxcar6():
    a = np.zeros((10,10))
    a[0,0] = 100
    a[5,5] = 1000
    a[9,9] = 10000
    result = boxcar(a, (3,3), mode="wrap").astype(np.longlong)
    test = np.array([[1122,   11,    0,    0,    0,    0,    0,    0, 1111, 1122],
           [  11,   11,    0,    0,    0,    0,    0,    0,    0,   11],
           [   0,    0,    0,    0,    0,    0,    0,    0,    0,    0],
           [   0,    0,    0,    0,    0,    0,    0,    0,    0,    0],
           [   0,    0,    0,    0,  111,  111,  111,    0,    0,    0],
           [   0,    0,    0,    0,  111,  111,  111,    0,    0,    0],
           [   0,    0,    0,    0,  111,  111,  111,    0,    0,    0],
           [   0,    0,    0,    0,    0,    0,    0,    0,    0,    0],
           [1111,    0,    0,    0,    0,    0,    0,    0, 1111, 1111],
           [1122,   11,    0,    0,    0,    0,    0,    0, 1111, 1122]], dtype=np.int64)
    assert_equal(result,test)

def test_boxcar7():    
    a = np.zeros((10,10))
    a[0,0] = 100
    a[5,5] = 1000
    a[9,9] = 10000
    result = boxcar(a, (3,3), mode="reflect").astype(np.longlong)
    test = np.array([[  44,   22,    0,    0,    0,    0,    0,    0,    0,    0],
           [  22,   11,    0,    0,    0,    0,    0,    0,    0,    0],
           [   0,    0,    0,    0,    0,    0,    0,    0,    0,    0],
           [   0,    0,    0,    0,    0,    0,    0,    0,    0,    0],
           [   0,    0,    0,    0,  111,  111,  111,    0,    0,    0],
           [   0,    0,    0,    0,  111,  111,  111,    0,    0,    0],
           [   0,    0,    0,    0,  111,  111,  111,    0,    0,    0],
           [   0,    0,    0,    0,    0,    0,    0,    0,    0,    0],
           [   0,    0,    0,    0,    0,    0,    0,    0, 1111, 2222],
           [   0,    0,    0,    0,    0,    0,    0,    0, 2222, 4444]], dtype=np.int64)
    assert_equal(result,test)

def test_boxcar8():
    a = np.zeros((10,10))
    a[0,0] = 100
    a[5,5] = 1000
    a[9,9] = 10000
    result = boxcar(a, (3,3), mode="constant").astype(np.longlong)
    test = np.array([[  11,   11,    0,    0,    0,    0,    0,    0,    0,    0],
          [  11,   11,    0,    0,    0,    0,    0,    0,    0,    0],
          [   0,    0,    0,    0,    0,    0,    0,    0,    0,    0],
          [   0,    0,    0,    0,    0,    0,    0,    0,    0,    0],
          [   0,    0,    0,    0,  111,  111,  111,    0,    0,    0],
          [   0,    0,    0,    0,  111,  111,  111,    0,    0,    0],
          [   0,    0,    0,    0,  111,  111,  111,    0,    0,    0],
          [   0,    0,    0,    0,    0,    0,    0,    0,    0,    0],
          [   0,    0,    0,    0,    0,    0,    0,    0, 1111, 1111],
          [   0,    0,    0,    0,    0,    0,    0,    0, 1111, 1111]], dtype=np.int64)
    assert_equal(result,test)

def test_boxcar9():
    a = np.zeros((10,10))
    a[3:6,3:6] = 111
    result = boxcar(a, (3,3)).astype(np.longlong)
    test = np.array([[  0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
           [  0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
           [  0,   0,  12,  24,  37,  24,  12,   0,   0,   0],
           [  0,   0,  24,  49,  74,  49,  24,   0,   0,   0],
           [  0,   0,  37,  74, 111,  74,  37,   0,   0,   0],
           [  0,   0,  24,  49,  74,  49,  24,   0,   0,   0],
           [  0,   0,  12,  24,  37,  24,  12,   0,   0,   0],
           [  0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
           [  0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
           [  0,   0,   0,   0,   0,   0,   0,   0,   0,   0]], dtype=np.int64)
    assert_equal(result,test)

if __name__ == "__main__":
    run_module_suite()
