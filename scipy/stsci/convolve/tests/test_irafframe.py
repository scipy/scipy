#!/usr/bin/env python
import numpy as np
import nose
from scipy.stsci.convolve import *
from scipy.stsci.convolve.iraf_frame import *
from numpy.testing import *

def test_frame_nearest():
    """
    frame_nearest creates an oversized copy of 'a' with new 'shape'
    and the contents of 'a' in the center.  The boundary pixels are
    copied from the nearest edge pixel in 'a'.
    """    
    a = np.arange(16)
    a.shape=(4,4)
    result = frame_nearest(a, (8,8))
    test = np.array([[ 0,  0,  0,  1,  2,  3,  3,  3],
           [ 0,  0,  0,  1,  2,  3,  3,  3],
           [ 0,  0,  0,  1,  2,  3,  3,  3],
           [ 4,  4,  4,  5,  6,  7,  7,  7],
           [ 8,  8,  8,  9, 10, 11, 11, 11],
           [12, 12, 12, 13, 14, 15, 15, 15],
           [12, 12, 12, 13, 14, 15, 15, 15],
           [12, 12, 12, 13, 14, 15, 15, 15]])
    assert_equal(result,test)

def test_frame_reflect():
    """
    frame_reflect creates an oversized copy of 'a' with new 'shape'
    and the contents of 'a' in the center.  The boundary pixels are
    reflected from the nearest edge pixels in 'a'.
    """
    a = np.arange(16)
    a.shape = (4,4)
    result = frame_reflect(a, (8,8))
    test = np.array([[ 5,  4,  4,  5,  6,  7,  7,  6],
           [ 1,  0,  0,  1,  2,  3,  3,  2],
           [ 1,  0,  0,  1,  2,  3,  3,  2],
           [ 5,  4,  4,  5,  6,  7,  7,  6],
           [ 9,  8,  8,  9, 10, 11, 11, 10],
           [13, 12, 12, 13, 14, 15, 15, 14],
           [13, 12, 12, 13, 14, 15, 15, 14],
           [ 9,  8,  8,  9, 10, 11, 11, 10]])
    assert_equal(result,test)

def test_frame_wrap():
    """
    frame_wrap creates an oversized copy of 'a' with new 'shape'
    and the contents of 'a' in the center.  The boundary pixels are
    wrapped around to the opposite edge pixels in 'a'.
    """
    a = np.arange(16)
    a.shape=(4,4)
    result = frame_wrap(a, (8,8))
    test=np.array([[10, 11,  8,  9, 10, 11,  8,  9],
           [14, 15, 12, 13, 14, 15, 12, 13],
           [ 2,  3,  0,  1,  2,  3,  0,  1],
           [ 6,  7,  4,  5,  6,  7,  4,  5],
           [10, 11,  8,  9, 10, 11,  8,  9],
           [14, 15, 12, 13, 14, 15, 12, 13],
           [ 2,  3,  0,  1,  2,  3,  0,  1],
           [ 6,  7,  4,  5,  6,  7,  4,  5]])
    assert_equal(result,test)

def test_frame_constant():
    """
    frame_nearest creates an oversized copy of 'a' with new 'shape'
    and the contents of 'a' in the center.  The boundary pixels are
    copied from the nearest edge pixel in 'a'.
    """
    a = np.arange(16)
    a.shape=(4,4)
    result = frame_constant(a, (8,8), cval=42)
    test = np.array([[42, 42, 42, 42, 42, 42, 42, 42],
           [42, 42, 42, 42, 42, 42, 42, 42],
           [42, 42,  0,  1,  2,  3, 42, 42],
           [42, 42,  4,  5,  6,  7, 42, 42],
           [42, 42,  8,  9, 10, 11, 42, 42],
           [42, 42, 12, 13, 14, 15, 42, 42],
           [42, 42, 42, 42, 42, 42, 42, 42],
           [42, 42, 42, 42, 42, 42, 42, 42]])
    assert_equal(result,test)
    
if __name__ == "__main__":
    run_module_suite()

