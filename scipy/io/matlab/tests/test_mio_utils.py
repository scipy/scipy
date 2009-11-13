""" Testing 

"""

import numpy as np

from nose.tools import assert_true, assert_false, \
     assert_equal, assert_raises

from numpy.testing import assert_array_equal, assert_array_almost_equal

from scipy.io.matlab.mio_utils import cproduct, process_element, \
    FileReadOpts, chars_to_strings


def test_cproduct():
    yield assert_equal, cproduct(()), 1
    yield assert_equal, cproduct((1,)), 1
    yield assert_equal, cproduct((1,3)), 3
    yield assert_equal, cproduct([1,3]), 3


def test_process_element():
    a = np.zeros((1,3))
    fro = FileReadOpts(0, 0, 0)
    pa = process_element(a, fro, None)
    yield assert_array_equal, a, pa
    # squeeze
    fro = FileReadOpts(1, 0, 0)
    yield assert_array_equal, a, process_element(a, fro, None)
    fro = FileReadOpts(1, 1, 0)
    yield assert_array_equal, a, process_element(a, fro, None)
    fro = FileReadOpts(1, 1, 1)
    yield (assert_array_equal,
           np.squeeze(a),
           process_element(a, fro, None))
    # mat_dtype
    fro = FileReadOpts(0, 0, 0) # not applying mat_dtype
    mdt = np.dtype('i2')
    pa = process_element(a, fro, mdt)
    yield assert_array_equal, a, pa
    yield assert_equal, pa.dtype, a.dtype
    fro = FileReadOpts(0, 1, 0) # applying mat_dtype
    pa = process_element(a, fro, mdt)
    yield assert_array_equal, a, pa
    yield assert_equal, pa.dtype, mdt
    # chars as strings
    strings = ['learn ', 'python', 'fast  ', 'here  ']
    str_arr = np.array(strings, dtype='U6') # shape (4,)
    chars = [list(s) for s in strings]
    char_arr = np.array(chars, dtype='U1') # shape (4,6)
    fro = FileReadOpts(0, 0, 0)
    pac = process_element(char_arr, fro, None)
    yield assert_array_equal, char_arr, pac
    fro = FileReadOpts(1, 0, 0)
    pac = process_element(char_arr, fro, None)
    yield assert_array_equal, str_arr, pac
    
    
def test_chars_strings():
    # chars as strings
    strings = ['learn ', 'python', 'fast  ', 'here  ']
    str_arr = np.array(strings, dtype='U6') # shape (4,)
    chars = [list(s) for s in strings]
    char_arr = np.array(chars, dtype='U1') # shape (4,6)
    yield assert_array_equal, chars_to_strings(char_arr), str_arr
    ca2d = char_arr.reshape((2,2,6))
    sa2d = str_arr.reshape((2,2))
    yield assert_array_equal, chars_to_strings(ca2d), sa2d
    ca3d = char_arr.reshape((1,2,2,6))
    sa3d = str_arr.reshape((1,2,2))
    yield assert_array_equal, chars_to_strings(ca3d), sa3d
    # Fortran ordered arrays
    char_arrf = np.array(chars, dtype='U1', order='F') # shape (4,6)
    yield assert_array_equal, chars_to_strings(char_arrf), str_arr
    
    
