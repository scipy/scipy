''' Tests for byteorder module '''

import sys

import numpy as np

from numpy.testing import assert_raises

import scipy.io.matlab.byteordercodes as sibc

def test_native():
    native_is_le = sys.byteorder == 'little'
    assert sibc.sys_is_le == native_is_le

def test_to_numpy():
    if sys.byteorder == 'little':
        assert sibc.to_numpy_code('native') == '<'
        assert sibc.to_numpy_code('swapped') == '>'
    else:
        assert sibc.to_numpy_code('native') == '>'
        assert sibc.to_numpy_code('swapped') == '<'
    assert sibc.to_numpy_code('native') == sibc.to_numpy_code('=')
    assert sibc.to_numpy_code('big') == '>'
    for code in ('little', '<', 'l', 'L', 'le'):
        assert sibc.to_numpy_code(code) == '<'
    for code in ('big', '>', 'b', 'B', 'be'):
        assert sibc.to_numpy_code(code) == '>'
    assert_raises(ValueError, sibc.to_numpy_code, 'silly string')
