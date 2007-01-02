"""Miscellaneous functions for testing masked arrays and subclasses

:author: Pierre Gerard-Marchant
:contact: pierregm_at_uga_dot_edu
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant ($Author$)"
__version__ = "1.0"
__revision__ = "$Revision$"
__date__ = "$Date$"


import numpy as N
from numpy.core.numerictypes import float_
import numpy.core.umath as umath
from numpy.testing import NumpyTest, NumpyTestCase
from numpy.testing.utils import build_err_msg, rand

import core
from core import mask_or, getmask, getmaskarray, masked_array, nomask
from core import filled, equal, less

#------------------------------------------------------------------------------
def approx (a, b, fill_value=1, rtol=1.e-5, atol=1.e-8):
    """Returns true if all components of a and b are equal subject to given tolerances.
    If fill_value is 1, masked values considered equal.
    If fill_value is 0, masked values considered unequal.
    The relative error rtol should be positive and << 1.0
    The absolute error atol comes into play for those elements of b that are 
    very small or zero; it says how small a must be also.
    """
    m = mask_or(getmask(a), getmask(b))
    d1 = filled(a)
    d2 = filled(b)
    if d1.dtype.char == "O" or d2.dtype.char == "O":
        return N.equal(d1,d2).ravel()
    x = filled(masked_array(d1, copy=False, mask=m), fill_value).astype(float_)
    y = filled(masked_array(d2, copy=False, mask=m), 1).astype(float_)
    d = N.less_equal(umath.absolute(x-y), atol + rtol * umath.absolute(y))
    return d.ravel()
#............................
def assert_equal(actual,desired,err_msg=''):
    """Asserts that two items are equal.
    """
    if isinstance(desired, dict):
        assert isinstance(actual, dict), repr(type(actual))
        assert_equal(len(actual),len(desired),err_msg)
        for k,i in desired.items():
            assert actual.has_key(k), repr(k)
            assert_equal(actual[k], desired[k], 'key=%r\n%s' % (k,err_msg))
        return
    if isinstance(desired, (list,tuple)) and isinstance(actual, (list,tuple)):
        assert_equal(len(actual),len(desired),err_msg)
        for k in range(len(desired)):
            assert_equal(actual[k], desired[k], 'item=%r\n%s' % (k,err_msg))
        return
    from numpy.core import ndarray
    if isinstance(actual, ndarray) or isinstance(desired, ndarray):
        return assert_array_equal(actual, desired, err_msg)
    msg = build_err_msg([actual, desired], err_msg,)
    assert desired == actual, msg
#.............................
def fail_if_equal(actual,desired,err_msg='',):
    """Raises an assertion error if two items are equal.
    """
    if isinstance(desired, dict):
        assert isinstance(actual, dict), repr(type(actual))
        fail_if_equal(len(actual),len(desired),err_msg)
        for k,i in desired.items():
            assert actual.has_key(k), repr(k)
            fail_if_equal(actual[k], desired[k], 'key=%r\n%s' % (k,err_msg))
        return
    if isinstance(desired, (list,tuple)) and isinstance(actual, (list,tuple)):
        fail_if_equal(len(actual),len(desired),err_msg)
        for k in range(len(desired)):
            fail_if_equal(actual[k], desired[k], 'item=%r\n%s' % (k,err_msg))
        return
    if isinstance(actual, N.ndarray) or isinstance(desired, N.ndarray):
        return fail_if_array_equal(actual, desired, err_msg)
    msg = build_err_msg([actual, desired], err_msg)
    assert desired != actual, msg
assert_not_equal = fail_if_equal
#............................
def assert_almost_equal(actual,desired,decimal=7,err_msg=''):
    """Asserts that two items are almost equal.
    The test is equivalent to abs(desired-actual) < 0.5 * 10**(-decimal)
    """
    if isinstance(actual, N.ndarray) or isinstance(desired, N.ndarray):
        return assert_array_almost_equal(actual, desired, decimal, err_msg)
    msg = build_err_msg([actual, desired], err_msg)
    assert round(abs(desired - actual),decimal) == 0, msg
#............................
def assert_array_compare(comparison, x, y, err_msg='', header='', 
                         fill_value=True):
    """Asserts that a comparison relation between two masked arrays is satisfied
    elementwise."""
    xf = filled(x)
    yf = filled(y)
    m = mask_or(getmask(x), getmask(y))
    
    x = filled(masked_array(xf, copy=False, mask=m), fill_value)
    y = filled(masked_array(yf, copy=False, mask=m), fill_value)
    if (x.dtype.char != "O"):
        x = x.astype(float_)
        if isinstance(x, N.ndarray) and x.size > 1:
            x[N.isnan(x)] = 0
        elif N.isnan(x):
            x = 0
    if (y.dtype.char != "O"):
        y = y.astype(float_)
        if isinstance(y, N.ndarray) and y.size > 1:
            y[N.isnan(y)] = 0
        elif N.isnan(y):
            y = 0
    try:
        cond = (x.shape==() or y.shape==()) or x.shape == y.shape
        if not cond:
            msg = build_err_msg([x, y],
                                err_msg
                                + '\n(shapes %s, %s mismatch)' % (x.shape,
                                                                  y.shape),
                                header=header,
                                names=('x', 'y'))
            assert cond, msg
        val = comparison(x,y)
        if m is not nomask and fill_value:
            val = masked_array(val, mask=m, copy=False)
        if isinstance(val, bool):
            cond = val
            reduced = [0]
        else:
            reduced = val.ravel()
            cond = reduced.all()
            reduced = reduced.tolist()
        if not cond:
            match = 100-100.0*reduced.count(1)/len(reduced)
            msg = build_err_msg([x, y],
                                err_msg
                                + '\n(mismatch %s%%)' % (match,),
                                header=header,
                                names=('x', 'y'))
            assert cond, msg
    except ValueError:
        msg = build_err_msg([x, y], err_msg, header=header, names=('x', 'y'))
        raise ValueError(msg)
#............................
def assert_array_equal(x, y, err_msg=''):
    """Checks the elementwise equality of two masked arrays."""
    assert_array_compare(equal, x, y, err_msg=err_msg,
                         header='Arrays are not equal')
##............................
def fail_if_array_equal(x, y, err_msg=''):
    """Raises an assertion error if two masked arrays are not equal 
    (elem by elem.)"""
    def compare(x,y):
        
        return (not N.alltrue(approx(x, y)))
    assert_array_compare(compare, x, y, err_msg=err_msg,
                         header='Arrays are not equal')
#............................
def assert_array_almost_equal(x, y, decimal=6, err_msg=''):
    """Checks the elementwise equality of two masked arrays, up to a given 
    number of decimals."""
    def compare(x, y):
        return approx(x,y)
    assert_array_compare(compare, x, y, err_msg=err_msg, 
                         header='Arrays are not almost equal')
#............................
def assert_array_less(x, y, err_msg=''):
    assert_array_compare(less, x, y, err_msg=err_msg,
                         header='Arrays are not less-ordered')
#............................
assert_close = assert_almost_equal
#............................
def assert_mask_equal(m1, m2):
    """Asserts the equality of two masks."""
    if m1 is nomask:
        assert(m2 is nomask)
    if m2 is nomask:
        assert(m1 is nomask)
    assert_array_equal(m1, m2)