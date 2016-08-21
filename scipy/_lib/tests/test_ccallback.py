from __future__ import division, print_function, absolute_import

from numpy.testing import assert_equal, assert_raises

from scipy._lib import _test_ccallback, _test_ccallback_cython


ERROR_VALUE = 2.0


def callback_python(a, user_data=None):
    if a == ERROR_VALUE:
        raise ValueError("bad value")

    if user_data is None:
        return a + 1
    else:
        return a + user_data


def test_callbacks():
    callers = {
        'simple': _test_ccallback.call_simple,
        'nodata': _test_ccallback.call_nodata,
        'nonlocal': _test_ccallback.call_nonlocal
    }

    funcs = {
        'capsule': _test_ccallback.get_plus1_capsule(),
        'cython': (_test_ccallback_cython, "plus1_cython"),
        'python': callback_python,
        'ctypes': _test_ccallback_cython.plus1_ctypes,
    }

    def check(caller, func):
        caller = callers[caller]
        func = funcs[func]

        # Test basic call
        assert_equal(caller(func, 1.0), 2.0)

        # Test 'bad' value resulting to an error
        assert_raises(ValueError, caller, func, ERROR_VALUE)

        # Test passing in user_data
        if isinstance(func, tuple):
            func2 = func + (2.0,)
        else:
            func2 = (func, 2.0)
        assert_equal(caller(func2, 1.0), 3.0)

    for caller in callers.keys():
        for func in funcs.keys():
            yield check, caller, func

