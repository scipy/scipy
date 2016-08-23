from __future__ import division, print_function, absolute_import

from numpy.testing import assert_equal, assert_raises

import nose
import ctypes
from scipy._lib import _test_ccallback, _test_ccallback_cython

try:
    import cffi
    HAVE_CFFI = True
except ImportError:
    HAVE_CFFI = False


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

    def _get_cffi_func():
        if not HAVE_CFFI:
            raise nose.SkipTest("cffi not installed")

        # Get function address
        voidp = ctypes.cast(_test_ccallback_cython.plus1_ctypes,
                            ctypes.c_void_p)
        address = voidp.value

        # Create corresponding cffi handle
        ffi = cffi.FFI()
        func = ffi.cast('double (*)(double, int *, void *)', address)
        return func

    funcs = {
        'capsule': lambda: _test_ccallback.get_plus1_capsule(),
        'cython': lambda: (_test_ccallback_cython, "plus1_cython"),
        'python': lambda: callback_python,
        'ctypes': lambda: _test_ccallback_cython.plus1_ctypes,
        'cffi': _get_cffi_func
    }

    def check(caller, func):
        caller = callers[caller]
        func = funcs[func]()

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

