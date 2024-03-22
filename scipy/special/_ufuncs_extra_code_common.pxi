cimport numpy as np
from numpy cimport (
    npy_float, npy_double, npy_longdouble,
    npy_cfloat, npy_cdouble, npy_clongdouble,
    npy_int, npy_long,
    NPY_FLOAT, NPY_DOUBLE, NPY_LONGDOUBLE,
    NPY_CFLOAT, NPY_CDOUBLE, NPY_CLONGDOUBLE,
    NPY_INT, NPY_LONG)

import warnings

ctypedef double complex double_complex

cdef extern from "numpy/ufuncobject.h":
    int PyUFunc_getfperr() nogil

cdef public int wrap_PyUFunc_getfperr() noexcept nogil:
    """
    Call PyUFunc_getfperr in a context where PyUFunc_API array is initialized;
    this avoids messing with the UNIQUE_SYMBOL #defines
    """
    return PyUFunc_getfperr()

cimport libc

from . cimport sf_error

np.import_array()
np.import_ufunc()

cdef void _set_action(sf_error.sf_error_t code,
                      sf_error.sf_action_t action) noexcept nogil:
    sf_error.set_action(code, action)

cdef void ufunc_warnorraise(sf_error.sf_action_t action, const char *msg) except *:
    from scipy.special import SpecialFunctionWarning, SpecialFunctionError

    if (action == sf_error.WARN):
        warnings.warn(msg, SpecialFunctionWarning)
    elif (action == sf_error.RAISE):
        raise SpecialFunctionError(msg)

sf_error.sf_error_set_callback(ufunc_warnorraise)

cdef int ufunc_getfpe() noexcept:
    return wrap_PyUFunc_getfperr()

sf_error.sf_error_set_callback_fpe(ufunc_getfpe)
