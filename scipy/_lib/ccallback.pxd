# -*-cython-*-


from libc.setjmp cimport jmp_buf
from cpython.object cimport PyObject

cdef extern from "ccallback.h":
    ctypedef struct ccallback_signature_t:
        char *signature
        int value

    ctypedef struct ccallback_t:
        void *c_function
        PyObject *py_function
        void *user_data
        jmp_buf error_buf
        ccallback_t *prev_callback
        ccallback_signature_t *signature

        # Unused variables that can be used by the thunk etc. code for any purpose
        long info
        void *info_p

    int CCALLBACK_DEFAULTS
    int CCALLBACK_OBTAIN
    int CCALLBACK_PARSE

    ccallback_t *ccallback_obtain() nogil
    int ccallback_prepare(ccallback_t *callback, ccallback_signature_t *sigs,
                          object func, int flags) except -1
    void ccallback_release(ccallback_t *callback)
