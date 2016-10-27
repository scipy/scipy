cimport cython
from cython cimport parallel

cimport numpy as np
import numpy as np
from numpy.testing import assert_equal

cimport cython_special as sc

DEF NERR = 10


def test_cython_special_error_cfunc():
    # C function
    sc.spence(-1.0)
    assert_equal(sc.get_errno(), sc.DOMAIN)
    assert_equal(sc.get_errno(), sc.OK)


def test_cython_special_error_cyfunc():    
    # Cython function
    sc.loggamma(0)
    assert_equal(sc.get_errno(), sc.SINGULAR)
    assert_equal(sc.get_errno(), sc.OK)


def test_cython_special_error_cppfunc():
    # C++ function
    sc.wrightomega(-1000)
    assert_equal(sc.get_errno(), sc.UNDERFLOW)
    assert_equal(sc.get_errno(), sc.OK)


@cython.boundscheck(False)
def test_cython_special_error_parallel():
    # Initialize to a value sf_error will never return
    errors = -np.ones(NERR, dtype=np.int64)

    cdef int n, tid
    cdef np.npy_int64 [:] errors_buf = errors

    for n in parallel.prange(NERR, nogil=True, num_threads=NERR):
        tid = parallel.threadid()
        sc._error_test_function(tid)
        errors_buf[tid] = <int>sc.get_errno()

    for n in range(NERR):
        assert_equal(errors_buf[n], n)


def test_cython_special_error_serial():
    cdef int n

    for n in range(10):
        sc._error_test_function(n)
        assert_equal(<int>sc.get_errno(), n)
