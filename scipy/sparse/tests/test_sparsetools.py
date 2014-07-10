from __future__ import division, print_function, absolute_import

import sys
import os
import gc
import re
import time
import threading

from nose import SkipTest
import numpy as np
from numpy.testing import assert_raises, assert_equal, dec, run_module_suite, assert_
from scipy.sparse import (_sparsetools, coo_matrix, csr_matrix, csc_matrix,
                          bsr_matrix, dia_matrix)
from scipy.sparse.sputils import supported_dtypes
from scipy.lib.decorator import decorator


@decorator
def xslow(func, *a, **kw):
    try:
        v = int(os.environ.get('SCIPY_XSLOW', '0'))
        if not v:
            raise ValueError()
    except ValueError:
        raise SkipTest("very slow test; set environment variable "
                       "SCIPY_XSLOW=1 to run it")
    return func(*a, **kw)


def test_exception():
    assert_raises(MemoryError, _sparsetools.test_throw_error)


def test_threads():
    # Smoke test for parallel threaded execution; doesn't actually
    # check that code runs in parallel, but just that it produces
    # expected results.
    nthreads = 10
    niter = 100

    n = 20
    a = csr_matrix(np.ones([n, n]))
    bres = []

    class Worker(threading.Thread):
        def run(self):
            b = a.copy()
            for j in range(niter):
                _sparsetools.csr_plus_csr(n, n,
                                          a.indptr, a.indices, a.data,
                                          a.indptr, a.indices, a.data,
                                          b.indptr, b.indices, b.data)
            bres.append(b)

    threads = [Worker() for _ in range(nthreads)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

    for b in bres:
        assert_(np.all(b.toarray() == 2))


def test_regression_std_vector_dtypes():
    # Regression test for gh-3780, checking the std::vector typemaps
    # in sparsetools.cxx are complete.
    for dtype in supported_dtypes:
        ad = np.matrix([[1, 2], [3, 4]]).astype(dtype)
        a = csr_matrix(ad, dtype=dtype)

        # getcol is one function using std::vector typemaps, and should not fail
        assert_equal(a.getcol(0).todense(), ad[:,0])


class TestInt32Overflow(object):
    """
    Some of the sparsetools routines use dense 2D matrices whose
    total size is not bounded by the nnz of the sparse matrix. These
    routines used to suffer from int32 wraparounds; here, we try to
    check that the wraparounds don't occur any more.
    """

    def __init__(self):
        # choose n large enough
        self.n = 50000
        assert self.n**2 > np.iinfo(np.int32).max

    @dec.skipif(not sys.platform.startswith('linux'), "test requires Linux")
    @dec.skipif(np.dtype(np.intp).itemsize < 8, "test requires 64-bit system")
    def setUp(self):
        check_free_memory(5000)

    def tearDown(self):
        gc.collect()

    def test_coo_todense(self):
        # Check *_todense routines (cf. gh-2179)
        #
        # All of them in the end call coo_matrix.todense

        n = self.n

        i = np.array([0, n-1])
        j = np.array([0, n-1])
        data = np.array([1, 2], dtype=np.int8)
        m = coo_matrix((data, (i, j)))

        r = m.todense()
        assert_equal(r[0,0], 1)
        assert_equal(r[-1,-1], 2)
        del r
        gc.collect()

    @dec.slow
    def test_matvecs(self):
        # Check *_matvecs routines
        n = self.n

        i = np.array([0, n-1])
        j = np.array([0, n-1])
        data = np.array([1, 2], dtype=np.int8)
        m = coo_matrix((data, (i, j)))

        b = np.ones((n, n), dtype=np.int8)
        for sptype in (csr_matrix, csc_matrix, bsr_matrix):
            m2 = sptype(m)
            r = m2.dot(b)
            assert_equal(r[0,0], 1)
            assert_equal(r[-1,-1], 2)
            del r
            gc.collect()

        del b
        gc.collect()

    @dec.slow
    def test_dia_matvec(self):
        # Check: huge dia_matrix _matvec
        n = self.n
        data = np.ones((n, n), dtype=np.int8)
        offsets = np.arange(n)
        m = dia_matrix((data, offsets), shape=(n, n))
        v = np.ones(m.shape[1], dtype=np.int8)
        r = m.dot(v)
        assert_equal(r[0], np.int8(n))
        del data, offsets, m, v, r
        gc.collect()

    @dec.slow
    def test_bsr_1_block(self):
        # Check: huge bsr_matrix (1-block)
        #
        # The point here is that indices inside a block may overflow.

        def check(op):
            n = self.n
            data = np.ones((1, n, n), dtype=np.int8)
            indptr = np.array([0, 1], dtype=np.int32)
            indices = np.array([0], dtype=np.int32)
            m = bsr_matrix((data, indices, indptr), blocksize=(n, n), copy=False)
            del data, indptr, indices
            getattr(self, "_check_bsr_" + op)(m)
            del m
            gc.collect()

        for op in ("matmat", "matvecs", "matvec", "diagonal",
                   "sort_indices", "transpose"):
            yield check, op

    @dec.slow
    def test_bsr_n_block(self):
        # Check: huge bsr_matrix (n-block)
        #
        # The point here is that while indices within a block don't
        # overflow, accumulators across many block may.

        def check(op):
            n = self.n
            data = np.ones((n, n, 1), dtype=np.int8)
            indptr = np.array([0, n], dtype=np.int32)
            indices = np.arange(n, dtype=np.int32)
            m = bsr_matrix((data, indices, indptr), blocksize=(n, 1), copy=False)
            del data, indptr, indices
            getattr(self, "_check_bsr_" + op)(m)
            del m
            gc.collect()

        for op in ("matmat", "matvecs", "matvec", "diagonal",
                   "sort_indices", "transpose"):
            yield check, op

    @xslow
    def _check_bsr_matvecs(self, m):
        n = self.n

        # _matvecs
        r = m.dot(np.ones((n, 2), dtype=np.int8))
        assert_equal(r[0,0], np.int8(n))
        del r
        gc.collect()

    def _check_bsr_matvec(self, m):
        n = self.n

        # _matvec
        r = m.dot(np.ones((n,), dtype=np.int8))
        assert_equal(r[0], np.int8(n))
        del r
        gc.collect()

    def _check_bsr_diagonal(self, m):
        n = self.n

        # _diagonal
        r = m.diagonal()
        assert_equal(r, np.ones(n))
        del r
        gc.collect()

    def _check_bsr_sort_indices(self, m):
        # _sort_indices
        m.sort_indices()

    @xslow
    def _check_bsr_transpose(self, m):
        # _transpose
        m.transpose()

    @xslow
    def _check_bsr_matmat(self, m):
        n = self.n

        # _bsr_matmat
        m2 = bsr_matrix(np.ones((n, 2), dtype=np.int8), blocksize=(m.blocksize[1], 2))
        m.dot(m2)  # shouldn't SIGSEGV

        # _bsr_matmat
        m2 = bsr_matrix(np.ones((2, n), dtype=np.int8), blocksize=(2, m.blocksize[0]))
        m2.dot(m)  # shouldn't SIGSEGV

    @dec.slow
    def test_csr_matmat(self):
        n = self.n

        # cf. gh-3212
        a = csr_matrix(np.ones((n, 1), dtype=np.int8))
        b = csr_matrix(np.ones((1, n), dtype=np.int8))
        assert_raises(RuntimeError, a.dot, b)


@dec.skipif(True, "64-bit indices in sparse matrices not available")
def test_csr_matmat_int64_overflow():
    n = 3037000500
    assert n**2 > np.iinfo(np.int64).max

    # the test would take crazy amounts of memory
    check_free_memory(n * (8*2 + 1) * 3 / 1e6)

    # int64 overflow
    data = np.ones((n,), dtype=np.int8)
    indptr = np.arange(n+1, dtype=np.int64)
    indices = np.zeros(n, dtype=np.int64)
    a = csr_matrix((data, indices, indptr))
    b = a.T

    assert_raises(RuntimeError, a.dot, b)


def check_free_memory(free_mb):
    meminfo = get_mem_info_linux()

    if meminfo['memfree'] + meminfo['cached'] < free_mb * 1e6:
        raise SkipTest("test requires %d MB of free memory" % (int(free_mb),))


def get_mem_info_linux():
    """
    Get information about available memory.

    Returns a dict of items in /proc/meminfo
    """
    info = {}
    with open('/proc/meminfo', 'r') as f:
        for line in f:
            p = line.split()
            info[p[0].strip(':').lower()] = float(p[1]) * 1e3
    return info


def get_own_memusage_linux():
    """
    Return the memory usage of the current process (in bytes)
    """
    with open('/proc/%d/status' % os.getpid(), 'r') as f:
        procdata = f.read()

    m = re.search('VmRSS:\s*(\d+)\s*kB', procdata, re.S | re.I)
    if m is not None:
        memusage = float(m.group(1)) * 1e3
        return memusage

    return np.nan


if __name__ == "__main__":
    run_module_suite()
