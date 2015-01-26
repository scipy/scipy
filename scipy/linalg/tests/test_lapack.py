#!/usr/bin/env python
#
# Created by: Pearu Peterson, September 2002
#

from __future__ import division, print_function, absolute_import

from numpy.testing import TestCase, run_module_suite, assert_equal, \
    assert_array_almost_equal, assert_, assert_raises, assert_allclose

import numpy as np

from scipy.linalg import _flapack as flapack
from scipy.linalg import inv
try:
    from scipy.linalg import _clapack as clapack
except ImportError:
    clapack = None
from scipy.linalg.lapack import get_lapack_funcs

REAL_DTYPES = [np.float32, np.float64]
COMPLEX_DTYPES = [np.complex64, np.complex128]
DTYPES = REAL_DTYPES + COMPLEX_DTYPES


class TestFlapackSimple(TestCase):

    def test_gebal(self):
        a = [[1,2,3],[4,5,6],[7,8,9]]
        a1 = [[1,0,0,3e-4],
              [4,0,0,2e-3],
              [7,1,0,0],
              [0,1,0,0]]
        for p in 'sdzc':
            f = getattr(flapack,p+'gebal',None)
            if f is None:
                continue
            ba,lo,hi,pivscale,info = f(a)
            assert_(not info,repr(info))
            assert_array_almost_equal(ba,a)
            assert_equal((lo,hi),(0,len(a[0])-1))
            assert_array_almost_equal(pivscale, np.ones(len(a)))

            ba,lo,hi,pivscale,info = f(a1,permute=1,scale=1)
            assert_(not info,repr(info))
            # print a1
            # print ba,lo,hi,pivscale

    def test_gehrd(self):
        a = [[-149, -50,-154],
             [537, 180, 546],
             [-27, -9, -25]]
        for p in 'd':
            f = getattr(flapack,p+'gehrd',None)
            if f is None:
                continue
            ht,tau,info = f(a)
            assert_(not info,repr(info))

    def test_trsyl(self):
        a = np.array([[1, 2], [0, 4]])
        b = np.array([[5, 6], [0, 8]])
        c = np.array([[9, 10], [11, 12]])
        trans = 'T'

        # Test single and double implementations, including most
        # of the options
        for dtype in 'fdFD':
            a1, b1, c1 = a.astype(dtype), b.astype(dtype), c.astype(dtype)
            trsyl, = get_lapack_funcs(('trsyl',), (a1,))
            if dtype.isupper():  # is complex dtype
                a1[0] += 1j
                trans = 'C'

            x, scale, info = trsyl(a1, b1, c1)
            assert_array_almost_equal(np.dot(a1, x) + np.dot(x, b1), scale * c1)

            x, scale, info = trsyl(a1, b1, c1, trana=trans, tranb=trans)
            assert_array_almost_equal(np.dot(a1.conjugate().T, x) + np.dot(x, b1.conjugate().T),
                scale * c1, decimal=4)

            x, scale, info = trsyl(a1, b1, c1, isgn=-1)
            assert_array_almost_equal(np.dot(a1, x) - np.dot(x, b1), scale * c1, decimal=4)


class TestLapack(TestCase):

    def test_flapack(self):
        if hasattr(flapack,'empty_module'):
            # flapack module is empty
            pass

    def test_clapack(self):
        if hasattr(clapack,'empty_module'):
            # clapack module is empty
            pass

class TestLeastSquaresSolvers(TestCase):

    def test_gelsd(self):
        for dtype in REAL_DTYPES:
            a1 = np.array([[1.0,2.0],
                          [4.0,5.0],
                          [7.0,8.0]], dtype=dtype)
            b1 = np.array([16.0, 17.0, 20.0], dtype=dtype)
            gelsd, = get_lapack_funcs(('gelsd',), (a1, b1))
            #x,s,rank,work,iwork,info = *gelsd(a,b,lwork,size_iwork,[cond,overwrite_a,overwrite_b])
            x, s, rank, work, iwork, info = gelsd(a1, b1, -1, 1)  # Request of sizes
            lwork = work[0].real.astype(np.int)
            iwork_size = iwork[0].real.astype(np.int)

            x, s, rank, work, iwork, info = gelsd(a1, b1, lwork, iwork_size, -1, False, False)
            assert_allclose(x[:-1], np.array([-14.333333333333323, 14.999999999999991], dtype=dtype), rtol=10*np.finfo(dtype).eps)
            assert_allclose(s, np.array([12.596017180511966, 0.583396253199685], dtype=dtype), rtol=10*np.finfo(dtype).eps)

        for dtype in COMPLEX_DTYPES:
            a1 = np.array([[1.0+4.0j,2.0],
                          [4.0+0.5j,5.0-3.0j],
                          [7.0-2.0j,8.0+0.7j]], dtype=dtype)
            b1 = np.array([16.0, 17.0+2.0j, 20.0-4.0j], dtype=dtype)
            gelsd, = get_lapack_funcs(('gelsd',), (a1, b1))
            #x,s,rank,work,rwork,iwork,info = cgelsd(a,b,lwork,size_rwork,size_iwork,[cond,overwrite_a,overwrite_b])
            x, s, rank, work, rwork, iwork, info = gelsd(a1, b1, -1, 1,1)  # Request of sizes
            lwork = work[0].real.astype(np.int)
            rwork_size = rwork[0].real.astype(np.int)
            iwork_size = iwork[0].real.astype(np.int)

            x, s, rank, work, rwork, iwork, info = gelsd(a1, b1, lwork, rwork_size, iwork_size, -1, False, False)
            assert_allclose(x[:-1], np.array([1.161753632288328-1.901075709391912j, 1.735882340522193+1.521240901196909j],
                            dtype=dtype), rtol=10*np.finfo(dtype).eps)
            assert_allclose(s, np.array([13.035514762572043, 4.337666985231382], dtype=dtype), rtol=10*np.finfo(dtype).eps)

    def test_gelss(self):

        for dtype in REAL_DTYPES:
            a1 = np.array([[1.0,2.0],
                          [4.0,5.0],
                          [7.0,8.0]], dtype=dtype)
            b1 = np.array([16.0, 17.0, 20.0], dtype=dtype)
            gelss, = get_lapack_funcs(('gelss',), (a1, b1))
            #v,x,s,rank,work,info = dgelss(a,b,[cond,lwork,overwrite_a,overwrite_b])
            v,x,s,rank,work,info = gelss(a1, b1,-1, -1, False, False)  # Request of sizes
            lwork = work[0].real.astype(np.int)

            v,x,s,rank,work,info = gelss(a1, b1,-1,lwork, False, False)
            assert_allclose(x[:-1], np.array([-14.333333333333323, 14.999999999999991], dtype=dtype),rtol=10*np.finfo(dtype).eps)
            assert_allclose(s, np.array([12.596017180511966, 0.583396253199685], dtype=dtype), rtol=10*np.finfo(dtype).eps)

        for dtype in COMPLEX_DTYPES:
            a1 = np.array([[1.0+4.0j,2.0],
                          [4.0+0.5j,5.0-3.0j],
                          [7.0-2.0j,8.0+0.7j]], dtype=dtype)
            b1 = np.array([16.0, 17.0+2.0j, 20.0-4.0j], dtype=dtype)
            gelss, = get_lapack_funcs(('gelss',), (a1, b1))
            #v,x,s,rank,work,info = cgelss(a,b,[cond,lwork,overwrite_a,overwrite_b])
            v, x, s, rank, work, info = gelss(a1, b1,-1, -1, False, False)  # Request of sizes
            lwork = work[0].real.astype(np.int)

            v,x,s,rank,work,info = gelss(a1, b1,-1,lwork, False, False)
            assert_allclose(x[:-1], np.array([1.161753632288328-1.901075709391912j, 1.735882340522193+1.521240901196909j],
                            dtype=dtype), rtol=10*np.finfo(dtype).eps)
            assert_allclose(s, np.array([13.035514762572043, 4.337666985231382], dtype=dtype), rtol=10*np.finfo(dtype).eps)

    def test_gelsy(self):

        for dtype in REAL_DTYPES:
            a1 = np.array([[1.0,2.0],
                          [4.0,5.0],
                          [7.0,8.0]], dtype=dtype)
            b1 = np.array([16.0, 17.0, 20.0], dtype=dtype)
            gelsy, = get_lapack_funcs(('gelsy',), (a1, b1))
            #x,s,rank,work,iwork,info = *gelsd(a,b,lwork,size_iwork,[cond,overwrite_a,overwrite_b])
            jptv = np.zeros((a1.shape[1],1), dtype=np.int32)
            v,x,j,rank,work,info = gelsy(a1, b1, jptv, np.finfo(dtype).eps,-1)  # Request of sizes
            lwork = work[0].real.astype(np.int)

            jptv = np.zeros((a1.shape[1],1), dtype=np.int32)
            v,x,j,rank,work,info = gelsy(a1, b1, jptv, np.finfo(dtype).eps, lwork, False, False)
            assert_allclose(x[:-1], np.array([-14.333333333333323, 14.999999999999991], dtype=dtype),rtol=10*np.finfo(dtype).eps)
            #assert_allclose(s, np.array([ 0.66666669,], dtype=dtype ) )

        for dtype in COMPLEX_DTYPES:
            a1 = np.array([[1.0+4.0j,2.0],
                          [4.0+0.5j,5.0-3.0j],
                          [7.0-2.0j,8.0+0.7j]], dtype=dtype)
            b1 = np.array([16.0, 17.0+2.0j, 20.0-4.0j], dtype=dtype)
            gelsy, = get_lapack_funcs(('gelsy',), (a1, b1))
            #x,s,rank,work,rwork,iwork,info = cgelsd(a,b,lwork,size_rwork,size_iwork,[cond,overwrite_a,overwrite_b])
            jptv = np.zeros((a1.shape[1],1), dtype=np.int32)
            v,x,j,rank,work,rwork,info = gelsy(a1, b1, jptv, np.finfo(dtype).eps,-1)  # Request of sizes
            lwork = work[0].real.astype(np.int)

            jptv = np.zeros((a1.shape[1],1), dtype=np.int32)
            v,x,j,rank,work,rwork,info = gelsy(a1, b1, jptv, np.finfo(dtype).eps, lwork, False, False)
            assert_allclose(x[:-1], np.array([1.161753632288328-1.901075709391912j, 1.735882340522193+1.521240901196909j],
                            dtype=dtype), rtol=10*np.finfo(dtype).eps)
            #assert_allclose( s, np.array([ 106.12267169], dtype=dtype ) )

class TestRegression(TestCase):

    def test_ticket_1645(self):
        # Check that RQ routines have correct lwork
        for dtype in DTYPES:
            a = np.zeros((300, 2), dtype=dtype)

            gerqf, = get_lapack_funcs(['gerqf'], [a])
            assert_raises(Exception, gerqf, a, lwork=2)
            rq, tau, work, info = gerqf(a)

            if dtype in REAL_DTYPES:
                orgrq, = get_lapack_funcs(['orgrq'], [a])
                assert_raises(Exception, orgrq, rq[-2:], tau, lwork=1)
                orgrq(rq[-2:], tau, lwork=2)
            elif dtype in COMPLEX_DTYPES:
                ungrq, = get_lapack_funcs(['ungrq'], [a])
                assert_raises(Exception, ungrq, rq[-2:], tau, lwork=1)
                ungrq(rq[-2:], tau, lwork=2)


class TestDpotr(TestCase):
    def test_gh_2691(self):
        # 'lower' argument of dportf/dpotri
        for lower in [True, False]:
            for clean in [True, False]:
                np.random.seed(42)
                x = np.random.normal(size=(3, 3))
                a = x.dot(x.T)

                dpotrf, dpotri = get_lapack_funcs(("potrf", "potri"), (a, ))

                c, info = dpotrf(a, lower, clean=clean)
                dpt = dpotri(c, lower)[0]

                if lower:
                    assert_allclose(np.tril(dpt), np.tril(inv(a)))
                else:
                    assert_allclose(np.triu(dpt), np.triu(inv(a)))

if __name__ == "__main__":
    run_module_suite()
