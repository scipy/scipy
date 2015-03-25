#!/usr/bin/env python
#
# Created by: Pearu Peterson, September 2002
#

from __future__ import division, print_function, absolute_import

from numpy.testing import TestCase, run_module_suite, assert_equal, \
    assert_array_almost_equal, assert_, assert_raises, assert_allclose, \
    assert_almost_equal

import numpy as np

from scipy.linalg import _flapack as flapack
from scipy.linalg import inv
from scipy.linalg import svd

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

    def test_lange(self):
        a = np.array([
            [-149, -50,-154],
            [537, 180, 546],
            [-27, -9, -25]])

        for dtype in 'fdFD':
            for norm in 'Mm1OoIiFfEe':
                a1 = a.astype(dtype)
                if dtype.isupper():
                    # is complex dtype
                    a1[0,0] += 1j

                lange, = get_lapack_funcs(('lange',), (a1,))
                value = lange(norm, a1)

                if norm in 'FfEe':
                    if dtype in 'Ff':
                        decimal = 3
                    else:
                        decimal = 7
                    ref = np.sqrt(np.sum(np.square(np.abs(a1))))
                    assert_almost_equal(value, ref, decimal)
                else:
                    if norm in 'Mm':
                        ref = np.max(np.abs(a1))
                    elif norm in '1Oo':
                        ref = np.max(np.sum(np.abs(a1), axis=0))
                    elif norm in 'Ii':
                        ref = np.max(np.sum(np.abs(a1), axis=1))

                    assert_equal(value, ref)


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
            gelsd, gelsd_lwork = get_lapack_funcs(('gelsd','gelsd_lwork'),
                                                  (a1, b1))

            m, n = a1.shape
            if len(b1.shape) == 2:
                nrhs = b1.shape[1]
            else:
                nrhs = 1

            # Request of sizes
            work,iwork,info = gelsd_lwork(m,n,nrhs,-1)
            lwork = np.int(np.real(work))
            iwork_size = iwork

            x, s, rank, info = gelsd(a1, b1, lwork, iwork_size,
                                    -1, False, False)
            assert_allclose(x[:-1], np.array([-14.333333333333323,
                                            14.999999999999991], dtype=dtype),
                                            rtol=25*np.finfo(dtype).eps)
            assert_allclose(s, np.array([12.596017180511966,
                                         0.583396253199685], dtype=dtype),
                                         rtol=25*np.finfo(dtype).eps)

        for dtype in COMPLEX_DTYPES:
            a1 = np.array([[1.0+4.0j,2.0],
                          [4.0+0.5j,5.0-3.0j],
                          [7.0-2.0j,8.0+0.7j]], dtype=dtype)
            b1 = np.array([16.0, 17.0+2.0j, 20.0-4.0j], dtype=dtype)
            gelsd, gelsd_lwork = get_lapack_funcs(('gelsd','gelsd_lwork'),
                                                  (a1, b1))

            m, n = a1.shape
            if len(b1.shape) == 2:
                nrhs = b1.shape[1]
            else:
                nrhs = 1

            # Request of sizes
            work, rwork, iwork, info = gelsd_lwork(m,n,nrhs,-1)
            lwork = np.int(np.real(work))
            rwork_size = np.int(rwork)
            iwork_size = iwork

            x, s, rank, info = gelsd(a1, b1, lwork, rwork_size, iwork_size,
                                     -1, False, False)
            assert_allclose(x[:-1],
                            np.array([1.161753632288328-1.901075709391912j,
                                      1.735882340522193+1.521240901196909j],
                            dtype=dtype), rtol=25*np.finfo(dtype).eps)
            assert_allclose(s,
                            np.array([13.035514762572043, 4.337666985231382],
                                     dtype=dtype), rtol=25*np.finfo(dtype).eps)

    def test_gelss(self):

        for dtype in REAL_DTYPES:
            a1 = np.array([[1.0,2.0],
                          [4.0,5.0],
                          [7.0,8.0]], dtype=dtype)
            b1 = np.array([16.0, 17.0, 20.0], dtype=dtype)
            gelss, gelss_lwork = get_lapack_funcs(('gelss','gelss_lwork'),
                                                  (a1, b1))

            m, n = a1.shape
            if len(b1.shape) == 2:
                nrhs = b1.shape[1]
            else:
                nrhs = 1

            # Request of sizes
            work,info = gelss_lwork(m,n,nrhs,-1)
            lwork = np.int(np.real(work))

            v,x,s,rank,work,info = gelss(a1, b1,-1,lwork, False, False)
            assert_allclose(x[:-1], np.array([-14.333333333333323,
                            14.999999999999991], dtype=dtype),
                            rtol=25*np.finfo(dtype).eps)
            assert_allclose(s, np.array([12.596017180511966,
                                         0.583396253199685], dtype=dtype),
                                         rtol=25*np.finfo(dtype).eps)

        for dtype in COMPLEX_DTYPES:
            a1 = np.array([[1.0+4.0j,2.0],
                          [4.0+0.5j,5.0-3.0j],
                          [7.0-2.0j,8.0+0.7j]], dtype=dtype)
            b1 = np.array([16.0, 17.0+2.0j, 20.0-4.0j], dtype=dtype)
            gelss, gelss_lwork = get_lapack_funcs(('gelss','gelss_lwork'),
                                                  (a1, b1))

            m, n = a1.shape
            if len(b1.shape) == 2:
                nrhs = b1.shape[1]
            else:
                nrhs = 1

            # Request of sizes
            work,info = gelss_lwork(m,n,nrhs,-1)
            lwork = np.int(np.real(work))

            v,x,s,rank,work,info = gelss(a1, b1,-1,lwork, False, False)
            assert_allclose(x[:-1],
                            np.array([1.161753632288328-1.901075709391912j,
                                      1.735882340522193+1.521240901196909j],
                            dtype=dtype), rtol=25*np.finfo(dtype).eps)
            assert_allclose(s, np.array([13.035514762572043,
                                         4.337666985231382], dtype=dtype),
                                         rtol=25*np.finfo(dtype).eps)

    def test_gelsy(self):

        for dtype in REAL_DTYPES:
            a1 = np.array([[1.0,2.0],
                          [4.0,5.0],
                          [7.0,8.0]], dtype=dtype)
            b1 = np.array([16.0, 17.0, 20.0], dtype=dtype)
            gelsy, gelsy_lwork = get_lapack_funcs(('gelsy','gelss_lwork'), (a1, b1))

            m, n = a1.shape
            if len(b1.shape) == 2:
                nrhs = b1.shape[1]
            else:
                nrhs = 1

            # Request of sizes
            work, info = gelsy_lwork(m,n,nrhs,10*np.finfo(dtype).eps)
            lwork = np.int(np.real(work))

            jptv = np.zeros((a1.shape[1],1), dtype=np.int32)
            v, x, j, rank, info = gelsy(a1, b1, jptv, np.finfo(dtype).eps,
                                        lwork, False, False)
            assert_allclose(x[:-1], np.array([-14.333333333333323,
                                            14.999999999999991], dtype=dtype),
                                            rtol=25*np.finfo(dtype).eps)

        for dtype in COMPLEX_DTYPES:
            a1 = np.array([[1.0+4.0j,2.0],
                          [4.0+0.5j,5.0-3.0j],
                          [7.0-2.0j,8.0+0.7j]], dtype=dtype)
            b1 = np.array([16.0, 17.0+2.0j, 20.0-4.0j], dtype=dtype)
            gelsy, gelsy_lwork = get_lapack_funcs(('gelsy','gelss_lwork'), (a1, b1))

            m, n = a1.shape
            if len(b1.shape) == 2:
                nrhs = b1.shape[1]
            else:
                nrhs = 1

            # Request of sizes
            work, info = gelsy_lwork(m,n,nrhs,10*np.finfo(dtype).eps)
            lwork = np.int(np.real(work))

            jptv = np.zeros((a1.shape[1],1), dtype=np.int32)
            v, x, j, rank, info = gelsy(a1, b1, jptv, np.finfo(dtype).eps,
                                        lwork, False, False)
            assert_allclose(x[:-1],
                            np.array([1.161753632288328-1.901075709391912j,
                                      1.735882340522193+1.521240901196909j],
                            dtype=dtype), rtol=25*np.finfo(dtype).eps)


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
                    
class TestDlasd4(TestCase):
    def test_sing_val_update(self):

        sigmas = np.array([4., 3., 2., 0])
        m_vec = np.array([3.12, 5.7, -4.8, -2.2])

        M = np.hstack((np.vstack((np.diag(sigmas[0:-1]),
                        np.zeros((1,len(m_vec) - 1)))), m_vec[:, np.newaxis]))
        SM = svd(M, full_matrices=False, compute_uv=False, overwrite_a=False,
                 check_finite=False)

        it_len = len(sigmas)
        sgm = np.concatenate((sigmas[::-1], (sigmas[0] +
                              it_len*np.sqrt(np.sum(np.power(m_vec,2))),)))
        mvc = np.concatenate((m_vec[::-1], (0,)))

        lasd4 = get_lapack_funcs('lasd4',(sigmas,))

        roots = []
        for i in range(0, it_len):
            res = lasd4(i, sgm, mvc)
            roots.append(res[1])

            assert_((res[3] <= 0),"LAPACK root finding dlasd4 failed to find \
                                    the singular value %i" % i)
        roots = np.array(roots)[::-1]

        assert_((not np.any(np.isnan(roots)),"There are NaN roots"))
        assert_allclose(SM, roots, atol=100*np.finfo(np.float64).eps,
                        rtol=100*np.finfo(np.float64).eps)


if __name__ == "__main__":
    run_module_suite()
