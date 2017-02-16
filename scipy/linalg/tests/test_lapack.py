#
# Created by: Pearu Peterson, September 2002
#

from __future__ import division, print_function, absolute_import

import sys
import subprocess
import time

from numpy.testing import TestCase, run_module_suite, assert_equal, \
    assert_array_almost_equal, assert_, assert_raises, assert_allclose, \
    assert_almost_equal

import numpy as np

from scipy.linalg import _flapack as flapack
from scipy.linalg import inv
from scipy.linalg import svd
from scipy._lib._testutils import xslow

try:
    from scipy.linalg import _clapack as clapack
except ImportError:
    clapack = None
from scipy.linalg.lapack import get_lapack_funcs
from scipy.linalg.blas import get_blas_funcs

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
            lwork = int(np.real(work))
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
            lwork = int(np.real(work))
            rwork_size = int(rwork)
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
            lwork = int(np.real(work))

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
            lwork = int(np.real(work))

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
            lwork = int(np.real(work))

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
            lwork = int(np.real(work))

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


def test_lartg():
    for dtype in 'fdFD':
        lartg = get_lapack_funcs('lartg', dtype=dtype)

        f = np.array(3, dtype)
        g = np.array(4, dtype)

        if np.iscomplexobj(g):
            g *= 1j

        cs, sn, r = lartg(f, g)

        assert_allclose(cs, 3.0/5.0)
        assert_allclose(r, 5.0)

        if np.iscomplexobj(g):
            assert_allclose(sn, -4.0j/5.0)
            assert_(type(r) == complex)
            assert_(type(cs) == float)
        else:
            assert_allclose(sn, 4.0/5.0)

def test_rot():
    # srot, drot from blas and crot and zrot from lapack.

    for dtype in 'fdFD':
        c = 0.6
        s = 0.8

        u = np.ones(4, dtype) * 3
        v = np.ones(4, dtype) * 4
        atol = 10**-(np.finfo(dtype).precision-1)

        if dtype in 'fd':
            rot = get_blas_funcs('rot', dtype=dtype)
            f = 4
        else:
            rot = get_lapack_funcs('rot', dtype=dtype)
            s *= -1j
            v *= 1j
            f = 4j

        assert_allclose(rot(u, v, c, s), [[5,5,5,5],[0,0,0,0]], atol=atol)
        assert_allclose(rot(u, v, c, s, n=2), [[5,5,3,3],[0,0,f,f]], atol=atol)
        assert_allclose(rot(u, v, c, s, offx=2,offy=2), [[3,3,5,5],[f,f,0,0]], atol=atol)
        assert_allclose(rot(u, v, c, s, incx=2, offy=2, n=2), [[5,3,5,3],[f,f,0,0]], atol=atol)
        assert_allclose(rot(u, v, c, s, offx=2, incy=2, n=2), [[3,3,5,5],[0,f,0,f]], atol=atol)
        assert_allclose(rot(u, v, c, s, offx=2, incx=2, offy=2, incy=2, n=1), [[3,3,5,3],[f,f,0,f]], atol=atol)
        assert_allclose(rot(u, v, c, s, incx=-2, incy=-2, n=2), [[5,3,5,3],[0,f,0,f]], atol=atol)
    
        a, b = rot(u, v, c, s, overwrite_x=1, overwrite_y=1)
        assert_(a is u)
        assert_(b is v)
        assert_allclose(a, [5,5,5,5], atol=atol)
        assert_allclose(b, [0,0,0,0], atol=atol)

def test_larfg_larf():
    np.random.seed(1234)
    a0 = np.random.random((4,4))
    a0 = a0.T.dot(a0)

    a0j = np.random.random((4,4)) + 1j*np.random.random((4,4))
    a0j = a0j.T.conj().dot(a0j)

    # our test here will be to do one step of reducing a hermetian matrix to
    # tridiagonal form using householder transforms.

    for dtype in 'fdFD':
        larfg, larf = get_lapack_funcs(['larfg', 'larf'], dtype=dtype)

        if dtype in 'FD':
            a = a0j.copy()
        else:
            a = a0.copy()

        # generate a householder transform to clear a[2:,0]
        alpha, x, tau = larfg(a.shape[0]-1, a[1,0], a[2:,0])

        # create expected output
        expected = np.zeros_like(a[:,0])
        expected[0] = a[0,0]
        expected[1] = alpha
        
        # assemble householder vector
        v = np.zeros_like(a[1:,0])
        v[0] = 1.0
        v[1:] = x

        # apply transform from the left
        a[1:,:] = larf(v, tau.conjugate(), a[1:,:], np.zeros(a.shape[1]))

        # apply transform from the right
        a[:,1:] = larf(v, tau, a[:,1:], np.zeros(a.shape[0]), side='R')
        
        assert_allclose(a[:,0], expected, atol=1e-5)
        assert_allclose(a[0,:], expected, atol=1e-5)


@xslow
def test_sgesdd_lwork_bug_workaround():
    # Test that SGESDD lwork is sufficiently large for LAPACK.
    #
    # This checks that workaround around an apparent LAPACK bug
    # actually works. cf. gh-5401
    #
    # xslow: requires 1GB+ of memory

    p = subprocess.Popen([sys.executable, '-c',
                          'import numpy as np; '
                          'from scipy.linalg import svd; '
                          'a = np.zeros([9537, 9537], dtype=np.float32); '
                          'svd(a)'],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)

    # Check if it an error occurred within 5 sec; the computation can
    # take substantially longer, and we will not wait for it to finish
    for j in range(50):
        time.sleep(0.1)
        if p.poll() is not None:
            returncode = p.returncode
            break
    else:
        # Didn't exit in time -- probably entered computation.  The
        # error is raised before entering computation, so things are
        # probably OK.
        returncode = 0
        p.terminate()

    assert_equal(returncode, 0,
                 "Code apparently failed: " + p.stdout.read())


if __name__ == "__main__":
    run_module_suite()
