#!/usr/bin/env python
#
# Created by: Pearu Peterson, April 2002
#
from __future__ import division, print_function, absolute_import


__usage__ = """
Build linalg:
  python setup.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.linalg.test()'
"""

import math

import numpy as np
from numpy.testing import (TestCase, run_module_suite, assert_equal,
    assert_almost_equal, assert_array_almost_equal, assert_raises, assert_,
    assert_allclose)

from scipy.linalg import _fblas as fblas, get_blas_funcs

try:
    from scipy.linalg import _cblas as cblas
except ImportError:
    cblas = None


def test_get_blas_funcs():
    # check that it returns Fortran code for arrays that are
    # fortran-ordered
    f1, f2, f3 = get_blas_funcs(
        ('axpy', 'axpy', 'axpy'),
        (np.empty((2,2), dtype=np.complex64, order='F'),
         np.empty((2,2), dtype=np.complex128, order='C'))
        )

    # get_blas_funcs will choose libraries depending on most generic
    # array
    assert_equal(f1.typecode, 'z')
    assert_equal(f2.typecode, 'z')
    if cblas is not None:
        assert_equal(f1.module_name, 'cblas')
        assert_equal(f2.module_name, 'cblas')

    # check defaults.
    f1 = get_blas_funcs('rotg')
    assert_equal(f1.typecode, 'd')

    # check also dtype interface
    f1 = get_blas_funcs('gemm', dtype=np.complex64)
    assert_equal(f1.typecode, 'c')
    f1 = get_blas_funcs('gemm', dtype='F')
    assert_equal(f1.typecode, 'c')

    # extended precision complex
    f1 = get_blas_funcs('gemm', dtype=np.longcomplex)
    assert_equal(f1.typecode, 'z')

    # check safe complex upcasting
    f1 = get_blas_funcs('axpy',
                        (np.empty((2,2), dtype=np.float64),
                         np.empty((2,2), dtype=np.complex64))
                        )
    assert_equal(f1.typecode, 'z')


def test_get_blas_funcs_alias():
    # check alias for get_blas_funcs
    f, g = get_blas_funcs(('nrm2', 'dot'), dtype=np.complex64)
    assert f.typecode == 'c'
    assert g.typecode == 'c'

    f, g, h = get_blas_funcs(('dot', 'dotc', 'dotu'), dtype=np.float64)
    assert f is g
    assert f is h


class TestCBLAS1Simple(TestCase):

    def test_axpy(self):
        for p in 'sd':
            f = getattr(cblas,p+'axpy',None)
            if f is None:
                continue
            assert_array_almost_equal(f([1,2,3],[2,-1,3],a=5),[7,9,18])
        for p in 'cz':
            f = getattr(cblas,p+'axpy',None)
            if f is None:
                continue
            assert_array_almost_equal(f([1,2j,3],[2,-1,3],a=5),[7,10j-1,18])


class TestFBLAS1Simple(TestCase):

    def test_axpy(self):
        for p in 'sd':
            f = getattr(fblas,p+'axpy',None)
            if f is None:
                continue
            assert_array_almost_equal(f([1,2,3],[2,-1,3],a=5),[7,9,18])
        for p in 'cz':
            f = getattr(fblas,p+'axpy',None)
            if f is None:
                continue
            assert_array_almost_equal(f([1,2j,3],[2,-1,3],a=5),[7,10j-1,18])

    def test_copy(self):
        for p in 'sd':
            f = getattr(fblas,p+'copy',None)
            if f is None:
                continue
            assert_array_almost_equal(f([3,4,5],[8]*3),[3,4,5])
        for p in 'cz':
            f = getattr(fblas,p+'copy',None)
            if f is None:
                continue
            assert_array_almost_equal(f([3,4j,5+3j],[8]*3),[3,4j,5+3j])

    def test_asum(self):
        for p in 'sd':
            f = getattr(fblas,p+'asum',None)
            if f is None:
                continue
            assert_almost_equal(f([3,-4,5]),12)
        for p in ['sc','dz']:
            f = getattr(fblas,p+'asum',None)
            if f is None:
                continue
            assert_almost_equal(f([3j,-4,3-4j]),14)

    def test_dot(self):
        for p in 'sd':
            f = getattr(fblas,p+'dot',None)
            if f is None:
                continue
            assert_almost_equal(f([3,-4,5],[2,5,1]),-9)

    def test_complex_dotu(self):
        for p in 'cz':
            f = getattr(fblas,p+'dotu',None)
            if f is None:
                continue
            assert_almost_equal(f([3j,-4,3-4j],[2,3,1]),-9+2j)

    def test_complex_dotc(self):
        for p in 'cz':
            f = getattr(fblas,p+'dotc',None)
            if f is None:
                continue
            assert_almost_equal(f([3j,-4,3-4j],[2,3j,1]),3-14j)

    def test_nrm2(self):
        for p in 'sd':
            f = getattr(fblas,p+'nrm2',None)
            if f is None:
                continue
            assert_almost_equal(f([3,-4,5]),math.sqrt(50))
        for p in ['c', 'z', 'sc','dz']:
            f = getattr(fblas,p+'nrm2',None)
            if f is None:
                continue
            assert_almost_equal(f([3j,-4,3-4j]),math.sqrt(50))

    def test_scal(self):
        for p in 'sd':
            f = getattr(fblas,p+'scal',None)
            if f is None:
                continue
            assert_array_almost_equal(f(2,[3,-4,5]),[6,-8,10])
        for p in 'cz':
            f = getattr(fblas,p+'scal',None)
            if f is None:
                continue
            assert_array_almost_equal(f(3j,[3j,-4,3-4j]),[-9,-12j,12+9j])
        for p in ['cs','zd']:
            f = getattr(fblas,p+'scal',None)
            if f is None:
                continue
            assert_array_almost_equal(f(3,[3j,-4,3-4j]),[9j,-12,9-12j])

    def test_swap(self):
        for p in 'sd':
            f = getattr(fblas,p+'swap',None)
            if f is None:
                continue
            x,y = [2,3,1],[-2,3,7]
            x1,y1 = f(x,y)
            assert_array_almost_equal(x1,y)
            assert_array_almost_equal(y1,x)
        for p in 'cz':
            f = getattr(fblas,p+'swap',None)
            if f is None:
                continue
            x,y = [2,3j,1],[-2,3,7-3j]
            x1,y1 = f(x,y)
            assert_array_almost_equal(x1,y)
            assert_array_almost_equal(y1,x)

    def test_amax(self):
        for p in 'sd':
            f = getattr(fblas,'i'+p+'amax')
            assert_equal(f([-2,4,3]),1)
        for p in 'cz':
            f = getattr(fblas,'i'+p+'amax')
            assert_equal(f([-5,4+3j,6]),1)
    #XXX: need tests for rot,rotm,rotg,rotmg


class TestFBLAS2Simple(TestCase):

    def test_gemv(self):
        for p in 'sd':
            f = getattr(fblas,p+'gemv',None)
            if f is None:
                continue
            assert_array_almost_equal(f(3,[[3]],[-4]),[-36])
            assert_array_almost_equal(f(3,[[3]],[-4],3,[5]),[-21])
        for p in 'cz':
            f = getattr(fblas,p+'gemv',None)
            if f is None:
                continue
            assert_array_almost_equal(f(3j,[[3-4j]],[-4]),[-48-36j])
            assert_array_almost_equal(f(3j,[[3-4j]],[-4],3,[5j]),[-48-21j])

    def test_ger(self):

        for p in 'sd':
            f = getattr(fblas,p+'ger',None)
            if f is None:
                continue
            assert_array_almost_equal(f(1,[1,
                                           2],[3,4]),[[3,4],[6,8]])
            assert_array_almost_equal(f(2,[1,
                                           2,
                                           3],[3,4]),[[6,8],[12,16],[18,24]])

            assert_array_almost_equal(f(1,[1,
                                           2],[3,4],
                                        a=[[1,2],[3,4]]
                                        ),[[4,6],[9,12]])

        for p in 'cz':
            f = getattr(fblas,p+'geru',None)
            if f is None:
                continue
            assert_array_almost_equal(f(1,[1j,
                                           2],[3,4]),[[3j,4j],[6,8]])
            assert_array_almost_equal(f(-2,[1j,
                                           2j,
                                           3j],[3j,4j]),[[6,8],[12,16],[18,24]])

        for p in 'cz':
            for name in ('ger', 'gerc'):
                f = getattr(fblas,p+name,None)
                if f is None:
                    continue
                assert_array_almost_equal(f(1,[1j,
                                               2],[3,4]),[[3j,4j],[6,8]])
                assert_array_almost_equal(f(2,[1j,
                                               2j,
                                               3j],[3j,4j]),[[6,8],[12,16],[18,24]])

    def test_syr_her(self):
        x = np.arange(1, 5, dtype='d')
        resx = np.triu(x[:, np.newaxis] * x)
        resx_reverse = np.triu(x[::-1, np.newaxis] * x[::-1])

        y = np.linspace(0,8.5,17,endpoint=False)

        z = np.arange(1, 9, dtype='d').view('D')
        resz = np.triu(z[:, np.newaxis] * z)
        resz_reverse = np.triu(z[::-1, np.newaxis] * z[::-1])
        rehz = np.triu(z[:, np.newaxis] * z.conj())
        rehz_reverse = np.triu(z[::-1, np.newaxis] * z[::-1].conj())

        w = np.c_[np.zeros(4), z, np.zeros(4)].ravel()
     
        for p, rtol in zip('sd',[1e-7,1e-14]):
            f = getattr(fblas, p+'syr', None)
            if f is None:
                continue
            assert_allclose(f(1.0, x), resx, rtol=rtol)
            assert_allclose(f(1.0, x, lower=True), resx.T, rtol=rtol)
            assert_allclose(f(1.0, y, incx=2, offx=2, n=4), resx, rtol=rtol)
            # negative increments imply reversed vectors in blas
            assert_allclose(f(1.0, y, incx=-2, offx=2, n=4),
                resx_reverse, rtol=rtol)

            a = np.zeros((4,4), 'f' if p == 's' else 'd', 'F')
            b = f(1.0, x, a=a, overwrite_a=True)
            assert_(a is b)
            assert_allclose(a, resx, rtol=rtol)

            b = f(2.0, x, a=a)
            assert_(a is not b)
            assert_allclose(b, 3*resx, rtol=rtol)

            assert_raises(Exception, f, 1.0, x, incx=0)
            assert_raises(Exception, f, 1.0, x, offx=5)
            assert_raises(Exception, f, 1.0, x, offx=-2)
            assert_raises(Exception, f, 1.0, x, n=-2)
            assert_raises(Exception, f, 1.0, x, n=5)
            assert_raises(Exception, f, 1.0, x, lower=2)
            assert_raises(Exception, f, 1.0, x, a=np.zeros((2,2), 'd', 'F'))
        
        for p, rtol in zip('cz',[1e-7,1e-14]):
            f = getattr(fblas, p+'syr', None)
            if f is None:
                continue
            assert_allclose(f(1.0, z), resz, rtol=rtol)
            assert_allclose(f(1.0, z, lower=True), resz.T, rtol=rtol)
            assert_allclose(f(1.0, w, incx=3, offx=1, n=4), resz, rtol=rtol)
            # negative increments imply reversed vectors in blas
            assert_allclose(f(1.0, w, incx=-3, offx=1, n=4),
                resz_reverse, rtol=rtol)

            a = np.zeros((4,4), 'F' if p == 'c' else 'D', 'F')
            b = f(1.0, z, a=a, overwrite_a=True)
            assert_(a is b)
            assert_allclose(a, resz, rtol=rtol)

            b = f(2.0, z, a=a)
            assert_(a is not b)
            assert_allclose(b, 3*resz, rtol=rtol)

            assert_raises(Exception, f, 1.0, x, incx=0)
            assert_raises(Exception, f, 1.0, x, offx=5)
            assert_raises(Exception, f, 1.0, x, offx=-2)
            assert_raises(Exception, f, 1.0, x, n=-2)
            assert_raises(Exception, f, 1.0, x, n=5)
            assert_raises(Exception, f, 1.0, x, lower=2)
            assert_raises(Exception, f, 1.0, x, a=np.zeros((2,2), 'd', 'F'))
 
        for p, rtol in zip('cz',[1e-7,1e-14]):
            f = getattr(fblas, p+'her', None)
            if f is None:
                continue
            assert_allclose(f(1.0, z), rehz, rtol=rtol)
            assert_allclose(f(1.0, z, lower=True), rehz.T.conj(), rtol=rtol)
            assert_allclose(f(1.0, w, incx=3, offx=1, n=4), rehz, rtol=rtol)
            # negative increments imply reversed vectors in blas
            assert_allclose(f(1.0, w, incx=-3, offx=1, n=4),
                rehz_reverse, rtol=rtol)

            a = np.zeros((4,4), 'F' if p == 'c' else 'D', 'F')
            b = f(1.0, z, a=a, overwrite_a=True)
            assert_(a is b)
            assert_allclose(a, rehz, rtol=rtol)

            b = f(2.0, z, a=a)
            assert_(a is not b)
            assert_allclose(b, 3*rehz, rtol=rtol)

            assert_raises(Exception, f, 1.0, x, incx=0)
            assert_raises(Exception, f, 1.0, x, offx=5)
            assert_raises(Exception, f, 1.0, x, offx=-2)
            assert_raises(Exception, f, 1.0, x, n=-2)
            assert_raises(Exception, f, 1.0, x, n=5)
            assert_raises(Exception, f, 1.0, x, lower=2)
            assert_raises(Exception, f, 1.0, x, a=np.zeros((2,2), 'd', 'F'))

    def test_syr2(self):
        x = np.arange(1, 5, dtype='d')
        y = np.arange(5, 9, dtype='d')
        resxy = np.triu(x[:, np.newaxis] * y + y[:, np.newaxis] * x)
        resxy_reverse = np.triu(x[::-1, np.newaxis] * y[::-1]
                                + y[::-1, np.newaxis] * x[::-1])

        q = np.linspace(0,8.5,17,endpoint=False)

        for p, rtol in zip('sd',[1e-7,1e-14]):
            f = getattr(fblas, p+'syr2', None)
            if f is None:
                continue
            assert_allclose(f(1.0, x, y), resxy, rtol=rtol)
            assert_allclose(f(1.0, x, y, n=3), resxy[:3,:3], rtol=rtol)
            assert_allclose(f(1.0, x, y, lower=True), resxy.T, rtol=rtol)

            assert_allclose(f(1.0, q, q, incx=2, offx=2, incy=2, offy=10),
                            resxy, rtol=rtol)
            assert_allclose(f(1.0, q, q, incx=2, offx=2, incy=2, offy=10, n=3),
                            resxy[:3,:3], rtol=rtol)
            # negative increments imply reversed vectors in blas
            assert_allclose(f(1.0, q, q, incx=-2, offx=2, incy=-2, offy=10),
                            resxy_reverse, rtol=rtol)

            a = np.zeros((4,4), 'f' if p == 's' else 'd', 'F')
            b = f(1.0, x, y, a=a, overwrite_a=True)
            assert_(a is b)
            assert_allclose(a, resxy, rtol=rtol)

            b = f(2.0, x, y, a=a)
            assert_(a is not b)
            assert_allclose(b, 3*resxy, rtol=rtol)

            assert_raises(Exception, f, 1.0, x, y, incx=0)
            assert_raises(Exception, f, 1.0, x, y, offx=5)
            assert_raises(Exception, f, 1.0, x, y, offx=-2)
            assert_raises(Exception, f, 1.0, x, y, incy=0)
            assert_raises(Exception, f, 1.0, x, y, offy=5)
            assert_raises(Exception, f, 1.0, x, y, offy=-2)
            assert_raises(Exception, f, 1.0, x, y, n=-2)
            assert_raises(Exception, f, 1.0, x, y, n=5)
            assert_raises(Exception, f, 1.0, x, y, lower=2)
            assert_raises(Exception, f, 1.0, x, y, a=np.zeros((2,2), 'd', 'F'))
  
    def test_her2(self):
        x = np.arange(1, 9, dtype='d').view('D')
        y = np.arange(9, 17, dtype='d').view('D')
        resxy = x[:, np.newaxis] * y.conj() + y[:, np.newaxis] * x.conj()
        resxy = np.triu(resxy)

        resxy_reverse = x[::-1, np.newaxis] * y[::-1].conj()
        resxy_reverse += y[::-1, np.newaxis] * x[::-1].conj()
        resxy_reverse = np.triu(resxy_reverse)

        u = np.c_[np.zeros(4), x, np.zeros(4)].ravel()
        v = np.c_[np.zeros(4), y, np.zeros(4)].ravel()

        for p, rtol in zip('cz',[1e-7,1e-14]):
            f = getattr(fblas, p+'her2', None)
            if f is None:
                continue
            assert_allclose(f(1.0, x, y), resxy, rtol=rtol)
            assert_allclose(f(1.0, x, y, n=3), resxy[:3,:3], rtol=rtol)
            assert_allclose(f(1.0, x, y, lower=True), resxy.T.conj(), rtol=rtol)

            assert_allclose(f(1.0, u, v, incx=3, offx=1, incy=3, offy=1),
                            resxy, rtol=rtol)
            assert_allclose(f(1.0, u, v, incx=3, offx=1, incy=3, offy=1, n=3),
                            resxy[:3,:3], rtol=rtol)
            # negative increments imply reversed vectors in blas
            assert_allclose(f(1.0, u, v, incx=-3, offx=1, incy=-3, offy=1),
                            resxy_reverse, rtol=rtol)

            a = np.zeros((4,4), 'F' if p == 'c' else 'D', 'F')
            b = f(1.0, x, y, a=a, overwrite_a=True)
            assert_(a is b)
            assert_allclose(a, resxy, rtol=rtol)

            b = f(2.0, x, y, a=a)
            assert_(a is not b)
            assert_allclose(b, 3*resxy, rtol=rtol)

            assert_raises(Exception, f, 1.0, x, y, incx=0)
            assert_raises(Exception, f, 1.0, x, y, offx=5)
            assert_raises(Exception, f, 1.0, x, y, offx=-2)
            assert_raises(Exception, f, 1.0, x, y, incy=0)
            assert_raises(Exception, f, 1.0, x, y, offy=5)
            assert_raises(Exception, f, 1.0, x, y, offy=-2)
            assert_raises(Exception, f, 1.0, x, y, n=-2)
            assert_raises(Exception, f, 1.0, x, y, n=5)
            assert_raises(Exception, f, 1.0, x, y, lower=2)
            assert_raises(Exception, f, 1.0, x, y, a=np.zeros((2,2), 'd', 'F'))
 

class TestFBLAS3Simple(TestCase):

    def test_gemm(self):
        for p in 'sd':
            f = getattr(fblas,p+'gemm',None)
            if f is None:
                continue
            assert_array_almost_equal(f(3,[3],[-4]),[[-36]])
            assert_array_almost_equal(f(3,[3],[-4],3,[5]),[-21])
        for p in 'cz':
            f = getattr(fblas,p+'gemm',None)
            if f is None:
                continue
            assert_array_almost_equal(f(3j,[3-4j],[-4]),[[-48-36j]])
            assert_array_almost_equal(f(3j,[3-4j],[-4],3,[5j]),[-48-21j])


def _get_func(func, ps='sdzc'):
    """Just a helper: return a specified BLAS function w/typecode."""
    for p in ps:
        f = getattr(fblas, p+func, None)
        if f is None:
            continue
        yield f


class TestBLAS3Symm(TestCase):

    def setUp(self):
        self.a = np.array([[1., 2.],
                           [0., 1.]])
        self.b = np.array([[1., 0., 3.],
                           [0., -1., 2.]])
        self.c = np.ones((2,3))
        self.t = np.array([[2., -1., 8.],
                           [3., 0., 9.]])

    def test_symm(self):
        for f in _get_func('symm'):
            res = f(a=self.a, b=self.b, c=self.c, alpha=1., beta=1.)
            assert_array_almost_equal(res, self.t)

            res = f(a=self.a.T, b=self.b, lower=1, c=self.c, alpha=1., beta=1.)
            assert_array_almost_equal(res, self.t)

            res = f(a=self.a, b=self.b.T, side=1, c=self.c.T, alpha=1., beta=1.)
            assert_array_almost_equal(res, self.t.T)

    def test_summ_wrong_side(self):
        f = getattr(fblas, 'dsymm', None)
        if f is not None:
            assert_raises(Exception, f, **{'a': self.a, 'b': self.b, 'alpha': 1,
                    'side': 1})
            # `side=1` means C <- B*A, hence shapes of A and B are to be
            #  compatible. Otherwise, f2py exception is raised

    def test_symm_wrong_uplo(self):
        """SYMM only considers the upper/lower part of A. Hence setting
        wrong value for `lower` (default is lower=0, meaning upper triangle)
        gives a wrong result.
        """
        f = getattr(fblas,'dsymm',None)
        if f is not None:
            res = f(a=self.a, b=self.b, c=self.c, alpha=1., beta=1.)
            assert np.allclose(res, self.t)

            res = f(a=self.a, b=self.b, lower=1, c=self.c, alpha=1., beta=1.)
            assert not np.allclose(res, self.t)


class TestBLAS3Syrk(TestCase):
    def setUp(self):
        self.a = np.array([[1., 0.],
                           [0., -2.],
                           [2., 3.]])
        self.t = np.array([[1., 0., 2.],
                           [0., 4., -6.],
                           [2., -6., 13.]])
        self.tt = np.array([[5., 6.],
                            [6., 13.]])

    def test_syrk(self):
        for f in _get_func('syrk'):
            c = f(a=self.a, alpha=1.)
            assert_array_almost_equal(np.triu(c), np.triu(self.t))

            c = f(a=self.a, alpha=1., lower=1)
            assert_array_almost_equal(np.tril(c), np.tril(self.t))

            c0 = np.ones(self.t.shape)
            c = f(a=self.a, alpha=1., beta=1., c=c0)
            assert_array_almost_equal(np.triu(c), np.triu(self.t+c0))

            c = f(a=self.a, alpha=1., trans=1)
            assert_array_almost_equal(np.triu(c), np.triu(self.tt))

    #prints '0-th dimension must be fixed to 3 but got 5', FIXME: suppress?
    # FIXME: how to catch the _fblas.error?
    def test_syrk_wrong_c(self):
        f = getattr(fblas, 'dsyrk', None)
        if f is not None:
            assert_raises(Exception, f, **{'a': self.a, 'alpha': 1.,
                    'c': np.ones((5, 8))})
        # if C is supplied, it must have compatible dimensions


class TestBLAS3Syr2k(TestCase):
    def setUp(self):
        self.a = np.array([[1., 0.],
                           [0., -2.],
                           [2., 3.]])
        self.b = np.array([[0., 1.],
                           [1., 0.],
                           [0, 1.]])
        self.t = np.array([[0., -1., 3.],
                           [-1., 0., 0.],
                           [3., 0., 6.]])
        self.tt = np.array([[0., 1.],
                            [1., 6]])

    def test_syr2k(self):
        for f in _get_func('syr2k'):
            c = f(a=self.a, b=self.b, alpha=1.)
            assert_array_almost_equal(np.triu(c), np.triu(self.t))

            c = f(a=self.a, b=self.b, alpha=1., lower=1)
            assert_array_almost_equal(np.tril(c), np.tril(self.t))

            c0 = np.ones(self.t.shape)
            c = f(a=self.a, b=self.b, alpha=1., beta=1., c=c0)
            assert_array_almost_equal(np.triu(c), np.triu(self.t+c0))

            c = f(a=self.a, b=self.b, alpha=1., trans=1)
            assert_array_almost_equal(np.triu(c), np.triu(self.tt))

    #prints '0-th dimension must be fixed to 3 but got 5', FIXME: suppress?
    def test_syr2k_wrong_c(self):
        f = getattr(fblas, 'dsyr2k', None)
        if f is not None:
            assert_raises(Exception, f, **{'a': self.a, 'b': self.b, 'alpha': 1.,
                    'c': np.zeros((15, 8))})
        # if C is supplied, it must have compatible dimensions


class TestSyHe(TestCase):
    """Quick and simple tests for (zc)-symm, syrk, syr2k."""
    def setUp(self):
        self.sigma_y = np.array([[0., -1.j],
                                 [1.j, 0.]])

    def test_symm_zc(self):
        for f in _get_func('symm', 'zc'):
            # NB: a is symmetric w/upper diag of ONLY
            res = f(a=self.sigma_y, b=self.sigma_y, alpha=1.)
            assert_array_almost_equal(np.triu(res), np.diag([1, -1]))

    def test_hemm_zc(self):
        for f in _get_func('hemm', 'zc'):
            # NB: a is hermitian w/upper diag of ONLY
            res = f(a=self.sigma_y, b=self.sigma_y, alpha=1.)
            assert_array_almost_equal(np.triu(res), np.diag([1, 1]))

    def test_syrk_zr(self):
        for f in _get_func('syrk', 'zc'):
            res = f(a=self.sigma_y, alpha=1.)
            assert_array_almost_equal(np.triu(res), np.diag([-1, -1]))

    def test_herk_zr(self):
        for f in _get_func('herk', 'zc'):
            res = f(a=self.sigma_y, alpha=1.)
            assert_array_almost_equal(np.triu(res), np.diag([1, 1]))

    def test_syr2k_zr(self):
        for f in _get_func('syr2k', 'zc'):
            res = f(a=self.sigma_y, b=self.sigma_y, alpha=1.)
            assert_array_almost_equal(np.triu(res), 2.*np.diag([-1, -1]))

    def test_her2k_zr(self):
        for f in _get_func('her2k', 'zc'):
            res = f(a=self.sigma_y, b=self.sigma_y, alpha=1.)
            assert_array_almost_equal(np.triu(res), 2.*np.diag([1, 1]))


class TestTRMM(TestCase):
    """Quick and simple tests for dtrmm."""
    def setUp(self):
        self.a = np.array([[1., 2., ],
                           [-2., 1.]])
        self.b = np.array([[3., 4., -1.],
                           [5., 6., -2.]])

    def test_ab(self):
        f = getattr(fblas, 'dtrmm', None)
        if f is not None:
            result = f(1., self.a, self.b)
            expected = np.array([[13., 16., -5.],
                                 [5., 6., -2.]])  # default a is upper triangular
            assert_array_almost_equal(result, expected)

    def test_ab_lower(self):
        f = getattr(fblas, 'dtrmm', None)
        if f is not None:
            result = f(1., self.a, self.b, lower=True)
            expected = np.array([[3., 4., -1.],
                                 [-1., -2., 0.]])  # now a is lower triangular
            assert_array_almost_equal(result, expected)

    def test_b_overwrites(self):
        # BLAS dtrmm modifies B argument in-place.
        # Here the default is to copy, but this can be overridden
        f = getattr(fblas, 'dtrmm', None)
        if f is not None:
            for overwr in [True, False]:
                bcopy = self.b.copy()
                result = f(1., self.a, bcopy, overwrite_b=overwr)
                # C-contiguous arrays are copied
                assert_(bcopy.flags.f_contiguous is False and
                        np.may_share_memory(bcopy, result) is False)
                assert_equal(bcopy, self.b)

            bcopy = np.asfortranarray(self.b.copy())  # or just transpose it
            result = f(1., self.a, bcopy, overwrite_b=True)
            assert_(bcopy.flags.f_contiguous is True and
                    np.may_share_memory(bcopy, result) is True)
            assert_array_almost_equal(bcopy, result)


if __name__ == "__main__":
    run_module_suite()
