from scipy.linalg cimport blas_pointers
import numpy as np
from numpy.testing import TestCase, assert_allclose

class test_dgemm_pointer(TestCase):
    
    def test_basic(self):
        
        cdef:
            double[:,:] a, b, c
            double alpha, beta
            int lda, ldb, ldc, m, n, k
        
        alpha, beta = 1., 0.
        lda = ldb = ldc = m = n = k = 2
        a = np.array([[1, 2], [3, 4]], float, order="F")
        b = np.array([[5, 6], [7, 8]], float, order="F")
        c = np.empty((2, 2), float, order="F")
        
        blas_pointers.dgemm("N", "N", &m, &n, &k, &alpha, &a[0,0],
                            &lda, &b[0,0], &ldb, &beta, &c[0,0], &ldc)
        
        assert_allclose(c, np.array([[19., 22.],
                                     [43., 50.]], float))

class test_wfunc_pointers(TestCase):
    """ Test the function pointers that are expected to fail on
    Mac OS X without the additional entry statement in their definitions
    in fblas_l1.pyf.src. """
    
    def test_complex_args(self):
        
        cdef:
            float complex[:] cx, cy
            float complex cout
            float sout
            int n, incx, incy
        
        cx = np.array([.5 + 1.j, .25 - .375j, 12.5 - 4.j], np.complex64)
        cy = np.array([.8 + 2.j, .875 - .625j, -1. + 2.j], np.complex64)
        
        n, incx, incy = 3, 1, 1
        cout = blas_pointers.cdotc(&n, &cx[0], &incx, &cy[0], &incy)
        assert_allclose(cout, -17.6468753815+21.3718757629j, 5)
        cout = blas_pointers.cdotu(&n, &cx[0], &incx, &cy[0], &incy)
        assert_allclose(cout, -6.11562538147+30.3156242371j, 5)
        sout = blas_pointers.scasum(&n, &cx[0], &incx)
        assert_allclose(sout, 18.625, 5)
        sout = blas_pointers.scnrm2(&n, &cx[0], &incx)
        assert_allclose(sout, 13.1796483994, 5)
        
        n = incx = incy = 2
        cout = blas_pointers.cdotc(&n, &cx[0], &incx, &cy[0], &incy)
        assert_allclose(cout, -18.1000003815+21.2000007629j, 5)
        cout = blas_pointers.cdotu(&n, &cx[0], &incx, &cy[0], &incy)
        assert_allclose(cout, -6.10000038147+30.7999992371j, 5)
        sout = blas_pointers.scasum(&n, &cx[0], &incx)
        assert_allclose(sout, 18., 5)
        sout = blas_pointers.scnrm2(&n, &cx[0], &incx)
        assert_allclose(sout, 13.1719398499, 5)
    
    def test_float_args(self):
        
        cdef:
            float[:] x, y
            float out
            int n, incx, incy
        
        x = np.array([5., -3, -.5], np.float32)
        y = np.array([2, 1, .5], np.float32)
        
        n, incx, incy = 3, 1, 1
        out = blas_pointers.sasum(&n, &x[0], &incx)
        assert_allclose(out, 8.5, 5)
        out = blas_pointers.sdot(&n, &x[0], &incx, &y[0], &incy)
        assert_allclose(out, 6.75, 5)
        out = blas_pointers.snrm2(&n, &x[0], &incx)
        assert_allclose(out, 5.85234975815, 5)
        
        n = incx = incy = 2
        out = blas_pointers.sasum(&n, &x[0], &incx)
        assert_allclose(out, 5.5, 5)
        out = blas_pointers.sdot(&n, &x[0], &incx, &y[0], &incy)
        assert_allclose(out, 9.75, 5)
        out = blas_pointers.snrm2(&n, &x[0], &incx)
        assert_allclose(out, 5.0249376297, 5)
    
    def test_double_complex_args(self):
        
        cdef:
            double complex[:] cx, cy
            double complex out
            int n, incx, incy
        
        cx = np.array([.5 + 1.j, .25 - .375j, 13. - 4.j], np.complex128)
        cy = np.array([.875 + 2.j, .875 - .625j, -1. + 2.j], np.complex128)
        
        n, incx, incy = 3, 1, 1
        out = blas_pointers.zdotc(&n, &cx[0], &incx, &cy[0], &incy)
        assert_allclose(out, -18.109375+22.296875j, 5)
        out = blas_pointers.zdotu(&n, &cx[0], &incx, &cy[0], &incy)
        assert_allclose(out, -6.578125+31.390625j, 5)
        
        n = incx = incy = 2
        out = blas_pointers.zdotc(&n, &cx[0], &incx, &cy[0], &incy)
        assert_allclose(out, -18.5625+22.125j, 5)
        out = blas_pointers.zdotu(&n, &cx[0], &incx, &cy[0], &incy)
        assert_allclose(out, -6.5625+31.875j, 5)
