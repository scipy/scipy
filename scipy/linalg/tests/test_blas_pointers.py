import numpy as np
from numpy.testing import TestCase, run_module_suite, assert_allclose
import scipy.linalg._cython_blas_wrappers as blas

class test_dgemm(TestCase):
    
    def test_transposes(self):

        a = np.arange(12, dtype='d').reshape((3, 4))[:2,:2]
        b = np.arange(1, 13, dtype='d').reshape((4, 3))[:2,:2]
        c = np.empty((2, 4))[:2,:2]

        blas.dgemm(1., a, b, 0., c)
        assert_allclose(c, a.dot(b))

        blas.dgemm(1., a.T, b, 0., c)
        assert_allclose(c, a.T.dot(b))

        blas.dgemm(1., a, b.T, 0., c)
        assert_allclose(c, a.dot(b.T))

        blas.dgemm(1., a.T, b.T, 0., c)
        assert_allclose(c, a.T.dot(b.T))

        blas.dgemm(1., a, b, 0., c.T)
        assert_allclose(c, a.dot(b).T)

        blas.dgemm(1., a.T, b, 0., c.T)
        assert_allclose(c, a.T.dot(b).T)

        blas.dgemm(1., a, b.T, 0., c.T)
        assert_allclose(c, a.dot(b.T).T)

        blas.dgemm(1., a.T, b.T, 0., c.T)
        assert_allclose(c, a.T.dot(b.T).T)
    
    def test_shapes(self):
        a = np.arange(6, dtype='d').reshape((3, 2))
        b = np.arange(-6, 2, dtype='d').reshape((2, 4))
        c = np.empty((3, 4))

        blas.dgemm(1., a, b, 0., c)
        assert_allclose(c, a.dot(b))

        blas.dgemm(1., b.T, a.T, 0., c.T)
        assert_allclose(c, b.T.dot(a.T).T)
        
class test_wfunc_pointers(TestCase):
    """ Test the function pointers that are expected to fail on
    Mac OS X without the additional entry statement in their definitions
    in fblas_l1.pyf.src. """

    def test_complex_args(self):

        cx = np.array([.5 + 1.j, .25 - .375j, 12.5 - 4.j], np.complex64)
        cy = np.array([.8 + 2.j, .875 - .625j, -1. + 2.j], np.complex64)

        assert_allclose(blas.cdotc(cx, cy), -17.6468753815+21.3718757629j, 5)
        assert_allclose(blas.cdotu(cx, cy), -6.11562538147+30.3156242371j, 5)
        assert_allclose(blas.scasum(cx), 18.625, 5)
        assert_allclose(blas.scnrm2(cx), 13.1796483994, 5)

        assert_allclose(blas.cdotc(cx[::2], cy[::2]),
                        -18.1000003815+21.2000007629j, 5)
        assert_allclose(blas.cdotu(cx[::2], cy[::2]),
                        -6.10000038147+30.7999992371j, 5)
        assert_allclose(blas.scasum(cx[::2]), 18., 5)
        assert_allclose(blas.scnrm2(cx[::2]), 13.1719398499, 5)

    def test_float_args(self):

        x = np.array([5., -3, -.5], np.float32)
        y = np.array([2, 1, .5], np.float32)

        assert_allclose(blas.sasum(x), 8.5, 5)
        assert_allclose(blas.sdot(x, y), 6.75, 5)
        assert_allclose(blas.snrm2(x), 5.85234975815, 5)

        assert_allclose(blas.sasum(x[::2]), 5.5, 5)
        assert_allclose(blas.sdot(x[::2], y[::2]), 9.75, 5)
        assert_allclose(blas.snrm2(x[::2]), 5.0249376297, 5)

    def test_double_complex_args(self):

        cx = np.array([.5 + 1.j, .25 - .375j, 13. - 4.j], np.complex128)
        cy = np.array([.875 + 2.j, .875 - .625j, -1. + 2.j], np.complex128)

        assert_allclose(blas.zdotc(cx, cy), -18.109375+22.296875j, 5)
        assert_allclose(blas.zdotu(cx, cy), -6.578125+31.390625j, 5)

        assert_allclose(blas.zdotc(cx[::2], cy[::2]), -18.5625+22.125j, 5)
        assert_allclose(blas.zdotu(cx[::2], cy[::2]), -6.5625+31.875j, 5)

if __name__ == "__main__":
    run_module_suite()
