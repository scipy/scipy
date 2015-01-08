from numpy.testing import TestCase, run_module_suite, assert_allclose
from scipy.linalg import _cython_lapack_wrappers as cython_lapack
from scipy.linalg import lapack


class test_lamch(TestCase):
    
    def test_slamch(self):
        for c in b'esbpnrmulo':
            assert_allclose(cython_lapack.slamch(c), lapack.slamch(c))
    
    def test_dlamch(self):
        for c in b'esbpnrmulo':
            assert_allclose(cython_lapack.slamch(c), lapack.slamch(c))


if __name__ == "__main__":
    run_module_suite()
