from numpy.testing import TestCase, run_module_suite, assert_allclose
from scipy.linalg import cython_lapack as cython_lapack
from scipy.linalg import lapack


class test_lamch(TestCase):

    def test_slamch(self):
        for c in [b'e', b's', b'b', b'p', b'n', b'r', b'm', b'u', b'l', b'o']:
            assert_allclose(cython_lapack._test_slamch(c),
                            lapack.slamch(c))

    def test_dlamch(self):
        for c in [b'e', b's', b'b', b'p', b'n', b'r', b'm', b'u', b'l', b'o']:
            assert_allclose(cython_lapack._test_dlamch(c),
                            lapack.dlamch(c))


if __name__ == "__main__":
    run_module_suite()
