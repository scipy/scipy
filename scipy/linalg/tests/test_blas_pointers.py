from numpy.testing import run_module_suite

from scipy.linalg._test_blas_pointers import (test_dgemm_pointer,
                                             test_wfunc_pointers)


if __name__ == "__main__":
    run_module_suite()
