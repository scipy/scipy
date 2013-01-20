#!/usr/bin/env python
#
# Created by: Pearu Peterson, March 2002
#
""" Test functions for scipy.linalg.matfuncs module

"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy import array, eye, dot, sqrt, double, exp, random
from numpy.testing import TestCase, run_module_suite, assert_array_almost_equal, \
     assert_array_almost_equal_nulp

from scipy.sparse import csc_matrix
from scipy.sparse.construct import eye as speye
from scipy.sparse.linalg import expm
from scipy.linalg import logm

class TestExpM(TestCase):
    def test_zero(self):
        a = array([[0.,0],[0,0]])
        assert_array_almost_equal(expm(a),[[1,0],[0,1]])

    def test_zero_sparse(self):
        a = csc_matrix([[0.,0],[0,0]])
        assert_array_almost_equal(expm(a).toarray(),[[1,0],[0,1]])

    def test_padecases_dtype(self):
        for dtype in [np.float32, np.float64, np.complex64, np.complex128]:
            # test double-precision cases
            for scale in [1e-2, 1e-1, 5e-1, 1, 10]:
                a = scale * eye(3, dtype=dtype)
                e = exp(scale) * eye(3, dtype=dtype)
                assert_array_almost_equal_nulp(expm(a), e, nulp=100)

    def test_padecases_dtype_sparse(self):
        # float32 and complex64 lead to errors in spsolve/UMFpack
        for dtype in [np.float64, np.complex128]:
            # test double-precision cases
            for scale in [1e-2, 1e-1, 5e-1, 1, 10]:
                a = scale * speye(3, 3, dtype=dtype, format='csc')
                e = exp(scale) * eye(3, dtype=dtype)
                assert_array_almost_equal_nulp(expm(a).toarray(), e, nulp=100)

    def test_logm_consistency(self):
        random.seed(1234)
        for dtype in [np.float32, np.float64, np.complex64, np.complex128]:
            for n in range(1, 10):
                for scale in [1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2]:
                    # make logm(a) be of a given scale
                    a = (eye(n) + random.rand(n, n) * scale).astype(dtype)
                    if np.iscomplexobj(a):
                        a = a + 1j * random.rand(n, n) * scale
                    assert_array_almost_equal(expm(logm(a)), a)

if __name__ == "__main__":
    run_module_suite()
