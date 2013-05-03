#!/usr/bin/env python
#
# Created by: Pearu Peterson, March 2002
#
""" Test functions for scipy.linalg.matfuncs module

"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy import array, eye, dot, sqrt, double, exp, random
from numpy.testing import (TestCase, run_module_suite,
        assert_array_almost_equal, assert_array_almost_equal_nulp,
        assert_allclose)

from scipy.sparse import csc_matrix
from scipy.sparse.construct import eye as speye
from scipy.sparse.linalg import expm
from scipy.sparse.linalg.matfuncs import expm_2009
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

    def test_overscaling_example(self):
        # See the blog post
        # http://blogs.mathworks.com/cleve/2012/07/23/a-balancing-act-for-the-matrix-exponential/
        a = 2e10
        b = 4e8/6.
        c = 200/3.
        d = 3
        e = 1e-8
        A = np.array([[0,e,0],[-(a+b), -d, a], [c, 0, -c]])

        # This answer is wrong, and it is caused by overscaling.
        wrong_solution = np.array([
            [1.7465684381715e+17, -923050477.783131, -1.73117355055901e+17],
            [-3.07408665108297e+25, 1.62463553675545e+17, 3.04699053651329e+25],
            [1.09189154376804e+17, -577057840.468934, -1.08226721572342e+17]])

        # This is the correct answer, to not great precision.
        correct_solution = np.array([
            [0.446849, 1.54044e-09, 0.462811],
            [-5743070, -0.015283, -4526540],
            [0.447723, 1.5427e-09, 0.463481]])

        # Assert that the Higham 2005 expm gives the wrong answer.
        assert_allclose(expm(A), wrong_solution)

        # Assert that the Higham 2009 expm gives the correct answer.
        assert_allclose(expm_2009(A), correct_solution, rtol=1e-4)


if __name__ == "__main__":
    run_module_suite()
