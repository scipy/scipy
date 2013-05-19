#!/usr/bin/env python
#
# Created by: Pearu Peterson, March 2002
#
""" Test functions for scipy.linalg.matfuncs module

"""
from __future__ import division, print_function, absolute_import

import math

import numpy as np
from numpy import array, eye, dot, sqrt, double, exp, random
from numpy.linalg import matrix_power
from numpy.testing import (TestCase, run_module_suite,
        assert_allclose, assert_, assert_array_almost_equal,
        assert_array_almost_equal_nulp)

from scipy.sparse import csc_matrix
from scipy.sparse.construct import eye as speye
from scipy.sparse.linalg.matfuncs import (expm_2005, expm_2009,
        ProductOperator, MatrixPowerOperator,
        _is_upper_triangular)
from scipy.linalg import logm
import scipy.sparse
import scipy.sparse.linalg


class TestExpM(TestCase):
    def test_zero(self):
        a = array([[0.,0],[0,0]])
        for expm in (expm_2005, expm_2009):
            assert_array_almost_equal(expm(a),[[1,0],[0,1]])

    def test_zero_sparse(self):
        a = csc_matrix([[0.,0],[0,0]])
        for expm in (expm_2005, expm_2009):
            assert_array_almost_equal(expm(a).toarray(),[[1,0],[0,1]])

    def test_bidiagonal_sparse(self):
        A = csc_matrix([
            [1, 3, 0],
            [0, 1, 5],
            [0, 0, 2]], dtype=float)
        e1 = math.exp(1)
        e2 = math.exp(2)
        expected = np.array([
            [e1, 3*e1, 15*(e2 - 2*e1)],
            [0, e1, 5*(e2 - e1)],
            [0, 0, e2]], dtype=float)
        for expm in (expm_2005, expm_2009):
            observed = expm(A).toarray()
            assert_array_almost_equal(observed, expected)

    def test_padecases_dtype_float(self):
        for expm in (expm_2005, expm_2009):
            for dtype in [np.float32, np.float64]:
                for scale in [1e-2, 1e-1, 5e-1, 1, 10]:
                    A = scale * eye(3, dtype=dtype)
                    observed = expm(A)
                    expected = exp(scale) * eye(3, dtype=dtype)
                    assert_array_almost_equal_nulp(observed, expected, nulp=100)

    def test_padecases_dtype_complex(self):
        for expm in (expm_2005, expm_2009):
            for dtype in [np.complex64, np.complex128]:
                for scale in [1e-2, 1e-1, 5e-1, 1, 10]:
                    a = scale * eye(3, dtype=dtype)
                    e = exp(scale) * eye(3, dtype=dtype)
                    assert_array_almost_equal_nulp(expm(a), e, nulp=100)

    def test_padecases_dtype_sparse_float(self):
        # float32 and complex64 lead to errors in spsolve/UMFpack
        dtype = np.float64
        for scale in [1e-2, 1e-1, 5e-1, 1, 10]:
            a = scale * speye(3, 3, dtype=dtype, format='csc')
            e = exp(scale) * eye(3, dtype=dtype)
            for expm in (expm_2005, expm_2009):
                assert_array_almost_equal_nulp(expm(a).toarray(), e, nulp=100)

    def test_padecases_dtype_sparse_complex(self):
        # float32 and complex64 lead to errors in spsolve/UMFpack
        dtype = np.complex128
        for scale in [1e-2, 1e-1, 5e-1, 1, 10]:
            a = scale * speye(3, 3, dtype=dtype, format='csc')
            e = exp(scale) * eye(3, dtype=dtype)
            for expm in (expm_2005, expm_2009):
                assert_array_almost_equal_nulp(expm(a).toarray(), e, nulp=100)

    def test_logm_consistency(self):
        random.seed(1234)
        for dtype in [np.float64, np.complex128]:
            for n in range(1, 10):
                for scale in [1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2]:
                    # make logm(A) be of a given scale
                    A = (eye(n) + random.rand(n, n) * scale).astype(dtype)
                    if np.iscomplexobj(A):
                        A = A + 1j * random.rand(n, n) * scale
                    for expm in (expm_2005, expm_2009):
                        assert_array_almost_equal(expm(logm(A)), A)


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

        # This is the correct answer.
        correct_solution = np.array([
            [0.446849468283175, 1.54044157383952e-09, 0.462811453558774],
            [-5743067.77947947, -0.0152830038686819, -4526542.71278401],
            [0.447722977849494, 1.54270484519591e-09, 0.463480648837651]])

        # Assert that the Higham 2005 expm gives the wrong answer.
        assert_allclose(expm_2005(A), wrong_solution)

        # Assert that the Higham 2009 expm gives the correct answer.
        assert_allclose(expm_2009(A), correct_solution)

    def test_expm_2009_random_upper_triangular(self):
        random.seed(1234)
        n = 10
        nsamples = 20
        for i in range(nsamples):
            A = np.triu(np.random.randn(n, n))
            assert_(_is_upper_triangular(A))
            assert_allclose(expm_2009(A), expm_2005(A))


class TestOperators(TestCase):

    def test_product_operator(self):
        random.seed(1234)
        n = 5
        k = 2
        nsamples = 10
        for i in range(nsamples):
            A = np.random.randn(n, n)
            B = np.random.randn(n, n)
            C = np.random.randn(n, n)
            D = np.random.randn(n, k)
            op = ProductOperator(A, B, C)
            assert_allclose(op.matmat(D), A.dot(B).dot(C).dot(D))
            assert_allclose(op.T.matmat(D), (A.dot(B).dot(C)).T.dot(D))

    def test_matrix_power_operator(self):
        random.seed(1234)
        n = 5
        k = 2
        p = 3
        nsamples = 10
        for i in range(nsamples):
            A = np.random.randn(n, n)
            B = np.random.randn(n, k)
            op = MatrixPowerOperator(A, p)
            assert_allclose(op.matmat(B), matrix_power(A, p).dot(B))
            assert_allclose(op.T.matmat(B), matrix_power(A, p).T.dot(B))


if __name__ == "__main__":
    run_module_suite()
