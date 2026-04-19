"""Tests for Part 1: sklearn-like wrappers with internal .pxd and int->blas_int."""

import numpy as np
from numpy.testing import assert_allclose
import scipy.linalg.cython_blas as cython_blas

from ilp64_test_package._consumer import (
    consumer_ddot,
    consumer_daxpy,
    consumer_dnrm2,
    consumer_dgemm,
    consumer_dgetrf,
    consumer_blas_int_size,
)


class TestBlasIntSize:
    def test_blas_int_matches_scipy(self):
        expected = cython_blas._blas_int_size()
        assert consumer_blas_int_size() == expected


class TestBLASLevel1:
    def test_ddot(self):
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([4.0, 5.0, 6.0])
        result = consumer_ddot(x, y)
        assert_allclose(result, np.dot(x, y))

    def test_daxpy(self):
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([4.0, 5.0, 6.0])
        consumer_daxpy(2.0, x, y)
        assert_allclose(y, [6.0, 9.0, 12.0])     

    def test_dnrm2_large_vector(self):
        x = np.zeros(2**31, dtype=np.float64)
        x[-1] = 1.0
        nrm2 = consumer_dnrm2(x)
        # LP64 answer due to blas_int->int casts in consumer_dnrm2
        # cf a similar test in test_direct.py
        assert_allclose(nrm2, 0.0, atol=1e-14)


class TestBLASLevel3:
    def test_dgemm_identity(self):
        a = np.eye(3, dtype=np.float64, order='F')
        b = np.arange(9, dtype=np.float64).reshape(3, 3, order='F')
        c = np.empty((3, 3), dtype=np.float64, order='F')
        consumer_dgemm(1.0, a, b, 0.0, c)
        assert_allclose(c, b)

    def test_dgemm_product(self):
        a = np.array([[1, 2], [3, 4]], dtype=np.float64, order='F')
        b = np.array([[5, 6], [7, 8]], dtype=np.float64, order='F')
        c = np.zeros((2, 2), dtype=np.float64, order='F')
        consumer_dgemm(1.0, a, b, 0.0, c)
        assert_allclose(c, a @ b)


class TestLAPACK:
    def test_dgetrf(self):
        a = np.array([[2, 1], [1, 3]], dtype=np.float64, order='F')
        ipiv = np.empty(2, dtype=np.intc)
        info = consumer_dgetrf(a, ipiv)
        assert info == 0

    def test_dgetrf_rectangular(self):
        a = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float64, order='F')
        ipiv = np.empty(2, dtype=np.intc)
        info = consumer_dgetrf(a, ipiv)
        assert info == 0
