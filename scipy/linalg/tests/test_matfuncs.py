#!/usr/bin/env python
#
# Created by: Pearu Peterson, March 2002
#
""" Test functions for linalg.matfuncs module

"""

from __future__ import division, print_function, absolute_import

import random

import numpy as np
from numpy import array, identity, dot, sqrt, double
from numpy.testing import (TestCase, run_module_suite,
        assert_array_almost_equal, assert_array_almost_equal_nulp,
        assert_allclose, assert_, assert_raises, decorators)

import scipy.linalg
from scipy.linalg import funm, signm, logm, sqrtm, expm, expm_frechet
from scipy.linalg import fractional_matrix_power
from scipy.linalg.matfuncs import expm2, expm3
from scipy.linalg import _matfuncs_inv_ssq
import scipy.linalg._expm_frechet


class TestSignM(TestCase):

    def test_nils(self):
        a = array([[29.2, -24.2, 69.5, 49.8, 7.],
                   [-9.2, 5.2, -18., -16.8, -2.],
                   [-10., 6., -20., -18., -2.],
                   [-9.6, 9.6, -25.5, -15.4, -2.],
                   [9.8, -4.8, 18., 18.2, 2.]])
        cr = array([[11.94933333,-2.24533333,15.31733333,21.65333333,-2.24533333],
                    [-3.84266667,0.49866667,-4.59066667,-7.18666667,0.49866667],
                    [-4.08,0.56,-4.92,-7.6,0.56],
                    [-4.03466667,1.04266667,-5.59866667,-7.02666667,1.04266667],
                    [4.15733333,-0.50133333,4.90933333,7.81333333,-0.50133333]])
        r = signm(a)
        assert_array_almost_equal(r,cr)

    def test_defective1(self):
        a = array([[0.0,1,0,0],[1,0,1,0],[0,0,0,1],[0,0,1,0]])
        r = signm(a, disp=False)
        #XXX: what would be the correct result?

    def test_defective2(self):
        a = array((
            [29.2,-24.2,69.5,49.8,7.0],
            [-9.2,5.2,-18.0,-16.8,-2.0],
            [-10.0,6.0,-20.0,-18.0,-2.0],
            [-9.6,9.6,-25.5,-15.4,-2.0],
            [9.8,-4.8,18.0,18.2,2.0]))
        r = signm(a, disp=False)
        #XXX: what would be the correct result?

    def test_defective3(self):
        a = array([[-2., 25., 0., 0., 0., 0., 0.],
                   [0., -3., 10., 3., 3., 3., 0.],
                   [0., 0., 2., 15., 3., 3., 0.],
                   [0., 0., 0., 0., 15., 3., 0.],
                   [0., 0., 0., 0., 3., 10., 0.],
                   [0., 0., 0., 0., 0., -2., 25.],
                   [0., 0., 0., 0., 0., 0., -3.]])
        r = signm(a, disp=False)
        #XXX: what would be the correct result?


class TestLogM(TestCase):

    def test_nils(self):
        a = array([[-2., 25., 0., 0., 0., 0., 0.],
                   [0., -3., 10., 3., 3., 3., 0.],
                   [0., 0., 2., 15., 3., 3., 0.],
                   [0., 0., 0., 0., 15., 3., 0.],
                   [0., 0., 0., 0., 3., 10., 0.],
                   [0., 0., 0., 0., 0., -2., 25.],
                   [0., 0., 0., 0., 0., 0., -3.]])
        m = (identity(7)*3.1+0j)-a
        logm(m, disp=False)
        #XXX: what would be the correct result?

    def _get_al_mohy_higham_2012_experiment_1(self):
        A = np.array([
            [3.2346e-1, 3e4, 3e4, 3e4],
            [0, 3.0089e-1, 3e4, 3e4],
            [0, 0, 3.221e-1, 3e4],
            [0, 0, 0, 3.0744e-1]], dtype=float)
        return A

    def test_al_mohy_higham_2012_experiment_1_logm(self):
        # The logm completes the round trip successfully.
        A = self._get_al_mohy_higham_2012_experiment_1()
        A_logm = logm(A)
        A_round_trip = expm(A_logm)
        assert_allclose(A_round_trip, A, rtol=1e-5)

    def test_al_mohy_higham_2012_experiment_1_funm_log(self):
        # The raw funm with np.log does not complete the round trip.
        A = self._get_al_mohy_higham_2012_experiment_1()
        A_funm_log = funm(A, np.log)
        A_funm_log_round_trip = expm(A_funm_log)
        assert_(not np.allclose(A_funm_log_round_trip, A, rtol=1e-5))

    def test_round_trip_random_float(self):
        np.random.seed(1234)
        for n in range(2, 6):
            M_unscaled = np.random.randn(n, n)
            for scale in np.logspace(-4, 4, 9):
                M = M_unscaled * scale
                W = np.linalg.eigvals(M)
                # We cannot take logm when an eigenvalue is real and negative.
                if any(w.real == w and w < 0 for w in W):
                    assert_raises(ValueError, logm, M)
                else:
                    M_round_trip = expm(logm(M))
                    assert_allclose(M_round_trip, M)

    def test_round_trip_random_complex(self):
        np.random.seed(1234)
        for n in range(2, 6):
            M_unscaled = np.random.randn(n, n) + 1j * np.random.randn(n, n)
            for scale in np.logspace(-4, 4, 9):
                M = M_unscaled * scale
                M_round_trip = expm(logm(M))
                assert_allclose(M_round_trip, M)


class TestSqrtM(TestCase):
    def test_bad(self):
        # See http://www.maths.man.ac.uk/~nareports/narep336.ps.gz
        e = 2**-5
        se = sqrt(e)
        a = array([[1.0,0,0,1],
                   [0,e,0,0],
                   [0,0,e,0],
                   [0,0,0,1]])
        sa = array([[1,0,0,0.5],
                    [0,se,0,0],
                    [0,0,se,0],
                    [0,0,0,1]])
        n = a.shape[0]
        assert_array_almost_equal(dot(sa,sa),a)
        # Check default sqrtm.
        esa = sqrtm(a, disp=False, blocksize=n)[0]
        assert_array_almost_equal(dot(esa,esa),a)
        # Check sqrtm with 2x2 blocks.
        esa = sqrtm(a, disp=False, blocksize=2)[0]
        assert_array_almost_equal(dot(esa,esa),a)

    def test_blocksizes(self):
        # Make sure I do not goof up the blocksizes when they do not divide n.
        np.random.seed(1234)
        for n in range(1, 8):
            A = np.random.rand(n, n) + 1j*np.random.randn(n, n)
            A_sqrtm_default, info = sqrtm(A, disp=False, blocksize=n)
            assert_allclose(A, np.linalg.matrix_power(A_sqrtm_default, 2))
            for blocksize in range(1, 10):
                A_sqrtm_new, info = sqrtm(A, disp=False, blocksize=blocksize)
                assert_allclose(A_sqrtm_default, A_sqrtm_new)


class TestFractionalMatrixPower(TestCase):
    def test_fractional_matrix_power_round_trip(self):
        np.random.seed(1234)
        for p in range(1, 5):
            for n in range(1, 5):
                M_unscaled = np.random.randn(n, n) + 1j * np.random.randn(n, n)
                for scale in np.logspace(-4, 4, 9):
                    M = M_unscaled * scale
                    M_root = fractional_matrix_power(M, 1/p)
                    M_round_trip = np.linalg.matrix_power(M_root, p)
                    assert_allclose(M, M_round_trip)

    def test_larger_abs_fractional_matrix_powers(self):
        np.random.seed(1234)
        for n in (2, 3, 5):
            for i in range(10):
                M = np.random.randn(n, n) + 1j * np.random.randn(n, n)
                M_one_fifth = fractional_matrix_power(M, 0.2)
                # Test the round trip.
                M_round_trip = np.linalg.matrix_power(M_one_fifth, 5)
                assert_allclose(M, M_round_trip)
                # Test a large abs fractional power.
                X = fractional_matrix_power(M, -5.4)
                Y = np.linalg.matrix_power(M_one_fifth, -27)
                assert_allclose(X, Y)
                # Test another large abs fractional power.
                X = fractional_matrix_power(M, 3.8)
                Y = np.linalg.matrix_power(M_one_fifth, 19)
                assert_allclose(X, Y)

    def test_briggs_helper_function(self):
        np.random.seed(1234)
        for a in np.random.randn(10) + 1j * np.random.randn(10):
            for k in range(5):
                x_observed = _matfuncs_inv_ssq._briggs_helper_function(a, k)
                x_expected = a ** np.exp2(-k) - 1
                assert_allclose(x_observed, x_expected)


class TestExpM(TestCase):
    def test_zero(self):
        a = array([[0.,0],[0,0]])
        assert_array_almost_equal(expm(a),[[1,0],[0,1]])
        assert_array_almost_equal(expm2(a),[[1,0],[0,1]])
        assert_array_almost_equal(expm3(a),[[1,0],[0,1]])

    def test_consistency(self):
        a = array([[0.,1],[-1,0]])
        assert_array_almost_equal(expm(a), expm2(a))
        assert_array_almost_equal(expm(a), expm3(a))

        a = array([[1j,1],[-1,-2j]])
        assert_array_almost_equal(expm(a), expm2(a))
        assert_array_almost_equal(expm(a), expm3(a))


class TestExpmFrechet(TestCase):

    def test_expm_frechet(self):
        # a test of the basic functionality
        M = np.array([
            [1, 2, 3, 4],
            [5, 6, 7, 8],
            [0, 0, 1, 2],
            [0, 0, 5, 6],
            ], dtype=float)
        A = np.array([
            [1, 2],
            [5, 6],
            ], dtype=float)
        E = np.array([
            [3, 4],
            [7, 8],
            ], dtype=float)
        expected_expm = scipy.linalg.expm(A)
        expected_frechet = scipy.linalg.expm(M)[:2, 2:]
        for kwargs in ({}, {'method':'SPS'}, {'method':'blockEnlarge'}):
            observed_expm, observed_frechet = expm_frechet(A, E, **kwargs)
            assert_allclose(expected_expm, observed_expm)
            assert_allclose(expected_frechet, observed_frechet)

    def test_small_norm_expm_frechet(self):
        # methodically test matrices with a range of norms, for better coverage
        M_original = np.array([
            [1, 2, 3, 4],
            [5, 6, 7, 8],
            [0, 0, 1, 2],
            [0, 0, 5, 6],
            ], dtype=float)
        A_original = np.array([
            [1, 2],
            [5, 6],
            ], dtype=float)
        E_original = np.array([
            [3, 4],
            [7, 8],
            ], dtype=float)
        A_original_norm_1 = scipy.linalg.norm(A_original, 1)
        selected_m_list = [1, 3, 5, 7, 9, 11, 13, 15]
        m_neighbor_pairs = zip(selected_m_list[:-1], selected_m_list[1:])
        for ma, mb in m_neighbor_pairs:
            ell_a = scipy.linalg._expm_frechet.ell_table_61[ma]
            ell_b = scipy.linalg._expm_frechet.ell_table_61[mb]
            target_norm_1 = 0.5 * (ell_a + ell_b)
            scale = target_norm_1 / A_original_norm_1
            M = scale * M_original
            A = scale * A_original
            E = scale * E_original
            expected_expm = scipy.linalg.expm(A)
            expected_frechet = scipy.linalg.expm(M)[:2, 2:]
            observed_expm, observed_frechet = expm_frechet(A, E)
            assert_allclose(expected_expm, observed_expm)
            assert_allclose(expected_frechet, observed_frechet)

    def test_fuzz(self):
        # try a bunch of crazy inputs
        rfuncs = (
                np.random.uniform,
                np.random.normal,
                np.random.standard_cauchy,
                np.random.exponential)
        ntests = 100
        for i in range(ntests):
            rfunc = random.choice(rfuncs)
            target_norm_1 = random.expovariate(1.0)
            n = random.randrange(2, 16)
            A_original = rfunc(size=(n,n))
            E_original = rfunc(size=(n,n))
            A_original_norm_1 = scipy.linalg.norm(A_original, 1)
            scale = target_norm_1 / A_original_norm_1
            A = scale * A_original
            E = scale * E_original
            M = np.vstack([
                np.hstack([A, E]),
                np.hstack([np.zeros_like(A), A])])
            expected_expm = scipy.linalg.expm(A)
            expected_frechet = scipy.linalg.expm(M)[:n, n:]
            observed_expm, observed_frechet = expm_frechet(A, E)
            assert_allclose(expected_expm, observed_expm)
            assert_allclose(expected_frechet, observed_frechet)

    def test_problematic_matrix(self):
        # this test case uncovered a bug which has since been fixed
        A = np.array([
                [1.50591997, 1.93537998],
                [0.41203263, 0.23443516],
                ], dtype=float)
        E = np.array([
                [1.87864034, 2.07055038],
                [1.34102727, 0.67341123],
                ], dtype=float)
        A_norm_1 = scipy.linalg.norm(A, 1)
        sps_expm, sps_frechet = expm_frechet(
                A, E, method='SPS')
        blockEnlarge_expm, blockEnlarge_frechet = expm_frechet(
                A, E, method='blockEnlarge')
        assert_allclose(sps_expm, blockEnlarge_expm)
        assert_allclose(sps_frechet, blockEnlarge_frechet)

    @decorators.slow
    @decorators.skipif(True, 'this test is deliberately slow')
    def test_medium_matrix(self):
        # profile this to see the speed difference
        n = 1000
        A = np.random.exponential(size=(n, n))
        E = np.random.exponential(size=(n, n))
        sps_expm, sps_frechet = expm_frechet(
                A, E, method='SPS')
        blockEnlarge_expm, blockEnlarge_frechet = expm_frechet(
                A, E, method='blockEnlarge')
        assert_allclose(sps_expm, blockEnlarge_expm)
        assert_allclose(sps_frechet, blockEnlarge_frechet)


if __name__ == "__main__":
    run_module_suite()
