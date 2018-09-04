"""
Tests for Cython Optimize API
=============================

``test_zeros`` - tests ``newton`` and ``bisect``.
"""

from __future__ import division, print_function, absolute_import
import numpy as np
from ..examples import (zeros_tuple_examples, zeros_struct_examples,
    zeros_struct_alt_examples, zeros_array_examples)


EXPECTED_NEWTON = [
    5.255320079106907,
    6.093781591553449,
    6.161364035402291,
    5.395937648109023,
    4.501196600138535,
    4.299785923137901,
    4.976896720121363,
    5.909959792054446,
    6.241136571144537,
    5.66596871302032
]


# test newton
def test_zeros_cython_newton():
    assert np.allclose(EXPECTED_NEWTON,
                       list(zeros_tuple_examples.test_cython_newton()))


def test_zeros_struct_cython_newton():
    assert np.allclose(EXPECTED_NEWTON,
                       list(zeros_struct_examples.test_cython_newton()))


def test_zeros_struct_alt_cython_newton():
    assert np.allclose(EXPECTED_NEWTON,
                       list(zeros_struct_alt_examples.test_cython_newton()))


def test_zeros_array_cython_newton():
    assert np.allclose(EXPECTED_NEWTON,
                       list(zeros_array_examples.test_cython_newton()))


# test secant
def test_zeros_cython_secant():
    assert np.allclose(EXPECTED_NEWTON,
                       list(zeros_tuple_examples.test_cython_secant()))


def test_zeros_struct_cython_secant():
    assert np.allclose(EXPECTED_NEWTON,
                       list(zeros_struct_examples.test_cython_secant()))


def test_zeros_array_cython_secant():
    assert np.allclose(EXPECTED_NEWTON,
                       list(zeros_array_examples.test_cython_secant()))


# test halley
def test_zeros_cython_halley():
    assert np.allclose(EXPECTED_NEWTON,
                       list(zeros_tuple_examples.test_cython_halley()))


def test_zeros_struct_cython_halley():
    assert np.allclose(EXPECTED_NEWTON,
                       list(zeros_struct_examples.test_cython_halley()))


def test_zeros_array_cython_halley():
    assert np.allclose(EXPECTED_NEWTON,
                       list(zeros_array_examples.test_cython_halley()))


def test_zeros_cython_newton_full_output():
    full_output = zeros_tuple_examples.test_newton_full_output()
    assert full_output['error_num'] == 0
    assert full_output['flag'] == b'Converged successfully'
    assert full_output['funcalls'] == 6
    assert full_output['iterations'] == 3
    assert np.isclose(full_output['root'], 5.255320079106907)
    full_output = zeros_tuple_examples.test_newton_full_output(tol=-2)
    assert full_output['error_num'] == -1
    assert full_output['flag'] == b'TOL and MAXITER must be positive integers'
    assert full_output['funcalls'] == 0
    assert full_output['iterations'] == 0
    assert full_output['root'] == 6.0
    full_output = zeros_tuple_examples.test_newton_full_output(maxiter=0)
    assert full_output['error_num'] == -1
    assert full_output['flag'] == b'TOL and MAXITER must be positive integers'
    assert full_output['funcalls'] == 0
    assert full_output['iterations'] == 0
    assert full_output['root'] == 6.0
    full_output = zeros_tuple_examples.test_newton_full_output(v=25.0)
    assert full_output['error_num'] == -2
    assert full_output['flag'] == b'Failed to converge'
    assert full_output['funcalls'] == 100
    assert full_output['iterations'] == 50
    assert np.isclose(full_output['root'], -3425.9998163580435)


def test_zeros_cython_secant_full_output():
    full_output = zeros_tuple_examples.test_secant_full_output()
    assert full_output['error_num'] == 0
    assert full_output['flag'] == b'Converged successfully'
    assert full_output['funcalls'] == 4
    assert full_output['iterations'] == 3
    assert np.isclose(full_output['root'], 5.255320079106908)
    full_output = zeros_tuple_examples.test_secant_full_output(tol=0)
    assert full_output['error_num'] == -1
    assert full_output['flag'] == b'TOL and MAXITER must be positive integers'
    assert full_output['funcalls'] == 0
    assert full_output['iterations'] == 0
    assert full_output['root'] == 6.0
    full_output = zeros_tuple_examples.test_secant_full_output(maxiter=-3)
    assert full_output['error_num'] == -1
    assert full_output['flag'] == b'TOL and MAXITER must be positive integers'
    assert full_output['funcalls'] == 0
    assert full_output['iterations'] == 0
    assert full_output['root'] == 6.0
    full_output = zeros_tuple_examples.test_secant_full_output(v=25.0)
    assert full_output['error_num'] == -2
    assert full_output['flag'] == b'Failed to converge'
    assert full_output['funcalls'] == 52
    assert full_output['iterations'] == 50
    assert np.isclose(full_output['root'], -2388.496387739789)


EXPECTED_BISECT = [
    5.2568359375,
    6.0908203125,
    6.1591796875,
    5.400390625,
    4.5048828125,
    4.2998046875,
    4.9833984375,
    5.9130859375,
    6.2412109375,
    5.6669921875
]


# test bisect
def test_zeros_cython_bisect():
    assert np.allclose(EXPECTED_BISECT,
                       list(zeros_tuple_examples.test_cython_bisect()))


def test_zeros_struct_cython_bisect():
    assert np.allclose(EXPECTED_BISECT,
                       list(zeros_struct_examples.test_cython_bisect()))


def test_zeros_struct_alt_cython_bisect():
    assert np.allclose(EXPECTED_BISECT,
                       list(zeros_struct_alt_examples.test_cython_bisect()))


def test_zeros_array_cython_bisect():
    assert np.allclose(EXPECTED_BISECT,
                       list(zeros_array_examples.test_cython_bisect()))


EXPECTED_RIDDER = [
    5.258445977956501,
    6.097326910150778,
    6.1649432076344155,
    5.399133820263375,
    4.503945969500776,
    4.302434803463735,
    4.979883542842023,
    5.913413069001644,
    6.244755713677152,
    5.669299903167427
]


# test ridder
def test_zeros_cython_ridder():
    assert np.allclose(EXPECTED_RIDDER,
                       list(zeros_tuple_examples.test_cython_ridder()))


def test_zeros_struct_cython_ridder():
    assert np.allclose(EXPECTED_RIDDER,
                       list(zeros_struct_examples.test_cython_ridder()))


def test_zeros_array_cython_ridder():
    assert np.allclose(EXPECTED_RIDDER,
                       list(zeros_array_examples.test_cython_ridder()))


EXPECTED_BRENT = [
    5.255112621677981,
    6.093656129793108,
    6.1612466036887055,
    5.3957416733746015,
    4.5009430414747245,
    4.2995244435825075,
    4.9766692097106,
    5.909813564693917,
    6.241028891976474,
    5.665797347132429
]


# test brenth
def test_zeros_cython_brenth():
    assert np.allclose(EXPECTED_BRENT,
                       list(zeros_tuple_examples.test_cython_brenth()))


def test_zeros_struct_cython_brenth():
    assert np.allclose(EXPECTED_BRENT,
                       list(zeros_struct_examples.test_cython_brenth()))


def test_zeros_array_cython_brenth():
    assert np.allclose(EXPECTED_BRENT,
                       list(zeros_array_examples.test_cython_brenth()))


# test brentq
def test_zeros_cython_brentq():
    assert np.allclose(EXPECTED_BRENT,
                       list(zeros_tuple_examples.test_cython_brentq()))


def test_zeros_struct_cython_brentq():
    assert np.allclose(EXPECTED_BRENT,
                       list(zeros_struct_examples.test_cython_brentq()))


def test_zeros_array_cython_brentq():
    assert np.allclose(EXPECTED_BRENT,
                       list(zeros_array_examples.test_cython_brentq()))


def test_zeros_cython_bisect_full_output():
    full_output = zeros_tuple_examples.test_bisect_full_output()
    assert full_output['error_num'] == -2
    assert full_output['flag'] == b'Failed to converge'
    assert full_output['funcalls'] == 12
    assert full_output['iterations'] == 10
    assert np.isclose(full_output['root'], 5.2568359375)
    full_output = zeros_tuple_examples.test_bisect_full_output(xa=0.0, xb=1.0)
    assert full_output['error_num'] == -1
    assert full_output['flag'] == b'F(XA) and F(XB) must have opposite signs'
    assert full_output['funcalls'] == 0
    assert full_output['iterations'] == 0
    assert full_output['root'] == 0.0
    full_output = zeros_tuple_examples.test_bisect_full_output(mitr=-1)
    assert full_output['error_num'] == -2
    assert full_output['flag'] == b'Failed to converge'
    assert full_output['funcalls'] == 2
    assert full_output['iterations'] == 0
    assert full_output['root'] == 7.0
    full_output = zeros_tuple_examples.test_bisect_full_output(mitr=15)
    assert full_output['error_num'] == 0
    assert full_output['flag'] == b'Converged successfully'
    assert full_output['funcalls'] == 13
    assert full_output['iterations'] == 11
    assert np.isclose(full_output['root'], 5.25341796875)


def test_zeros_cython_brentq_full_output():
    full_output = zeros_tuple_examples.test_brentq_full_output()
    assert full_output['error_num'] == 0
    assert full_output['flag'] == b'Converged successfully'
    assert full_output['funcalls'] == 4
    assert full_output['iterations'] == 3
    assert np.isclose(full_output['root'], 5.255112621677981)
    full_output = zeros_tuple_examples.test_brentq_full_output(xa=8.0, xb=10.0)
    assert full_output['error_num'] == -1
    assert full_output['flag'] == b'F(XA) and F(XB) must have opposite signs'
    assert full_output['funcalls'] == 0
    assert full_output['iterations'] == 0
    assert full_output['root'] == 0.0
    full_output = zeros_tuple_examples.test_brentq_full_output(v=9.0, xa=-3000.0, xb=0.0)
    assert full_output['error_num'] == -2
    assert full_output['flag'] == b'Failed to converge'
    assert full_output['funcalls'] == 12
    assert full_output['iterations'] == 10
    assert full_output['root'] == -413.02279389470186
    full_output = zeros_tuple_examples.test_brentq_full_output(xtol=-1, rtol=-1)
    assert full_output['error_num'] == 0
    assert full_output['flag'] == b'Converged successfully'
    assert full_output['funcalls'] == 5
    assert full_output['iterations'] == 4
    assert np.isclose(full_output['root'], 5.255320079106907)
