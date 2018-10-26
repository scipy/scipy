"""
Tests for Cython Optimize API
=============================

``test_zeros`` - tests ``newton`` and ``bisect``.
"""

from __future__ import division, print_function, absolute_import
import numpy as np
from . import zeros_examples, zeros_alt_examples


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
                       list(zeros_examples.test_cython_bisect()))


def test_zeros_alt_cython_bisect():
    assert np.allclose(EXPECTED_BISECT,
                       list(zeros_alt_examples.test_cython_bisect()))


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
                       list(zeros_examples.test_cython_ridder()))


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
                       list(zeros_examples.test_cython_brenth()))


# test brentq
def test_zeros_cython_brentq():
    assert np.allclose(EXPECTED_BRENT,
                       list(zeros_examples.test_cython_brentq()))
