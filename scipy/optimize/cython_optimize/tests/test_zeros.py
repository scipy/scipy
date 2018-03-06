"""
Tests for Cython Optimize API
=============================

``test_zeros`` - tests ``newton`` and ``bisect``.
"""

import numpy as np
from scipy.optimize.cython_optimize.tests import zeros_examples


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


def test_zeros_cython_newton():
    assert np.allclose(EXPECTED_NEWTON, zeros_examples.test_cython_newton())


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

def test_zeros_cython_bisect():
    assert np.allclose(EXPECTED_BISECT, zeros_examples.test_cython_bisect())
