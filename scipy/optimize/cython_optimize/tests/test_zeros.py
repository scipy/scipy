"""
test cython optimize zeros api functions
"""

from __future__ import division, print_function, absolute_import
import numpy as np
from . import zeros_examples


EXPECTED_BISECT = [
    1.259765625,
    1.279296875,
    1.298828125,
    1.318359375,
    1.337890625,
    1.357421875,
    1.376953125,
    1.392578125,
    1.408203125,
    1.427734375]


# test bisect
def test_bisect():
    assert np.allclose(EXPECTED_BISECT,
                       list(zeros_examples.loop_example('bisect')))


EXPECTED_RIDDER = [
    1.2588478785767947,
    1.2795040615075954,
    1.299514441316524,
    1.318927124420269,
    1.3377847289304623,
    1.356125211864335,
    1.373982543571637,
    1.3913872624129802,
    1.408366934614972,
    1.4249465383291897]


# test ridder
def test_ridder():
    assert np.allclose(EXPECTED_RIDDER,
                       list(zeros_examples.loop_example('ridder')))


EXPECTED_BRENTH = [
    1.2599013632116787,
    1.2805183764606847,
    1.300478153836079,
    1.319836062211477,
    1.3386396735339707,
    1.3569302110890944,
    1.3747436816264305,
    1.392111769620191,
    1.4090625496526121,
    1.4256210584575042]


# test brenth
def test_brenth():
    assert np.allclose(EXPECTED_BRENTH,
                       list(zeros_examples.loop_example('brenth')))


EXPECTED_BRENTQ = [
    1.259872799904563,
    1.28042866862737,
    1.3003083443276644,
    1.3195775152132223,
    1.3382912269178606,
    1.3576742112124918,
    1.375418396530433,
    1.3927268676453197,
    1.4096289546536032,
    1.4261502005552444]


# test brentq
def test_brentq():
    assert np.allclose(EXPECTED_BRENTQ,
                       list(zeros_examples.loop_example('brentq')))
