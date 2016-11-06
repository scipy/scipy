from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (assert_array_equal, assert_almost_equal,
                           assert_array_almost_equal, assert_equal, assert_)

from scipy.misc import pade, logsumexp, face, ascent
from scipy.special import logsumexp as sc_logsumexp

def test_logsumexp():
    # make sure logsumexp can be imported from either scipy.misc or
    # scipy.special
    assert_(logsumexp is sc_logsumexp)


def test_pade_trivial():
    nump, denomp = pade([1.0], 0)
    assert_array_equal(nump.c, [1.0])
    assert_array_equal(denomp.c, [1.0])


def test_pade_4term_exp():
    # First four Taylor coefficients of exp(x).
    # Unlike poly1d, the first array element is the zero-order term.
    an = [1.0, 1.0, 0.5, 1.0/6]

    nump, denomp = pade(an, 0)
    assert_array_almost_equal(nump.c, [1.0/6, 0.5, 1.0, 1.0])
    assert_array_almost_equal(denomp.c, [1.0])

    nump, denomp = pade(an, 1)
    assert_array_almost_equal(nump.c, [1.0/6, 2.0/3, 1.0])
    assert_array_almost_equal(denomp.c, [-1.0/3, 1.0])

    nump, denomp = pade(an, 2)
    assert_array_almost_equal(nump.c, [1.0/3, 1.0])
    assert_array_almost_equal(denomp.c, [1.0/6, -2.0/3, 1.0])

    nump, denomp = pade(an, 3)
    assert_array_almost_equal(nump.c, [1.0])
    assert_array_almost_equal(denomp.c, [-1.0/6, 0.5, -1.0, 1.0])


def test_face():
    assert_equal(face().shape, (768, 1024, 3))


def test_ascent():
    assert_equal(ascent().shape, (512, 512))
