from __future__ import division, print_function, absolute_import

from numpy.testing import assert_equal, assert_

from scipy.misc import pade, logsumexp, face, ascent
from scipy.special import logsumexp as sc_logsumexp
from scipy.interpolate import pade as i_pade


def test_logsumexp():
    # make sure logsumexp can be imported from either scipy.misc or
    # scipy.special
    assert_(logsumexp is sc_logsumexp)


def test_pade():
    assert_(pade is i_pade)


def test_face():
    assert_equal(face().shape, (768, 1024, 3))


def test_ascent():
    assert_equal(ascent().shape, (512, 512))
