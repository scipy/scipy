from __future__ import division, print_function, absolute_import

import warnings

from numpy.testing import assert_equal, assert_allclose

from scipy.misc import pade, logsumexp, face, ascent
from scipy.special import logsumexp as sc_logsumexp


def test_logsumexp():
    # make sure logsumexp can be imported from either scipy.misc or
    # scipy.special
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DeprecationWarning)
        assert_allclose(logsumexp([0, 1]),
                        sc_logsumexp([0, 1]), atol=1e-16)


def test_pade():
    # make sure scipy.misc.pade exists
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DeprecationWarning)
        pade([1, 2], 1)


def test_face():
    assert_equal(face().shape, (768, 1024, 3))


def test_ascent():
    assert_equal(ascent().shape, (512, 512))
