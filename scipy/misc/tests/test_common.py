from __future__ import division, print_function, absolute_import

import pytest
from numpy.testing import assert_equal, assert_allclose, assert_almost_equal
from scipy._lib._numpy_compat import suppress_warnings

from scipy.misc import pade, logsumexp, face, ascent
from scipy.special import logsumexp as sc_logsumexp


def test_logsumexp():
    # make sure logsumexp can be imported from either scipy.misc or
    # scipy.special
    with suppress_warnings() as sup:
        sup.filter(DeprecationWarning, "`logsumexp` is deprecated")
        assert_allclose(logsumexp([0, 1]), sc_logsumexp([0, 1]), atol=1e-16)


def test_pade():
    # make sure scipy.misc.pade exists
    with suppress_warnings() as sup:
        sup.filter(DeprecationWarning, "`pade` is deprecated")
        pade([1, 2], 1)


def test_face():
    with suppress_warnings() as sup:
        sup.filter(DeprecationWarning, "`face` is deprecated")
        assert_equal(face().shape, (768, 1024, 3))


def test_ascent():
    with suppress_warnings() as sup:
        sup.filter(DeprecationWarning, "`ascent` is deprecated")
        assert_equal(ascent().shape, (512, 512))
