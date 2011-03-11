import numpy as np
from numpy.testing import assert_array_equal, assert_almost_equal, \
                          assert_array_almost_equal

from scipy.misc import pade, logsumexp


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
    assert_array_almost_equal(denomp.c, [-1.0/6, 0.5,  -1.0, 1.0])

def test_logsumexp():
    """Test whether logsumexp() function correctly handles large inputs."""
    a = np.arange(200)
    desired = np.log(np.sum(np.exp(a)))
    assert_almost_equal(logsumexp(a), desired)

    # Now test with large numbers
    b = [1000,1000]
    desired = 1000.0 + np.log(2.0)
    assert_almost_equal(logsumexp(b), desired)

    n = 1000
    b = np.ones(n)*10000
    desired = 10000.0 + np.log(n)
    assert_almost_equal(logsumexp(b), desired)
