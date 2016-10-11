from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import dec

from scipy.special import _test_round


@dec.skipif(not _test_round.have_fenv())
def test_add_round_up():
    np.random.seed(1234)
    _test_round.test_add_round(10**5, 'up')


@dec.skipif(not _test_round.have_fenv())
def test_add_round_down():
    np.random.seed(1234)
    _test_round.test_add_round(10**5, 'down')
