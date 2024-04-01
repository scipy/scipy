import numpy as np
import pytest  # type: ignore[import]

from scipy.special import _test_internal  # type: ignore[import]


@pytest.mark.skipif(not _test_internal.have_fenv(), reason="no fenv()")
def test_add_round_up():
    np.random.seed(1234)
    _test_internal.test_add_round(10**5, 'up')


@pytest.mark.skipif(not _test_internal.have_fenv(), reason="no fenv()")
def test_add_round_down():
    np.random.seed(1234)
    _test_internal.test_add_round(10**5, 'down')
