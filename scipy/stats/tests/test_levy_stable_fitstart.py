import pytest

# Test scipy.stats.levy_stable._fitstart().

from numpy.testing import assert_allclose
from scipy.stats import levy_stable
from scipy.stats import rv_continuous as rvc

class TestLevy_stable:
    def test_neg_beta(self):
        # Confirm that _fitstart() gets same estimate of alpha for +beta and -beta.
        rvc.random_state = 2021
        draw = levy_stable.rvs(alpha=1.5, beta=1, size=1001)
        # Fit stable distribution to positive and negative mirror image data, examine alpha.
        alpha_pos, alpha_neg = [levy_stable._fitstart(data)[0] for data in [draw, -draw]]
        
        # assert np.isclose(alpha_pos, alpha_neg)
        assert_allclose(alpha_pos, alpha_neg)
        
    def test_neg_delta(self):
        # Confirm that _fitstart() returns symmetric +delta and -delta if beta flips.
        rvc.random_state = 2021
        draw = levy_stable.rvs(alpha=1.5, beta=1, size=1001)
        # Fit stable distribution to positive and negative mirror image data, examine location.
        delta_pos, delta_neg = [levy_stable._fitstart(data)[2] for data in [draw, -draw]]
        
        # assert np.isclose(delta_pos, -delta_neg)
        assert_allclose(delta_pos, -delta_neg)
        
    def test_delta_shift(self):
        # Confirm that _fitstart() returns delta that slides up and down if data shifts.
        rvc.random_state = 2021
        SHIFT = 1
        
        draw = levy_stable.rvs(alpha=1.5, beta=1, loc=-2, size=1001)
        # Fit stable distribution to data and shifted data, examine location parameter.
        delta0, delta_up = [levy_stable._fitstart(draw + x)[2] for x in [0, SHIFT]]
        
        # assert np.isclose(delta0, delta_up-SHIFT)
        assert_allclose(delta0, delta_up-SHIFT)
