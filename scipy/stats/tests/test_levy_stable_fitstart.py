# Test scipy.stats.levy_stable._fitstart().

import numpy as np
import pytest

from scipy.stats import levy_stable
from scipy.stats import rv_continuous as rvc

class TestLevy_stable:
    def test_neg_beta(self):
        # Confirm that _fitstart() gets same estimate of alpha for +beta and -beta.
        rvc.random_state = 2021
        NPOINTS = 1001
        ALPHA = 1.5
        BETA_POS = 1
        
        draw = levy_stable.rvs(alpha=ALPHA, beta=BETA_POS, size=NPOINTS)
        # Fit stable distribution to positive and negative mirror image data, examine alpha.
        alpha_pos, alpha_neg = [levy_stable._fitstart(data)[0] for data in [draw, -draw]]
        
        assert np.isclose(alpha_pos, alpha_neg)
        
    def test_neg_delta(self):
        # Confirm that _fitstart() returns symmetric +delta and -delta if beta flips.
        rvc.random_state = 2021
        NPOINTS = 1001
        ALPHA = 1.5
        BETA_POS = 1
        
        draw = levy_stable.rvs(alpha=ALPHA, beta=BETA_POS, size=NPOINTS)
        # Fit stable distribution to positive and negative mirror image data, examine location.
        delta_pos, delta_neg = [levy_stable._fitstart(data)[2] for data in [draw, -draw]]
        
        assert np.isclose(delta_pos, -delta_neg)
        
    def test_delta_shift(self):
        # Confirm that _fitstart() returns delta that slides up and down if data shifts.
        rvc.random_state = 2021
        NPOINTS = 1001
        ALPHA = 1.5
        BETA_POS = 1
        LOC = -2
        SHIFT = 1
        
        draw = levy_stable.rvs(alpha=ALPHA, beta=BETA_POS, loc=LOC, size=NPOINTS)
        # Fit stable distribution to data and shifted data, examine location parameter.
        delta0, delta_up = [levy_stable._fitstart(draw + x)[2] for x in [0, SHIFT]]
        
        assert np.isclose(delta0, delta_up-SHIFT)
        
    def test_alpha1(self):
        # Confirm that stable distribution can fit alpha==1 without blowing up.
        SHIFT = 0
        SCALE = 1
                
        # Make a dataset whose quintiles lead to alpha==1.
        x5 = x25 = x50 = 5  # Percentiles 5, 10, 25 -> 5.
        x95 = 100 - x5  # Percentile 95.
        # Calculate value for percentile 75 that will interpolate alpha lookup to 1.
        #   We assume beta==1 so interpolation is between rows 4 and 5 of last column
        #   of _fitstart() alpha_table.
        x75 = x25 + (x95 - x5) / (4 + (1.150 - 1) / (1.150 - .973))
        #assert x5 <= x25 <= x50 <= x75 <= x95
        
        xx = np.array([x5, x25, x50, x75, x95])
        # np.percentile(xx, [5, 25, 50, 75, 95]) -> [5, 5, 5, ..., 95]
        
        # Place target percentiles either side.
        data = np.array([x5, x5,  # Percentiles 0, 10.
                         x25, x25,  # 20, 30.
                         x50, x50, x50,  # 40, 50, 60.
                         x75, x75,  # 70, 80.
                         x95, x95  # 90, 100.
                         ])
        #assert np.array_equal(np.percentile(data, [5, 25, 50, 75, 95]), xx)
        
        alpha, beta, loc, scale = fitstart_fix._fitstart_fix(SCALE*data + SHIFT)
        
        # Not sure what value of delta should be returned when alpha==1 so
        #  for now just verify that alpha==1.
        assert alpha == 1
        # Might want to check loc too--what value should it have when alpha==1?
        
#%%

if __name__ == "__main__":
    t = TestLevy_stable()
    t.test_neg_beta()
    t.test_neg_delta()
    t.test_delta_shift()
    t.test_alpha1()
