# This is added as a separate file to demonstrate how users can run standard
# distribution tests on custom distributions; it relies only on `DistributionsTest`,
# no other definitions in `test_continuous.py`.
import pytest
import numpy as np
from scipy import stats, special
from .test_continuous import DistributionsTest


class MyNormal:
    __make_distribution_version__ = "1.16.0"
    parameters = {'u': {'endpoints': (-np.inf, np.inf), 'typical': (-1, 1)},
                    's': {'endpoints': (0, np.inf), 'typical': (0.5, 2)}}
    support = {'endpoints': (-np.inf, np.inf), 'typical': (-3, 3)}

    def pdf(self, x, u, s):
        return 1 / np.sqrt(2*np.pi) / s * np.exp(-((x-u)/s)**2/2)

    def cdf(self, x, u, s):
        return special.ndtr((x-u)/s)

    def logcdf(self, x, u, s):
        return special.log_ndtr((x-u)/s)

    def icdf(self, p, u, s):
        return special.ndtri(p)*s + u

    def entropy(self, u, s):
        return np.full_like(u, fill_value=np.log(2*np.pi*np.e*s**2)/2)


class TestMyNormal(DistributionsTest):
    family = stats.make_distribution(MyNormal())
    seed =7694871135

    @pytest.mark.xslow
    def test_lmoment(self, case):
        return super().test_lmoment()
