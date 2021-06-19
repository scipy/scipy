import numpy as np
from .common import Benchmark, safe_import

with safe_import():
    from scipy import stats
with safe_import():
    from scipy.stats._distr_params import distdiscrete


# Simple Custom Distribution
class contdist1:
    def pdf(self, x):
        return 3/4 * (1-x*x) if abs(x) <= 1 else 0
    def dpdf(self, x):
        return 3/4 * (-2*x) if abs(x) <= 1 else 0
    def cdf(self, x):
        return 3/4 * (x - x**3/3 + 2/3) if abs(x) <= 1 else 1*(x >= 1)


# Standard Normal Distribution
class contdist2:
    def pdf(self, x):
        return stats.norm._pdf(x)
    def dpdf(self, x):
        return -x * stats.norm._pdf(x)
    def cdf(self, x):
        return stats.norm._cdf(x)


# pdf with piecewise linear function as transformed density with T = -1/sqrt
# Taken from UNU.RAN test suite (from file t_tdr_ps.c)
class contdist3:
    def pdf(self, x):
        y = 1. / (abs(x) + 1.)
        return y * y
    def dpdf(self, x):
        y = 1. / (abs(x) + 1.)
        y = 2. * y * y * y
        return y if (x < 0.) else -y
    def cdf(self, x):
        if x <= 0.:
            return 0.5 / (1. - x)
        else:
            return 1. - 0.5 / (1. + x)


allcontdists = [contdist1(), contdist2(), contdist3()]


class TransformedDensityRejection(Benchmark):

    param_names = ['dist']

    params = [allcontdists]

    def setup(self, dist):
        self.urng = np.random.default_rng(0xfaad7df1c89e050200dbe258636b3265)
        self.rng = stats.TransformedDensityRejection(dist, seed=self.urng)

    def time_tdr_setup(self, dist):
        rng = stats.TransformedDensityRejection(dist, seed=self.urng)

    def time_tdr_rvs(self, dist):
        rvs = self.rng.rvs(100000)


class DiscreteAliasUrn(Benchmark):

    params = ['distribution']

    params = [distdiscrete]

    def setup(self, distribution):
        distname, params = distribution
        if not isinstance(distname, str):
            dist, params = distname, params
        else:
            dist = getattr(stats, distname)
        domain = dist.support(*params)
        if not np.isfinite(domain[1] - domain[0]):
            raise NotImplementedError("Skipped")
        self.urng = np.random.default_rng(0x2fc9eb71cd5120352fa31b7a048aa867)
        x = np.arange(domain[0], domain[1] + 1)
        self.pv = dist.pmf(x, *params)
        self.rng = stats.DiscreteAliasUrn(self.pv, seed=self.urng)

    def time_dau_setup(self, distribution):
        rng = stats.DiscreteAliasUrn(self.pv, seed=self.urng)

    def time_dau_rvs(self, distribution):
        rvs = self.rng.rvs(100000)
