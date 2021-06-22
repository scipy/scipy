import numpy as np
from .common import Benchmark, safe_import

with safe_import():
    from scipy import stats
with safe_import():
    from scipy.stats._distr_params import distdiscrete


# Beta distribution with a = 2, b = 3
class contdist1:
    def pdf(self, x):
        if 0 < x < 1:
            return x * (1-x)**2
        return 0
    def dpdf(self, x):
        if 0 < x < 1:
            return (1-x)**2 - 2*x*(1-x)
        return 0
    def cdf(self, x):
        if x < 0:
            return 0
        if x > 1:
            return 1
        return stats.beta._cdf(x, 2, 3)
    def __repr__(self):
        # asv prints this.
        return 'beta(2, 3)'


# Standard Normal Distribution
class contdist2:
    def pdf(self, x):
        return stats.norm._pdf(x)
    def dpdf(self, x):
        return -x * stats.norm._pdf(x)
    def cdf(self, x):
        return stats.norm._cdf(x)
    def __repr__(self):
        return 'norm(0, 1)'


# pdf with piecewise linear function as transformed density with T = -1/sqrt
# Taken from UNU.RAN test suite (from file t_tdr_ps.c)
class contdist3:
    def __init__(self, shift=0.):
        self.shift = shift
    def pdf(self, x):
        x -= self.shift
        y = 1. / (abs(x) + 1.)
        return y * y
    def dpdf(self, x):
        x -= self.shift
        y = 1. / (abs(x) + 1.)
        y = 2. * y * y * y
        return y if (x < 0.) else -y
    def cdf(self, x):
        x -= self.shift
        if x <= 0.:
            return 0.5 / (1. - x)
        return 1. - 0.5 / (1. + x)
    def __repr__(self):
        return f'sqrtlinshft({self.shift})'


allcontdists = [contdist1(), contdist2(), contdist3(), contdist3(10000.)]


class TransformedDensityRejection(Benchmark):

    param_names = ['dist', 'c', 'cpoints']

    params = [allcontdists, [0., -0.5], [10, 20, 30, 50]]

    def setup(self, dist, c, cpoints):
        self.urng = np.random.default_rng(0xfaad7df1c89e050200dbe258636b3265)
        with np.testing.suppress_warnings() as sup:
            sup.filter(UserWarning)
            try:
                self.rng = stats.TransformedDensityRejection(dist, c=c,
                                                             cpoints=cpoints,
                                                             seed=self.urng)
            except RuntimeError:
                # contdist3 is not T-concave for c=0. So, skip such test-cases
                raise NotImplementedError(f"{dist} not T-concave for c={c}")

    def time_tdr_setup(self, dist, c, cpoints):
        with np.testing.suppress_warnings() as sup:
            sup.filter(UserWarning)
            rng = stats.TransformedDensityRejection(dist, c=c,
                                                    cpoints=cpoints,
                                                    seed=self.urng)

    def time_tdr_rvs(self, dist, c, cpoints):
        rvs = self.rng.rvs(100000)


class DiscreteAliasUrn(Benchmark):

    param_names = ['distribution']

    params = [
        # a subset of discrete distributions with finite domain.
        [['nhypergeom', (20, 7, 1)],
         ['hypergeom', (30, 12, 6)],
         ['nchypergeom_wallenius', (140, 80, 60, 0.5)],
         ['binom', (5, 0.4)]]
    ]

    def setup(self, distribution):
        distname, params = distribution
        dist = getattr(stats, distname)
        domain = dist.support(*params)
        self.urng = np.random.default_rng(0x2fc9eb71cd5120352fa31b7a048aa867)
        x = np.arange(domain[0], domain[1] + 1)
        self.pv = dist.pmf(x, *params)
        self.rng = stats.DiscreteAliasUrn(self.pv, seed=self.urng)

    def time_dau_setup(self, distribution):
        rng = stats.DiscreteAliasUrn(self.pv, seed=self.urng)

    def time_dau_rvs(self, distribution):
        rvs = self.rng.rvs(100000)
