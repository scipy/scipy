import numpy as np
from .common import Benchmark, safe_import

with safe_import():
    from scipy import stats
with safe_import():
    from scipy.stats import sampling
with safe_import():
    from scipy import special


# Beta distribution with a = 2, b = 3
class contdist1:
    def __init__(self):
        self.mode = 1/3

    def pdf(self, x):
        return 12 * x * (1-x)**2

    def dpdf(self, x):
        return 12 * ((1-x)**2 - 2*x*(1-x))

    def cdf(self, x):
        return 12 * (x**2/2 - x**3/3 + x**4/4)

    def support(self):
        return 0, 1

    def __repr__(self):
        # asv prints this.
        return 'beta(2, 3)'


# Standard Normal Distribution
class contdist2:
    def __init__(self):
        self.mode = 0

    def pdf(self, x):
        return 1./np.sqrt(2*np.pi) * np.exp(-0.5 * x*x)

    def dpdf(self, x):
        return 1./np.sqrt(2*np.pi) * -x * np.exp(-0.5 * x*x)

    def cdf(self, x):
        return special.ndtr(x)

    def __repr__(self):
        return 'norm(0, 1)'


# pdf with piecewise linear function as transformed density with T = -1/sqrt
# Taken from UNU.RAN test suite (from file t_tdr_ps.c)
class contdist3:
    def __init__(self, shift=0.):
        self.shift = shift
        self.mode = shift

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


# Sin 2 distribution
#          /  0.05 + 0.45*(1 +sin(2 Pi x))  if |x| <= 1
#  f(x) = <
#          \  0        otherwise
# Taken from UNU.RAN test suite (from file t_pinv.c)
class contdist4:
    def __init__(self):
        self.mode = 0

    def pdf(self, x):
        return 0.05 + 0.45 * (1 + np.sin(2*np.pi*x))

    def dpdf(self, x):
        return 0.2 * 0.45 * (2*np.pi) * np.cos(2*np.pi*x)

    def cdf(self, x):
        return (0.05*(x + 1) +
                0.9*(1. + 2.*np.pi*(1 + x) - np.cos(2.*np.pi*x)) /
                (4.*np.pi))

    def support(self):
        return -1, 1

    def __repr__(self):
        return 'sin2'


# Sin 10 distribution
#          /  0.05 + 0.45*(1 +sin(2 Pi x))  if |x| <= 5
#  f(x) = <
#          \  0        otherwise
# Taken from UNU.RAN test suite (from file t_pinv.c)
class contdist5:
    def __init__(self):
        self.mode = 0

    def pdf(self, x):
        return 0.2 * (0.05 + 0.45 * (1 + np.sin(2*np.pi*x)))

    def dpdf(self, x):
        return 0.2 * 0.45 * (2*np.pi) * np.cos(2*np.pi*x)

    def cdf(self, x):
        return x/10. + 0.5 + 0.09/(2*np.pi) * (np.cos(10*np.pi) -
                                               np.cos(2*np.pi*x))

    def support(self):
        return -5, 5

    def __repr__(self):
        return 'sin10'


allcontdists = [contdist1(), contdist2(), contdist3(), contdist3(10000.),
                contdist4(), contdist5()]


class TransformedDensityRejection(Benchmark):

    param_names = ['dist', 'c']

    params = [allcontdists, [0., -0.5]]

    def setup(self, dist, c):
        self.urng = np.random.default_rng(0xfaad7df1c89e050200dbe258636b3265)
        with np.testing.suppress_warnings() as sup:
            sup.filter(RuntimeWarning)
            try:
                self.rng = sampling.TransformedDensityRejection(
                    dist, c=c, random_state=self.urng
                )
            except sampling.UNURANError:
                # contdist3 is not T-concave for c=0. So, skip such test-cases
                raise NotImplementedError(f"{dist} not T-concave for c={c}")

    def time_tdr_setup(self, dist, c):
        with np.testing.suppress_warnings() as sup:
            sup.filter(RuntimeWarning)
            sampling.TransformedDensityRejection(
                dist, c=c, random_state=self.urng
            )

    def time_tdr_rvs(self, dist, c):
        self.rng.rvs(100000)


class SimpleRatioUniforms(Benchmark):

    param_names = ['dist', 'cdf_at_mode']

    params = [allcontdists, [0, 1]]

    def setup(self, dist, cdf_at_mode):
        self.urng = np.random.default_rng(0xfaad7df1c89e050200dbe258636b3265)
        try:
            if cdf_at_mode:
                cdf_at_mode = dist.cdf(dist.mode)
            else:
                cdf_at_mode = None
            self.rng = sampling.SimpleRatioUniforms(
                dist, mode=dist.mode,
                cdf_at_mode=cdf_at_mode,
                random_state=self.urng
            )
        except sampling.UNURANError:
            raise NotImplementedError(f"{dist} not T-concave")

    def time_srou_setup(self, dist, cdf_at_mode):
        if cdf_at_mode:
            cdf_at_mode = dist.cdf(dist.mode)
        else:
            cdf_at_mode = None
        sampling.SimpleRatioUniforms(
            dist, mode=dist.mode,
            cdf_at_mode=cdf_at_mode,
            random_state=self.urng
        )

    def time_srou_rvs(self, dist, cdf_at_mode):
        self.rng.rvs(100000)


class NumericalInversePolynomial(Benchmark):

    param_names = ['dist']

    params = [allcontdists]

    def setup(self, dist):
        self.urng = np.random.default_rng(0xb235b58c1f616c59c18d8568f77d44d1)
        with np.testing.suppress_warnings() as sup:
            sup.filter(RuntimeWarning)
            try:
                self.rng = sampling.NumericalInversePolynomial(
                    dist, random_state=self.urng
                )
            except sampling.UNURANError:
                raise NotImplementedError(f"setup failed for {dist}")

    def time_pinv_setup(self, dist):
        with np.testing.suppress_warnings() as sup:
            sup.filter(RuntimeWarning)
            sampling.NumericalInversePolynomial(
                dist, random_state=self.urng
            )

    def time_pinv_rvs(self, dist):
        self.rng.rvs(100000)


class NumericalInverseHermite(Benchmark):

    param_names = ['dist', 'order']
    params = [allcontdists, [3, 5]]

    def setup(self, dist, order):
        self.urng = np.random.default_rng(0xb235b58c1f616c59c18d8568f77d44d1)
        with np.testing.suppress_warnings() as sup:
            sup.filter(RuntimeWarning)
            try:
                self.rng = sampling.NumericalInverseHermite(
                    dist, order=order, random_state=self.urng
                )
            except sampling.UNURANError:
                raise NotImplementedError(f"setup failed for {dist}")

    def time_hinv_setup(self, dist, order):
        with np.testing.suppress_warnings() as sup:
            sup.filter(RuntimeWarning)
            sampling.NumericalInverseHermite(
                dist, order=order, random_state=self.urng
            )

    def time_hinv_rvs(self, dist, order):
        self.rng.rvs(100000)


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
        self.rng = sampling.DiscreteAliasUrn(self.pv, random_state=self.urng)

    def time_dau_setup(self, distribution):
        sampling.DiscreteAliasUrn(self.pv, random_state=self.urng)

    def time_dau_rvs(self, distribution):
        self.rng.rvs(100000)


class DiscreteGuideTable(Benchmark):

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
        self.rng = sampling.DiscreteGuideTable(self.pv, random_state=self.urng)

    def time_dgt_setup(self, distribution):
        sampling.DiscreteGuideTable(self.pv, random_state=self.urng)

    def time_dgt_rvs(self, distribution):
        self.rng.rvs(100000)
