import warnings

import numpy as np
from .common import Benchmark, safe_import, is_xslow

with safe_import():
    import scipy.stats as stats
with safe_import():
    from scipy.stats._distr_params import distcont, distdiscrete

try:  # builtin lib
    from itertools import compress
except ImportError:
    pass


class Anderson_KSamp(Benchmark):
    def setup(self, *args):
        self.rand = [np.random.normal(loc=i, size=1000) for i in range(3)]

    def time_anderson_ksamp(self):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            stats.anderson_ksamp(self.rand)


class CorrelationFunctions(Benchmark):
    param_names = ['alternative']
    params = [
        ['two-sided', 'less', 'greater']
    ]

    def setup(self, mode):
        a = np.random.rand(2,2) * 10
        self.a = a

    def time_fisher_exact(self, alternative):
        oddsratio, pvalue = stats.fisher_exact(self.a, alternative=alternative)

    def time_barnard_exact(self, alternative):
        resBarnard = stats.barnard_exact(self.a, alternative=alternative)

    def time_boschloo_exact(self, alternative):
        resBoschloo = stats.boschloo_exact(self.a, alternative=alternative)


class ANOVAFunction(Benchmark):
    def setup(self):
        rng = np.random.default_rng(12345678)
        self.a = rng.random((6,3)) * 10
        self.b = rng.random((6,3)) * 10
        self.c = rng.random((6,3)) * 10

    def time_f_oneway(self):
        statistic, pvalue = stats.f_oneway(self.a, self.b, self.c)
        statistic, pvalue = stats.f_oneway(self.a, self.b, self.c, axis=1)


class Kendalltau(Benchmark):
    param_names = ['nan_policy','method','variant']
    params = [
        ['propagate', 'raise', 'omit'],
        ['auto', 'asymptotic', 'exact'],
        ['b', 'c']
    ]

    def setup(self, nan_policy, method, variant):
        rng = np.random.default_rng(12345678)
        a = np.arange(200)
        rng.shuffle(a)
        b = np.arange(200)
        rng.shuffle(b)
        self.a = a
        self.b = b

    def time_kendalltau(self, nan_policy, method, variant):
        tau, p_value = stats.kendalltau(self.a, self.b, nan_policy=nan_policy, method=method, variant=variant)


class InferentialStats(Benchmark):
    def setup(self):
        rng = np.random.default_rng(12345678)
        self.a = stats.norm.rvs(loc=5, scale=10, size=500, random_state=rng)
        self.b = stats.norm.rvs(loc=8, scale=10, size=20, random_state=rng)
        self.c = stats.norm.rvs(loc=8, scale=20, size=20, random_state=rng)

    def time_ttest_ind_same_var(self):
        # test different sized sample with variances
        stats.ttest_ind(self.a, self.b)
        stats.ttest_ind(self.a, self.b, equal_var=False)

    def time_ttest_ind_diff_var(self):
        # test different sized sample with different variances
        stats.ttest_ind(self.a, self.c)
        stats.ttest_ind(self.a, self.c, equal_var=False)


class DistributionsAll(Benchmark):
    # all distributions are in this list. A conversion to a set is used to
    # remove duplicates that appear more than once in either `distcont` or
    # `distdiscrete`.
    dists = sorted(list(set([d[0] for d in distcont + distdiscrete])))

    param_names = ['dist_name', 'method']
    params = [
        dists, ['pdf/pmf', 'logpdf/logpmf', 'cdf', 'logcdf', 'rvs', 'fit',
                'sf', 'logsf', 'ppf', 'isf', 'moment', 'stats_s', 'stats_v',
                'stats_m', 'stats_k', 'stats_mvsk', 'entropy']
    ]
    # stats_mvsk is tested separately because of gh-11742
    # `moment` tests a higher moment (order 5)

    dist_data = dict(distcont + distdiscrete)
    # custom shape values can be provided for any distribution in the format
    # `dist_name`: [shape1, shape2, ...]
    custom_input = {}

    # these are the distributions that are the slowest
    slow_dists = ['nct', 'ncx2', 'argus', 'cosine', 'foldnorm', 'gausshyper',
                  'kappa4', 'invgauss', 'wald', 'vonmises_line', 'ksone',
                  'genexpon', 'exponnorm', 'recipinvgauss', 'vonmises',
                  'foldcauchy', 'kstwo', 'levy_stable', 'skewnorm']
    slow_methods = ['moment']

    def setup(self, dist_name, method):
        if not is_xslow() and (dist_name in self.slow_dists
                               or method in self.slow_methods):
            raise NotImplementedError("Skipped")

        self.dist = getattr(stats, dist_name)

        dist_shapes = self.dist_data[dist_name]

        if isinstance(self.dist, stats.rv_discrete):
            # discrete distributions only use location
            self.isCont = False
            kwds = {'loc': 4}
        else:
            # continuous distributions use location and scale
            self.isCont = True
            kwds = {'loc': 4, 'scale': 10}

        bounds = self.dist.interval(.99, *dist_shapes, **kwds)
        x = np.linspace(*bounds, 100)
        args = [x, *self.custom_input.get(dist_name, dist_shapes)]
        self.args = args
        self.kwds = kwds
        if method == 'fit':
            # there are no fit methods for discrete distributions
            if isinstance(self.dist, stats.rv_discrete):
                raise NotImplementedError("This attribute is not a member "
                                          "of the distribution")
            # the only positional argument is the data to be fitted
            self.args = [self.dist.rvs(*dist_shapes, size=100, random_state=0, **kwds)]
        elif method == 'rvs':
            # add size keyword argument for data creation
            kwds['size'] = 1000
            kwds['random_state'] = 0
            # keep shapes as positional arguments, omit linearly spaced data
            self.args = args[1:]
        elif method == 'pdf/pmf':
            method = ('pmf' if isinstance(self.dist, stats.rv_discrete)
                      else 'pdf')
        elif method == 'logpdf/logpmf':
            method = ('logpmf' if isinstance(self.dist, stats.rv_discrete)
                      else 'logpdf')
        elif method in ['ppf', 'isf']:
            self.args = [np.linspace((0, 1), 100), *args[1:]]
        elif method == 'moment':
            # the first four moments may be optimized, so compute the fifth
            self.args = [5, *args[1:]]
        elif method.startswith('stats_'):
            kwds['moments'] = method[6:]
            method = 'stats'
            self.args = args[1:]
        elif method == 'entropy':
            self.args = args[1:]

        self.method = getattr(self.dist, method)

    def time_distribution(self, dist_name, method):
        self.method(*self.args, **self.kwds)


class Distribution(Benchmark):
    # though there is a new version of this benchmark that runs all the
    # distributions, at the time of writing there was odd behavior on
    # the asv for this benchmark, so it is retained.
    # https://pv.github.io/scipy-bench/#stats.Distribution.time_distribution

    param_names = ['distribution', 'properties']
    params = [
        ['cauchy', 'gamma', 'beta'],
        ['pdf', 'cdf', 'rvs', 'fit']
    ]

    def setup(self, distribution, properties):
        rng = np.random.default_rng(12345678)
        self.x = rng.random(100)

    def time_distribution(self, distribution, properties):
        if distribution == 'gamma':
            if properties == 'pdf':
                stats.gamma.pdf(self.x, a=5, loc=4, scale=10)
            elif properties == 'cdf':
                stats.gamma.cdf(self.x, a=5, loc=4, scale=10)
            elif properties == 'rvs':
                stats.gamma.rvs(size=1000, a=5, loc=4, scale=10)
            elif properties == 'fit':
                stats.gamma.fit(self.x, loc=4, scale=10)
        elif distribution == 'cauchy':
            if properties == 'pdf':
                stats.cauchy.pdf(self.x, loc=4, scale=10)
            elif properties == 'cdf':
                stats.cauchy.cdf(self.x, loc=4, scale=10)
            elif properties == 'rvs':
                stats.cauchy.rvs(size=1000, loc=4, scale=10)
            elif properties == 'fit':
                stats.cauchy.fit(self.x, loc=4, scale=10)
        elif distribution == 'beta':
            if properties == 'pdf':
                stats.beta.pdf(self.x, a=5, b=3, loc=4, scale=10)
            elif properties == 'cdf':
                stats.beta.cdf(self.x, a=5, b=3, loc=4, scale=10)
            elif properties == 'rvs':
                stats.beta.rvs(size=1000, a=5, b=3, loc=4, scale=10)
            elif properties == 'fit':
                stats.beta.fit(self.x, loc=4, scale=10)

    # Retain old benchmark results (remove this if changing the benchmark)
    time_distribution.version = "fb22ae5386501008d945783921fe44aef3f82c1dafc40cddfaccaeec38b792b0"


class DescriptiveStats(Benchmark):
    param_names = ['n_levels']
    params = [
        [10, 1000]
    ]

    def setup(self, n_levels):
        rng = np.random.default_rng(12345678)
        self.levels = rng.integers(n_levels, size=(1000, 10))

    def time_mode(self, n_levels):
        stats.mode(self.levels, axis=0)


class GaussianKDE(Benchmark):
    def setup(self):
        rng = np.random.default_rng(12345678)
        n = 2000
        m1 = rng.normal(size=n)
        m2 = rng.normal(scale=0.5, size=n)

        xmin = m1.min()
        xmax = m1.max()
        ymin = m2.min()
        ymax = m2.max()

        X, Y = np.mgrid[xmin:xmax:200j, ymin:ymax:200j]
        self.positions = np.vstack([X.ravel(), Y.ravel()])
        values = np.vstack([m1, m2])
        self.kernel = stats.gaussian_kde(values)

    def time_gaussian_kde_evaluate_few_points(self):
        # test gaussian_kde evaluate on a small number of points
        self.kernel(self.positions[:, :10])

    def time_gaussian_kde_evaluate_many_points(self):
        # test gaussian_kde evaluate on many points
        self.kernel(self.positions)


class GroupSampling(Benchmark):
    param_names = ['dim']
    params = [[3, 10, 50, 200]]

    def setup(self, dim):
        self.rng = np.random.default_rng(12345678)

    def time_unitary_group(self, dim):
        stats.unitary_group.rvs(dim, random_state=self.rng)

    def time_ortho_group(self, dim):
        stats.ortho_group.rvs(dim, random_state=self.rng)

    def time_special_ortho_group(self, dim):
        stats.special_ortho_group.rvs(dim, random_state=self.rng)


class BinnedStatisticDD(Benchmark):

    params = ["count", "sum", "mean", "min", "max", "median", "std", np.std]

    def setup(self, statistic):
        rng = np.random.default_rng(12345678)
        self.inp = rng.random(9999).reshape(3, 3333) * 200
        self.subbin_x_edges = np.arange(0, 200, dtype=np.float32)
        self.subbin_y_edges = np.arange(0, 200, dtype=np.float64)
        self.ret = stats.binned_statistic_dd(
            [self.inp[0], self.inp[1]], self.inp[2], statistic=statistic,
            bins=[self.subbin_x_edges, self.subbin_y_edges])

    def time_binned_statistic_dd(self, statistic):
        stats.binned_statistic_dd(
            [self.inp[0], self.inp[1]], self.inp[2], statistic=statistic,
            bins=[self.subbin_x_edges, self.subbin_y_edges])

    def time_binned_statistic_dd_reuse_bin(self, statistic):
        stats.binned_statistic_dd(
            [self.inp[0], self.inp[1]], self.inp[2], statistic=statistic,
            binned_statistic_result=self.ret)


class ContinuousFitAnalyticalMLEOverride(Benchmark):
    # list of distributions to time
    dists = ["pareto", "laplace", "rayleigh",
             "invgauss", "gumbel_r", "gumbel_l"]
    # add custom values for rvs and fit, if desired, for any distribution:
    # key should match name in dists and value should be list of loc, scale,
    # and shapes
    custom_input = {}
    fnames = ['floc', 'fscale', 'f0', 'f1', 'f2']
    fixed = {}
    distcont = dict(distcont)

    param_names = ["distribution", "loc_fixed", "scale_fixed",
                   "shape1_fixed", "shape2_fixed", "shape3_fixed"]
    params = [dists, * [[True, False]] * 5]

    def setup(self, dist_name, loc_fixed, scale_fixed, shape1_fixed,
              shape2_fixed, shape3_fixed):
        self.distn = eval("stats." + dist_name)

        # default `loc` and `scale` are .834 and 4.342, and shapes are from
        # `_distr_params.py`
        default_shapes = self.distcont[dist_name]
        param_values = self.custom_input.get(dist_name, [.834, 4.342,
                                                         *default_shapes])
        # separate relevant and non-relevant parameters for this distribution
        # based on the number of shapes
        nparam = len(param_values)
        all_parameters = [loc_fixed, scale_fixed, shape1_fixed, shape2_fixed,
                          shape3_fixed]
        relevant_parameters = all_parameters[:nparam]
        nonrelevant_parameters = all_parameters[nparam:]

        # skip if all parameters are fixed or if non relevant parameters are
        # not all false
        if True in nonrelevant_parameters or False not in relevant_parameters:
            raise NotImplementedError("skip non-relevant case")

        # add fixed values if fixed in relevant_parameters to self.fixed
        # with keys from self.fnames and values from parameter_values
        self.fixed = dict(zip(compress(self.fnames, relevant_parameters),
                          compress(param_values, relevant_parameters)))
        self.data = self.distn.rvs(*param_values, size=1000)

    def time_fit(self, dist_name, loc_fixed, scale_fixed, shape1_fixed,
                 shape2_fixed, shape3_fixed):
        self.distn.fit(self.data, **self.fixed)


class BenchMoment(Benchmark):
    params = [
        [1, 2, 3, 8],
        [100, 1000, 10000],
    ]
    param_names = ["order", "size"]

    def setup(self, order, size):
        np.random.random(1234)
        self.x = np.random.random(size)

    def time_moment(self, order, size):
        stats.moment(self.x, order)


class BenchSkewKurtosis(Benchmark):
    params = [
        [1, 2, 3, 8],
        [100, 1000, 10000],
        [False, True]
    ]
    param_names = ["order", "size", "bias"]

    def setup(self, order, size, bias):
        np.random.random(1234)
        self.x = np.random.random(size)

    def time_skew(self, order, size, bias):
        stats.skew(self.x, bias=bias)

    def time_kurtosis(self, order, size, bias):
        stats.kurtosis(self.x, bias=bias)


class BenchQMCDiscrepancy(Benchmark):
    param_names = ['method']
    params = [
        ["CD", "WD", "MD", "L2-star",]
    ]

    def setup(self, method):
        rng = np.random.default_rng(1234)
        sample = rng.random((1000, 10))
        self.sample = sample

    def time_discrepancy(self, method):
        disc = stats.qmc.discrepancy(self.sample, method=method)


class NumericalInverseHermite(Benchmark):

    param_names = ['distribution']
    params = [distcont]

    def setup(self, *args):
        self.rand = [np.random.normal(loc=i, size=1000) for i in range(3)]

    def time_fni(self, distcase):
        distname, shapes = distcase
        slow_dists = {'ksone', 'kstwo', 'levy_stable', 'skewnorm'}
        fail_dists = {'beta', 'gausshyper', 'geninvgauss', 'ncf', 'nct',
                      'norminvgauss', 'genhyperbolic', 'studentized_range'}

        if distname in slow_dists or distname in fail_dists:
            raise NotImplementedError("skipped")

        dist = getattr(stats, distname)(*shapes)

        with np.testing.suppress_warnings() as sup:
            sup.filter(RuntimeWarning, "overflow encountered")
            sup.filter(RuntimeWarning, "divide by zero")
            sup.filter(RuntimeWarning, "invalid value encountered")
            stats.NumericalInverseHermite(dist)


class DistanceFunctions(Benchmark):
    param_names = ['n_size']
    params = [
        [10, 4000]
    ]

    def setup(self, n_size):
        rng = np.random.default_rng(12345678)
        self.u_values = rng.random(n_size) * 10
        self.u_weights = rng.random(n_size) * 10
        self.v_values = rng.random(n_size // 2) * 10
        self.v_weights = rng.random(n_size // 2) * 10

    def time_energy_distance(self, n_size):
        distance = stats.energy_distance(
                 self.u_values, self.v_values,
                 self.u_weights, self.v_weights)

    def time_wasserstein_distance(self, n_size):
        distance = stats.wasserstein_distance(
                 self.u_values, self.v_values,
                 self.u_weights, self.v_weights)


class Somersd(Benchmark):
    param_names = ['n_size']
    params = [
        [10, 100]
    ]

    def setup(self, n_size):
        rng = np.random.default_rng(12345678)
        self.x = rng.choice(n_size, size=n_size)
        self.y = rng.choice(n_size, size=n_size)

    def time_somersd(self, n_size):
        res = stats.somersd(self.x, self.y)
