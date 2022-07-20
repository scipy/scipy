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
        stats.kendalltau(self.a, self.b, nan_policy=nan_policy,
                         method=method, variant=variant)


class KS(Benchmark):
    param_names = ['alternative', 'mode']
    params = [
        ['two-sided', 'less', 'greater'],
        ['auto', 'exact', 'asymp'],
    ]

    def setup(self, alternative, mode):
        rng = np.random.default_rng(0x2e7c964ff9a5cd6be22014c09f1dbba9)
        self.a = stats.norm.rvs(loc=5, scale=10, size=500, random_state=rng)
        self.b = stats.norm.rvs(loc=8, scale=10, size=500, random_state=rng)

    def time_ks_1samp(self, alternative, mode):
        stats.ks_1samp(self.a, stats.norm.cdf,
                       alternative=alternative, mode=mode)

    def time_ks_2samp(self, alternative, mode):
        stats.ks_2samp(self.a, self.b, alternative=alternative, mode=mode)


class RankSums(Benchmark):
    param_names = ['alternative']
    params = [
        ['two-sided', 'less', 'greater']
    ]

    def setup(self, alternative):
        rng = np.random.default_rng(0xb6acd7192d6e5da0f68b5d8ab8ce7af2)
        self.u1 = rng.uniform(-1, 1, 200)
        self.u2 = rng.uniform(-0.5, 1.5, 300)

    def time_ranksums(self, alternative):
        stats.ranksums(self.u1, self.u2, alternative=alternative)


class BrunnerMunzel(Benchmark):
    param_names = ['alternative', 'nan_policy', 'distribution']
    params = [
        ['two-sided', 'less', 'greater'],
        ['propagate', 'raise', 'omit'],
        ['t', 'normal']
    ]

    def setup(self, alternative, nan_policy, distribution):
        rng = np.random.default_rng(0xb82c4db22b2818bdbc5dbe15ad7528fe)
        self.u1 = rng.uniform(-1, 1, 200)
        self.u2 = rng.uniform(-0.5, 1.5, 300)

    def time_brunnermunzel(self, alternative, nan_policy, distribution):
        stats.brunnermunzel(self.u1, self.u2, alternative=alternative,
                            distribution=distribution, nan_policy=nan_policy)


class InferentialStats(Benchmark):
    def setup(self):
        rng = np.random.default_rng(0x13d756fadb635ae7f5a8d39bbfb0c931)
        self.a = stats.norm.rvs(loc=5, scale=10, size=500, random_state=rng)
        self.b = stats.norm.rvs(loc=8, scale=10, size=500, random_state=rng)
        self.c = stats.norm.rvs(loc=8, scale=20, size=500, random_state=rng)
        self.chisq = rng.integers(1, 20, 500)

    def time_ttest_ind_same_var(self):
        # test different sized sample with variances
        stats.ttest_ind(self.a, self.b)
        stats.ttest_ind(self.a, self.b, equal_var=False)

    def time_ttest_ind_diff_var(self):
        # test different sized sample with different variances
        stats.ttest_ind(self.a, self.c)
        stats.ttest_ind(self.a, self.c, equal_var=False)

    def time_chisqure(self):
        stats.chisquare(self.chisq)

    def time_friedmanchisquare(self):
        stats.friedmanchisquare(self.a, self.b, self.c)

    def time_epps_singleton_2samp(self):
        stats.epps_singleton_2samp(self.a, self.b)

    def time_kruskal(self):
        stats.mstats.kruskal(self.a, self.b)


# Benchmark data for the truncnorm stats() method.
# The data in each row is:
#   a, b, mean, variance, skewness, excess kurtosis. Generated using
# https://gist.github.com/WarrenWeckesser/636b537ee889679227d53543d333a720
truncnorm_cases = [[-20, -19, -19.052343945976656, 0.002725073018195613,
                   -1.9838693623377885, 5.871801893091683],
                   [-30, -29, -29.034401237736176, 0.0011806604886186853,
                    -1.9929615171469608, 5.943905539773037],
                   [-40, -39, -39.02560741993011, 0.0006548827702932775,
                    -1.9960847672775606, 5.968744357649675],
                   [39, 40, 39.02560741993011, 0.0006548827702932775,
                    1.9960847672775606, 5.968744357649675]]
truncnorm_cases = np.array(truncnorm_cases)


class TruncnormStats(Benchmark):
    param_names = ['case', 'moment']
    params = [list(range(len(truncnorm_cases))), ['m', 'v', 's', 'k']]

    def track_truncnorm_stats_error(self, case, moment):
        result_indices = dict(zip(['m', 'v', 's', 'k'], range(2, 6)))
        ref = truncnorm_cases[case, result_indices[moment]]
        a, b = truncnorm_cases[case, 0:2]
        res = stats.truncnorm(a, b).stats(moments=moment)
        return np.abs((res - ref)/ref)


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


class TrackContinuousRoundtrip(Benchmark):
    # Benchmarks that track a value for every distribution can go here
    param_names = ['dist_name']
    params = list(dict(distcont).keys())
    dist_data = dict(distcont)

    def setup(self, dist_name):
        # Distribution setup follows `DistributionsAll` benchmark.
        # This focuses on ppf, so the code for handling other functions is
        # removed for simplicity.
        self.dist = getattr(stats, dist_name)
        self.shape_args = self.dist_data[dist_name]

    def track_distribution_ppf_roundtrip(self, dist_name):
        # Tracks the worst relative error of a
        # couple of round-trip ppf -> cdf calculations.
        vals = [0.001, 0.5, 0.999]

        ppf = self.dist.ppf(vals, *self.shape_args)
        round_trip = self.dist.cdf(ppf, *self.shape_args)

        err_rel = np.abs(vals - round_trip) / vals
        return np.max(err_rel)

    def track_distribution_ppf_roundtrip_extrema(self, dist_name):
        # Tracks the absolute error of an "extreme" round-trip
        # ppf -> cdf calculation.
        v = 1e-6
        ppf = self.dist.ppf(v, *self.shape_args)
        round_trip = self.dist.cdf(ppf, *self.shape_args)

        err_abs = np.abs(v - round_trip)
        return err_abs

    def track_distribution_isf_roundtrip(self, dist_name):
        # Tracks the worst relative error of a
        # couple of round-trip isf -> sf calculations.
        vals = [0.001, 0.5, 0.999]

        isf = self.dist.isf(vals, *self.shape_args)
        round_trip = self.dist.sf(isf, *self.shape_args)

        err_rel = np.abs(vals - round_trip) / vals
        return np.max(err_rel)

    def track_distribution_isf_roundtrip_extrema(self, dist_name):
        # Tracks the absolute error of an "extreme" round-trip
        # isf -> sf calculation.
        v = 1e-6
        ppf = self.dist.isf(v, *self.shape_args)
        round_trip = self.dist.sf(ppf, *self.shape_args)

        err_abs = np.abs(v - round_trip)
        return err_abs


class PDFPeakMemory(Benchmark):
    # Tracks peak memory when a distribution is given a large array to process
    # See gh-14095

    # Run for up to 30 min - some dists are quite slow.
    timeout = 1800.0

    x = np.arange(1e6)

    param_names = ['dist_name']
    params = list(dict(distcont).keys())
    dist_data = dict(distcont)

    # So slow that 30min isn't enough time to finish.
    slow_dists = ["levy_stable"]

    def setup(self, dist_name):
        # This benchmark is demanding. Skip it if the env isn't xslow.
        if not is_xslow():
            raise NotImplementedError("skipped - enviroment is not xslow. "
                                      "To enable this benchamark, set the "
                                      "enviroment variable SCIPY_XSLOW=1")

        if dist_name in self.slow_dists:
            raise NotImplementedError("skipped - dist is too slow.")

        self.dist = getattr(stats, dist_name)
        self.shape_args = self.dist_data[dist_name]

    def peakmem_bigarr_pdf(self, dist_name):
        self.dist.pdf(self.x, *self.shape_args)


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


class BenchQMCHalton(Benchmark):
    param_names = ['d', 'scramble', 'n', 'workers']
    params = [
        [1, 10],
        [True, False],
        [10, 1_000, 100_000],
        [1, 4]
    ]

    def setup(self, d, scramble, n, workers):
        self.rng = np.random.default_rng(1234)

    def time_halton(self, d, scramble, n, workers):
        seq = stats.qmc.Halton(d, scramble=scramble, seed=self.rng)
        seq.random(n, workers=workers)


class BenchQMCSobol(Benchmark):
    param_names = ['d', 'base2']
    params = [
        [1, 50, 100],
        [3, 10, 11, 12],
    ]

    def setup(self, d, base2):
        self.rng = np.random.default_rng(168525179735951991038384544)
        stats.qmc.Sobol(1, bits=32)  # make it load direction numbers

    def time_sobol(self, d, base2):
        # scrambling is happening at init only, not worth checking
        seq = stats.qmc.Sobol(d, scramble=False, bits=32, seed=self.rng)
        seq.random_base2(base2)


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
