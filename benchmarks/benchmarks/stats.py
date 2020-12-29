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


class InferentialStats(Benchmark):
    def setup(self):
        np.random.seed(12345678)
        self.a = stats.norm.rvs(loc=5, scale=10, size=500)
        self.b = stats.norm.rvs(loc=8, scale=10, size=20)
        self.c = stats.norm.rvs(loc=8, scale=20, size=20)

    def time_ttest_ind_same_var(self):
        # test different sized sample with variances
        stats.ttest_ind(self.a, self.b)
        stats.ttest_ind(self.a, self.b, equal_var=False)

    def time_ttest_ind_diff_var(self):
        # test different sized sample with different variances
        stats.ttest_ind(self.a, self.c)
        stats.ttest_ind(self.a, self.c, equal_var=False)


class Distribution(Benchmark):
    # add distributions here
    dists = ([d[0] for d in distcont + distdiscrete] +
             ['frechet_l', 'frechet_r', 'trapz', 'laplace_asymmetric',
              'nhypergeom', 'reciprocal', 'trapezoid'])

    param_names = ['distribution', 'properties']
    params = [
        dists, ['pdf/pmf', 'logpdf/logpmf', 'cdf', 'logcdf', 'rvs', 'fit',
                'sf', 'logsf', 'ppf', 'isf', 'moment', 'stats_s', 'stats_v',
                'stats_m', 'stats_k', 'stats_mvsk', 'entropy', 'median',
                'mean', 'var', 'std', 'interval', 'expect']
    ]
    distcont = dict(distcont + distdiscrete)
    # maintain previous benchmarks' values
    custom_input = {'gamma': [5], 'beta': [5, 3]}

    def setup(self, distribution, properties):
        if not is_xslow():
            raise NotImplementedError("Skipped")

        self.dist = getattr(stats, distribution)

        shapes = self.distcont[distribution]

        if isinstance(self.dist, stats.rv_discrete):
            self.isCont = False
            kwds = {'loc': 4}
        else:
            self.isCont = True
            kwds = {'loc': 4, 'scale': 10}

        rng = self.dist.interval(.99, *shapes, **kwds)
        x = np.linspace(*rng, 100)
        args = [x, *self.custom_input.get(distribution, shapes)]

        if properties == 'fit':
            # provide only the data to fit in args
            self.args = args[:1]
        elif properties == 'rvs':
            # add size for creation of data
            kwds['size'] = 1000
            self.args = args[1:]
        elif properties == 'pdf/pmf':
            # picking arbitrary value for this
            self.args = [.99, *args[1:]]
            properties = 'pmf' if isinstance(self.dist, stats.rv_discrete) else 'pdf' 
        elif properties == 'logpdf/logpmf':
            properties = 'logpmf' if isinstance(self.dist, stats.rv_discrete) else 'logpdf'
            self.args = args
        elif properties == 'isf':
            self.args = [.99, *args[1:]]
        elif properties == 'moment':
            self.args = [10, *args[1:]]
        elif properties.startswith('stats_'):
            properties = 'stats'
            self.args = args[1:]
            kwds['moments'] = properties[6:]
            self.kwds = kwds
        elif properties in ['entropy', 'median', 'mean', 'var', 'std']:
            self.args = args[1:]
        elif properties == 'expect':
            self.args = []
            kwds['args'] = tuple(args[1:])
        elif properties == 'interval':
            self.args = [.99, *args[1:]]
        else:
            self.args = args
        self.kwds = kwds

        try:
            self.method = getattr(self.dist, properties)
        except AttributeError:
            raise NotImplementedError("This attribute is not a member "
                                      "of the distribution")

    def time_distribution(self, distribution, properties):
        self.method(*self.args, **self.kwds)

    # Retain old benchmark results (remove this if changing the benchmark)
    time_distribution.version = "fb22ae5386501008d945783921fe44aef3f82c1dafc40cddfaccaeec38b792b0"


class DescriptiveStats(Benchmark):
    param_names = ['n_levels']
    params = [
        [10, 1000]
    ]

    def setup(self, n_levels):
        np.random.seed(12345678)
        self.levels = np.random.randint(n_levels, size=(1000, 10))

    def time_mode(self, n_levels):
        stats.mode(self.levels, axis=0)


class GaussianKDE(Benchmark):
    def setup(self):
        np.random.seed(12345678)
        n = 2000
        m1 = np.random.normal(size=n)
        m2 = np.random.normal(scale=0.5, size=n)

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
        np.random.seed(12345678)

    def time_unitary_group(self, dim):
        stats.unitary_group.rvs(dim)

    def time_ortho_group(self, dim):
        stats.ortho_group.rvs(dim)

    def time_special_ortho_group(self, dim):
        stats.special_ortho_group.rvs(dim)


class BinnedStatisticDD(Benchmark):

    params = ["count", "sum", "mean", "min", "max", "median", "std", np.std]

    def setup(self, statistic):
        np.random.seed(12345678)
        self.inp = np.random.rand(9999).reshape(3, 3333) * 200
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
