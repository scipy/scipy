from functools import partial
import threading
import pickle
import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal, suppress_warnings
from numpy.lib import NumpyVersion
from scipy.stats import (
    TransformedDensityRejection,
    DiscreteAliasUrn,
    NumericalInversePolynomial
)
from scipy.stats import UNURANError
from scipy import stats
from scipy import special
from scipy.stats import chisquare, cramervonmises
from scipy.stats._distr_params import distdiscrete, distcont


# common test data: this data can be shared between all the tests.


# Normal distribution shared between all the continuous methods
class StandardNormal:
    def pdf(self, x):
        return np.exp(-0.5 * x*x)

    def dpdf(self, x):
        return -x * np.exp(-0.5 * x*x)

    def cdf(self, x):
        return special.ndtr(x)


# A binomial distribution to share between all the discrete methods
class Binomial:
    def __init__(self, n, p):
        self.n = n
        self.p = p

    def pmf(self, k):
        return self.p**k * (1-self.p)**(self.n-k)

    def cdf(self, k):
        k = np.asarray(k)
        return stats.binom._cdf(k, self.n, self.p)

    def support(self):
        return 0, self.n


all_methods = [
    ("TransformedDensityRejection", {"dist": StandardNormal()}),
    ("DiscreteAliasUrn", {"dist": [0.02, 0.18, 0.8]}),
    ("NumericalInversePolynomial", {"dist": StandardNormal()})
]


# Make sure an internal error occurs in UNU.RAN when invalid callbacks are
# passed. Moreover, different generators throw different error messages.
# So, in case of an `UNURANError`, we do not validate the error message.
bad_pdfs_common = [
    # Negative PDF
    (lambda x: -x, UNURANError, r"..."),
    # Returning wrong type
    (lambda x: [], TypeError, r"must be real number, not list"),
    # Undefined name inside the function
    (lambda x: foo, NameError, r"name 'foo' is not defined"),  # type: ignore[name-defined]  # noqa
    # Infinite value returned => Overflow error.
    (lambda x: np.inf, UNURANError, r"..."),
    # NaN value => internal error in UNU.RAN
    (lambda x: np.nan, UNURANError, r"..."),
    # signature of PDF wrong
    (lambda: 1.0, TypeError, r"takes 0 positional arguments but 1 was given")
]

# Make sure an internal error occurs in UNU.RAN when invalid callbacks are
# passed. Moreover, different generators throw different error messages.
# So, in case of an `UNURANError`, we do not validate the messages.
bad_dpdf_common = [
    # Infinite value returned.
    (lambda x: np.inf, UNURANError, r"..."),
    # NaN value => internal error in UNU.RAN
    (lambda x: np.nan, UNURANError, r"..."),
    # Returning wrong type
    (lambda x: [], TypeError, r"must be real number, not list"),
    # Undefined name inside the function
    (lambda x: foo, NameError, r"name 'foo' is not defined"),  # type: ignore[name-defined]  # noqa
    # signature of dPDF wrong
    (lambda: 1.0, TypeError, r"takes 0 positional arguments but 1 was given")
]


bad_pv_common = [
    ([], r"must contain at least one element"),
    ([[1.0, 0.0]], r"wrong number of dimensions \(expected 1, got 2\)"),
    ([0.2, 0.4, np.nan, 0.8], r"must contain only finite / non-nan values"),
    ([0.2, 0.4, np.inf, 0.8], r"must contain only finite / non-nan values"),
    ([0.0, 0.0], r"must contain at least one non-zero value"),
]


# size of the domains is incorrect
bad_sized_domains = [
    # > 2 elements in the domain
    ((1, 2, 3), ValueError, r"must be a length 2 tuple"),
    # empty domain
    ((), ValueError, r"must be a length 2 tuple")
]

# domain values are incorrect
bad_domains = [
    ((2, 1), UNURANError, r"left >= right"),
    ((1, 1), UNURANError, r"left >= right"),
]

# infinite and nan values present in domain.
inf_nan_domains = [
    # left >= right
    ((10, 10), UNURANError, r"left >= right"),
    ((np.inf, np.inf), UNURANError, r"left >= right"),
    ((-np.inf, -np.inf), UNURANError, r"left >= right"),
    ((np.inf, -np.inf), UNURANError, r"left >= right"),
    # Also include nans in some of the domains.
    ((-np.inf, np.nan), ValueError, r"only non-nan values"),
    ((np.nan, np.inf), ValueError, r"only non-nan values")
]

# `nan` values present in domain. Some distributions don't support
# infinite tails, so don't mix the nan values with infinities.
nan_domains = [
    ((0, np.nan), ValueError, r"only non-nan values"),
    ((np.nan, np.nan), ValueError, r"only non-nan values")
]


# all the methods should throw errors for nan, bad sized, and bad valued
# domains.
@pytest.mark.parametrize("domain, err, msg",
                         bad_domains + bad_sized_domains +
                         nan_domains)  # type: ignore[operator]
@pytest.mark.parametrize("method, kwargs", all_methods)
def test_bad_domain(domain, err, msg, method, kwargs):
    Method = getattr(stats, method)
    with pytest.raises(err, match=msg):
        Method(**kwargs, domain=domain)


@pytest.mark.parametrize("method, kwargs", all_methods)
def test_random_state(method, kwargs):
    Method = getattr(stats, method)

    # simple seed that works for any version of NumPy
    seed = 123
    rng1 = Method(**kwargs, random_state=seed)
    rng2 = Method(**kwargs, random_state=seed)
    assert_equal(rng1.rvs(100), rng2.rvs(100))

    # global seed
    np.random.seed(123)
    rng1 = Method(**kwargs)
    rvs1 = rng1.rvs(100)
    np.random.seed(None)
    rng2 = Method(**kwargs, random_state=123)
    rvs2 = rng2.rvs(100)
    assert_equal(rvs1, rvs2)

    # RandomState seed for old numpy
    if NumpyVersion(np.__version__) < '1.19.0':
        seed1 = np.random.RandomState(123)
        seed2 = 123
        rng1 = Method(**kwargs, random_state=seed1)
        rng2 = Method(**kwargs, random_state=seed2)
        assert_equal(rng1.rvs(100), rng2.rvs(100))
        rvs11 = rng1.rvs(550)
        rvs12 = rng1.rvs(50)
        rvs2 = rng2.rvs(600)
        assert_equal(rvs11, rvs2[:550])
        assert_equal(rvs12, rvs2[550:])
    else:  # Generator seed for new NumPy
        # when a RandomState is given, it should take the bitgen_t
        # member of the class and create a Generator instance.
        seed1 = np.random.RandomState(np.random.MT19937(123))
        seed2 = np.random.Generator(np.random.MT19937(123))
        rng1 = Method(**kwargs, random_state=seed1)
        rng2 = Method(**kwargs, random_state=seed2)
        assert_equal(rng1.rvs(100), rng2.rvs(100))


def test_set_random_state():
    rng1 = TransformedDensityRejection(StandardNormal(), random_state=123)
    rng2 = TransformedDensityRejection(StandardNormal())
    rng2.set_random_state(123)
    assert_equal(rng1.rvs(100), rng2.rvs(100))
    rng = TransformedDensityRejection(StandardNormal(), random_state=123)
    rvs1 = rng.rvs(100)
    rng.set_random_state(123)
    rvs2 = rng.rvs(100)
    assert_equal(rvs1, rvs2)


def test_threading_behaviour():
    # Test if the API is thread-safe.
    # This verifies if the lock mechanism and the use of `PyErr_Occurred`
    # is correct.
    errors = {"err1": None, "err2": None}

    class Distribution:
        def __init__(self, pdf_msg):
            self.pdf_msg = pdf_msg

        def pdf(self, x):
            if 49.9 < x < 50.0:
                raise ValueError(self.pdf_msg)
            return x

        def dpdf(self, x):
            return 1

    def func1():
        dist = Distribution('foo')
        rng = TransformedDensityRejection(dist, domain=(10, 100),
                                          random_state=12)
        try:
            rng.rvs(100000)
        except ValueError as e:
            errors['err1'] = e.args[0]

    def func2():
        dist = Distribution('bar')
        rng = TransformedDensityRejection(dist, domain=(10, 100),
                                          random_state=2)
        try:
            rng.rvs(100000)
        except ValueError as e:
            errors['err2'] = e.args[0]

    t1 = threading.Thread(target=func1)
    t2 = threading.Thread(target=func2)

    t1.start()
    t2.start()

    t1.join()
    t2.join()

    assert errors['err1'] == 'foo'
    assert errors['err2'] == 'bar'


@pytest.mark.parametrize("method, kwargs", all_methods)
def test_pickle(method, kwargs):
    Method = getattr(stats, method)
    rng1 = Method(**kwargs, random_state=123)
    obj = pickle.dumps(rng1)
    rng2 = pickle.loads(obj)
    assert_equal(rng1.rvs(100), rng2.rvs(100))


@pytest.mark.parametrize("size", [None, 0, (0, ), 1, (10, 3), (2, 3, 4, 5),
                                  (0, 0), (0, 1)])
def test_rvs_size(size):
    # As the `rvs` method is present in the base class and shared between
    # all the classes, we can just test with one of the methods.
    rng = TransformedDensityRejection(StandardNormal())
    if size is None:
        assert np.isscalar(rng.rvs(size))
    else:
        if np.isscalar(size):
            size = (size, )
        assert rng.rvs(size).shape == size


def check_cont_samples(rng, dist, mv_ex):
    rvs = rng.rvs(100000)
    mv = rvs.mean(), rvs.var()
    # test the moments only if the variance is finite
    if np.isfinite(mv_ex[1]):
        assert_allclose(mv, mv_ex, rtol=1e-7, atol=1e-1)
    # Cramer Von Mises test for goodness-of-fit
    rvs = rng.rvs(500)
    dist.cdf = np.vectorize(dist.cdf)
    pval = cramervonmises(rvs, dist.cdf).pvalue
    assert pval > 0.1


def check_discr_samples(rng, pv, mv_ex):
    rvs = rng.rvs(100000)
    # test if the first few moments match
    mv = rvs.mean(), rvs.var()
    assert_allclose(mv, mv_ex, rtol=1e-3, atol=1e-1)
    # normalize
    pv = pv / pv.sum()
    # chi-squared test for goodness-of-fit
    obs_freqs = np.zeros_like(pv)
    _, freqs = np.unique(rvs, return_counts=True)
    freqs = freqs / freqs.sum()
    obs_freqs[:freqs.size] = freqs
    pval = chisquare(obs_freqs, pv).pvalue
    assert pval > 0.1


class TestTransformedDensityRejection:
    # Simple Custom Distribution
    class dist0:
        def pdf(self, x):
            return 3/4 * (1-x*x)

        def dpdf(self, x):
            return 3/4 * (-2*x)

        def cdf(self, x):
            return 3/4 * (x - x**3/3 + 2/3)

        def support(self):
            return -1, 1

    # Standard Normal Distribution
    class dist1:
        def pdf(self, x):
            return stats.norm._pdf(x / 0.1)

        def dpdf(self, x):
            return -x / 0.01 * stats.norm._pdf(x / 0.1)

        def cdf(self, x):
            return stats.norm._cdf(x / 0.1)

    # pdf with piecewise linear function as transformed density
    # with T = -1/sqrt with shift. Taken from UNU.RAN test suite
    # (from file t_tdr_ps.c)
    class dist2:
        def __init__(self, shift):
            self.shift = shift

        def pdf(self, x):
            x -= self.shift
            y = 1. / (abs(x) + 1.)
            return 0.5 * y * y

        def dpdf(self, x):
            x -= self.shift
            y = 1. / (abs(x) + 1.)
            y = y * y * y
            return y if (x < 0.) else -y

        def cdf(self, x):
            x -= self.shift
            if x <= 0.:
                return 0.5 / (1. - x)
            else:
                return 1. - 0.5 / (1. + x)

    dists = [dist0(), dist1(), dist2(0.), dist2(10000.)]

    # exact mean and variance of the distributions in the list dists
    mv0 = [0., 4./15.]
    mv1 = [0., 0.01]
    mv2 = [0., np.inf]
    mv3 = [10000., np.inf]
    mvs = [mv0, mv1, mv2, mv3]

    @pytest.mark.parametrize("dist, mv_ex",
                             zip(dists, mvs))
    def test_basic(self, dist, mv_ex):
        with suppress_warnings() as sup:
            # filter the warnings thrown by UNU.RAN
            sup.filter(RuntimeWarning)
            rng = TransformedDensityRejection(dist, random_state=42)
        check_cont_samples(rng, dist, mv_ex)

    # PDF 0 everywhere => bad construction points
    bad_pdfs = [(lambda x: 0, UNURANError, r"50 : bad construction points.")]
    bad_pdfs += bad_pdfs_common  # type: ignore[arg-type]

    @pytest.mark.parametrize("pdf, err, msg", bad_pdfs)
    def test_bad_pdf(self, pdf, err, msg):
        class dist:
            pass
        dist.pdf = pdf
        dist.dpdf = lambda x: 1  # an arbitrary dPDF
        with pytest.raises(err, match=msg):
            TransformedDensityRejection(dist)

    @pytest.mark.parametrize("dpdf, err, msg", bad_dpdf_common)
    def test_bad_dpdf(self, dpdf, err, msg):
        class dist:
            pass
        dist.pdf = lambda x: x
        dist.dpdf = dpdf
        with pytest.raises(err, match=msg):
            TransformedDensityRejection(dist, domain=(1, 10))

    # test domains with inf + nan in them. need to write a custom test for
    # this because not all methods support infinite tails.
    @pytest.mark.parametrize("domain, err, msg", inf_nan_domains)
    def test_inf_nan_domains(self, domain, err, msg):
        with pytest.raises(err, match=msg):
            TransformedDensityRejection(StandardNormal(), domain=domain)

    @pytest.mark.parametrize("construction_points", [-1, 0, 0.1])
    def test_bad_construction_points_scalar(self, construction_points):
        with pytest.raises(ValueError, match=r"`construction_points` must be "
                                             r"a positive integer."):
            TransformedDensityRejection(
                StandardNormal(), construction_points=construction_points
            )

    def test_bad_construction_points_array(self):
        # empty array
        construction_points = []
        with pytest.raises(ValueError, match=r"`construction_points` must "
                                             r"either be a "
                                             r"scalar or a non-empty array."):
            TransformedDensityRejection(
                StandardNormal(), construction_points=construction_points
            )

        # construction_points not monotonically increasing
        construction_points = [1, 1, 1, 1, 1, 1]
        with pytest.warns(RuntimeWarning, match=r"33 : starting points not "
                                                r"strictly monotonically "
                                                r"increasing"):
            TransformedDensityRejection(
                StandardNormal(), construction_points=construction_points
            )

        # construction_points containing nans
        construction_points = [np.nan, np.nan, np.nan]
        with pytest.raises(UNURANError, match=r"50 : bad construction "
                                              r"points."):
            TransformedDensityRejection(
                StandardNormal(), construction_points=construction_points
            )

        # construction_points out of domain
        construction_points = [-10, 10]
        with pytest.warns(RuntimeWarning, match=r"50 : starting point out of "
                                                r"domain"):
            TransformedDensityRejection(
                StandardNormal(), domain=(-3, 3),
                construction_points=construction_points
            )

    @pytest.mark.parametrize("c", [-1., np.nan, np.inf, 0.1, 1.])
    def test_bad_c(self, c):
        msg = r"`c` must either be -0.5 or 0."
        with pytest.raises(ValueError, match=msg):
            TransformedDensityRejection(StandardNormal(), c=-1.)

    def test_bad_variant(self):
        msg = r"Invalid option for the `variant`"
        with pytest.raises(ValueError, match=msg):
            TransformedDensityRejection(StandardNormal(), variant='foo')

    u = [np.linspace(0, 1, num=1000), [], [[]], [np.nan],
         [-np.inf, np.nan, np.inf], 0,
         [[np.nan, 0.5, 0.1], [0.2, 0.4, np.inf], [-2, 3, 4]]]

    @pytest.mark.parametrize("u", u)
    def test_ppf_hat(self, u):
        # Increase the `max_squeeze_hat_ratio` so the ppf_hat is more
        # accurate.
        rng = TransformedDensityRejection(StandardNormal(),
                                          max_squeeze_hat_ratio=0.9999,
                                          max_intervals=10000)
        # Older versions of NumPy throw RuntimeWarnings for comparisons
        # with nan.
        with suppress_warnings() as sup:
            sup.filter(RuntimeWarning, "invalid value encountered in greater")
            sup.filter(RuntimeWarning, "invalid value encountered in "
                                       "greater_equal")
            sup.filter(RuntimeWarning, "invalid value encountered in less")
            sup.filter(RuntimeWarning, "invalid value encountered in "
                                       "less_equal")
            res = rng.ppf_hat(u)
            expected = stats.norm.ppf(u)
        assert_allclose(res, expected, rtol=1e-3, atol=1e-5)
        assert res.shape == expected.shape

    def test_bad_dist(self):
        # Empty distribution
        class dist:
            ...

        msg = r"`pdf` required but not found."
        with pytest.raises(ValueError, match=msg):
            TransformedDensityRejection(dist)

        # dPDF not present in dist
        class dist:
            pdf = lambda x: 1-x*x  # noqa: E731

        msg = r"`dpdf` required but not found."
        with pytest.raises(ValueError, match=msg):
            TransformedDensityRejection(dist)


class TestDiscreteAliasUrn:
    # DAU fails on these probably because of large domains and small
    # computation errors in PMF. Mean/SD match but chi-squared test fails.
    basic_fail_dists = {
        'nchypergeom_fisher',  # numerical erros on tails
        'nchypergeom_wallenius',  # numerical erros on tails
        'randint'  # fails on 32-bit ubuntu
    }

    @pytest.mark.parametrize("distname, params", distdiscrete)
    def test_basic(self, distname, params):
        if distname in self.basic_fail_dists:
            msg = ("DAU fails on these probably because of large domains "
                   "and small computation errors in PMF.")
            pytest.skip(msg)
        if not isinstance(distname, str):
            dist = distname
        else:
            dist = getattr(stats, distname)
        dist = dist(*params)
        domain = dist.support()
        if not np.isfinite(domain[1] - domain[0]):
            # DAU only works with finite domain. So, skip the distributions
            # with infinite tails.
            pytest.skip("DAU only works with a finite domain.")
        k = np.arange(domain[0], domain[1]+1)
        pv = dist.pmf(k)
        mv_ex = dist.stats('mv')
        rng = DiscreteAliasUrn(dist, random_state=42)
        check_discr_samples(rng, pv, mv_ex)

    # Can't use bad_pmf_common here as we evaluate PMF early on to avoid
    # unhelpful errors from UNU.RAN.
    bad_pmf = [
        # inf returned
        (lambda x: np.inf, ValueError,
         r"must contain only finite / non-nan values"),
        # nan returned
        (lambda x: np.nan, ValueError,
         r"must contain only finite / non-nan values"),
        # all zeros
        (lambda x: 0.0, ValueError,
         r"must contain at least one non-zero value"),
        # Undefined name inside the function
        (lambda x: foo, NameError,  # type: ignore[name-defined]  # noqa
         r"name 'foo' is not defined"),
        # Returning wrong type.
        (lambda x: [], ValueError,
         r"setting an array element with a sequence."),
        # probabilities < 0
        (lambda x: -x, UNURANError,
         r"50 : probability < 0"),
        # signature of PMF wrong
        (lambda: 1.0, TypeError,
         r"takes 0 positional arguments but 1 was given")
    ]

    @pytest.mark.parametrize("pmf, err, msg", bad_pmf)
    def test_bad_pmf(self, pmf, err, msg):
        class dist:
            pass
        dist.pmf = pmf
        with pytest.raises(err, match=msg):
            DiscreteAliasUrn(dist, domain=(1, 10))

    @pytest.mark.parametrize("pv", [[0.18, 0.02, 0.8],
                                    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]])
    def test_sampling_with_pv(self, pv):
        pv = np.asarray(pv, dtype=np.float64)
        rng = DiscreteAliasUrn(pv, random_state=123)
        rvs = rng.rvs(100_000)
        pv = pv / pv.sum()
        variates = np.arange(0, len(pv))
        # test if the first few moments match
        m_expected = np.average(variates, weights=pv)
        v_expected = np.average((variates - m_expected) ** 2, weights=pv)
        mv_expected = m_expected, v_expected
        check_discr_samples(rng, pv, mv_expected)

    @pytest.mark.parametrize("pv, msg", bad_pv_common)
    def test_bad_pv(self, pv, msg):
        with pytest.raises(ValueError, match=msg):
            DiscreteAliasUrn(pv)

    # DAU doesn't support infinite tails. So, it should throw an error when
    # inf is present in the domain.
    inf_domain = [(-np.inf, np.inf), (np.inf, np.inf), (-np.inf, -np.inf),
                  (0, np.inf), (-np.inf, 0)]

    @pytest.mark.parametrize("domain", inf_domain)
    def test_inf_domain(self, domain):
        with pytest.raises(ValueError, match=r"must be finite"):
            DiscreteAliasUrn(Binomial(10, 0.2), domain=domain)

    def test_bad_urn_factor(self):
        with pytest.warns(RuntimeWarning, match=r"relative urn size < 1."):
            DiscreteAliasUrn([0.5, 0.5], urn_factor=-1)

    def test_bad_args(self):
        msg = (r"`domain` must be provided when the "
               r"probability vector is not available.")

        class dist:
            def pmf(self, x):
                return x

        with pytest.raises(ValueError, match=msg):
            DiscreteAliasUrn(dist)


class TestNumericalInversePolynomial:
    # Simple Custom Distribution
    class dist0:
        def pdf(self, x):
            return 3/4 * (1-x*x)

        def cdf(self, x):
            return 3/4 * (x - x**3/3 + 2/3)

        def support(self):
            return -1, 1

    # Standard Normal Distribution
    class dist1:
        def pdf(self, x):
            return stats.norm._pdf(x / 0.1)

        def cdf(self, x):
            return stats.norm._cdf(x / 0.1)

    # Sin 2 distribution
    #          /  0.05 + 0.45*(1 +sin(2 Pi x))  if |x| <= 1
    #  f(x) = <
    #          \  0        otherwise
    # Taken from UNU.RAN test suite (from file t_pinv.c)
    class dist2:
        def pdf(self, x):
            return 0.05 + 0.45 * (1 + np.sin(2*np.pi*x))

        def cdf(self, x):
            return (0.05*(x + 1) +
                    0.9*(1. + 2.*np.pi*(1 + x) - np.cos(2.*np.pi*x)) /
                    (4.*np.pi))

        def support(self):
            return -1, 1

    # Sin 10 distribution
    #          /  0.05 + 0.45*(1 +sin(2 Pi x))  if |x| <= 5
    #  f(x) = <
    #          \  0        otherwise
    # Taken from UNU.RAN test suite (from file t_pinv.c)
    class dist3:
        def pdf(self, x):
            return 0.2 * (0.05 + 0.45 * (1 + np.sin(2*np.pi*x)))

        def cdf(self, x):
            return x/10. + 0.5 + 0.09/(2*np.pi) * (np.cos(10*np.pi) -
                                                   np.cos(2*np.pi*x))

        def support(self):
            return -5, 5

    dists = [dist0(), dist1(), dist2(), dist3()]

    # exact mean and variance of the distributions in the list dists
    mv0 = [0., 4./15.]
    mv1 = [0., 0.01]
    mv2 = [-0.45/np.pi, 2/3*0.5 - 0.45**2/np.pi**2]
    mv3 = [-0.45/np.pi, 0.2 * 250/3 * 0.5 - 0.45**2/np.pi**2]
    mvs = [mv0, mv1, mv2, mv3]

    @pytest.mark.parametrize("dist, mv_ex",
                             zip(dists, mvs))
    def test_basic(self, dist, mv_ex):
        rng = NumericalInversePolynomial(dist, random_state=42)
        check_cont_samples(rng, dist, mv_ex)

    very_slow_dists = ['studentized_range', 'trapezoid', 'triang', 'vonmises',
                       'levy_stable', 'kappa4', 'ksone', 'kstwo', 'levy_l',
                       'gausshyper', 'anglit']
    # for some reason, UNU.RAN segmentation faults for the uniform.
    fatal_fail_dists = ['uniform']
    # fails for unbounded PDFs
    unbounded_pdf_fail_dists = ['beta']
    # for these distributions, some assertions fail due to minor
    # numerical differences. They can be avoided either by changing
    # the seed or by increasing the u_resolution.
    fail_dists = ['ncf', 'pareto', 'chi2', 'fatiguelife', 'halfgennorm',
                  'gilbrat', 'lognorm', 'ncx2', 't']

    @pytest.mark.xslow
    @pytest.mark.parametrize("distname, params", distcont)
    def test_basic_all_scipy_dists(self, distname, params):
        if distname in self.very_slow_dists:
            pytest.skip(f"PINV too slow for {distname}")
        if distname in self.fail_dists:
            pytest.skip(f"PINV fails for {distname}")
        if distname in self.unbounded_pdf_fail_dists:
            pytest.skip("PINV fails for unbounded PDFs.")
        if distname in self.fatal_fail_dists:
            pytest.xfail(f"PINV segmentation faults for {distname}")
        dist = (getattr(stats, distname)
                if isinstance(distname, str)
                else distname)
        dist = dist(*params)
        with suppress_warnings() as sup:
            sup.filter(RuntimeWarning)
            rng = NumericalInversePolynomial(dist, random_state=42)
        check_cont_samples(rng, dist, [dist.mean(), dist.var()])

    @pytest.mark.parametrize("pdf, err, msg", bad_pdfs_common)
    def test_bad_pdf(self, pdf, err, msg):
        class dist:
            pass
        dist.pdf = pdf
        with pytest.raises(err, match=msg):
            NumericalInversePolynomial(dist, domain=[0, 5])

    # test domains with inf + nan in them. need to write a custom test for
    # this because not all methods support infinite tails.
    @pytest.mark.parametrize("domain, err, msg", inf_nan_domains)
    def test_inf_nan_domains(self, domain, err, msg):
        with pytest.raises(err, match=msg):
            NumericalInversePolynomial(StandardNormal(), domain=domain)

    u = [
        # test if quantile 0 and 1 return -inf and inf respectively and check
        # the correctness of the PPF for equidistant points between 0 and 1.
        np.linspace(0, 1, num=10000),
        # test the PPF method for empty arrays
        [], [[]],
        # test if nans and infs return nan result.
        [np.nan], [-np.inf, np.nan, np.inf],
        # test if a scalar is returned for a scalar input.
        0,
        # test for arrays with nans, values greater than 1 and less than 0,
        # and some valid values.
        [[np.nan, 0.5, 0.1], [0.2, 0.4, np.inf], [-2, 3, 4]]
    ]

    @pytest.mark.parametrize("u", u)
    def test_ppf(self, u):
        dist = StandardNormal()
        rng = NumericalInversePolynomial(dist, u_resolution=1e-14)
        # Older versions of NumPy throw RuntimeWarnings for comparisons
        # with nan.
        with suppress_warnings() as sup:
            sup.filter(RuntimeWarning, "invalid value encountered in greater")
            sup.filter(RuntimeWarning, "invalid value encountered in "
                                       "greater_equal")
            sup.filter(RuntimeWarning, "invalid value encountered in less")
            sup.filter(RuntimeWarning, "invalid value encountered in "
                                       "less_equal")
            res = rng.ppf(u)
            expected = stats.norm.ppf(u)
        assert_allclose(res, expected, rtol=1e-11, atol=1e-11)
        assert res.shape == expected.shape

    x = [np.linspace(-10, 10, num=10000), [], [[]], [np.nan],
         [-np.inf, np.nan, np.inf], 0,
         [[np.nan, 0.5, 0.1], [0.2, 0.4, np.inf], [-np.inf, 3, 4]]]

    @pytest.mark.parametrize("x", x)
    def test_cdf(self, x):
        dist = StandardNormal()
        rng = NumericalInversePolynomial(dist, keep_cdf=True,
                                         u_resolution=1e-14)
        # Older versions of NumPy throw RuntimeWarnings for comparisons
        # with nan.
        with suppress_warnings() as sup:
            sup.filter(RuntimeWarning, "invalid value encountered in greater")
            sup.filter(RuntimeWarning, "invalid value encountered in "
                                       "greater_equal")
            sup.filter(RuntimeWarning, "invalid value encountered in less")
            sup.filter(RuntimeWarning, "invalid value encountered in "
                                       "less_equal")
            res = rng.cdf(x)
            expected = stats.norm.cdf(x)
        assert_allclose(res, expected, rtol=1e-11, atol=1e-11)
        assert res.shape == expected.shape

    def test_u_error(self):
        dist = StandardNormal()
        rng = NumericalInversePolynomial(dist, u_resolution=1e-10)
        max_error, mae = rng.u_error()
        assert max_error < 1e-10
        assert mae <= max_error
        rng = NumericalInversePolynomial(dist, u_resolution=1e-14)
        max_error, mae = rng.u_error()
        assert max_error < 1e-14
        assert mae <= max_error

    bad_orders = [1, 4.5, 20, np.inf, np.nan]
    bad_u_resolution = [1e-20, 1e-1, np.inf, np.nan]
    bad_max_intervals = [10, 10000000, 1000.5, np.inf, np.nan]

    @pytest.mark.parametrize("order", bad_orders)
    def test_bad_orders(self, order):
        dist = StandardNormal()

        msg = r"`order` must be an integer in the range \[3, 17\]."
        with pytest.raises(ValueError, match=msg):
            NumericalInversePolynomial(dist, order=order)

    @pytest.mark.parametrize("u_resolution", bad_u_resolution)
    def test_bad_u_resolution(self, u_resolution):
        msg = r"`u_resolution` must be between 1e-15 and 1e-5."
        with pytest.raises(ValueError, match=msg):
            NumericalInversePolynomial(StandardNormal(),
                                       u_resolution=u_resolution)

    @pytest.mark.parametrize("max_intervals", bad_max_intervals)
    def test_bad_max_intervals(self, max_intervals):
        msg = (r"`max_intervals` must be an integer in the range "
               r"\[100, 1000000\].")
        with pytest.raises(ValueError, match=msg):
            NumericalInversePolynomial(StandardNormal(),
                                       max_intervals=max_intervals)

    def test_bad_args(self):
        dist = StandardNormal()
        rng = NumericalInversePolynomial(dist)
        msg = r"CDF is not available."
        with pytest.raises(ValueError, match=msg):
            rng.cdf([1, 2, 3])
        msg = r"`sample_size` must be greater than or equal to 1000."
        with pytest.raises(ValueError, match=msg):
            rng.u_error(10)

        class Distribution:
            def pdf(self, x):
                return np.exp(-0.5 * x*x)

        dist = Distribution()
        rng = NumericalInversePolynomial(dist)
        msg = r"Exact CDF required but not found."
        with pytest.raises(ValueError, match=msg):
            rng.u_error()
