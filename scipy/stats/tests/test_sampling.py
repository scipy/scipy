from functools import partial
import threading
import pickle
import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal, suppress_warnings
from numpy.lib import NumpyVersion
from scipy.stats import TransformedDensityRejection, DiscreteAliasUrn
from scipy.stats import UNURANError
from scipy import stats
from scipy.stats import chisquare, cramervonmises
from scipy.stats._distr_params import distdiscrete


# common test data: this data can be shared between all the tests.


# Normal distribution shared between all the continuous methods
class common_cont_dist:
    # Private methods are used for performance. Public methods do a lot of
    # validations and expensive numpy operations (like broadcasting), hence,
    # slowing down some of the methods significantly. We can use private
    # methods as we know that the data is guaranteed to be inside the domain
    # and a scalar (so, no masking bad values and broadcasting required).
    pdf = stats.norm._pdf
    dpdf = lambda x: -x * stats.norm._pdf(x)  # noqa: E731
    cdf = stats.norm._cdf


# A binomial distribution to share between all the discrete methods
class common_discr_dist:
    params = (10, 0.2)
    pmf = partial(stats.binom._pmf, *params)
    cdf = partial(stats.binom._cdf, *params)


all_methods = [
    ("TransformedDensityRejection", {"dist": common_cont_dist}),
    ("DiscreteAliasUrn", {"dist": [0.02, 0.18, 0.8]})
]


bad_pdfs_common = [
    # Negative PDF
    (lambda x: -x, UNURANError, r"50 : PDF\(x\) < 0.!"),
    # Returning wrong type
    (lambda x: [], TypeError, r"must be real number, not list"),
    # Undefined name inside the function
    (lambda x: foo, NameError, r"name 'foo' is not defined"),  # type: ignore[name-defined]  # noqa
    # Infinite value returned => Overflow error.
    (lambda x: np.inf, UNURANError, r"50 : PDF\(x\) overflow"),
    # NaN value => internal error in UNU.RAN
    # Currently, UNU.RAN just returns a "cannot create bounded hat!"
    # error which is not very useful. So, instead of testing the error
    # message, we just ensure an error is raised.
    (lambda x: np.nan, UNURANError, r"..."),
    # signature of PDF wrong
    (lambda: 1.0, TypeError, r"takes 0 positional arguments but 1 was given")
]

bad_dpdf_common = [
    # Infinite value returned. For dPDF, UNU.RAN complains that it cannot
    # create a bounded hat instead of an overflow error.
    (lambda x: np.inf, UNURANError, r"51 : cannot create bounded hat"),
    # NaN value => internal error in UNU.RAN
    # Currently, UNU.RAN just returns a "cannot create bounded hat!"
    # error which is not very useful. So, instead of testing the error
    # message, we just ensure an error is raised.
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
    rng1 = TransformedDensityRejection(common_cont_dist, random_state=123)
    rng2 = TransformedDensityRejection(common_cont_dist)
    rng2.set_random_state(123)
    assert_equal(rng1.rvs(100), rng2.rvs(100))
    rng = TransformedDensityRejection(common_cont_dist, random_state=123)
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
    rng = TransformedDensityRejection(common_cont_dist)
    if size is None:
        assert np.isscalar(rng.rvs(size))
    else:
        if np.isscalar(size):
            size = (size, )
        assert rng.rvs(size).shape == size


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
            TransformedDensityRejection(common_cont_dist, domain=domain)

    @pytest.mark.parametrize("construction_points", [-1, 0, 0.1])
    def test_bad_construction_points_scalar(self, construction_points):
        with pytest.raises(ValueError, match=r"`construction_points` must be "
                                             r"a positive integer."):
            TransformedDensityRejection(
                common_cont_dist, construction_points=construction_points
            )

    def test_bad_construction_points_array(self):
        # empty array
        construction_points = []
        with pytest.raises(ValueError, match=r"`construction_points` must "
                                             r"either be a "
                                             r"scalar or a non-empty array."):
            TransformedDensityRejection(
                common_cont_dist, construction_points=construction_points
            )

        # construction_points not monotonically increasing
        construction_points = [1, 1, 1, 1, 1, 1]
        with pytest.warns(RuntimeWarning, match=r"33 : starting points not "
                                                r"strictly monotonically "
                                                r"increasing"):
            TransformedDensityRejection(
                common_cont_dist, construction_points=construction_points
            )

        # construction_points containing nans
        construction_points = [np.nan, np.nan, np.nan]
        with pytest.raises(UNURANError, match=r"50 : bad construction "
                                              r"points."):
            TransformedDensityRejection(
                common_cont_dist, construction_points=construction_points
            )

        # construction_points out of domain
        construction_points = [-10, 10]
        with pytest.warns(RuntimeWarning, match=r"50 : starting point out of "
                                                r"domain"):
            TransformedDensityRejection(
                common_cont_dist, domain=(-3, 3),
                construction_points=construction_points
            )

    @pytest.mark.parametrize("c", [-1., np.nan, np.inf, 0.1, 1.])
    def test_bad_c(self, c):
        msg = r"`c` must either be -0.5 or 0."
        with pytest.raises(ValueError, match=msg):
            TransformedDensityRejection(common_cont_dist, c=-1.)

    def test_bad_variant(self):
        msg = r"Invalid option for the `variant`"
        with pytest.raises(ValueError, match=msg):
            TransformedDensityRejection(common_cont_dist, variant='foo')

    u = [np.linspace(0, 1, num=1000), [], [[]], [np.nan],
         [-np.inf, np.nan, np.inf], 0,
         [[np.nan, 0.5, 0.1], [0.2, 0.4, np.inf], [-2, 3, 4]]]

    @pytest.mark.parametrize("u", u)
    def test_ppf_hat(self, u):
        # Increase the `max_squeeze_hat_ratio` so the ppf_hat is more
        # accurate.
        rng = TransformedDensityRejection(common_cont_dist,
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
        mv = rvs.mean(), rvs.var()
        assert_allclose(mv, mv_expected, atol=1e-2)
        # chi-squared test for goodness-of-fit
        _, freqs = np.unique(rvs, return_counts=True)
        freqs = freqs / freqs.sum()
        obs_freqs = np.zeros_like(pv)
        obs_freqs[: freqs.size] = freqs
        pval = chisquare(obs_freqs, pv).pvalue
        assert pval > 0.1

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
            DiscreteAliasUrn(common_discr_dist, domain=domain)

    def test_bad_urn_factor(self):
        with pytest.warns(RuntimeWarning, match=r"relative urn size < 1."):
            DiscreteAliasUrn([0.5, 0.5], urn_factor=-1)

    def test_bad_args(self):
        msg = (r"`domain` must be provided when the "
               r"probability vector is not available.")
        with pytest.raises(ValueError, match=msg):
            DiscreteAliasUrn(common_discr_dist)
