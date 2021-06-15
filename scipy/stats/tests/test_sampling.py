import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal
from numpy.lib import NumpyVersion
from scipy.stats import TransformedDensityRejection, DiscreteAliasUrn
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
    dpdf = lambda x: -x * stats.norm._pdf(x)
    cdf = stats.norm._cdf


# A binomial distribution to share between all the discrete methods
class common_discr_dist:
    params = (10, 0.2)
    pmf = stats.binom._pmf
    cdf = stats.binom._cdf


all_methods = [
    ("TransformedDensityRejection", {"dist": common_cont_dist}),
    ("DiscreteAliasUrn", {"pv": [0.02, 0.18, 0.8]})
]


bad_pdfs_common = [
    # Negative PDF
    (lambda x: -x, RuntimeError, r"50 : PDF\(x\) < 0.!"),
    # Returning wrong type
    (lambda x: [], TypeError, r"must be real number, not list"),
    # Undefined name inside the function
    (lambda x: foo, NameError, r"name 'foo' is not defined"),  # type: ignore[name-defined]
    # Infinite value returned => Overflow error.
    (lambda x: np.inf, RuntimeError, r"50 : PDF\(x\) overflow"),
    # NaN value => internal error in UNU.RAN
    # Currently, UNU.RAN just returns a "cannot create bounded hat!"
    # error which is not very useful. So, instead of testing the error
    # message, we just ensure an error is raised.
    (lambda x: np.nan, RuntimeError, r"..."),
    # signature of PDF wrong
    (lambda: 1.0, TypeError, r"takes 0 positional arguments but 1 was given")
]

bad_dpdf_common = [
    # Infinite value returned. For dPDF, UNU.RAN complains that it cannot
    # create a bounded hat instead of an overflow error.
    (lambda x: np.inf, RuntimeError, r"51 : cannot create bounded hat"),
    # NaN value => internal error in UNU.RAN
    # Currently, UNU.RAN just returns a "cannot create bounded hat!"
    # error which is not very useful. So, instead of testing the error
    # message, we just ensure an error is raised.
    (lambda x: np.nan, RuntimeError, r"..."),
    # Returning wrong type
    (lambda x: [], TypeError, r"must be real number, not list"),
    # Undefined name inside the function
    (lambda x: foo, NameError, r"name 'foo' is not defined"),  # type: ignore[name-defined]
    # signature of dPDF wrong
    (lambda: 1.0, TypeError, r"takes 0 positional arguments but 1 was given")
]

bad_pmf_common = [
    # UNU.RAN fails to validate float inf, nan values returned
    # by the PMF and throws an unhelpful "unknown error". One
    # helpful thing to do here is to calculate the PV ourselves
    # and check for inf, nan, if the domain is known and distribution
    # doesn't have infinite tails.
    (lambda x: np.inf, RuntimeError, r"240 : unknown error"),
    (lambda x: np.nan, RuntimeError, r"240 : unknown error"),
    (lambda x: 0.0, RuntimeError, r"240 : unknown error"),
    # Undefined name inside the function
    (lambda x: foo, NameError, r"name 'foo' is not defined"),  # type: ignore[name-defined]
    # Returning wrong type
    (lambda x: [], TypeError, r"must be real number, not list"),
    # probabilities < 0
    (lambda x: -x, RuntimeError, r"50 : probability < 0"),
    # signature of PMF wrong
    (lambda: 1.0, TypeError, r"takes 0 positional arguments but 1 was given")
]

bad_pv_common = [
    ([], r"must have at least one element"),
    ([[1.0, 0.0]], r"wrong number of dimensions \(expected 1, got 2\)"),
    ([0.2, 0.4, np.nan, 0.8], r"must contain only non-nan values"),
    ([0.2, 0.4, np.inf, 0.8], r"must contain only finite values"),
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
    ((2, 1), RuntimeError, r"left >= right"),
    ((1, 1), RuntimeError, r"left >= right"),
]

# infinite and nan values present in domain.
inf_nan_domains = [
    # left >= right
    ((10, 10), RuntimeError, r"left >= right"),
    ((np.inf, np.inf), RuntimeError, r"left >= right"),
    ((-np.inf, -np.inf), RuntimeError, r"left >= right"),
    ((np.inf, -np.inf), RuntimeError, r"left >= right"),
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
                         bad_domains + bad_sized_domains + nan_domains)  # type: ignore[operator]
@pytest.mark.parametrize("method, kwargs", all_methods)
def test_bad_domain(domain, err, msg, method, kwargs):
    Method = getattr(stats, method)
    with pytest.raises(err, match=msg):
        Method(**kwargs, domain=domain)


@pytest.mark.parametrize("method, kwargs", all_methods)
def test_seed(method, kwargs):
    Method = getattr(stats, method)

    # simple seed that works for any version of NumPy
    seed = 123
    rng1 = Method(**kwargs, seed=seed)
    rng2 = Method(**kwargs, seed=seed)
    assert_equal(rng1.rvs(100), rng2.rvs(100))

    # RandomState seed for old numpy
    if NumpyVersion(np.__version__) < '1.19.0':
        seed1 = np.random.RandomState(123)
        seed2 = 123
        rng1 = Method(**kwargs, seed=seed1)
        rng2 = Method(**kwargs, seed=seed2)
        assert_equal(rng1.rvs(100), rng2.rvs(100))
    else:  # Generator seed for new NumPy
        seed1 = np.random.default_rng(123)
        seed2 = np.random.PCG64(123)
        rng1 = Method(**kwargs, seed=seed1)
        rng2 = Method(**kwargs, seed=seed2)
        assert_equal(rng1.rvs(100), rng2.rvs(100))

        # when a RandomState is given, it should take the bitgen_t
        # member of the class and create a Generator instance.
        seed1 = np.random.RandomState(np.random.MT19937(123))
        seed2 = np.random.Generator(np.random.MT19937(123))
        rng1 = Method(**kwargs, seed=seed1)
        rng2 = Method(**kwargs, seed=seed2)
        assert_equal(rng1.rvs(100), rng2.rvs(100))

    # testing with seed sequence
    seed = [1, 2, 3]
    rng1 = Method(**kwargs, seed=seed)
    rng2 = Method(**kwargs, seed=seed)
    assert_equal(rng1.rvs(100), rng2.rvs(100))


class TestTransformedDensityRejection:
    pdf1 = lambda x: 1-x*x if abs(x) <= 1 else 0
    pdf2 = lambda x: stats.norm._pdf(x / 0.1)  # test case in gh-13051
    pdfs = [pdf1, pdf2]

    dpdf1 = lambda x: -2*x if abs(x) <= 1 else 0
    dpdf2 = lambda x: -x / 0.01 * stats.norm._pdf(x / 0.1)
    dpdfs = [dpdf1, dpdf2]

    cdf1 = lambda x: 3/4 * (x - x**3/3 + 2/3) if abs(x) <= 1 else 1*(x >= 1)
    cdf1 = np.vectorize(cdf1)
    cdf2 = stats.norm(0, 0.1).cdf
    cdfs = [cdf1, cdf2]

    mv1 = [0., 4./15.]
    mv2 = [0., 0.01]
    mvs = [mv1, mv2]

    @pytest.mark.parametrize("pdf, dpdf, mv_ex, cdf",
                             zip(pdfs, dpdfs, mvs, cdfs))
    def test_basic(self, pdf, dpdf, mv_ex, cdf):
        class dist:
            pass
        dist.pdf = pdf
        dist.dpdf = dpdf
        rng = TransformedDensityRejection(dist, seed=42)
        rvs = rng.rvs(100000)
        mv = rvs.mean(), rvs.var()
        assert_allclose(mv, mv_ex, rtol=1e-7, atol=1e-1)
        # Cramer Von Mises test for goodness-of-fit
        rvs = rng.rvs(500)
        pval = cramervonmises(rvs, cdf).pvalue
        assert pval > 0.1

    # PDF 0 everywhere => bad construction points
    bad_pdfs = [(lambda x: 0, RuntimeError, r"50 : bad construction points.")]
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


    # TODO: for cpoints < 0, UNU.RAN throws a warning and sets cpoints
    #       to a default value. This is not consistent with other invalid
    #       values of cpoints for which it throws errors. Hence, We should
    #       validate this parameter ourself.
    @pytest.mark.parametrize("cpoints", [0, 0.1])
    def test_bad_cpoints_scalar(self, cpoints):
        with pytest.raises(RuntimeError, match=r"50 : bad construction "
                                               r"points."):
            TransformedDensityRejection(common_cont_dist, cpoints=cpoints)

    def test_bad_cpoints_array(self):
        # empty array
        cpoints = []
        with pytest.raises(ValueError, match=r"`cpoints` must either be a "
                                             r"scalar or a non-empty array."):
            TransformedDensityRejection(common_cont_dist, cpoints=cpoints)

        # cpoints not monotonically increasing
        cpoints = [1, 1, 1, 1, 1, 1]
        with pytest.warns(UserWarning, match=r"33 : starting points not "
                                             r"strictly monotonically "
                                             r"increasing"):
            TransformedDensityRejection(common_cont_dist, cpoints=cpoints)

        # cpoints containing nans
        cpoints = [np.nan, np.nan, np.nan]
        with pytest.raises(RuntimeError, match=r"50 : bad construction "
                                               r"points."):
            TransformedDensityRejection(common_cont_dist, cpoints=cpoints)

        # cpoints out of domain
        cpoints = [-10, 10]
        with pytest.warns(UserWarning, match=r"50 : starting point out of "
                                             r"domain"):
            TransformedDensityRejection(common_cont_dist, domain=(-3, 3),
                                        cpoints=cpoints)

    def test_bad_c(self):
        # c < -0.5 => Not Implemented Error
        msg = r"33 : c < -0.5 not implemented yet"
        with pytest.raises(RuntimeError, match=msg):
            TransformedDensityRejection(common_cont_dist, c=-1.)
        # -0.5 < c < 0. => Warning: Not recommended. Using default
        msg = (r"33 : -0.5 < c < 0 not recommended. using c = -0.5 "
                r"instead.")
        with pytest.warns(UserWarning, match=msg):
            TransformedDensityRejection(common_cont_dist, c=-0.1)
        # c > 0. => Warning: Using default
        msg = r"33 : c > 0"
        with pytest.warns(UserWarning, match=msg):
            TransformedDensityRejection(common_cont_dist, c=1.)
        # nan c
        msg = r"`c` must be a non-nan value."
        with pytest.raises(ValueError, match=msg):
            TransformedDensityRejection(common_cont_dist, c=np.nan)

    def test_bad_variant(self):
        msg = r"Invalid option for the `variant`"
        with pytest.raises(ValueError, match=msg):
            TransformedDensityRejection(common_cont_dist, variant='foo')


class TestDiscreteAliasUrn:
    # DAU fails on these probably because of large domains and small
    # computation errors in PMF. Mean/SD match but chi-squared test fails.
    basic_fail_dists = {'nchypergeom_fisher', 'nchypergeom_wallenius'}

    @pytest.mark.parametrize("distname, params", distdiscrete)
    def test_basic(self, distname, params):
        if distname in self.basic_fail_dists:
            msg = ("DAU fails on these probably because of large domains "
                   "and small computation errors in PMF.")
            pytest.skip(msg)
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
        rng = DiscreteAliasUrn(pv, domain=domain, seed=42)
        rvs = rng.rvs(100000)
        # test if the first few moments match
        mv = rvs.mean(), rvs.var()
        assert_allclose(mv, mv_ex, rtol=1e-3, atol=1e-1)
        # correct for some numerical errors
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
         r"must contain only finite values"),
        # nan returned
        (lambda x: np.nan, ValueError,
         r"must contain only non-nan values"),
        # all zeros
        (lambda x: 0.0, ValueError,
         r"must contain at least one non-zero value"),
        # Undefined name inside the function
        (lambda x: foo, NameError,  # type: ignore[name-defined]
         r"name 'foo' is not defined"),
        # Returning wrong type.
        (lambda x: [], ValueError,
         r"setting an array element with a sequence."),
        # probabilities < 0
        (lambda x: -x, RuntimeError,
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
            DiscreteAliasUrn(dist=dist, domain=(1, 10))

    @pytest.mark.parametrize("pv", [[0.18, 0.02, 0.8],
                                    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]])
    def test_sampling_with_pv(self, pv):
        pv = np.asarray(pv, dtype=np.float64)
        rng = DiscreteAliasUrn(pv, seed=123)
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
            DiscreteAliasUrn(dist=common_discr_dist, domain=domain,
                             params=common_discr_dist.params)

    def test_bad_urn_factor(self):
        with pytest.warns(UserWarning, match=r"relative urn size < 1."):
            DiscreteAliasUrn([0.5, 0.5], urn_factor=-1)

    def test_bad_args(self):
        msg = (r"Either a `pv` or a `dist` object with a PMF method "
               r"required but none given.")
        with pytest.raises(ValueError, match=msg):
            DiscreteAliasUrn()

        msg = r"`domain` must be provided if `pv` is not available"
        with pytest.raises(ValueError, match=msg):
            DiscreteAliasUrn(dist=common_discr_dist,
                             params=common_discr_dist.params)
