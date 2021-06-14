import pytest
import numpy as np
from numpy.testing import assert_allclose, suppress_warnings
from scipy.stats import TransformedDensityRejection, DiscreteAliasUrn
import scipy.stats as stats
from scipy.stats import chisquare


class TestTransformedDensityRejection:
    @pytest.mark.parametrize(
        "dist, dpdf, params",
        [
            (  # test case in gh-13051
                stats.norm,
                lambda x, loc, scale: (
                    -(x - loc)
                    / (scale ** 2)
                    * stats.norm._pdf((x - loc) / scale)
                ),
                (0.0, 0.1),
            ),
            (stats.expon, lambda x: -stats.expon._pdf(x), ()),
            (
                stats.laplace,
                lambda x: 0
                if x == 0
                else (
                    -stats.laplace._pdf(x)
                    if x > 0
                    else stats.laplace._pdf(x)
                ),
                (),
            ),
        ],
    )
    def test_sampling(self, dist, dpdf, params):
        domain = dist.support(*params)
        # call the private method to avoid validations and expensive
        # numpy operations (like broadcasting).
        pdf = lambda x, loc=0, scale=1: dist._pdf((x - loc) / scale)
        dist.pdf = pdf
        dist.dpdf = dpdf
        rng = TransformedDensityRejection(
            dist, params=params, domain=domain, seed=42
        )
        rvs = rng.rvs(100_000)
        # test if the first few moments match.
        mv_expected = dist.stats(*params, moments="mv")
        mv = rvs.mean(), rvs.var()
        assert_allclose(mv, mv_expected, atol=1e-1)

    @pytest.mark.parametrize(
        "bad_pdf, msg",
        [
            (lambda x: foo, r"name 'foo' is not defined"),
            (
                lambda x, a, b: x + a + b,
                r"missing 2 required positional arguments: 'a' and 'b'",
            ),
            (lambda x: -x, r"50 : PDF\(x\) < 0.!"),
            (lambda x: x * x, r"51 : hat\(x\) < PDF\(x\)"),
            (lambda x: None, r"must be real number, not NoneType"),
            (lambda x: np.inf, r"50 : PDF\(x\) overflow"),
            (lambda x: np.nan, r"51 : cannot create bounded hat!"),
            (
                lambda x: x * x / 2 + 2,
                r"51 : dTfx0 < dTfx1 \(x0<x1\). PDF not T-concave!",
            ),
            (lambda: 1.0, r"takes 0 positional arguments but 1 was given")
        ],
    )
    def test_bad_pdf(self, bad_pdf, msg):
        class dist:
            pdf = bad_pdf
            dpdf = lambda x: x

        with pytest.raises(Exception, match=msg):
            TransformedDensityRejection(dist, domain=(1, 10))

    @pytest.mark.parametrize(
        "bad_dpdf, msg",
        [
            (lambda x: foo, r"name 'foo' is not defined"),
            (lambda x: x * x, r"51 : hat\(x\) < PDF\(x\)"),
            (lambda x: None, r"must be real number, not NoneType"),
            (lambda x: np.inf, r"51 : cannot create bounded hat!"),
            (lambda x: np.nan, r"51 : cannot create bounded hat!"),
            (lambda: 1.0, r"takes 0 positional arguments but 1 was given")
        ],
    )
    def test_bad_dpdf(self, bad_dpdf, msg):
        class dist:
            pdf = lambda x: x * x / 2
            dpdf = bad_dpdf

        with pytest.raises(Exception, match=msg):
            TransformedDensityRejection(dist, domain=(1, 10))

    @pytest.mark.parametrize(
        "domain", [(0, 0), (1, 0), (np.inf, np.inf), (-np.inf, -np.inf)]
    )
    def test_bad_domain(self, domain):
        class dist:
            pdf = lambda x: x
            dpdf = lambda x: 1
        with pytest.raises(RuntimeError, match=r"left >= right"):
            TransformedDensityRejection(dist, domain=domain)

    def test_bad_sized_domain(self):
        class dist:
            pdf = lambda x: x
            dpdf = lambda x: 1
        with pytest.raises(ValueError, match=r"must be a length 2 tuple"):
            TransformedDensityRejection(dist, domain=(1, 2, 3))

    @pytest.mark.parametrize(
        "domain",
        [
            (np.nan, np.nan),
            (np.inf, np.nan),
            (np.nan, -np.inf),
            (np.nan, 0),
            (-1, np.nan),
            (0, float("nan")),
        ],
    )
    def test_nan_domain(self, domain):
        class dist:
            pdf = lambda x: x
            dpdf = lambda x: 1

        with pytest.raises(ValueError, match=r"only non-nan values"):
            TransformedDensityRejection(dist, domain=domain)

    def test_bad_cpoints(self):
        class dist:
            pdf = lambda x: x
            dpdf = lambda x: 1

        with pytest.warns(
            UserWarning, match=r"number of starting points < 0"
        ):
            TransformedDensityRejection(dist, domain=(0, 10), cpoints=-10)

        dist.pdf = lambda x: 1 - x * x
        dist.dpdf = lambda x: -2 * x
        with pytest.warns(UserWarning, match=r"hat/squeeze ratio too small"):
            TransformedDensityRejection(dist, domain=(-1, 1), cpoints=1)

    def test_bad_cpoints_array(self):
        class dist:
            pdf = lambda x: 1 - x * x
            dpdf = lambda x: -2 * x

        with pytest.warns(
            UserWarning,
            match=r"starting points not strictly "
            r"monotonically increasing",
        ):
            rng = TransformedDensityRejection(dist, domain=(-1, 1),
                                              cpoints=[1, 1, 1, 1, 1, 1])

        with pytest.raises(RuntimeError, match=r"bad construction points"):
            with pytest.warns(
                UserWarning, match=r"starting point out of " r"domain"
            ):
                cpoints = [1, 2, 3, 4, 5, 6]
                rng = TransformedDensityRejection(dist, domain=(-1, 1),
                                                  cpoints=cpoints)

        with pytest.raises(RuntimeError, match=r"bad construction points"):
            with pytest.warns(
                UserWarning, match=r"starting point out of " r"domain"
            ):
                cpoints = [np.nan, np.inf, np.nan]
                rng = TransformedDensityRejection(dist, domain=(-1, 1),
                                                  cpoints=cpoints)

    def test_bad_c(self):
        class dist:
            pdf = lambda x: x
            dpdf = lambda x: 1

        # c < -0.5
        with pytest.raises(
            RuntimeError, match=r"c < -0.5 not implemented yet"
        ):
            TransformedDensityRejection(
                dist, domain=(0, 10), c=-1.0
            )
        with pytest.raises(
            RuntimeError, match=r"c < -0.5 not implemented yet"
        ):
            TransformedDensityRejection(
                dist, domain=(0, 10), c=-np.inf
            )
        #  c > 0
        with pytest.warns(UserWarning, match=r"c > 0"):
            TransformedDensityRejection(
                dist, domain=(0, 10), c=10.0
            )
        # c = nan
        with pytest.raises(ValueError, match=r"must be a non-nan value"):
            TransformedDensityRejection(
                dist, domain=(0, 10), c=np.nan
            )

    def test_bad_variant(self):
        class dist:
            pdf = lambda x: x
            dpdf = lambda x: 1

        with pytest.raises(
            ValueError, match=r"Invalid option for the `variant`"
        ):
            TransformedDensityRejection(
                dist, domain=(0, 10), variant="foo"
            )

    # TODO: test other parameters


class TestDiscreteAliasUrn:
    @pytest.mark.parametrize(
        "dist, params",
        [  # discrete distributions with finite support.
            (stats.hypergeom, (20, 7, 12)),
            (stats.nhypergeom, (20, 7, 12)),
            (stats.binom, (20, 0.3)),
        ],
    )
    def test_sampling_with_pmf(self, dist, params):
        domain = dist.support(*params)
        with suppress_warnings() as sup:
            sup.filter(UserWarning)
            rng = DiscreteAliasUrn(
                dist=dist, domain=domain, params=params, seed=123
            )
        rvs = rng.rvs(100_000)
        # test if the first few moments match
        mv = rvs.mean(), rvs.var()
        mv_expected = dist.stats(*params, moments="mv")
        assert_allclose(mv, mv_expected, atol=1e-1)

        pv = dist.pmf(np.arange(domain[0], domain[1] + 1), *params)
        # correct for some numerical errors
        pv = pv / pv.sum()
        obs_freqs = np.zeros_like(pv)
        # chi-squared test for goodness-of-fit
        _, freqs = np.unique(rvs, return_counts=True)
        freqs = freqs / freqs.sum()
        obs_freqs[:freqs.size] = freqs
        with suppress_warnings() as sup:
            sup.filter(RuntimeWarning, "divide by zero encountered in "
                                       "true_divide")
            sup.filter(RuntimeWarning, "invalid value encountered in "
                                       "true_divide")
            pval = chisquare(obs_freqs, pv).pvalue
        assert_allclose(pval, 1.0, atol=1e-2)

    @pytest.mark.parametrize(
        "bad_pmf, msg",
        [
            (lambda x: foo, r"name 'foo' is not defined"),
            (
                lambda x, a, b: x + a + b,
                r"missing 2 required positional arguments: 'a' and 'b'",
            ),
            (lambda x: -x, r"50 : probability < 0"),
            (lambda x: None, r"must be real number, not NoneType"),
            (lambda x: np.inf, r"240 : unknown error"),
            (lambda x: np.nan, r"240 : unknown error"),
            (lambda x: 0.0, r"240 : unknown error"),
            (lambda: 1.0, r"takes 0 positional arguments but 1 was given"),
        ],
    )
    def test_bad_pmf(self, bad_pmf, msg):
        class dist:
            pmf = bad_pmf

        with pytest.raises(Exception, match=msg):
            with suppress_warnings() as sup:
                sup.filter(UserWarning)
                DiscreteAliasUrn(dist=dist, domain=(1, 10))

    @pytest.mark.parametrize(
        "pv", [[0.18, 0.02, 0.8], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]]
    )
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
        assert_allclose(pval, 1.0, atol=1e-3)

    @pytest.mark.parametrize(
        "bad_pv, msg",
        [
            ([], r"must have at least one element"),
            ([[1.0, 0.0]], r"must be a one-dimensional vector"),
            ([0.2, 0.4, np.nan, 0.8], r"must contain only non-nan values"),
            ([0.2, 0.4, np.inf, 0.8], r"must contain only finite values"),
            ([0.0, 0.0], r"must contain at least one non-zero value"),
        ],
    )
    def test_bad_pv(self, bad_pv, msg):
        with pytest.raises(ValueError, match=msg):
            DiscreteAliasUrn(bad_pv, domain=(1, 10))

    @pytest.mark.parametrize("domain", [(0, 0), (1, 0)])
    def test_bad_domain(self, domain):
        with pytest.raises(RuntimeError, match=r"left >= right"):
            class dist:
                pmf = lambda x: x

            with suppress_warnings() as sup:
                sup.filter(UserWarning)
                DiscreteAliasUrn(dist=dist, domain=domain)

    def test_bad_sized_domain(self):
        with pytest.raises(ValueError, match=r"must be a length 2 tuple"):
            DiscreteAliasUrn([0.02, 0.18, 0.8], domain=(1, 2, 3))

    @pytest.mark.parametrize(
        "domain",
        [
            (np.nan, np.nan),
            (np.inf, np.nan),
            (np.nan, -np.inf),
            (np.nan, 0),
            (-1, np.nan),
            (0, float("nan")),
        ],
    )
    def test_nan_domain(self, domain):
        with pytest.raises(
            ValueError, match=r"must contain only non-nan values"
        ):
            class dist:
                pmf = lambda x: x

            with suppress_warnings() as sup:
                sup.filter(UserWarning)
                DiscreteAliasUrn(dist=dist, domain=domain)

    @pytest.mark.parametrize(
        "domain",
        [
            (-np.inf, np.inf),
            (np.inf, np.inf),
            (-np.inf, -np.inf),
            (0, np.inf),
            (-np.inf, 0),
        ],
    )
    def test_inf_domain(self, domain):
        class dist:
                pmf = lambda x: x

        with pytest.raises(ValueError, match=r"must be finite"):
            with suppress_warnings() as sup:
                sup.filter(UserWarning)
                DiscreteAliasUrn(dist=dist, domain=domain)

    @pytest.mark.parametrize(
        "pv",
        [
            [0.0],
            [],
            [0.0, 0.0, 0.0],
            [0.0, np.inf],
            [-np.inf, 0.0],
            [np.inf, np.inf],
            [-np.inf, -np.inf],
            [np.nan],
            [np.nan, np.inf, -np.inf],
            [[1.0, 0.0], [0.5, 0.5]],
        ],
    )
    def test_bad_pv(self, pv):
        with pytest.raises(ValueError):
            DiscreteAliasUrn(pv)

    def test_bad_urn_factor(self):
        with pytest.warns(UserWarning, match=r"relative urn size < 1."):
            DiscreteAliasUrn([0.5, 0.5], urn_factor=-1)
