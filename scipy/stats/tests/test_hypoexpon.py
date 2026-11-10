import numpy as np
from numpy.testing import assert_allclose, assert_equal
import pytest

from scipy import stats
from scipy.integrate import quad
from scipy.stats import hypoexpon, erlang, expon
from scipy.stats._hypoexpon import (_closed_form, _phase_type_series,
                                    _phase_type_squaring)

from .common_tests import check_random_state_property


# Reference values computed with mpmath at 500 significant digits from the
# closed form sum; the repeated-scale case [1, 1, 2] was evaluated with the
# middle scale perturbed by 1e-80.
_REFERENCE = {
    (1., 2., 4.): [
        (0.01, 6.21365534712053e-6, 2.0742414944748814e-8, 0.99999997925758506),
        (0.5, 0.011707371889203208, 0.0020996060130108543, 0.99790039398698915),
        (3., 0.10837656446820129, 0.17002049019819912, 0.82997950980180088),
        (12., 0.03071467480602738, 0.87218994063491109, 0.12781005936508891),
        (40., 3.0264558688035545e-5, 0.99987893764294062, 0.00012106235705938281),
    ],
    (1., 1.00000001, 1.00000002): [
        (0.001, 4.9950023493667535e-7, 1.6654171165777837e-10, 0.99999999983345829),
        (1., 0.1839397169069268, 0.080301395231997022, 0.91969860476800298),
        (3., 0.2240418076553877, 0.57680991215190229, 0.42319008784809771),
        (20., 4.1223079456694079e-7, 0.99999954448496699, 4.5551503300507387e-7),
    ],
    (2., 2.000000000001, 2.000000000002, 2.000000000003): [
        (0.5, 0.001014063519621373, 0.00013336965051368211, 0.99986663034948632),
        (8., 0.097683407406582295, 0.56652987963270497, 0.43347012036729503),
        (30., 8.6035027641848292e-5, 0.9997886214965313, 0.00021137850346869741),
    ],
    (0.001, 0.05, 1., 30.): [
        (1e-5, 1.1082805727130982e-13, 2.7721148147928049e-19, 1.0),
        (0.01, 2.5650029475944531e-5, 7.9845174893724297e-8, 0.99999992015482511),
        (1., 0.020042550616842662, 0.011094335436676217, 0.98890566456332378),
        (60., 0.0046746806020297552, 0.85975958193910734, 0.14024041806089266),
        (400., 5.5943265695882786e-8, 0.99999832170202912, 1.6782979708764836e-6),
    ],
    (1., 1., 2.): [
        (0.1, 0.0023002711259129145, 7.8297908618640443e-5, 0.99992170209138136),
        (1., 0.10942299591093988, 0.045395125835235592, 0.95460487416476441),
        (4., 0.1607767331408203, 0.58686833927468849, 0.41313166072531151),
        (15., 0.0011009684008471361, 0.9977931687611777, 0.0022068312388223015),
    ],
}


class TestHypoexpon:

    @pytest.mark.parametrize('scales, rows', _REFERENCE.items())
    def test_against_reference(self, scales, rows):
        x = np.array([row[0] for row in rows])
        ref = np.array([row[1:] for row in rows])
        assert_allclose(hypoexpon.pdf(x, scales), ref[:, 0], rtol=1e-12)
        assert_allclose(hypoexpon.cdf(x, scales), ref[:, 1], rtol=1e-12)
        assert_allclose(hypoexpon.sf(x, scales), ref[:, 2], rtol=1e-12)
        assert_allclose(hypoexpon.logpdf(x, scales), np.log(ref[:, 0]),
                        rtol=1e-12)
        # scalar input gives scalar output with the same values
        for xi, pdf_i, cdf_i, sf_i in rows:
            assert np.isscalar(hypoexpon.pdf(xi, scales))
            assert_allclose(hypoexpon.pdf(xi, scales), pdf_i, rtol=1e-12)
            assert_allclose(hypoexpon.cdf(xi, scales), cdf_i, rtol=1e-12)
            assert_allclose(hypoexpon.sf(xi, scales), sf_i, rtol=1e-12)

    def test_cancellation_in_closed_form_is_detected(self):
        # gh-25639: the closed form gives a negative density here
        scales = [1., 1.00000001, 1.00000002]
        value, rel_err = _closed_form(np.array([1.]), np.array(scales), 'pdf')
        assert value[0] < 0 and rel_err[0] > 1
        assert_allclose(hypoexpon.pdf(1., scales), 0.1839397169069268,
                        rtol=1e-13)

    def test_single_scale_is_exponential(self):
        x = np.array([0., 0.1, 1., 5., 50.])
        assert_allclose(hypoexpon.pdf(x, [2.]), expon.pdf(x, scale=2.),
                        rtol=1e-14)
        assert_allclose(hypoexpon.cdf(x, [2.]), expon.cdf(x, scale=2.),
                        rtol=1e-14)
        assert_allclose(hypoexpon.sf(x, [2.]), expon.sf(x, scale=2.),
                        rtol=1e-14)

    @pytest.mark.parametrize('n', [2, 3, 7, 40])
    def test_equal_scales_is_erlang(self, n):
        scale = 0.7
        x = np.linspace(0, 10 * n * scale, 25)
        rv = hypoexpon([scale] * n)
        assert_allclose(rv.pdf(x), erlang.pdf(x, n, scale=scale),
                        rtol=1e-12, atol=1e-300)
        assert_allclose(rv.cdf(x), erlang.cdf(x, n, scale=scale),
                        rtol=1e-12, atol=1e-300)
        assert_allclose(rv.sf(x), erlang.sf(x, n, scale=scale),
                        rtol=1e-12, atol=1e-300)

    def test_series_and_squaring_agree(self):
        rng = np.random.default_rng(1836459265409)
        scales = rng.uniform(0.5, 2, size=6)
        x = np.geomspace(1e-6, 1e3, 30)
        for kind in ['pdf', 'cdf', 'sf']:
            series = _phase_type_series(x, scales, kind)
            squaring = _phase_type_squaring(x, scales, kind)
            assert_allclose(series, squaring, rtol=1e-12, atol=1e-300)

    def test_far_tail(self):
        # x mu large: the squaring variant handles these points
        scales = [1e-3, 1.]
        x = np.array([10., 20., 50.])
        # the smallest scale contributes a negligible factor here
        assert_allclose(hypoexpon.sf(x, scales), np.exp(-x) / (1 - 1e-3),
                        rtol=1e-10)
        assert_allclose(hypoexpon.pdf(x, scales), np.exp(-x) / (1 - 1e-3),
                        rtol=1e-10)

    def test_consistency(self):
        scales = np.array([0.3, 1., 1., 2.5])
        x = np.array([0.2, 1., 3., 8.])
        assert_allclose(hypoexpon.cdf(x, scales) + hypoexpon.sf(x, scales), 1,
                        rtol=1e-14)
        for xi in x:
            integral, _ = quad(hypoexpon.pdf, 0, xi, args=(scales,))
            assert_allclose(integral, hypoexpon.cdf(xi, scales), rtol=1e-10)
        assert_allclose(hypoexpon.mean(scales), scales.sum())
        assert_allclose(hypoexpon.var(scales), (scales**2).sum())

    def test_edge_inputs(self):
        scales = [1., 2.]
        x = np.array([[-1., 0.], [np.inf, np.nan]])
        assert_equal(hypoexpon.pdf(x, scales), [[0., 0.], [0., np.nan]])
        assert_equal(hypoexpon.cdf(x, scales), [[0., 0.], [1., np.nan]])
        assert_equal(hypoexpon.sf(x, scales), [[1., 1.], [0., np.nan]])
        assert_equal(hypoexpon.logpdf(x, scales),
                     [[-np.inf, -np.inf], [-np.inf, np.nan]])
        assert_equal(hypoexpon.pdf(0., [2.]), 0.5)
        assert hypoexpon.pdf(np.array([1.]), scales).shape == ()
        assert hypoexpon.pdf(np.ones((2, 3)), scales).shape == (2, 3)

    @pytest.mark.parametrize('scales', [[], [[1., 2.]], [1., -1.], [1., 0.],
                                        [1., np.inf], [np.nan]])
    def test_invalid_scales(self, scales):
        message = "Parameter vector `scales`|All entries of `scales`"
        with pytest.raises(ValueError, match=message):
            hypoexpon.pdf(1., scales)
        with pytest.raises(ValueError, match=message):
            hypoexpon(scales)

    def test_frozen(self):
        rng = np.random.default_rng(2846)
        scales = rng.uniform(0.1, 10, size=5)
        x = rng.uniform(0, 30, size=10)
        rv = hypoexpon(scales)
        assert_equal(rv.pdf(x), hypoexpon.pdf(x, scales))
        assert_equal(rv.logpdf(x), hypoexpon.logpdf(x, scales))
        assert_equal(rv.cdf(x), hypoexpon.cdf(x, scales))
        assert_equal(rv.sf(x), hypoexpon.sf(x, scales))
        assert_equal(rv.mean(), hypoexpon.mean(scales))
        assert_equal(rv.var(), hypoexpon.var(scales))
        assert_equal(rv.rvs(size=3, random_state=1),
                     hypoexpon.rvs(scales, size=3, random_state=1))

    @pytest.mark.parametrize('scales', [[2.], [1., 2., 4.], [1., 1., 1.],
                                        [0.1, 0.1, 5.]])
    def test_rvs(self, scales):
        rng = np.random.default_rng(6148309432765)
        rv = hypoexpon(scales)
        x = rv.rvs(size=50000, random_state=rng)
        assert x.shape == (50000,)
        assert_allclose(x.mean(), rv.mean(), rtol=0.02)
        assert_allclose(x.var(), rv.var(), rtol=0.05)
        assert stats.kstest(x, rv.cdf).pvalue > 0.01

    def test_rvs_shape(self):
        rng = np.random.default_rng(19460293467)
        assert hypoexpon.rvs([1., 2.], random_state=rng).shape == (1,)
        assert hypoexpon.rvs([1., 2.], size=4, random_state=rng).shape == (4,)
        assert hypoexpon.rvs([1., 2.], size=(2, 3),
                             random_state=rng).shape == (2, 3)

    def test_random_state_property(self):
        check_random_state_property(hypoexpon, ([1., 2.],))
