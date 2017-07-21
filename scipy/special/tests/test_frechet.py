from __future__ import division, print_function, absolute_import

from numpy.testing import assert_allclose, run_module_suite

from scipy.special._ufuncs import (_frechet_pdf, _frechet_logpdf, _frechet_cdf,
                                   _frechet_ppf, _frechet_sf, _frechet_isf)


def test_frechet():
    # "Exact" values were computed using Wolfram Alpha online.
    # E.g. PDF[FrechetDistribution[1/2, 1], 5/2] for _frechet_pdf(2.5, 0.5)
    assert_allclose(_frechet_pdf(2.5, 0.5), 0.067202904517205343834,
                    rtol=1e-12)
    assert_allclose(_frechet_logpdf(2.5, 0.5), -2.7000388104048537736,
                    rtol=1e-12)
    assert_allclose(_frechet_cdf(2.5, 0.5), 0.53128560913296781152,
                    rtol=1e-12)
    assert_allclose(_frechet_cdf(1e-3, 0.5), 1.84672666240969310047e-14,
                    rtol=1e-12)
    assert_allclose(_frechet_sf(1e-3, 0.5), 0.99999999999998153273,
                    rtol=1e-12)
    assert_allclose(_frechet_ppf(0.975, 0.5), 1560.0833306626037620,
                    rtol=1e-12)
    assert_allclose(_frechet_ppf(0.5, 0.5), 2.0813689810056077979,
                    rtol=1e-12)
    assert_allclose(_frechet_isf(0.01, 0.5), 9900.0833329124631421,
                    rtol=1e-12)

    assert_allclose(_frechet_pdf(3, 13), 2.7179753508900831584e-6,
                    rtol=1e-12)
    assert_allclose(_frechet_logpdf(3, 13), -12.815623311117473330,
                    rtol=1e-12)
    assert_allclose(_frechet_cdf(3, 13), 0.99999937277472231955,
                    rtol=1e-12)
    assert_allclose(_frechet_sf(3, 13), 6.2722527768045018197e-7,
                    rtol=1e-12)
    assert_allclose(_frechet_ppf(0.975, 13), 1.3268241778612848659,
                    rtol=1e-12)
    assert_allclose(_frechet_ppf(0.5, 13), 1.0285944941498834248,
                    rtol=1e-12)
    assert_allclose(_frechet_ppf(0.01, 13), 0.88916242414515439302,
                    rtol=1e-12)
    assert_allclose(_frechet_ppf(1e-14, 13), 0.76554999794716080164,
                    rtol=1e-12)
    assert_allclose(_frechet_isf(0.995, 13), 0.87962401826455398956,
                    rtol=1e-12)


if __name__ == "__main__":
    run_module_suite()
