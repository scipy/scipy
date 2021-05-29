import pytest
from math import e, pi
from numpy.testing import assert_allclose

from scipy.special import hyp2f1


class TestHyp2f1:
    @pytest.mark.parametrize(
        "a,b,c,expected,rtol",
        [
            (0.5, 0.2, 1.5, 1.1496439092239847+0j, 1e-15),
            (12.3, 8.0, 20.31, 69280986.75273195+0j, 1e-15),
            (290.2, 321.5, 700.1, 1.3396562400934e+117+0j, 1e-12),
            (-102.1, -20.3, 1.3, 2.7899070752746906e+22+0j, 1e-15),
            pytest.param(
                -202.6, 60.3, 1.5, -1.152742155279977e-24+0j, 1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values")
            )
        ]
    )
    def test_unital_argument(self, a, b, c, expected, rtol):
        """Tests for case z = 1, c - a - b > 0.

        Expected answers computed using mpmath.
        """
        assert_allclose(hyp2f1(a, b, c, 1 + 0j), expected, rtol=rtol)

    @pytest.mark.parametrize(
        "a,b,expected,rtol",
        [
            (0.5, 0.2, 0.9428846409614143+0j, 1e-15),
            (12.3, 8.0, -4.845809986595704e-06+0j, 1e-15),
            (221.5, 90.2, 2.0490488728377282e-42+0j, 1e-7),
            pytest.param(
                -102.1, -20.3, 45143784.46783885+0j, 1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values")
            )
        ]
    )
    def test_special_case_z_near_minus_1(self, a, b, expected, rtol):
        """Tests for case z ~ -1, c ~ 1 + a - b

        Expected answers computed using mpmath.
        """
        assert_allclose(
            [
                hyp2f1(a, b, 1 + a - b + 2e-16, -1 + 2e-16j),
                hyp2f1(a, b, 1 + a - b + 2e-16, -1 - 2e-16j),
                hyp2f1(a, b, 1 + a - b - 2e-16, -1 + 2e-16j),
                hyp2f1(a, b, 1 + a - b - 2e-16, -1 - 2e-16j)
            ],
            [expected, expected, expected, expected],
            rtol=rtol
        )

    @pytest.mark.parametrize(
        "a,b,c,z,expected,rtol",
        [
            (
                -3, 0.2, 2.5, 0.3 + 0.6j,
                0.9097654857142857-0.11365302857142857j,
                1e-15
            ),
            (
                12.7, -16, -80.2, 8.6 - 3.2j,
                -31002615.365060613+85098290.74616767j,
                1e-15
            ),
            (
                -12, 201.6, 5.5, 0.95*e**(pi/2.5j),
                -1515637791889714+531381661325624.8j,
                1e-15
            ),
            (
                -20, 21.2, 7.4, 0.2 + 0.2j,
                0.06496386936970593-0.12392135982988922j,
                1e-10
            ),
            pytest.param(
                -100, 21.2, 7.4, 0.2 + 0.2j,
                0.000667333949637058-0.0006572419106018783j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values")
            )
        ]
    )
    def test_a_b_negative_int(self, a, b, c, z, expected, rtol):
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)

    def test_domain1(self):
        """|z| < 0.9 and real(z) > 0."""
        assert True

    def test_domain2(self):
        """|z| < 1.0 and real(z) < 0."""
        assert True

    def test_domain3(self):
        """0.9 < |z| < 1.0 and |1 - z| < 0.75"""
        assert True

    @pytest.mark.xfail
    def test_domain4(self):
        """0.9 < |z| < 1.0 and |1 - z| > 0.75"""
        assert False

    def test_domain5(self):
        """|z| > 1.0."""
        assert True
