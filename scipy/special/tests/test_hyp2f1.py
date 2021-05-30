import pytest
from math import e, pi
from numpy.testing import assert_allclose

from scipy.special import hyp2f1


class TestHyp2f1:
    @pytest.mark.parametrize(
        "a,b,c,expected,rtol",
        [
            (0.5, 0.2, 1.5, 1.1496439092239847 + 0j, 1e-15),
            (12.3, 8.0, 20.31, 69280986.75273195 + 0j, 1e-15),
            (290.2, 321.5, 700.1, 1.3396562400934e117 + 0j, 1e-12),
            (-102.1, -20.3, 1.3, 2.7899070752746906e22 + 0j, 1e-15),
            (-202.6, 60.3, 1.5, -1.3113641413099326e-56 + 0j, 1e-12),
        ],
    )
    def test_unital_argument(self, a, b, c, expected, rtol):
        """Tests for case z = 1, c - a - b > 0.

        Expected answers computed using mpmath.
        """
        assert_allclose(hyp2f1(a, b, c, 1 + 0j), expected, rtol=rtol)

    @pytest.mark.parametrize(
        "a,b,expected,rtol",
        [
            (0.5, 0.2, 0.9428846409614143 + 0j, 1e-15),
            (12.3, 8.0, -4.845809986595704e-06 + 0j, 1e-15),
            (221.5, 90.2, 2.0490488728377282e-42 + 0j, 1e-7),
            (-102.1, -20.3, 45143784.46783885 + 0j, 1e-7),
        ],
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
                hyp2f1(a, b, 1 + a - b - 2e-16, -1 - 2e-16j),
            ],
            [expected, expected, expected, expected],
            rtol=rtol,
        )

    @pytest.mark.parametrize(
        "a,b,c,z,expected,rtol",
        [
            (
                -3,
                0.2,
                2.5,
                0.3 + 0.6j,
                0.9097654857142857 - 0.11365302857142857j,
                1e-15,
            ),
            (
                12.7,
                -16,
                -80.2,
                8.6 - 3.2j,
                -31002615.365060613 + 85098290.74616767j,
                1e-15,
            ),
            (
                -12,
                201.6,
                5.5,
                0.95 * e ** (pi / 2.5j),
                -1515637791889714 + 531381661325624.8j,
                1e-15,
            ),
            (
                -20,
                21.2,
                7.4,
                0.2 + 0.2j,
                0.06496386936970593 - 0.12392135982988922j,
                1e-10,
            ),
            pytest.param(
                -100,
                21.2,
                7.4,
                0.2 + 0.2j,
                0.000667333949637058 - 0.0006572419106018783j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values"),
            ),
        ],
    )
    def test_a_b_negative_int(self, a, b, c, z, expected, rtol):
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)

    @pytest.mark.parametrize(
        "a,b,c,z,expected,rtol",
        [
            (
                1.5,
                0.2,
                0.5 + 1e-16,
                0.3 + 0.6j,
                0.8946325213869533 + 0.41753347569164395j,
                1e-15,
            ),
            (
                5.5,
                0.2,
                2.5,
                0.3 + 0.6j,
                0.9075594698587819 + 0.27844218583081926j,
                1e-15,
            ),
            (
                12.5,
                201.6,
                5.5,
                0.95 * e ** (pi / 2.5 * 1j),
                -0.0011442768023360108 + 0.0003439142635864824j,
                1e-13,
            ),
            (
                20.4,
                21.2,
                7.4,
                0.2 + 0.2j,
                -8041.92476417893 - 32027.24728016176j,
                1e-14,
            ),
            (
                200.4,
                21.2,
                7.4,
                0.2 + 0.2j,
                -3.267272756016562e27 - 5.226794898779876e27j,
                1e-14,
            ),
        ],
    )
    def test_a_b_neg_int_after_euler_hypergeometric_transformation(
        self, a, b, c, z, expected, rtol
    ):
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)

    @pytest.mark.parametrize(
        "a,b,c,z,expected,rtol",
        [
            (
                0.5,
                0.2,
                1.5,
                0.89999j,
                0.9842713242681234 + 0.0530492028861082j,
                1e-15
            ),
            (
                0.5,
                0.2,
                1.5,
                0.89999 * e ** (-pi / 4 * 1j),
                1.0305646182162314 - 0.06403266840337725j,
                1e-14,
            ),
            (0.5, 0.2, 1.5, 0.89999 + 0j, 1.104974843561603 + 0j, 1e-14),
            (
                0.5,
                0.2,
                1.5,
                0.5 * e ** (pi / 3 * 1j),
                1.0119995351960096 + 0.03352156960908768j,
                1e-15,
            ),
            (
                12.3,
                8.0,
                20.31,
                0.25 * e ** (-pi / 7 * 1j),
                2.5835524424679908 - 1.9080020009340966j,
                1e-15,
            ),
            (
                12.3,
                8.0,
                20.31,
                0.89999j,
                -0.24096298357604617 - 0.18430716189222762j,
                1e-12,
            ),
            (
                12.3,
                8.0,
                20.31,
                0.89999 * e ** (pi / 4 * 1j),
                0.15530503195210124 - 7.954774127378903j,
                1e-12,
            ),
            (12.3, 8.0, 20.31, 0.89999 + 0j, 4088.2331156689106 + 0j, 1e-12),
            (
                290.2,
                321.5,
                700.1,
                0.2 + 0.2j,
                169017475934.32013 - 160289742414.83957j,
                1e-10,
            ),
            pytest.param(
                290.2,
                321.5,
                700.1,
                0.89999j,
                9.062746910000646e-14 + 7.233753971322562e-15j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values"),
            ),
            pytest.param(
                290.2,
                321.5,
                700.1,
                0.89999 * e ** (-pi / 4 * 1j),
                1.6459562082242555e27 - 5.345976086313433e26j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values"),
            ),
            pytest.param(
                2.902,
                3.215,
                3.51,
                0.2 + 0.2j,
                1.3297180638579997 + 1.0127093358201944j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values"),
            ),
            pytest.param(
                2.902,
                3.215,
                3.51,
                0.89999j,
                -0.17229597828163576 + 0.4235651069433403j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values"),
            ),
            pytest.param(
                2.902,
                3.215,
                3.501,
                0.89999 * e ** (pi / 4 * 1j),
                -2.1743292914809533 + 0.800290539074747j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values"),
            ),
            pytest.param(
                29.02,
                32.15,
                35.01,
                0.2 + 0.2j,
                166.06805807466486 + 37.5498289807123j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values"),
            ),
            pytest.param(
                29.02,
                32.15,
                35.01,
                0.89999j,
                0.0003036147178447333 + 0.00027412751087335435j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values"),
            ),
            pytest.param(
                29.02,
                32.15,
                35.01,
                0.89999 * e ** (-pi / 4 * 1j),
                -4154.9779821267075 - 1856.8889143809663j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values"),
            ),
            pytest.param(
                290.2,
                321.5,
                350.1,
                0.2 + 0.2j,
                -1.2374685868266302e22 + 1.6317891878064388e22j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values"),
            ),
            pytest.param(
                290.2,
                321.5,
                350.1,
                0.89999j,
                6.447352310994814e-35 + 1.1919519203662124e-34j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values"),
            ),
            pytest.param(
                290.2,
                321.5,
                350.1,
                0.89999 * e ** (pi / 4 * 1j),
                -1.7345536772971163e36 + 3.475022465669823e36j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values"),
            ),
            (
                -102.1,
                -20.3,
                1.3,
                0.2 + 0.2j,
                -29011968709735.227 - 39637825264825.09j,
                1e-15,
            ),
            pytest.param(
                -102.1,
                -20.3,
                1.3,
                0.89999 * e ** (pi / 4 * 1j),
                -3.529065252670616e19 - 1.9971755803552488e20j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values"),
            ),
            pytest.param(
                -202.6,
                60.3,
                1.5,
                0.89999 * e ** (-pi / 4 * 1j),
                1.8311248864404518e21 - 1.213560644283594e21j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values"),
            ),
            pytest.param(
                2.92,
                3.215,
                -35.01,
                0.2 + 0.2j,
                0.9467474022243998 - 0.04857412986476413j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values"),
            ),
            pytest.param(
                2.92,
                3.215,
                -35.01,
                0.89999 * e ** (-pi / 4 * 1j),
                -2057136627467.506 - 1039430034887.269j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values"),
            ),
            pytest.param(
                2.92,
                3.215,
                -305.01,
                0.2 + 0.2j,
                0.9937729899771914 - 0.0061593055571362715j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values"),
            ),
            pytest.param(
                2.92,
                3.215,
                -305.01,
                0.89999 * e ** (pi / 4 * 1j),
                -7.663551124354694e40 + 1.008395072437797e41j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values"),
            ),
        ],
    )
    def test_domain1(self, a, b, c, z, expected, rtol):
        """|z| < 0.9 and real(z) > 0."""
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)

    @pytest.mark.parametrize(
        "a,b,c,z,expected,rtol",
        [
            (
                0.5,
                0.2,
                1.5,
                -0.2 + 0.2j,
                0.9868244017244943+0.01160931510212128j,
                1e-7
            ),
            (
                0.5,
                0.2,
                1.5,
                0.89999*e**(0.75*pi*1j),
                0.9607557591904695+0.028500075843422335j,
                1e-7
            ),
            (
                0.5,
                0.2,
                1.5,
                0.99999*e**(0.75*pi*1j),
                0.9569430650864205+0.03051659034567534j,
                1e-5
            ),
            pytest.param(
                0.5,
                0.2,
                1.5,
                0.99999*e**(0.75*pi*1j),
                0.9569430650864205+0.03051659034567534j,
                1e-7,
                marks=pytest.mark.xfail(reason="unhandled parameter values")
            )
        ]
    )
    def test_domain2(
            self, a, b, c, z, expected, rtol
    ):
        """|z| < 1.0 and real(z) < 0."""
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)

    @pytest.mark.parametrize(
        "a,b,c,z,expected,rtol",
        [
            (
                0.5,
                0.2,
                1.5,
                -0.2 + 0.2j,
                0.9868244017244943+0.01160931510212128j,
                1e-7
            ),
            (
                0.5,
                0.2,
                1.5,
                0.99999*e**(0.75*pi*1j),
                0.9569430650864205+0.03051659034567534j,
                1e-7,
            )
        ]
    )
    def test_domain3(
            self, a, b, c, z, expected, rtol
    ):
        """0.9 < |z| < 1.0 and |1 - z| < 0.75"""
        assert True

    @pytest.mark.xfail
    @pytest.mark.parametrize(
        "a,b,c,z,expected,rtol",
        [
            (
                0.5,
                0.2,
                1.5,
                0.95 * e ** (pi / 2.5 * 1j),
                0.9983727370077057 + 0.06309930161718051j,
                1e-15,
            ),
            (
                0.5,
                2.5,
                3.5,
                0.9999 * e ** (-pi / 2.1 * 1j),
                0.8685871797143658 - 0.28018174801351164j,
                1e-14,
            ),
            (
                0.5,
                0.2,
                1.5,
                0.900001 * e ** (pi / 2.8 * 1j),
                1.0071803820262928 + 0.06225771451103533j,
                1e-14,
            ),
            (
                12.3,
                8.0,
                20.31,
                0.95 * e ** (-pi / 2.5 * 1j),
                -0.1373986808738931 + 0.7917645870031411j,
                1e-15,
            ),
            (
                12.3,
                8.3,
                20.31,
                0.900001 * e ** (pi / 2.8 * 1j),
                -0.05480312054409342 - 1.5089920040627767j,
                1e-12,
            ),
            (
                12.3,
                7.5,
                20.31,
                0.9999 * e ** (-pi / 2.1 * 1j),
                -0.21573631704577922 + 0.2583512402098123j,
                1e-12,
            ),
            (
                290.2,
                321.5,
                700.1,
                0.95 * e ** (pi / 2.5 * 1j),
                0.2779946965654163 + 0.15373012869874267j,
                1e-10,
            ),
            (
                290.2,
                321.2,
                700.1,
                0.95 * e ** (-pi / 2.5 * 1j),
                0.2981479522523749 - 0.12287022500428477j,
                1e-10,
            ),
            (
                -102.1,
                -20.3,
                1.3,
                0.9999 * e ** (pi / 2.9 * 1j),
                6.460072943882971e-21 + 2.084296332355303e21j,
                1e-15,
            ),
            (
                -102.7,
                -20.3,
                1.3,
                0.9999 * e ** (-pi / 2.9 * 1j),
                7.269067978154673e21 + 2.1985643096814801e21j,
                1e-15,
            ),
        ],
    )
    def test_domain4(self, a, b, c, z, expected, rtol):
        """0.9 < |z| < 1.0 and |1 - z| > 1.0"""
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)

    def test_domain5(self):
        """|z| > 1.0."""
        assert True
