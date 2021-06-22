import pytest
import numpy as np
from numpy.testing import assert_allclose

from scipy.special import hyp2f1


class TestHyp2f1:
    """Tests for hyp2f1 for complex values.

    Expected values for test cases computed using mpmath. See
    scipy.special._precompute.hyp2f1_data. The pytest.mark.parametrize
    style of specifying test cases is used instead of FuncData from
    scipy.special._testutils to make it easier to mark individual cases as
    expected to fail. Expected failures are used to highlight cases where
    improvements are needed.
    """
    @pytest.mark.parametrize(
        "a,b,c,z,expected,rtol",
        [
            (0.5, 0.2, -10, 0.2 + 0.2j, np.inf + 0j, 0),
            pytest.param(
                0.5,
                0.2,
                -10,
                0 + 0j,
                1 + 0j,
                0,
                marks=pytest.mark.xfail(reason="gh-7340")
            ),
            pytest.param(
                0.5,
                0,
                -10,
                0.2 + 0.2j,
                1 + 0j,
                0,
                marks=pytest.mark.xfail(reason="gh-7340"),
            ),
            (
                0.5,
                0,
                0,
                0.2 + 0.2j,
                1 + 0j,
                0,
            ),
            pytest.param(
                0.5,
                0.2,
                0,
                0.2 + 0.2j,
                np.inf + 0j,
                0,
                marks=pytest.mark.xfail(reason="gh-7340"),
            ),
            pytest.param(
                0.5,
                0.2,
                0,
                0 + 0j,
                np.nan + 0j,
                0,
                marks=pytest.mark.xfail(reason="gh-7340"),
            ),
            pytest.param(
                0.5,
                -5,
                -10,
                0.2 + 0.2j,
                (1.0495404166666666+0.05708208333333334j),
                1e-15,
                marks=pytest.mark.xfail(reason="gh-7340"),
            ),
            pytest.param(
                0.5,
                -10,
                -10,
                0.2 + 0.2j,
                (1.092966013125+0.13455014673750001j),
                1e-15,
                marks=pytest.mark.xfail(reason="gh-7340"),
            ),
            pytest.param(
                -10,
                -20,
                -10,
                0.2 + 0.2j,
                (-0.07712512000000005+0.12752814080000005j),
                1e-13,
                marks=pytest.mark.xfail(reason="gh-7340"),
            ),
            pytest.param(
                -1,
                3.2,
                -1,
                0.2 + 0.2j,
                (1.6400000000000001+0.6400000000000001j),
                1e-13,
                marks=pytest.mark.xfail(reason="gh-7340"),
            ),
        ],
    )
    def test_c_non_positive_int(self, a, b, c, z, expected, rtol):
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)

    @pytest.mark.parametrize(
        "a,b,c,expected,rtol",
        [
            (0.5, 0.2, 1.5, 1.1496439092239847 + 0j, 1e-15),
            (12.3, 8.0, 20.31, 69280986.75273195 + 0j, 1e-15),
            pytest.param(
                290.2,
                321.5,
                700.1,
                1.3396562400934e117 + 0j,
                1e-12,
                marks=pytest.mark.xfail(reason="overflow"),
            ),
            (-102.1, -20.3, 1.3, 2.7899070752746906e22 + 0j, 1e-15),
            pytest.param(
                -202.6,
                60.3,
                1.5,
                -1.3113641413099326e-56 + 0j,
                1e-12,
                marks=pytest.mark.xfail(reason="underflow"),
            ),
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
            pytest.param(
                -102.1,
                -20.3,
                45143784.46783885 + 0j,
                1e-7,
                marks=pytest.mark.xfail,
            ),
        ],
    )
    def test_special_case_z_near_minus_1(self, a, b, expected, rtol):
        """Tests for case z ~ -1, c ~ 1 + a - b

        Expected answers computed using mpmath.
        """
        assert_allclose(
            [
                hyp2f1(a, b, 1 + a - b, -1 + 0j),
                hyp2f1(a, b, 1 + a - b, -1 + 2e-16j),
                hyp2f1(a, b, 1 + a - b, -1 - 2e-16j),
            ],
            [expected, expected, expected],
            rtol=rtol,
        )

    @pytest.mark.parametrize(
        "a,b,c,z,expected,rtol",
        [
            (
                -4,
                2.02764642551431,
                1.0561196186065624,
                (0.9473684210526314-0.10526315789473695j),
                (0.0031961077109535375-0.0011313924606557173j),
                1e-12,
            ),
            (
                -8,
                -7.937789122896016,
                -15.964218273004214,
                (2-0.10526315789473695j),
                (0.005543763196412503-0.0025948879065698306j),
                5e-13,
            ),
            (
                -8,
                8.095813935368371,
                4.0013768449590685,
                (0.9473684210526314-0.10526315789473695j),
                (-0.0003054674127221263-9.261359291755414e-05j),
                5e-11,
            ),
            (
                -4,
                -3.956227226099288,
                -3.9316537064827854,
                (1.1578947368421053-0.3157894736842106j),
                (-0.0020809502580892937-0.0041877333232365095j),
                1e-13,
            ),
            (
                2.02764642551431,
                -4,
                2.050308316530781,
                (0.9473684210526314-0.10526315789473695j),
                (0.0011282435590058734+0.0002027062303465851j),
                1e-13,
            ),
            (
                -7.937789122896016,
                -8,
                -15.964218273004214,
                (1.3684210526315788+0.10526315789473673j),
                (-9.134907719238265e-05-0.00040219233987390723j),
                5e-12,
            ),
            (
                4.080187217753502,
                -4,
                4.0013768449590685,
                (0.9473684210526314-0.10526315789473695j),
                (-0.000519013062087489-0.0005855883076830948j),
                5e-12,
            ),
        ]
    )
    def test_a_b_negative_int(self, a, b, c, z, expected, rtol):
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)

    @pytest.mark.parametrize(
        "a,b,c,z,expected,rtol",
        [
            (
                -0.5,
                -0.9629749245209605,
                -15.5,
                (1.1578947368421053-1.1578947368421053j),
                (0.9778506962676361+0.044083801141231616j),
                1e-12,
            ),
            (
                8.5,
                -3.9316537064827854,
                1.5,
                (0.9473684210526314-0.10526315789473695j),
                (4.0793167523167675-10.11694246310966j),
                5e-12,
            ),
            (
                8.5,
                -0.9629749245209605,
                2.5,
                (1.1578947368421053-0.10526315789473695j),
                (-2.9692999501916915+0.6394599899845594j),
                5e-12,
            ),
            (
                -0.5,
                -0.9629749245209605,
                -15.5,
                (1.5789473684210522-1.1578947368421053j),
                (0.9493076367106102-0.04316852977183447j),
                5e-12,
            ),
            (
                -0.9220024191881196,
                -0.5,
                -15.5,
                (0.5263157894736841+0.10526315789473673j),
                (0.9844377175631795-0.003120587561483841j),
                1e-10,
            ),
            (
                -15.980848054962111,
                4.5,
                8.5,
                (0.5263157894736841-0.5263157894736843j),
                (-0.008821385280750719+0.0026165583173886836j),
                1e-13,
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
                -0.9220024191881196,
                -0.9629749245209605,
                -15.963511401609862,
                (0.10526315789473673-0.3157894736842106j),
                (0.9941449585778349+0.01756335047931358j),
                1e-14,
            ),
            (
                1.0272592605282642,
                -0.9629749245209605,
                -15.963511401609862,
                (0.5263157894736841+0.5263157894736841j),
                (1.0388722293372104-0.09549450380041416j),
                1e-11,
            ),
            (
                2.02764642551431,
                1.0561196186065624,
                -7.93846038215665,
                (0.10526315789473673+0.7368421052631575j),
                (2.1948378809826434+24.934157235172222j),
                5e-15,
            ),
            (
                2.02764642551431,
                16.088264119063613,
                8.031683612216888,
                (0.3157894736842106-0.736842105263158j),
                (-0.4075277891264672-0.06819344579666956j),
                1e-12,
            ),
            (
                4.080187217753502,
                2.050308316530781,
                8.031683612216888,
                (0.7368421052631575-0.10526315789473695j),
                (2.833535530740603-0.6925373701408158j),
                5e-15,
            ),
            (
                2.02764642551431,
                2.050308316530781,
                4.078873014294075,
                (0.10526315789473673-0.3157894736842106j),
                (1.005347176329683-0.3580736009337313j),
                5e-16,
            ),
            (
                -0.9220024191881196,
                -0.9629749245209605,
                -15.963511401609862,
                (0.3157894736842106-0.5263157894736843j),
                (0.9824353641135369+0.029271018868990268j),
                5e-13,
            ),
            pytest.param(
                -0.9220024191881196,
                -0.9629749245209605,
                -159.63511401609862,
                (0.3157894736842106-0.5263157894736843j),
                (0.9982436200365834+0.002927268199671111j),
                1e-7,
                marks=pytest.mark.xfail(reason="Poor convergence.")
            ),
            (
                2.02764642551431,
                16.088264119063613,
                8.031683612216888,
                (0.5263157894736841-0.5263157894736843j),
                (-0.6906825165778091+0.8176575137504892j),
                5e-13,
            ),
        ]
    )
    def test_region1(self, a, b, c, z, expected, rtol):
        """|z| < 0.9 and real(z) > 0."""
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)

    @pytest.mark.parametrize(
        "a,b,c,z,expected,rtol",
        [
            (
                2.02764642551431,
                1.0561196186065624,
                4.078873014294075,
                (-0.3157894736842106+0.7368421052631575j),
                (0.7751915029081136+0.24068493258607315j),
                1e-15,
            ),
        ]
    )
    def test_region2(self, a, b, c, z, expected, rtol):
        """|z| < 1 and real(z) < 0."""
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)

    @pytest.mark.xfail(reason="gh-8054")
    @pytest.mark.parametrize(
        "a,b,c,z,expected,rtol",
        [
            (
                1.0272592605282642,
                1.0561196186065624,
                -0.906685989801748,
                (0.10526315789473673-0.9473684210526316j),
                (-3.9995506969395778-8.179533155337996j),
                1e-12,
            ),
            (
                8.5,
                4.5,
                8.077282662161238,
                (0.3157894736842106+0.9473684210526314j),
                (-0.11307039404123598-0.443195310438102j),
                1e-12,
            ),
        ]
    )
    def test_region4(self, a, b, c, z, expected, rtol):
        """0.9 <= |z| <= 1 and |1 - z| >= 1.

        This region is unhandled by of the standard transformations and
        needs special care.
        """
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)
