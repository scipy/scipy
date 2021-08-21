"""Tests for hyp2f1 for complex values.

Author: Albert Steppi, with credit to Adam Kullberg (FormerPhycisist) for
the implementation of mp_hyp2f1 below, which modifies mpmath's hyp2f1 to
return the same branch as scipy's on the standard branch cut.
"""

import pytest
import numpy as np
from typing import NamedTuple
from numpy.testing import assert_allclose

from scipy.special import hyp2f1
from scipy.special._testutils import check_version, MissingModule


try:
    import mpmath
except ImportError:
    mpmath = MissingModule("mpmath")


def mp_hyp2f1(a, b, c, z):
    """Return mpmath hyp2f1 calculated on same branch as scipy hyp2f1.

    For most values of a,b,c mpmath returns the x - 0j branch of hyp2f1 on the
    branch cut x=(1,inf) whereas scipy's hyp2f1 calculates the x + 0j branch.
    Thus, to generate the right comparison values on the branch cut, we
    evaluate mpmath.hyp2f1 at x + 1e-15*j.

    The exception to this occurs when c-a=-m in which case both mpmath and
    scipy calculate the x + 0j branch on the branch cut. When this happens
    mpmath.hyp2f1 will be evaluated at the original z point.
    """
    on_branch_cut = z.real > 1.0 and abs(z.imag) < 1.0e-15
    cond1 = abs(c - a - round(c - a)) < 1.0e-15 and round(c - a) <= 0
    cond2 = abs(c - b - round(c - b)) < 1.0e-15 and round(c - b) <= 0
    # Make sure imaginary part is *exactly* zero
    if on_branch_cut:
        z = z.real + 0.0j
    if on_branch_cut and not (cond1 or cond2):
        z_mpmath = z.real + 1.0e-15j
    else:
        z_mpmath = z
    return complex(mpmath.hyp2f1(a, b, c, z_mpmath))


class Hyp2f1TestCase(NamedTuple):
    a: float
    b: float
    c: float
    z: complex
    expected: complex
    rtol: float


class TestHyp2f1:
    """Tests for hyp2f1 for complex values.

    Expected values for test cases were computed using mpmath. See
    `scipy.special._precompute.hyp2f1_data`. The verbose style of specifying
    test cases is used for readability and to make it easier to mark individual
    cases as expected to fail. Expected failures are used to highlight cases
    where improvements are needed. See
    `scipy.special._precompute.hyp2f1_data.make_hyp2f1_test_cases` for a
    function to generate the boilerplate for the test cases.

    Assertions have been added to each test to ensure that the test cases match
    the situations that are intended. A final test `test_test_hyp2f1` checks
    that the expected values in the test cases actually match what is computed
    by mpmath. This test is marked slow even though it isn't particularly slow
    so that it won't run by default on continuous integration builds.
    """
    @pytest.mark.parametrize(
        "hyp2f1_test_case",
        [
            pytest.param(
                Hyp2f1TestCase(
                    a=0.5,
                    b=0.2,
                    c=-10,
                    z=0.2 + 0.2j,
                    expected=np.inf + 0j,
                    rtol=0
                )
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=0.5,
                    b=0.2,
                    c=-10,
                    z=0 + 0j,
                    expected=1 + 0j,
                    rtol=0
                ),
                marks=pytest.mark.xfail(reason="gh-7340")
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=0.5,
                    b=0,
                    c=-10,
                    z=0.2 + 0.2j,
                    expected=1 + 0j,
                    rtol=0
                ),
                marks=pytest.mark.xfail(reason="gh-7340"),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=0.5,
                    b=0,
                    c=0,
                    z=0.2 + 0.2j,
                    expected=1 + 0j,
                    rtol=0,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=0.5,
                    b=0.2,
                    c=0,
                    z=0.2 + 0.2j,
                    expected=np.inf + 0j,
                    rtol=0,
                ),
                marks=pytest.mark.xfail(reason="gh-7340"),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=0.5,
                    b=0.2,
                    c=0,
                    z=0 + 0j,
                    expected=np.nan + 0j,
                    rtol=0,
                ),
                marks=pytest.mark.xfail(reason="gh-7340"),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=0.5,
                    b=-5,
                    c=-10,
                    z=0.2 + 0.2j,
                    expected=(1.0495404166666666+0.05708208333333334j),
                    rtol=1e-15,
                ),
                marks=pytest.mark.xfail(reason="gh-7340"),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=0.5,
                    b=-10,
                    c=-10,
                    z=0.2 + 0.2j,
                    expected=(1.092966013125+0.13455014673750001j),
                    rtol=1e-15,
                ),
                marks=pytest.mark.xfail(reason="gh-7340"),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=-10,
                    b=-20,
                    c=-10,
                    z=0.2 + 0.2j,
                    expected=(-0.07712512000000005+0.12752814080000005j),
                    rtol=1e-13,
                ),
                marks=pytest.mark.xfail(reason="gh-7340"),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=-1,
                    b=3.2,
                    c=-1,
                    z=0.2 + 0.2j,
                    expected=(1.6400000000000001+0.6400000000000001j),
                    rtol=1e-13,
                ),
                marks=pytest.mark.xfail(reason="gh-7340"),
            ),
        ]
    )
    def test_c_non_positive_int(self, hyp2f1_test_case):
        a, b, c, z, expected, rtol = hyp2f1_test_case
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)

    @pytest.mark.parametrize(
        "hyp2f1_test_case",
        [
            pytest.param(
                Hyp2f1TestCase(
                    a=0.5,
                    b=0.2,
                    c=1.5,
                    z=1 + 0j,
                    expected=1.1496439092239847 + 0j,
                    rtol=1e-15
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=12.3,
                    b=8.0,
                    c=20.31,
                    z=1 + 0j,
                    expected=69280986.75273195 + 0j,
                    rtol=1e-15
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=290.2,
                    b=321.5,
                    c=700.1,
                    z=1 + 0j,
                    expected=1.3396562400934e117 + 0j,
                    rtol=1e-12,
                ),
                marks=pytest.mark.xfail(reason="overflow"),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=-102.1,
                    b=-20.3,
                    c=1.3,
                    z=1 + 0j,
                    expected=2.7899070752746906e22 + 0j,
                    rtol=1e-15
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=-202.6,
                    b=60.3,
                    c=1.5,
                    z=1 + 0j,
                    expected=-1.3113641413099326e-56 + 0j,
                    rtol=1e-12,
                ),
                marks=pytest.mark.xfail(reason="underflow"),
            ),
        ],
    )
    def test_unital_argument(self, hyp2f1_test_case):
        """Tests for case z = 1, c - a - b > 0.

        Expected answers computed using mpmath.
        """
        a, b, c, z, expected, rtol = hyp2f1_test_case
        assert z == 1 and c - a - b > 0  # Tests the test
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)

    @pytest.mark.parametrize(
        "hyp2f1_test_case",
        [
            pytest.param(
                Hyp2f1TestCase(
                    a=0.5,
                    b=0.2,
                    c=1.3,
                    z=-1 + 0j,
                    expected=0.9428846409614143 + 0j,
                    rtol=1e-15),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=12.3,
                    b=8.0,
                    c=5.300000000000001,
                    z=-1 + 0j,
                    expected=-4.845809986595704e-06 + 0j,
                    rtol=1e-15
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=221.5,
                    b=90.2,
                    c=132.3,
                    z=-1 + 0j,
                    expected=2.0490488728377282e-42 + 0j,
                    rtol=1e-7,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=-102.1,
                    b=-20.3,
                    c=-80.8,
                    z=-1 + 0j,
                    expected=45143784.46783885 + 0j,
                    rtol=1e-7,
                ),
                marks=pytest.mark.xfail,
            ),
        ],
    )
    def test_special_case_z_near_minus_1(self, hyp2f1_test_case):
        """Tests for case z ~ -1, c ~ 1 + a - b

        Expected answers computed using mpmath.
        """
        a, b, c, z, expected, rtol = hyp2f1_test_case
        assert abs(1 + a - b - c) < 1e-15 and abs(z + 1) < 1e-15
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)

    @pytest.mark.parametrize(
        "hyp2f1_test_case",
        [
            pytest.param(
                Hyp2f1TestCase(
                    a=-4,
                    b=2.02764642551431,
                    c=1.0561196186065624,
                    z=(0.9473684210526314-0.10526315789473695j),
                    expected=(0.0031961077109535375-0.0011313924606557173j),
                    rtol=1e-12,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=-8,
                    b=-7.937789122896016,
                    c=-15.964218273004214,
                    z=(2-0.10526315789473695j),
                    expected=(0.005543763196412503-0.0025948879065698306j),
                    rtol=5e-13,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=-8,
                    b=8.095813935368371,
                    c=4.0013768449590685,
                    z=(0.9473684210526314-0.10526315789473695j),
                    expected=(-0.0003054674127221263-9.261359291755414e-05j),
                    rtol=9e-11,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=-4,
                    b=-3.956227226099288,
                    c=-3.9316537064827854,
                    z=(1.1578947368421053-0.3157894736842106j),
                    expected=(-0.0020809502580892937-0.0041877333232365095j),
                    rtol=2e-13,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=2.02764642551431,
                    b=-4,
                    c=2.050308316530781,
                    z=(0.9473684210526314-0.10526315789473695j),
                    expected=(0.0011282435590058734+0.0002027062303465851j),
                    rtol=5e-13,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=-7.937789122896016,
                    b=-8,
                    c=-15.964218273004214,
                    z=(1.3684210526315788+0.10526315789473673j),
                    expected=(-9.134907719238265e-05-0.00040219233987390723j),
                    rtol=5e-12,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=4.080187217753502,
                    b=-4,
                    c=4.0013768449590685,
                    z=(0.9473684210526314-0.10526315789473695j),
                    expected=(-0.000519013062087489-0.0005855883076830948j),
                    rtol=5e-12,
                ),
            ),
        ]
    )
    def test_a_b_negative_int(self, hyp2f1_test_case):
        a, b, c, z, expected, rtol = hyp2f1_test_case
        assert a == int(a) and a < 0 or b == int(b) and b < 0  # Tests the test
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)

    @pytest.mark.parametrize(
        "hyp2f1_test_case",
        [
            pytest.param(
                Hyp2f1TestCase(
                    a=-0.5,
                    b=-0.9629749245209605,
                    c=-15.5,
                    z=(1.1578947368421053-1.1578947368421053j),
                    expected=(0.9778506962676361+0.044083801141231616j),
                    rtol=1e-12,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=8.5,
                    b=-3.9316537064827854,
                    c=1.5,
                    z=(0.9473684210526314-0.10526315789473695j),
                    expected=(4.0793167523167675-10.11694246310966j),
                    rtol=6e-12,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=8.5,
                    b=-0.9629749245209605,
                    c=2.5,
                    z=(1.1578947368421053-0.10526315789473695j),
                    expected=(-2.9692999501916915+0.6394599899845594j),
                    rtol=1e-11,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=-0.5,
                    b=-0.9629749245209605,
                    c=-15.5,
                    z=(1.5789473684210522-1.1578947368421053j),
                    expected=(0.9493076367106102-0.04316852977183447j),
                    rtol=1e-11,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=-0.9220024191881196,
                    b=-0.5,
                    c=-15.5,
                    z=(0.5263157894736841+0.10526315789473673j),
                    expected=(0.9844377175631795-0.003120587561483841j),
                    rtol=1e-10,
                ),
            ),
        ],
    )
    def test_a_b_neg_int_after_euler_hypergeometric_transformation(
        self, hyp2f1_test_case
    ):
        a, b, c, z, expected, rtol = hyp2f1_test_case
        assert (  # Tests the test
            (abs(c - a - int(c - a)) < 1e-15 and c - a < 0) or
            (abs(c - b - int(c - b)) < 1e-15 and c - b < 0)
        )
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)

    @pytest.mark.parametrize(
        "hyp2f1_test_case",
        [
            pytest.param(
                Hyp2f1TestCase(
                    a=-0.9220024191881196,
                    b=-0.9629749245209605,
                    c=-15.963511401609862,
                    z=(0.10526315789473673-0.3157894736842106j),
                    expected=(0.9941449585778349+0.01756335047931358j),
                    rtol=1e-14,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=1.0272592605282642,
                    b=-0.9629749245209605,
                    c=-15.963511401609862,
                    z=(0.5263157894736841+0.5263157894736841j),
                    expected=(1.0388722293372104-0.09549450380041416j),
                    rtol=5e-11,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=2.02764642551431,
                    b=1.0561196186065624,
                    c=-7.93846038215665,
                    z=(0.10526315789473673+0.7368421052631575j),
                    expected=(2.1948378809826434+24.934157235172222j),
                    rtol=5e-15,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=2.02764642551431,
                    b=16.088264119063613,
                    c=8.031683612216888,
                    z=(0.3157894736842106-0.736842105263158j),
                    expected=(-0.4075277891264672-0.06819344579666956j),
                    rtol=2e-12,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=4.080187217753502,
                    b=2.050308316530781,
                    c=8.031683612216888,
                    z=(0.7368421052631575-0.10526315789473695j),
                    expected=(2.833535530740603-0.6925373701408158j),
                    rtol=5e-15,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=2.02764642551431,
                    b=2.050308316530781,
                    c=4.078873014294075,
                    z=(0.10526315789473673-0.3157894736842106j),
                    expected=(1.005347176329683-0.3580736009337313j),
                    rtol=5e-16,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=-0.9220024191881196,
                    b=-0.9629749245209605,
                    c=-15.963511401609862,
                    z=(0.3157894736842106-0.5263157894736843j),
                    expected=(0.9824353641135369+0.029271018868990268j),
                    rtol=5e-13,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=-0.9220024191881196,
                    b=-0.9629749245209605,
                    c=-159.63511401609862,
                    z=(0.3157894736842106-0.5263157894736843j),
                    expected=(0.9982436200365834+0.002927268199671111j),
                    rtol=1e-7,
                ),
                marks=pytest.mark.xfail(reason="Poor convergence.")
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=2.02764642551431,
                    b=16.088264119063613,
                    c=8.031683612216888,
                    z=(0.5263157894736841-0.5263157894736843j),
                    expected=(-0.6906825165778091+0.8176575137504892j),
                    rtol=5e-13,
                ),
            ),
        ]
    )
    def test_region1(self, hyp2f1_test_case):
        """|z| < 0.9 and real(z) >= 0."""
        a, b, c, z, expected, rtol = hyp2f1_test_case
        assert abs(z) < 0.9 and z.real >= 0  # Tests the test
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)

    @pytest.mark.parametrize(
        "hyp2f1_test_case",
        [
            pytest.param(
                Hyp2f1TestCase(
                    a=2.02764642551431,
                    b=1.0561196186065624,
                    c=4.078873014294075,
                    z=(-0.3157894736842106+0.7368421052631575j),
                    expected=(0.7751915029081136+0.24068493258607315j),
                    rtol=5e-15,
                ),
            ),
        ]
    )
    def test_region2(self, hyp2f1_test_case):
        """|z| < 1 and real(z) < 0."""
        a, b, c, z, expected, rtol = hyp2f1_test_case
        assert abs(z) < 1 and z.real < 0  # Tests the test
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)

    @pytest.mark.xfail(reason="gh-8054")
    @pytest.mark.parametrize(
        "hyp2f1_test_case",
        [
            pytest.param(
                Hyp2f1TestCase(
                    a=1.0272592605282642,
                    b=1.0561196186065624,
                    c=-0.906685989801748,
                    z=(0.10526315789473673-0.9473684210526316j),
                    expected=(-3.9995506969395858-8.179533155338005j),
                    rtol=1e-12,
                ),
            ),
            pytest.param(
                Hyp2f1TestCase(
                    a=8.5,
                    b=4.5,
                    c=8.077282662161238,
                    z=(0.3157894736842106+0.9473684210526314j),
                    expected=(-0.11307039404123598-0.443195310438102j),
                    rtol=1e-12,
                ),
            ),
        ]
    )
    def test_region4(self, hyp2f1_test_case):
        """0.9 <= |z| <= 1 and |1 - z| >= 1.

        This region is unhandled by of the standard transformations and
        needs special care.
        """
        a, b, c, z, expected, rtol = hyp2f1_test_case
        assert 0.9 <= abs(z) <= 1 and abs(1 - z) >= 1  # Tests the test
        assert_allclose(hyp2f1(a, b, c, z), expected, rtol=rtol)

    # Marked as slow so it won't run by default. This test is not slow.
    # Including it only increases the running time of the entire suite by
    # a handful of hundreths of seconds. This test could become slow in the
    # future if enough test cases are added.
    @pytest.mark.slow
    @check_version(mpmath, "1.0.0")
    def test_test_hyp2f1(self):
        """Test that expected values match what is computed by mpmath.

        This gathers the parameters for the test cases out of the pytest marks.
        The parameters are a, b, c, z, expected, rtol, where expected should
        be the value of hyp2f1(a, b, c, z) computed with mpmath. The test
        recomputes hyp2f1(a, b, c, z) using mpmath and verifies that expected
        actually is the correct value. This allows the data for the tests to
        live within the test code instead of an external datafile, while
        avoiding having to compute the results with mpmath during the test,
        except for when slow tests are being run.
        """
        test_methods = [
            test_method for test_method in dir(self)
            if test_method.startswith('test') and
            # Filter properties and attributes (futureproofing).
            callable(getattr(self, test_method)) and
            # Filter out this test
            test_method != 'test_test_hyp2f1'
        ]
        for test_method in test_methods:
            params = self._get_test_parameters(getattr(self, test_method))
            for a, b, c, z, expected, _ in params:
                assert_allclose(mp_hyp2f1(a, b, c, z), expected, rtol=2.25e-16)

    def _get_test_parameters(self, test_method):
        """Get pytest.mark parameters for a test in this class."""
        return [
            case.values[0] for mark in test_method.pytestmark
            if mark.name == 'parametrize'
            for case in mark.args[1]
        ]
