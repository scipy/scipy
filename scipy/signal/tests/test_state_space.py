"""Unit tests for `ltisys.py`. """
import numpy as np
import pytest

from scipy._lib._array_api import make_xp_test_case, xp_assert_close, xp_assert_equal
from scipy.signal import freqs
from scipy.signal._state_space import (_calcSSTF_helper, _geomspace, _parse_om_spec,
                                       freqs_ss)


def _diag(x, xp):
    """Portable version of `numpy.diag` that works with all array API backends. """
    return xp.asarray(np.diag(np.asarray(x)))


@make_xp_test_case(_calcSSTF_helper)
class TestCalcSSTFHelper:
    """Tests for function `_calcSSTF_helper`.  """

    @pytest.mark.parametrize('B, C', [
        (np.eye(3), np.eye(3)), (np.ones((3, 2)), np.ones((1, 3)))])
    def test_basic(self, B, C, xp):
        """Basic verification.

        Verifying the results by hand is trivial when B, C are identity matrices.
        """
        B, C = xp.asarray(B), xp.asarray(C)
        s, T = xp.asarray([0., 2., 4., 8.]), _diag([2., 4., 8.], xp)
        G = _calcSSTF_helper(T, B, C, s, rtol=1e-12, atol=0, overwrite_b=False,
                             xp=xp, device=None)

        inf = 1e64  # placeholder for ∞ to propagate through C @ M @ B correctly
        M = xp.asarray(([[-1/2, 0., 0.], [0., -1/4, 0.], [0., 0., -1/8]],
                        [[ inf, 0., 0.], [0., -1/2, 0.], [0., 0., -1/6]],
                        [[ 1/2, 0., 0.], [0.,  inf, 0.], [0., 0., -1/4]],
                        [[ 1/6, 0., 0.], [0.,  1/4, 0.], [0., 0.,  inf]]))
        G_ref = C @ M @ B
        G_ref[abs(G_ref) > inf/32] = xp.inf

        xp_assert_close(G, G_ref)

    @pytest.mark.parametrize('rtol, atol, v', [
        (0., 0., 0), (1., 0., np.inf), (0., 0.5, np.inf)])
    def test_tol(self, rtol, atol, v, xp):
        """Verify that the `tol` parameter works. """
        T, B, C = _diag([1, 0.8], xp), xp.eye(2), xp.eye(2)
        s = xp.asarray([0.8, 1, 1.2])
        G = _calcSSTF_helper(T, B, C, s, rtol=rtol, atol=atol, overwrite_b=False,
                             xp=xp, device=None)
        G_ref = xp.asarray([[[-5.0+v, 0.0], [0.0, xp.inf]],
                            [[xp.inf, 0.0], [0.0,  5.0+v]],
                            [[ 5.0+v, 0.0], [0.0,  2.5+v]]])
        xp_assert_close(G, G_ref)


@make_xp_test_case(_parse_om_spec)
class TestParseOmSpec:
    """Tests for function `_parse_om_spec`. """

    def test_invalid_f_scale(self, xp):
        """Raise exceptions in `_parse_om_spec`. """
        with pytest.raises(ValueError, match="^Parameter f_scale='INVALID'"):
            # noinspection bad-argument-type
            _parse_om_spec(10, 'INVALID', 0.1, 1., xp=xp, device=None)

    @pytest.mark.parametrize('om_spec', [
        (0, 1, 'INVALID'), (0, 2, -5), (0, 1j, 5), (3, 2, 5)])
    def test_invalid_om_spec_tuple(self, om_spec, xp):
        """Test for invalid `om_spec` tuple values. """
        with pytest.raises(ValueError, match="^Parameter om_spec="):
            _parse_om_spec(om_spec, 'linear', 0.1, 1., xp=xp, device=None)

    @pytest.mark.parametrize('om_spec', ['', [[1, 2]], [], [3+4j]])
    def test_invalid_om_spec_array(self, om_spec, xp):
        """Test for invalid `om_spec` values, when not a tuple. """
        with pytest.raises(ValueError, match="^Parameter `om_spec` must either be"):
            _parse_om_spec(om_spec, 'linear', 0.1, 1, xp=xp, device=None)

    @pytest.mark.parametrize('f_scale', ['linear', 'log'])
    def test_basic_int(self, f_scale, xp):
        """Basic test for `om_spec` being an integer. """
        N, om_min, om_max = 10, 2, 16

        om = _parse_om_spec(N, f_scale, om_min, om_max, xp=xp, device=None)

        om0_ref = 0 if f_scale == 'linear' else om_min/5
        om1_ref = 1.5*om_max if f_scale == 'linear' else 5*om_max
        om_ref =(xp.linspace(om0_ref, om1_ref, N) if f_scale == 'linear' else
                 _geomspace(om0_ref, om1_ref, N, xp=xp))
        xp_assert_close(om, om_ref)

    @pytest.mark.parametrize('f_scale', ['linear', 'log'])
    def test_basic_tuple(self, f_scale, xp):
        """Basic test for `om_spec` being a tuple. """
        om0, om1, N = 0.1, 11, 10

        om = _parse_om_spec((om0, om1, N), f_scale, 0.1, 1, xp=xp, device=None)
        om_ref = (xp.linspace(om0, om1, N) if f_scale == 'linear' else
                  _geomspace(om0, om1, N, xp=xp))
        xp_assert_close(om, om_ref)

    def test_basic_array(self, xp):
        """Basic test for `om_spec` being an array. """
        om_ref = xp.arange(5, 10)
        om = _parse_om_spec(om_ref, 'linear', 0.1, 1., xp=xp, device=None)
        xp_assert_close(om, om_ref)


@make_xp_test_case(freqs_ss)
class TestFreqsSS:
    """Tests for function `freqs_ss`.

    Since `_calcSSTF_helper` and `_parse_om_spec` have their own tests, `freq_ss` needs
    only rudimentary tests.
    """

    @pytest.mark.parametrize('B, C', [
        (np.zeros((3, 2)), np.ones((2, 3))),
        (np.ones((3, 2)), np.zeros((2, 3))),])
    def test_constant_tf(self, B, C, xp):
        """Test for C or B being zero, which produces a constant transfer function. """
        A, B = xp.eye(3), xp.asarray(B)
        C, D = xp.asarray(C), xp.asarray([[1., 2.], [3., 4.]])
        G_ref = xp.moveaxis(xp.tile(D, (5, 1, 1)), 0, -1)
        om, G = freqs_ss(A, B, C, D, om_spec=(0, 10, 5))
        xp_assert_close(om, xp.linspace(0, 10, 5))
        xp_assert_close(G, G_ref)

    def test_transfer_functions(self, xp):
        """Compare with transfer functions.

         Consult the docstring of `freqs_ss` for details on the comparison of the
         transfer functions and the state space system.
         """
        om0, om1, d0, d1, gamma  = 5, 15, 0.1, 0.2, 1.
        AA = xp.asarray([[           0,         1,       0,         0],
                         [     -om0**2, -2*d0*om0,       0,         0],
                         [           0,         0,       0,         1],
                         [gamma*om1**2,         0, -om1**2, -2*d1*om1]])
        BB = xp.asarray([[0, om0**2, 0, 0], [0, 0, 0, om1**2]]).T
        CC, DD = xp.asarray([[1, 0, 0, 0], [0, 0, 1, 0]]), xp.zeros((2, 2))

        om, b_ = xp.arange(20, dtype=AA.dtype), xp.asarray([1.,])
        _, H0_ref = freqs(b_, xp.asarray([1/om0**2, 2*d0/om0, 1.]), worN=om)
        _, H1_ref = freqs(b_, xp.asarray([1/om1**2, 2*d1/om1, 1.]), worN=om)
        om_GG, GG = freqs_ss(AA, BB, CC, DD, om_spec=om)  # the response

        xp_assert_equal(om, om_GG)
        xp_assert_close(GG[0, 0, :], H0_ref)
        xp_assert_close(GG[0, 1, :], xp.zeros_like(H0_ref), atol=1e-12)
        xp_assert_close(GG[1, 1, :], H1_ref)
        xp_assert_close(GG[1, 0, :], H1_ref * gamma * H0_ref)

    @pytest.mark.parametrize('N', [1, 3])
    @pytest.mark.parametrize('p', [0, 1, 2])
    @pytest.mark.parametrize('q', [0, 1, 2])
    @pytest.mark.parametrize('n', [0, 1, 2])
    def test_degenerate_system(self, n, q, p, N, xp):
        """Test degenerate systems with empty matrices.. """
        ABCD = (xp.eye(n), xp.ones((n, p)), xp.ones((q, n)), xp.ones((q, p)))
        om, GG = freqs_ss(*ABCD, om_spec=(0, 1, N))
        xp_assert_close(GG.shape, (q, p, N))

    @pytest.mark.parametrize('f_scale, om_min, om_max', [
        ('linear', 0., 30.), ('log', 0.4, 100.)])
    def test_automatic_frequency_range(self, f_scale, om_min, om_max, xp):
        ABCD = _diag([2j, 1, -20j], xp), xp.eye(3), xp.eye(3), xp.eye(3)
        om, _ = freqs_ss(*ABCD, om_spec=2, f_scale=f_scale)
        xp_assert_close(om, xp.asarray((om_min, om_max)))

    @pytest.mark.parametrize('f_scale, om_min, om_max', [
        ('linear', 0., 1.5), ('log', 0.2, 5.)])
    def test_automatic_frequency_range_no_imag(self, f_scale, om_min, om_max, xp):
        """Test for A's eigenvalues being real. """
        ABCD = xp.eye(3), xp.eye(3), xp.eye(3), xp.eye(3)
        om, _ = freqs_ss(*ABCD, om_spec=2, f_scale=f_scale)
        xp_assert_close(om, (om_min, om_max))

    @pytest.mark.parametrize('f_scale, om_min, om_max', [
        ('linear', 0., 1.5), ('log', 0.2, 5.)])
    def test_automatic_frequency_range_one_ev(self, f_scale, om_min, om_max, xp):
        """Test for A's only having one eigenvalue. """
        ABCD = [xp.ones((1, 1)) for _ in range(4)]
        om, _ = freqs_ss(*ABCD, om_spec=2, f_scale=f_scale)
        xp_assert_close(om, xp.asarray((om_min, om_max)))

    @pytest.mark.parametrize('f_scale, om_min, om_max', [
        ('linear', 0., 1.5), ('log', 0.2, 5.)])
    def test_automatic_frequency_range_no_A(self, f_scale, om_min, om_max, xp):
        """Test for A's eigenvalues  not producing a range."""
        A, B, C = xp.zeros((0, 0)), xp.zeros((0, 1)), xp.zeros((2, 0))
        D = xp.ones((2, 1))

        om, _ = freqs_ss(A, B, C, D, om_spec=2, f_scale=f_scale)
        xp_assert_close(om, xp.asarray((om_min, om_max)))

