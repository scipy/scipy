"""Unit tests for file `_state_space.py`. """
from typing import Literal

import numpy as np
import pytest

from scipy import signal
from scipy._lib._array_api import Array, array_namespace, make_xp_test_case, \
    xp_assert_close, xp_assert_equal
from scipy.signal._state_space import _calcSSTF_helper, _diagsvd, _svd_rank, \
    _tran_zeros_reduce, _transmission_zeros, conjT


def _abcd_normalize(A: Array|None = None, B: Array|None = None, C: Array|None = None,
                    D: Array|None = None, *, xp=None) \
    -> tuple[Array, Array, Array, Array]:
    """Wrapper for `abcd_normalize()` which allows forcing an Array API namespace `xp`.

    This helper adds an optional `xp` parameter to `signal.abcd_normalize()`. This
    ensures that the following expression produces the correct arrays in all
    namespaces::

        A, B, C, D = _abcd_normalize(1, 2, 3, 4, xp=xp)

    This wrapper make the following unit tests more concise.
    """
    xp = xp or array_namespace(A, B, C, D)
    A_ = xp.asarray(A) if A is not None else None
    B_ = xp.asarray(B) if B is not None else None
    C_ = xp.asarray(C) if C is not None else None
    D_ = xp.asarray(D) if D is not None else None
    return signal.abcd_normalize(A_, B_, C_, D_)


@make_xp_test_case(conjT)
@pytest.mark.parametrize('dim', [(2, 3), (4, 3, 2), (0, 0)])
def test_conjT(dim, xp):
    """Test `conjT()`. """
    x = xp.reshape(xp.arange(int(np.prod(dim)), dtype=xp.complex128), dim) * (2.+4j)
    y, y_ref = conjT(x), xp.conj(x).mT
    assert xp.all(y == y_ref)


@make_xp_test_case(_svd_rank)
class TestSvdRank:
    """Tests for `_svd_rank()`. """

    def test_exceptions(self, xp):
        """Test exception handling in `_svd_rank()`. """
        A = xp.eye(3)
        with pytest.raises(ValueError, match="^Parameter m_size='INVALID' not in "):
            _svd_rank(A, 'INVALID', rtol=1e-12, atol=0.)

    @pytest.mark.parametrize('m_size', ['full', 'shape', 'rank'])
    @pytest.mark.parametrize('atol, rtol', [(1e-12, 0.), (1e-8, 0.), (0., 1e-8)])
    def test_basic(self, atol, rtol, m_size: Literal['full', 'shape', 'rank'], xp):
        """Test `_svd_rank()` function. """
        s_in = np.asarray([2., 1e-8, 1.])
        A = np.zeros((3, 4), dtype=s_in.dtype)
        for i_, v_ in enumerate(s_in):
            A[i_, i_] = v_

        s_ref = np.array([v_ if v_ > rtol*np.max(s_in) + atol else 0. for v_ in s_in])
        rank_ref = int(np.sum(s_ref > 0.))
        A_ref = _diagsvd(s_ref, *A.shape, xp=xp)

        U, s, Vh, rank = _svd_rank(xp.asarray(A), m_size, rtol=rtol, atol=atol)
        S = _diagsvd(s, U.shape[1], Vh.shape[0])
        xp_assert_close(U @ S @ Vh, A_ref)
        assert rank == rank_ref

    @pytest.mark.parametrize('m_size', ['full', 'shape', 'rank'])
    def test_empty(self, m_size: Literal['full', 'shape', 'rank'], xp):
        """Test `_svd_rank()` with an empty matrix. """
        A = xp.zeros((0, 0))
        U, s, Vh, rank = _svd_rank(A, m_size, rtol=1e-12, atol=0.)
        assert U.shape == (0, 0)
        assert s.shape == (0,)
        assert Vh.shape == (0, 0)
        assert rank == 0


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
        s, T = xp.asarray([0., 2., 4., 8.]), xp.asarray(np.diag([2., 4., 8.+0j]))
        G = _calcSSTF_helper(T, B, C, s, rtol=1e-12, atol=0, overwrite_b=False,
                             xp=xp, device=None)

        inf = 1e64  # placeholder for ∞ to propagate through C @ M @ B correctly
        M = xp.asarray(([[-1/2, 0., 0.], [0., -1/4, 0.], [0., 0., -1/8]],
                        [[ inf, 0., 0.], [0., -1/2, 0.], [0., 0., -1/6]],
                        [[ 1/2, 0., 0.], [0.,  inf, 0.], [0., 0., -1/4]],
                        [[ 1/6, 0., 0.], [0.,  1/4, 0.], [0., 0.,  inf]]))
        G_ref = xp.astype(C @ M @ B, T.dtype)
        G_ref[abs(G_ref) > inf/32] = xp.inf * (1+1j)

        xp_assert_close(G, G_ref)

    @pytest.mark.parametrize('rtol, atol, v', [
        (0., 0., 0), (1., 0., float('inf')), (0., 0.5, float('inf'))])
    def test_tol(self, rtol, atol, v, xp):
        """Verify that tolerance parameters (`rtol`, `atol`) are applied correctly. """
        T, B, C = xp.asarray(np.diag([1, 0.8+0j])), xp.eye(2), xp.eye(2)
        s = xp.asarray([0.8, 1, 1.2])
        G = _calcSSTF_helper(T, B, C, s, rtol=rtol, atol=atol, overwrite_b=False,
                             xp=xp, device=None)

        G_ref = xp.asarray([[[ -5.0+v, 0.0], [0.0, xp.inf]],
                            [[ xp.inf, 0.0], [0.0,  5.0+v]],
                            [[  5.0+v, 0.0], [0.0,  2.5+v]]], dtype=xp.complex128)
        G_ref[xp.isinf(G_ref)] = xp.inf * (1+1j)  # set real and imaginary part to ∞
        xp_assert_close(G, G_ref)


@make_xp_test_case(_tran_zeros_reduce)
class Test_tran_zeros_reduce:
    """Tests for `_tran_zeros_reduce()`.

    The tests for `_transmission_zeros()` verify the numeric correctness of
    _tran_zeros_reduce()`. The test here just cover special cases and provide coverage.
    """

    def test_basic(self, xp):
        """Basic test for `_tran_zeros_reduce()`. """
        ABCD = xp.eye(4), xp.eye(4), xp.eye(4), _diagsvd([1, 1, 0, 0], 4, 4, xp=xp)
        A_ref, B_ref = xp.eye(2), xp.asarray([[0., -1, 0, 0], [-1, 0, 0, 0]])
        C_ref = xp.asarray([[0., 0], [0, 0], [0, -1], [-1, 0]])
        D_ref = xp.asarray([[0., 0, 0, -1], [0, 0, -1, 0], [1, 0, 0, 0], [0, 1, 0, 0]])

        A, B, C, D = _tran_zeros_reduce(*ABCD, rtol=1e-13, atol=0.)
        xp_assert_close(A, A_ref)
        xp_assert_close(B, B_ref)
        xp_assert_close(C, C_ref)
        xp_assert_close(D, D_ref)

    def test_basic2(self, xp):
        """Test complete reduction of `C_tilde` in `_tran_zeros_reduce()`. """
        A0, B0, C0, D0 = _abcd_normalize(C=xp.eye(4), D=np.diag([1., 2, 0, 0]), xp=xp)
        C_ref = xp.asarray([[0., 0], [0, 0], [0, -1], [-1, 0]])
        D_ref = xp.asarray([[0., 0, 0, 0], [0, 0, 0, 0], [0, 2, 0, 0], [1, 0, 0, 0]])
        A_ref, B_ref, C_ref, D_ref = _abcd_normalize(C=C_ref, D=D_ref)

        A, B, C, D = _tran_zeros_reduce(A0, B0, C0, D0, rtol=1e-13, atol=0.)
        xp_assert_close(A, A_ref)
        xp_assert_close(B, B_ref)
        xp_assert_close(C, C_ref)
        xp_assert_close(D, D_ref)

    @pytest.mark.parametrize("B, D", [(1, np.ones((3, 1))),
        (np.ones((1, 3)), np.diag([1, 1, 0])), (0, np.zeros((3, 0)))])
    def test_no_system_zeros(self, B, D, xp):
        """Test (degenerated) systems with no zeros. """
        ABCD = _abcd_normalize(A=1, B=xp.asarray(B), C=xp.ones((3, 1)), D=xp.asarray(D))

        A, B, C, D = _tran_zeros_reduce(*ABCD, rtol=1e-12, atol=0.)
        xp_assert_close(A, xp.zeros((0, 0)))
        xp_assert_close(B, xp.zeros((0, 0)))
        xp_assert_close(C, xp.zeros((0, 0)))
        xp_assert_close(D, xp.zeros((0, 0)))

    def test_full_rankD(self, xp):
        """Test a system that is already reduced. """
        A0, B0, C0, D0 = xp.eye(4), xp.eye(4), xp.eye(4), xp.eye(4)

        A, B, C, D = _tran_zeros_reduce(A0, B0, C0, D0, rtol=1e-13, atol=0.)
        xp_assert_close(A, A0)
        xp_assert_close(B, B0)
        xp_assert_close(C, C0)
        xp_assert_close(D, D0)


@make_xp_test_case(_transmission_zeros)
class TestTransmissionZeros:
    """Tests for function `_transmission_zeros`.

    The examples are taken from:

    Abbas Emami-Naeini, Paul M. van Dooren, "Computation of zeros of linear
    multivariable systems" Automatica 18.4 (1982): 0005-1098.
    :doi:`10.1016/0005-1098(82)90070-X`
    """

    def test_scalar_system(self, xp):
        """A scalar SISO system should yield its single transmission zero. """
        ABCD = _abcd_normalize(0., 1., 1., 2., xp=xp)
        z = _transmission_zeros(*ABCD, rtol=1e-12, atol=0.)
        xp_assert_close(z, xp.asarray([-0.5 + 0j]))

    def test_static_system(self, xp):
        """Test a static system where `A = []`. """
        ABCD = _abcd_normalize(A=xp.asarray([[]]), D=xp.ones((3, 4)))
        z = _transmission_zeros(*ABCD, rtol=1e-12, atol=0.)
        xp_assert_close(z, xp.asarray([], dtype=z.dtype))

    def test_static_system2(self, xp):
        """Test for A = 0, i.e, zeros at ∞. """
        n, q, p = 3, 4, 8
        A, B, C, D = _abcd_normalize(A=xp.zeros((n, n)),
                                    B=_diagsvd(xp.ones(n), n, p),
                                    C=_diagsvd(xp.ones(n), q, n),
                                    D=_diagsvd(xp.ones(min(q, p)), q, p))
        z = _transmission_zeros(A, B, C, D, rtol=1e-12, atol=0.)
        xp_assert_close(z, xp.asarray([xp.inf+0j]*n))

    def test_full_rank_feedthrough(self, xp):
        """Test for no zeros if D has full row rank but not full column rank.

        This test is needed to achieve complete coverage of `_transmission_zeros()`.
        """
        A, B, C, D = _abcd_normalize([[-1.]], [[0, 1]], [[1]], [[1, 0]], xp=xp)

        z = _transmission_zeros(A, B, C, D)
        xp_assert_equal(z, xp.asarray([], dtype=z.dtype))

    def test_rank_deficient_feedthrough(self, xp):
        """A rank-deficient feedthrough matrix should reduce the system rank. """
        ABCD = _abcd_normalize(A=[[0., 0.], [0., -3.]], B=[[1., 0.], [0., 1.]],
                              C=[[1., 0.], [0.,  1.]], D=[[1., 0.], [0., 0.]], xp=xp)
        z = _transmission_zeros(*ABCD, rtol=1e-12, atol=0.)
        xp_assert_close(z, xp.asarray([-1.+0j]))

    def test_zero_feedthrough(self, xp):
        """A zero feedthrough matrix has no transmission zeros. """
        A, B, C, D = _abcd_normalize(A=0., B=1., C=1., D=0., xp=xp)
        z = _transmission_zeros(A, B, C, D, rtol=1e-12, atol=0.)
        xp_assert_close(z, xp.asarray([], dtype=z.dtype))

    def test_rectangular_D(self, xp):
        """Test a non-square feedthrough matrix. """
        A, B, C, D = _abcd_normalize(A=xp.eye(2), B=xp.ones((2, 3)), C=xp.ones((2, 2)))
        z = _transmission_zeros(A, B, C, D, rtol=1e-12, atol=0.)
        xp_assert_close(z, xp.asarray([1.+0.j], dtype=z.dtype))

    @pytest.mark.parametrize('alpha', [0, 1])
    def test_BC_zero(self, alpha, xp):
        """Test cases where input matrix `B` or output matrix `C` is zero. """
        ABCD = _abcd_normalize(A=xp.eye(2), B=alpha*xp.eye(2), C=(1-alpha)*xp.eye(2))
        z = _transmission_zeros(*ABCD, rtol=1e-12, atol=0.)
        xp_assert_close(z, xp.asarray([], dtype=z.dtype))

    def test_inf_zeros_parameter(self, xp):
        """Verify that parameter `inf_zeros` works as expected. """
        ABCD = _abcd_normalize(A=[[-1, 1, 0], [0, 0, 0], [1, 0, 0]],
                              B=[[0], [0], [1]], C=[[0, -1, 0]], xp=xp)
        z0 = _transmission_zeros(*ABCD, inf_zeros=True,  rtol=1e-12, atol=0.)
        z1 = _transmission_zeros(*ABCD, inf_zeros=False, rtol=1e-12, atol=0.)

        z0 = z0 if xp.isfinite(z0[0]) else xp.flip(z0, axis=0)  # sort z0
        xp_assert_close(z0, xp.asarray([-1+0j, xp.inf+0j]))
        xp_assert_close(z1, xp.asarray([-1 + 0j]))

    def test_harold_example(self, xp):
        """Test adapted from `harold` unit test `test_model_zeros()`. """
        A = xp.asarray([[0, 0, 1, 0, 0, 0],
                        [2, 0, 0, 3, 4, 0],
                        [0, 0, 5, 0, 0, 6],
                        [0, 7, 0, 0, 0, 0],
                        [0, 0, 0, 8, 9, 0],
                        [0, 0, 0, 0, 0, 0]])
        B = xp.asarray([[0, 0, 0, 0, 10, 0], [0, 0, 0, 0, 0, 11]]).T
        C = xp.asarray([[0, 12, 0, 0, 13, 0],
                        [14, 0,  0, 0, 0, 0],
                        [15, 0, 16, 0, 0, 0]])
        D = xp.zeros((3, 2))
        z_ref = xp.asarray([-6.78662791+0.j, 3.09432022+0.j])
        z = _transmission_zeros(A, B, C, D, rtol=0, atol=1e-12)
        z = z[xp.argsort(xp.real(z))]  # sort z into a suitable order
        xp_assert_close(z, z_ref)

    def test_example3(self, xp):
        """Example 3 from the Emami-Naeini and Van Dooren paper. """
        A, B, C, D = _abcd_normalize(A=xp.zeros((6,6)), D=[[1, 0], [1, 0]], xp=xp)
        A[0, 1] = A[1, 2] = A[3, 4] = A[4, 5] = 1
        B[2, 0] = B[5, 1] = C[0, 0] = C[0, 1] = C[1, 3] = 1
        C[1, 4] = -1
        z_ref = xp.asarray([-0.6823278038280190, 0.3411639019140096-1.161541399997251j,
                            0.3411639019140096+1.161541399997251j, 0.9999999999999997])
        z = _transmission_zeros(A, B, C, D, rtol=0, atol=5e-16)
        z = z[xp.argsort(xp.real(z) + xp.imag(z)/16)]  # sort z into a suitable order
        xp_assert_close(z, z_ref)

    def test_example6(self, xp):
        """Example 6 from the Emami-Naeini and Van Dooren paper. """
        A, B, C, D = _abcd_normalize(B=xp.asarray([[0, 0, 1]]).T, C=[[0, -1, 0]], xp=xp)
        A[0,0], A[0, 1], A[2, 0] = -2, 1, 1

        z = _transmission_zeros(A, B, C, D, rtol=0, atol=1e-12)
        z = z if xp.isfinite(z[0]) else xp.flip(z, axis=0)  # sort z
        xp_assert_close(z, xp.asarray([-2.+0j, xp.inf]))

    def test_example8(self, xp):
        """Example 8 from the Emami-Naeini and Van Dooren paper. """
        A = xp.asarray([[-2,  1,  0,  0,  0,  0],
                        [ 1, -2,  1,  0,  1, -1],
                        [ 0,  1, -2,  1,  0,  0],
                        [ 0,  0,  1, -1,  0,  1],
                        [ 0, -1,  0,  0,  0,  0],
                        [ 0,  1,  0, -1,  0,  0]])
        B = xp.asarray([[ 1,  0,  0,  0,  1,  0]]).T
        C = xp.asarray([[ 0,  0,  0,  1,  0,  0]])
        A, B, C, D = _abcd_normalize(A, B, C)

        z = _transmission_zeros(A, B, C, D, rtol=0, atol=1e-12)
        xp_assert_close(z, xp.asarray([-1+0j, -1+0j]))


@make_xp_test_case(signal.freqs_ss)
class TestFreqsSS:
    """Tests for function `freqs_ss`.

    Since `_calcSSTF_helper` has its own tests, `freqs_ss` needs only rudimentary tests.
    """

    def test_invalid_f_scale_parameter(self, xp):
        """Verify that `f_scale` parameter raises an exception. """
        with pytest.raises(ValueError, match="^Parameter f_scale='INVALID'"):
            signal.freqs_ss(A=1, B=1, C=1, D=1, om_spec=1, f_scale='INVALID')

    @pytest.mark.parametrize('om_spec', [
        (0, 1, 'INVALID'), (0, 2, -5), (0, 1j, 5), (3, 2, 5)])
    def test_invalid_om_spec_tuple(self, om_spec, xp):
        """Test for invalid `om_spec` tuple values. """
        with pytest.raises(ValueError, match="^Parameter om_spec="):
            signal.freqs_ss(A=1, B=1, C=1, D=1, om_spec=om_spec)

    @pytest.mark.parametrize('om_spec', ['INVALID', [[1, 2]], [], [3+4j]])
    def test_invalid_om_spec_array(self, om_spec, xp):
        """Test invalid `om_spec` values when `om_spec` is not a tuple. """
        with pytest.raises(ValueError, match="^Parameter `om_spec` must either be"):
            signal.freqs_ss(A=1, B=1, C=1, D=1, om_spec=om_spec)

    @pytest.mark.parametrize('B, C', [
        (np.zeros((3, 2)), np.ones((2, 3))),
        (np.ones((3, 2)), np.zeros((2, 3))),])
    def test_constant_tf(self, B, C, xp):
        """Test that `B == 0` or `C == 0` produces a constant transfer function. """
        A, B = xp.eye(3), xp.asarray(B)
        C, D = xp.asarray(C), xp.asarray([[1., 2.], [3., 4.]])
        G_ref = xp.moveaxis(xp.tile(D, (5, 1, 1)), 0, -1)
        om, G = signal.freqs_ss(A, B, C, D, om_spec=(0, 10, 5))
        xp_assert_close(om, xp.linspace(0, 10, 5))
        xp_assert_close(G, G_ref)

    def test_transfer_functions(self, xp):
        """Compare with transfer functions.

        Consult the docstring of `freqs_ss` for details on the comparison between
        transfer functions and the state-space system.
        """
        om0, om1, d0, d1, gamma  = 5, 15, 0.1, 0.2, 1.
        AA = xp.asarray([[           0,         1,       0,         0],
                         [     -om0**2, -2*d0*om0,       0,         0],
                         [           0,         0,       0,         1],
                         [gamma*om1**2,         0, -om1**2, -2*d1*om1]])
        BB = xp.asarray([[0, om0**2, 0, 0], [0, 0, 0, om1**2]]).T
        CC, DD = xp.asarray([[1, 0, 0, 0], [0, 0, 1, 0]]), xp.zeros((2, 2))

        om, b_ = xp.arange(20, dtype=AA.dtype), xp.asarray([1.,])
        _, H0_ref = signal.freqs(b_, xp.asarray([1/om0**2, 2*d0/om0, 1.]), worN=om)
        _, H1_ref = signal.freqs(b_, xp.asarray([1/om1**2, 2*d1/om1, 1.]), worN=om)
        om_GG, GG = signal.freqs_ss(AA, BB, CC, DD, om_spec=om)  # the response

        xp_assert_equal(om, om_GG)
        xp_assert_close(GG[0, 0, :], H0_ref)
        # pyrefly: ignore [bad-argument-type]
        xp_assert_close(GG[0, 1, :], xp.zeros_like(H0_ref), atol=1e-12)
        xp_assert_close(GG[1, 1, :], H1_ref)
        xp_assert_close(GG[1, 0, :], H1_ref * gamma * H0_ref)

    @pytest.mark.parametrize('N', [1, 3])
    @pytest.mark.parametrize('p', [0, 1, 2])
    @pytest.mark.parametrize('q', [0, 1, 2])
    @pytest.mark.parametrize('n', [0, 1, 2])
    def test_degenerate_system(self, n, q, p, N, xp):
        """Test degenerate systems with empty matrices. """
        ABCD = (xp.eye(n), xp.ones((n, p)), xp.ones((q, n)), xp.ones((q, p)))
        _, GG = signal.freqs_ss(*ABCD, om_spec=(0, 1, N))
        xp_assert_close(GG.shape, (q, p, N))

    @pytest.mark.parametrize('f_scale, om_min, om_max', [
        ('linear', 0., 30.), ('log', 4., 100.)])
    def test_automatic_frequency_range(self, f_scale, om_min, om_max, xp):
        ABCD = xp.asarray(np.diag([2j, 1, -20j])), xp.eye(3), xp.eye(3), xp.eye(3)
        om, _ = signal.freqs_ss(*ABCD, om_spec=2, f_scale=f_scale)
        xp_assert_close(om, xp.asarray((om_min, om_max)))

    @pytest.mark.parametrize('f_scale, om_min, om_max', [
        ('linear', 0., 15), ('log', 2., 50.)])
    def test_automatic_frequency_range_no_imag(self, f_scale, om_min, om_max, xp):
        """Test when `A` has only real eigenvalues. """
        ABCD = xp.eye(3), xp.eye(3), xp.eye(3), xp.eye(3)
        om, _ = signal.freqs_ss(*ABCD, om_spec=2, f_scale=f_scale)
        xp_assert_close(om, (om_min, om_max))

    @pytest.mark.parametrize('f_scale, om_min, om_max', [
        ('linear', 0., 15), ('log', 2., 50.)])
    def test_automatic_frequency_range_one_ev(self, f_scale, om_min, om_max, xp):
        """Test when `A` has only one eigenvalue. """
        ABCD = [xp.ones((1, 1)) for _ in range(4)]
        om, _ = signal.freqs_ss(*ABCD, om_spec=2, f_scale=f_scale)
        xp_assert_close(om, xp.asarray((om_min, om_max)))

    @pytest.mark.parametrize('f_scale, om_min, om_max', [
        ('linear', 0., 15.), ('log', 2., 50.)])
    def test_automatic_frequency_range_no_A(self, f_scale, om_min, om_max, xp):
        """Test when `A` has no eigenvalues to infer a frequency range. """
        A, B, C = xp.zeros((0, 0)), xp.zeros((0, 1)), xp.zeros((2, 0))
        D = xp.ones((2, 1))

        om, _ = signal.freqs_ss(A, B, C, D, om_spec=2, f_scale=f_scale)
        xp_assert_close(om, xp.asarray((om_min, om_max)))

