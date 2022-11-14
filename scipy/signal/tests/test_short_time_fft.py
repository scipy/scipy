from itertools import product
from typing import get_args

import numpy as np
import pytest
from numpy.testing import assert_allclose

from scipy.signal import ShortTimeFFT
from scipy.signal import get_window
from scipy.signal._short_time_fft import FFT_TYP_TYPE, _calc_dual_canonical_window
from scipy.signal.windows import gaussian


def test__calc_dual_canonical_window_roundtrip():
    """Test dual window calculation with a round trip to verify duality.

    Note that this works only for canonical window pairs (having minimal
    energy) like a Gaussian.

    The window is the same as in the example of `from ShortTimeFFT.from_dual`.
    """
    win = gaussian(51, std=10, sym=True)
    d_win = _calc_dual_canonical_window(win, 10)
    win2 = _calc_dual_canonical_window(d_win, 10)
    assert_allclose(win2, win)


def test__calc_dual_canonical_window_exceptions():
    """Raise all exceptions in `_calc_dual_canonical_window`."""
    # Verify that calculation can fail:
    with pytest.raises(ValueError, match="hop=5 is larger than window len.*"):
        _calc_dual_canonical_window(np.ones(4), 5)
    with pytest.raises(ValueError, match=".* Transform not invertible!"):
        _calc_dual_canonical_window(np.array([.1, .2, .3, 0]), 4)

    # Verify that parameter `win` may not be integers:
    with pytest.raises(ValueError, match="Parameter 'win' cannot be of int.*"):
        _calc_dual_canonical_window(np.ones(4, dtype=int), 1)


@pytest.mark.parametrize('n', [8, 9])
def test_fft_func_roundtrip(n: int):
    """Test roundtrip `ifft_func(fft_func(x)) == x` for all permutations of
    relevant parameters. """
    np.random.seed(2394795)
    x0 = np.random.rand(n)
    w, h_n = np.ones(n), 4

    pp = dict(
        fft_type=get_args(FFT_TYP_TYPE),
        mfft=[None, n, n+1, n+2],
        scaling=[None, 'magnitude', 'psd'],
        phase_shift=[None, -n+1, 0, n // 2, n-1])
    for f_typ, mfft, scaling, phase_shift in product(*pp.values()):
        if f_typ == 'onesided2X' and scaling is None:
            continue  # this combination is forbidden
        SFT = ShortTimeFFT(w, h_n, T=1/n, fft_typ=f_typ, mfft=mfft,
                           scale_to=scaling, phase_shift=phase_shift)
        X0 = SFT._fft_func(x0)
        x1 = SFT._ifft_func(X0)
        assert_allclose(x0, x1, err_msg="_fft_func() roundtrip failed for " +
                        f"{f_typ=}, {mfft=}, {scaling=}, {phase_shift=}")

    SFT = ShortTimeFFT(w, h_n, T=1)
    SFT._fft_typ = 'invalid_fft'  # type: ignore
    with pytest.raises(RuntimeError):
        SFT._fft_func(x0)
    with pytest.raises(RuntimeError):
        SFT._ifft_func(x0)  # type: ignore


def test_impulse_roundtrip():
    n = 19
    w, h_n = np.ones(8), 3
    x_in = np.zeros(n)
    x_in[0] = 1

    SFT = ShortTimeFFT(w, hop=h_n, T=1, scale_to=None, phase_shift=None)
    for i in range(0, n):
        x = np.roll(x_in, i)
        Sx = SFT.stft(x)
        # test slicing the input signal into two parts:
        n_q = SFT.nearest_k_p(n // 2)
        Sx0 = SFT.stft(x[:n_q], padding='zeros')
        Sx1 = SFT.stft(x[n_q:], padding='zeros')
        q0_ub = SFT.upper_border_begin(n_q)[1] - SFT.p_min
        q1_le = SFT.lower_border_end[1] - SFT.p_min
        assert_allclose(Sx0[:, :q0_ub], Sx[:, :q0_ub], err_msg=f"{i=}")
        assert_allclose(Sx1[:, q1_le:], Sx[:, q1_le-Sx1.shape[1]:],
                        err_msg=f"{i=}")

        Sx01 = np.hstack((Sx0[:, :q0_ub],
                          Sx0[:, q0_ub:] + Sx1[:, :q1_le],
                          Sx1[:, q1_le:]))
        assert_allclose(Sx, Sx01, atol=1e-8, err_msg=f"{i=}")

        y = SFT.istft(Sx, 0, n)
        assert_allclose(y, x, atol=1e-8, err_msg=f"{i=}")
        y0 = SFT.istft(Sx, 0, n//2)
        assert_allclose(x[:n//2], y0, atol=1e-8, err_msg=f"{i=}")
        y1 = SFT.istft(Sx, n // 2, n)
        assert_allclose(x[n // 2:], y1, atol=1e-8, err_msg=f"{i=}")


@pytest.mark.parametrize('hop', [1, 7, 8])
def test_asymmetric_window_roundtrip(hop: int):
    """An asymmetric window could uncover indexing problems. """
    np.random.seed(23371)

    w = np.arange(16) / 8  # must be of type float
    w[len(w)//2:] = 1
    SFT = ShortTimeFFT(w, hop, T=1)

    x = 10 * np.random.randn(64)
    Sx = SFT.stft(x)
    x1 = SFT.istft(Sx, k1=len(x))
    assert_allclose(x1, x1, err_msg="Roundtrip for asymmetric window with " +
                                    f" {hop=} failed!")


def test_tutorial_stft_sliding_win():
    """Verify example in "Sliding Windows" subsection from "User Guide".

    In :ref:`tutorial_stft_sliding_win` (file ``signal.rst``) of the
    :ref:`user_guide` the behavior the border behavior of
    ``ShortTimeFFT(np.ones(6), 2, T=1)`` with a 50 sample signal is discussed.
    This test verifies the presented indexes.
    """
    SFT = ShortTimeFFT(np.ones(6), 2, T=1)

    # Lower border:
    assert SFT.m_num_mid == 3, f"Slice middle is not 3 but {SFT.m_num_mid=}"
    assert SFT.p_min == -1, f"Lowest slice {SFT.p_min=} is not -1"
    assert SFT.k_min == -5, f"Lowest slice sample {SFT.p_min=} is not -5"
    k_lb, p_lb = SFT.lower_border_end
    assert p_lb == 2, f"First unaffected slice {p_lb=} is not 2"
    assert k_lb == 5, f"First unaffected sample {k_lb=} is not 5"

    n = 50  # upper signal border
    assert (p_max := SFT.p_max(n)) == 27, f"Last slice {p_max=} must be 27"
    assert (k_max := SFT.k_max(n)) == 55, f"Last sample {k_max=} must be 55"
    k_ub, p_ub = SFT.upper_border_begin(n)
    assert p_ub == 24, f"First upper border slice {p_ub=} must be 24"
    assert k_ub == 45, f"First upper border slice {k_ub=} must be 45"


@pytest.mark.parametrize('window, n, nperseg, noverlap',
                         [('boxcar', 100, 10, 0),     # Test no overlap
                          ('boxcar', 100, 10, 9),     # Test high overlap
                          ('bartlett', 101, 51, 26),  # Test odd nperseg
                          ('hann', 1024, 256, 128),   # Test defaults
                          (('tukey', 0.5), 1152, 256, 64),  # Test Tukey
                          ('hann', 1024, 256, 255),   # Test overlapped hann
                          ('boxcar', 100, 10, 3),     # NOLA True, COLA False
                          ('bartlett', 101, 51, 37),  # NOLA True, COLA False
                          ('hann', 1024, 256, 127),   # NOLA True, COLA False
                          (('tukey', 0.5), 1152, 256, 14),  # NOLA True, COLA False
                          ('hann', 1024, 256, 5)])    # NOLA True, COLA False
def test_roundtrip_windows(window, n: int, nperseg: int, noverlap: int):
    """Roundtrip test adapted from `test_spectral.TestSTFT`.

    The parameters are taken from the methods test_roundtrip_real(),
    test_roundtrip_nola_not_cola(), test_roundtrip_float32(),
    test_roundtrip_complex().
    """
    np.random.seed(2394655)

    w = get_window(window, nperseg)
    SFT = ShortTimeFFT(w, nperseg-noverlap, T=1, fft_typ='twosided',
                       phase_shift=None)

    z = 10 * np.random.randn(n) + 10j * np.random.randn(n)
    Sz = SFT.stft(z)
    z1 = SFT.istft(Sz, k1=len(z))
    assert_allclose(z, z1, err_msg="Roundtrip for complex values failed")

    x = 10 * np.random.randn(n)
    Sx = SFT.stft(x)
    x1 = SFT.istft(Sx, k1=len(z))
    assert_allclose(x, x1, err_msg="Roundtrip for float values failed")

    x32 = x.astype(np.float32)
    Sx32 = SFT.stft(x32)
    x32_1 = SFT.istft(Sx32, k1=len(x32))
    assert_allclose(x32, x32_1,
                    err_msg="Roundtrip for 32 Bit float values failed")
