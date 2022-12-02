import math
from itertools import product
from typing import cast, get_args, Literal

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_equal

from scipy.signal import get_window
from scipy.signal import ShortTimeFFT
from scipy.signal._short_time_fft import FFT_TYP_TYPE, \
    _calc_dual_canonical_window, PAD_TYPE
from scipy.signal.windows import gaussian


def chk_VE(match):
    """Assert for a ValueError matching regexp `match`.

    This little wrapper allows a more concise code layout.
    """
    return pytest.raises(ValueError, match=match)


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


def test_invalid_initializer_parameters():
    """Verify that exceptions get raised on invalid parameters when
    instantiating ShortTimeFFT. """
    with chk_VE(r"Parameter win must be 1d, but win.shape=\(2, 2\)!"):
        ShortTimeFFT(np.ones((2, 2)), hop=4, fs=1)
    with chk_VE("Parameter win must have finite entries"):
        ShortTimeFFT(np.array([1, np.inf, 2, 3]), hop=4, fs=1)
    with chk_VE("Parameter hop=0 is not an integer >= 1!"):
        ShortTimeFFT(np.ones(4), hop=0, fs=1)
    with chk_VE("Parameter hop=2.0 is not an integer >= 1!"):
        # noinspection PyTypeChecker
        ShortTimeFFT(np.ones(4), hop=2.0, fs=1)
    with chk_VE(r"dual_win.shape=\(5,\) must equal win.shape=\(4,\)!"):
        ShortTimeFFT(np.ones(4), hop=2, fs=1, dual_win=np.ones(5))
    with chk_VE("Parameter dual_win must be a finite array!"):
        ShortTimeFFT(np.ones(3), hop=2, fs=1, dual_win=np.array([np.nan, 2, 3]))


def test_exceptions_properties_methods():
    """Verify that exceptions get raised when setting properties or calling
    method of ShortTimeFFT to/with invalid values."""
    SFT = ShortTimeFFT(np.ones(8), hop=4, fs=1)
    with chk_VE("Sampling interval T=-1 must be positive!"):
        SFT.T = -1
    with chk_VE("Sampling frequency fs=-1 must be positive!"):
        SFT.fs = -1
    with chk_VE("fft_typ='invalid_typ' not in " +
                r"\('twosided', 'centered', 'onesided', 'onesided2X'\)!"):
        SFT.fft_typ = 'invalid_typ'
    with chk_VE("For scaling is None, fft_typ='onesided2X' is invalid.*"):
        SFT.fft_typ = 'onesided2X'
    with chk_VE("Attribute mfft=7 needs to be at least the window length.*"):
        SFT.mfft = 7
    with chk_VE("scaling='invalid' not in.*"):
        # noinspection PyTypeChecker
        SFT.scale_to('invalid')
    with chk_VE("phase_shift=3.0 has the unit samples.*"):
        SFT.phase_shift = 3.0
    with chk_VE("-mfft < phase_shift < mfft does not hold.*"):
        SFT.phase_shift = 2*SFT.mfft
    with chk_VE("Parameter padding='invalid' not in.*"):
        # noinspection PyTypeChecker
        g = SFT._x_slices(np.zeros(16), k_off=0, p0=0, p1=1, padding='invalid')
        next(g)  # execute generator
    with chk_VE("Trend type must be 'linear' or 'constant'"):
        # noinspection PyTypeChecker
        SFT.stft_detrend(np.zeros(16), detr='invalid')
    with chk_VE("Parameter detr=nan is not a str, function or None!"):
        # noinspection PyTypeChecker
        SFT.stft_detrend(np.zeros(16), detr=np.nan)
    with chk_VE("Invalid Parameter p0=0, p1=200.*"):
        SFT.p_range(100, 0, 200)

    with chk_VE("f_axis=0 may not be equal to t_axis=0!"):
        SFT.istft(np.zeros((SFT.f_pts, 2)), t_axis=0, f_axis=0)
    with chk_VE(r"S.shape\[f_axis\]=2 must be equal to self.f_pts=5.*"):
        SFT.istft(np.zeros((2, 2)))
    with chk_VE(r"S.shape\[t_axis\]=1 needs to have at least 2 slices.*"):
        SFT.istft(np.zeros((SFT.f_pts, 1)))
    with chk_VE(r".*\(k1=100\) <= \(k_max=12\) is false!$"):
        SFT.istft(np.zeros((SFT.f_pts, 3)), k1=100)
    with chk_VE(r"\(k1=1\) - \(k0=0\) = 1 has to be at least.* length 4!"):
        SFT.istft(np.zeros((SFT.f_pts, 3)), k0=0, k1=1)

    with chk_VE(r"Parameter axes_seq='invalid' not in \['tf', 'ft'\]!"):
        # noinspection PyTypeChecker
        SFT.extent(n=100, axes_seq='invalid')
    with chk_VE("Attribute fft_typ=twosided must be in.*"):
        SFT.fft_typ = 'twosided'
        SFT.extent(n=100)


def test_invalid_fft_typ_RuntimeError():
    """Ensure exception gets raised when property `fft_typ` is invalid. """
    SFT = ShortTimeFFT(np.ones(8), hop=4, fs=1)
    SFT._fft_typ = 'invalid_typ'

    with pytest.raises(RuntimeError):
        _ = SFT.f
    with pytest.raises(RuntimeError):
        SFT._fft_func(np.ones(8))
    with pytest.raises(RuntimeError):
        SFT._ifft_func(np.ones(8))


@pytest.mark.parametrize('win_params, Nx', [(('gaussian', 2.), 9),  # in docstr
                                            ('triang', 7),
                                            (('kaiser', 4.0), 9),
                                            (('exponential', None, 1.), 9),
                                            (4.0, 9)])
def test_from_window(win_params, Nx: int):
    """Verify that `from_window()` handels parameters correctly.

    The window parameterizations where take from the `get_window` docstring.
    """
    w, fs = get_window(win_params, Nx, fftbins=False), 16.
    SFT0 = ShortTimeFFT(w, hop=3, fs=fs, fft_typ='twosided', scale_to='psd',
                        phase_shift=1)
    nperseg = len(w)
    noverlap = nperseg - SFT0.hop
    SFT1 = ShortTimeFFT.from_window(win_params, fs, nperseg, noverlap,
                                    fft_typ='twosided', scale_to='psd',
                                    phase_shift=1)
    # Be informative when comparing instances:
    assert_equal(SFT1.win, SFT0.win)
    for n_ in ('hop', 'T', 'fft_typ', 'mfft', 'scaling', 'phase_shift'):
        v0, v1 = getattr(SFT0, n_), getattr(SFT1, n_)
        assert v1 == v0, f"SFT1.{n_}={v1} does not equal SFT0.{n_}={v0}"


def test_dual_win_roundtrip():
    """Verify the duality of `win` and `dual_win`.

    Note that this test does not work for any window, since dual windows are
    not unique. It always works if the windows do not overlap.
    """
    SFT0 = ShortTimeFFT(np.ones(4), hop=4, fs=1)
    SFT1 = ShortTimeFFT.from_dual(SFT0.dual_win, hop=4, fs=1)
    assert_allclose(SFT1.dual_win, SFT0.win)


@pytest.mark.parametrize('scale_to, fac_psd, fac_mag',
                         [(None, 0.25, 0.125),
                          ('magnitude', 2.0, 1),
                          ('psd', 1, 0.5)])
def test_scaling_init(scale_to: Literal['magnitude', 'psd'], fac_psd, fac_mag):
    """Verify scaling calculations when passing `scale_to` to ``__init__(). """
    SFT = ShortTimeFFT(np.ones(4) * 2, hop=4, fs=1, scale_to=scale_to)
    assert SFT.fac_psd == fac_psd
    assert SFT.fac_magnitude == fac_mag


def test_scaling():
    """Verify `scale_to()` method."""
    SFT = ShortTimeFFT(np.ones(4) * 2, hop=4, fs=1, scale_to=None)

    SFT.scale_to('magnitude')
    assert SFT.scaling == 'magnitude'
    assert SFT.fac_psd == 2.0
    assert SFT.fac_magnitude == 1

    SFT.scale_to('psd')
    assert SFT.scaling == 'psd'
    assert SFT.fac_psd == 1
    assert SFT.fac_magnitude == 0.5

    SFT.scale_to('psd')  # needed for coverage

    for scale, s_fac in zip(('magnitude', 'psd'), (8, 4)):
        SFT = ShortTimeFFT(np.ones(4) * 2, hop=4, fs=1, scale_to=None)
        dual_win = SFT.dual_win.copy()

        SFT.scale_to(cast(Literal['magnitude', 'psd'], scale))
        assert_allclose(SFT.dual_win, dual_win * s_fac)


def test_x_slices_padding():
    """Verify padding.

    The reference arrays were taken from  the docstrings of `zero_ext`,
    `const_ext`, `odd_ext()`, and `even_ext()` from the _array_tools module.
    """
    SFT = ShortTimeFFT(np.ones(5), hop=4, fs=1)
    x = np.array([[1, 2, 3, 4, 5], [0, 1, 4, 9, 16]], dtype=float)
    d = {'zeros': [[[0, 0, 1, 2, 3], [0, 0, 0, 1, 4]],
                   [[3, 4, 5, 0, 0], [4, 9, 16, 0, 0]]],
         'edge': [[[1, 1, 1, 2, 3], [0, 0, 0, 1, 4]],
                  [[3, 4, 5, 5, 5], [4, 9, 16, 16, 16]]],
         'even': [[[3, 2, 1, 2, 3], [4, 1, 0, 1, 4]],
                  [[3, 4, 5, 4, 3], [4, 9, 16, 9, 4]]],
         'odd': [[[-1, 0, 1, 2, 3], [-4, -1, 0, 1, 4]],
                 [[3, 4, 5, 6, 7], [4, 9, 16, 23, 28]]]}
    for p_, xx in d.items():
        gen = SFT._x_slices(np.array(x), 0, 0, 2, padding=cast(PAD_TYPE, p_))
        yy = np.array([y_.copy() for y_ in gen])  # due to inplace copying
        assert_equal(yy, xx, err_msg=f"Failed '{p_}' padding.")


def test_invertible():
    """Verify `invertible` property. """
    SFT = ShortTimeFFT(np.ones(8), hop=4, fs=1)
    assert SFT.invertible
    SFT = ShortTimeFFT(np.ones(8), hop=9, fs=1)
    assert not SFT.invertible


def test_border_values():
    """Ensure that minimum and maximum values of slices are correct."""
    SFT = ShortTimeFFT(np.ones(8), hop=4, fs=1)
    assert SFT.p_min == 0
    assert SFT.k_min == -4
    assert SFT.lower_border_end == (4, 1)
    assert SFT.p_max(10) == 4
    assert SFT.k_max(10) == 16
    assert SFT.upper_border_begin(10) == (4, 2)


def test_border_values_exotic():
    """Ensure that the border calculations are correct for windows with
    zeros. """
    w = np.array([0, 0, 0, 0, 0, 0, 0, 1.])
    SFT = ShortTimeFFT(w, hop=1, fs=1)
    assert SFT.lower_border_end == (0, 0)

    SFT = ShortTimeFFT(np.flip(w), hop=20, fs=1)
    assert SFT.upper_border_begin(4) == (0, 0)

    SFT._hop = -1  # provoke unreachable line
    with pytest.raises(RuntimeError):
        _ = SFT.k_max(4)
    with pytest.raises(RuntimeError):
        _ = SFT.k_min


def test_t():
    """Verify that the times of the slices are correct. """
    SFT = ShortTimeFFT(np.ones(8), hop=4, fs=2)
    assert SFT.T == 1/2
    assert SFT.fs == 2.
    assert SFT.delta_t == 4 * 1/2
    t_stft = np.arange(0, SFT.p_max(10)) * SFT.delta_t
    assert_equal(SFT.t(10), t_stft)
    assert_equal(SFT.t(10, 1, 3), t_stft[1:3])
    SFT.T = 1/4
    assert SFT.T == 1/4
    assert SFT.fs == 4
    SFT.fs = 1/8
    assert SFT.fs == 1/8
    assert SFT.T == 8


@pytest.mark.parametrize('fft_typ, f',
                         [('onesided', [0., 1., 2.]),
                          ('onesided2X', [0., 1., 2.]),
                          ('twosided', [0., 1., 2., -2., -1.]),
                          ('centered', [-2., -1., 0., 1., 2.])])
def test_f(fft_typ: FFT_TYP_TYPE, f):
    """Verify the frequency values property `f`."""
    SFT = ShortTimeFFT(np.ones(5), hop=4, fs=5, fft_typ=fft_typ,
                       scale_to='psd')
    assert_equal(SFT.f, f)


def test_extent():
    """Ensure that the `extent()` method is correct. """
    SFT = ShortTimeFFT(np.ones(32), hop=4, fs=32, fft_typ='onesided')
    assert SFT.extent(100, 'tf', False) == (-0.375, 3.625, 0.0, 17.0)
    assert SFT.extent(100, 'ft', False) == (0.0, 17.0, -0.375, 3.625)
    assert SFT.extent(100, 'tf', True) == (-0.4375, 3.5625, -0.5, 16.5)
    assert SFT.extent(100, 'ft', True) == (-0.5, 16.5, -0.4375, 3.5625)

    SFT = ShortTimeFFT(np.ones(32), hop=4, fs=32, fft_typ='centered')
    assert SFT.extent(100, 'tf', False) == (-0.375, 3.625, -16.0, 15.0)


def test_spectrogram():
    """Verify spectrogram and cross-spectrogram methods. """
    SFT = ShortTimeFFT(np.ones(8), hop=4, fs=1)
    x, y = np.ones(10), np.arange(10)
    X, Y = SFT.stft(x), SFT.stft(y)
    assert_allclose(SFT.spectrogram(x), X.real**2+X.imag**2)
    assert_allclose(SFT.spectrogram(x, y), X * Y.conj())


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
        SFT = ShortTimeFFT(w, h_n, fs=n, fft_typ=f_typ, mfft=mfft,
                           scale_to=scaling, phase_shift=phase_shift)
        X0 = SFT._fft_func(x0)
        x1 = SFT._ifft_func(X0)
        assert_allclose(x0, x1, err_msg="_fft_func() roundtrip failed for " +
                        f"{f_typ=}, {mfft=}, {scaling=}, {phase_shift=}")

    SFT = ShortTimeFFT(w, h_n, fs=1)
    SFT._fft_typ = 'invalid_fft'  # type: ignore
    with pytest.raises(RuntimeError):
        SFT._fft_func(x0)
    with pytest.raises(RuntimeError):
        SFT._ifft_func(x0)


def test_impulse_roundtrip():
    n = 19
    w, h_n = np.ones(8), 3
    x_in = np.zeros(n)
    x_in[0] = 1

    SFT = ShortTimeFFT(w, hop=h_n, fs=1, scale_to=None, phase_shift=None)
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
    SFT = ShortTimeFFT(w, hop, fs=1)

    x = 10 * np.random.randn(64)
    Sx = SFT.stft(x)
    x1 = SFT.istft(Sx, k1=len(x))
    assert_allclose(x1, x1, err_msg="Roundtrip for asymmetric window with " +
                                    f" {hop=} failed!")


@pytest.mark.parametrize('m_num', [6, 7])
def test_minimal_length_signal(m_num):
    """Verify that the shortest allowed signal works. """
    SFT = ShortTimeFFT(np.ones(m_num), m_num//2, fs=1)
    n = math.ceil(m_num/2)
    x = np.ones(n)
    Sx = SFT.stft(x)
    x1 = SFT.istft(Sx, k1=n)
    assert_allclose(x1, x1, err_msg="fRoundtrip minimal length signal ({n=})" +
                                    f" for {m_num} sample window failed!")
    with pytest.raises(ValueError, match=rf"len\(x\)={n-1} must be >= ceil.*"):
        SFT.stft(x[:-1])
    with pytest.raises(ValueError, match=rf"S.shape\[t_axis\]={Sx.shape[1]-1}"
                       f" needs to have at least {Sx.shape[1]} slices"):
        SFT.istft(Sx[:, :-1], k1=n)


def test_tutorial_stft_sliding_win():
    """Verify example in "Sliding Windows" subsection from "User Guide".

    In :ref:`tutorial_stft_sliding_win` (file ``signal.rst``) of the
    :ref:`user_guide` the behavior the border behavior of
    ``ShortTimeFFT(np.ones(6), 2, fs=1)`` with a 50 sample signal is discussed.
    This test verifies the presented indexes.
    """
    SFT = ShortTimeFFT(np.ones(6), 2, fs=1)

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


def test_permute_axes():
    """Verify correctness of four-dimensional signal by permuting its
    shape. """
    n = 25
    SFT = ShortTimeFFT(np.ones(8)/8, hop=3, fs=n)
    x0 = np.arange(n)
    Sx0 = SFT.stft(x0)
    Sx0 = Sx0.reshape((Sx0.shape[0], 1, 1, 1, Sx0.shape[-1]))
    SxT = np.moveaxis(Sx0, (0, -1), (-1, 0))

    atol = 2 * np.finfo(SFT.win.dtype).resolution
    for i in range(4):
        y = np.reshape(x0, np.roll((n, 1, 1, 1), i))
        Sy = SFT.stft(y, axis=i)
        assert_allclose(Sy, np.moveaxis(Sx0, 0, i))

        yb = SFT.istft(Sy, k1=n, f_axis=i)
        assert_allclose(yb, y, atol=atol)

        SyT = np.moveaxis(Sy, (i, -1), (-1, i))
        assert_allclose(SyT, np.moveaxis(SxT, 0, i))

        ybT = SFT.istft(SyT, k1=n, t_axis=i, f_axis=-1)
        assert_allclose(ybT, y, atol=atol)


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
                          # NOLA True, COLA False:
                          (('tukey', 0.5), 1152, 256, 14),
                          ('hann', 1024, 256, 5)])    # NOLA True, COLA False
def test_roundtrip_windows(window, n: int, nperseg: int, noverlap: int):
    """Roundtrip test adapted from `test_spectral.TestSTFT`.

    The parameters are taken from the methods test_roundtrip_real(),
    test_roundtrip_nola_not_cola(), test_roundtrip_float32(),
    test_roundtrip_complex().
    """
    np.random.seed(2394655)

    w = get_window(window, nperseg)
    SFT = ShortTimeFFT(w, nperseg-noverlap, fs=1, fft_typ='twosided',
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


@pytest.mark.parametrize('N_x', (128, 129, 255, 256, 1337))  # signal length
@pytest.mark.parametrize('w_size', (128, 256))  # window length
@pytest.mark.parametrize('t_step', (4, 64))  # SFT time hop
@pytest.mark.parametrize('f_c', (7., 23.))  # frequency of input sine
def test_energy_conservation(N_x: int, w_size: int, t_step: int, f_c: float):
    """Test if a `psd`-scaled STFT conserves the L2 norm.

    This test is adapted from MNE-Python [1]_. Besides being battle-tested,
    this test has the benefit of using non-standard window including
    non-positive values and a 2d input signal.

    Since `ShortTimeFFT` requires the signal length `N_x` to be at least the
    window length `w_size`, the parameter `N_x` was changed from
    ``(127, 128, 255, 256, 1337)`` to ``(128, 129, 255, 256, 1337)`` to be
    more useful.

    .. [1] File ``test_stft.py`` of MNE-Python
        https://github.com/mne-tools/mne-python/blob/main/mne/time_frequency/tests/test_stft.py
    """
    window = np.sin(np.arange(.5, w_size + .5) / w_size * np.pi)
    SFT = ShortTimeFFT(window, t_step, fs=1000, fft_typ='onesided2X',
                       scale_to='psd')
    atol = 2*np.finfo(window.dtype).resolution
    N_x = max(N_x, w_size)  # minimal sing
    # Test with low frequency signal
    t = np.arange(N_x).astype(np.float64)
    x = np.sin(2 * np.pi * f_c * t * SFT.T)
    x = np.array([x, x + 1.])
    X = SFT.stft(x)
    xp = SFT.istft(X, k1=N_x)

    max_freq = SFT.f[np.argmax(np.sum(np.abs(X[0]) ** 2, axis=1))]

    assert X.shape[1] == SFT.f_pts
    assert np.all(SFT.f >= 0.)
    assert np.abs(max_freq - f_c) < 1.
    assert_allclose(x, xp, atol=atol)

    # check L2-norm squared (i.e. energy) conservation:
    E_x = np.sum(x**2, axis=-1) * SFT.T  # numerical integration
    aX2 = X.real**2 + X.imag.real**2
    E_X = np.sum(np.sum(aX2, axis=-1) * SFT.delta_t, axis=-1) * SFT.delta_f
    assert_allclose(E_X, E_x, atol=atol)

    # Test with random signal
    np.random.seed(2392795)
    x = np.random.randn(2, N_x)
    X = SFT.stft(x)
    xp = SFT.istft(X, k1=N_x)

    assert X.shape[1] == SFT.f_pts
    assert np.all(SFT.f >= 0.)
    assert np.abs(max_freq - f_c) < 1.
    assert_allclose(x, xp, atol=atol)

    # check L2-norm squared (i.e. energy) conservation:
    E_x = np.sum(x**2, axis=-1) * SFT.T  # numeric integration
    aX2 = X.real ** 2 + X.imag.real ** 2
    E_X = np.sum(np.sum(aX2, axis=-1) * SFT.delta_t, axis=-1) * SFT.delta_f
    assert_allclose(E_X, E_x, atol=atol)

    # Try with empty array
    x = np.zeros((0, N_x))
    X = SFT.stft(x)
    xp = SFT.istft(X, k1=N_x)
    assert xp.shape == x.shape
