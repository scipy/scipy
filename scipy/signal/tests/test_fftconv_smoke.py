import numpy as np
import pytest

from scipy import signal


@pytest.mark.parametrize("x", [-0.3, -0.1, -0.01, -1e-4])
def test_lambertw_wm1_inverse(x):
    # W_{-1}(x) satisfies W*exp(W) = x on the lower branch (W <= -1).
    from scipy.signal import _fftconv
    from scipy.special import lambertw
    w = _fftconv._lambertw_wm1(x)
    assert w <= -1.0
    np.testing.assert_allclose(w * np.exp(w), x, rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(w, lambertw(x, k=-1).real, rtol=1e-12)


def test_lambertw_wm1_boundary():
    # At x = -1/e the W_{-1} and W_0 branches meet: W_{-1}(-1/e) = -1 exactly.
    from scipy.signal import _fftconv
    w = _fftconv._lambertw_wm1(-1.0 / np.e)
    assert np.isfinite(w)
    np.testing.assert_allclose(w, -1.0, atol=1e-9)


@pytest.mark.parametrize("s1,s2", [
    (500, 500), (4000, 4000),      # equal sizes -> fallback
    (2000, 1), (1, 2000),          # unit length -> fallback
    (2000, 1500),                  # s2 >= s1/2 -> fallback
    (2000, 65), (100000, 65), (100000, 512),
    (65, 2000), (512, 100000),     # swapped (kernel-first) orientation
])
def test_calc_oa_lens_matches_python(s1, s2):
    from scipy.signal import _fftconv
    from scipy.signal._signaltools import _calc_oa_lens
    assert _fftconv._calc_oa_lens(s1, s2) == _calc_oa_lens(s1, s2)


def test_calc_oa_lens_known_values():
    from scipy.signal import _fftconv
    # (block_size, overlap, in1_step, in2_step)
    assert _fftconv._calc_oa_lens(2000, 65) == (512, 64, 448, 65)
    assert _fftconv._calc_oa_lens(100000, 65) == (512, 64, 448, 65)
    assert _fftconv._calc_oa_lens(100000, 512) == (5250, 511, 4739, 512)
    assert _fftconv._calc_oa_lens(500, 500)[1] is None   # fallback


@pytest.mark.parametrize("mode", ["full", "same", "valid"])
@pytest.mark.parametrize("dtype", [np.float64, np.float32])
def test_fftconv_oaconvolve_1d_matches_reference(mode, dtype):
    from scipy.signal import _fftconv
    rng = np.random.default_rng(0)
    a = rng.standard_normal(5000).astype(dtype)
    b = rng.standard_normal(65).astype(dtype)
    got = _fftconv.oaconvolve(a, b, mode, None)
    assert got is not None
    ref = signal.fftconvolve(a, b, mode=mode)
    rtol = 1e-5 if dtype == np.float32 else 1e-10
    np.testing.assert_allclose(got, ref, rtol=rtol, atol=rtol)


def test_fftconv_oaconvolve_non_contiguous_input():
    from scipy.signal import _fftconv
    rng = np.random.default_rng(0)
    a = rng.standard_normal(10000)[::2]   # non-contiguous view (stride 2)
    b = rng.standard_normal(65)
    assert not a.flags['C_CONTIGUOUS']
    got = _fftconv.oaconvolve(a, b, "full", None)
    assert got is not None   # handled (copied), not declined
    np.testing.assert_allclose(got, signal.fftconvolve(a, b, mode="full"),
                               rtol=1e-10, atol=1e-10)


def test_fftconv_oaconvolve_list_axes():
    from scipy.signal import _fftconv
    rng = np.random.default_rng(0)
    a = rng.standard_normal(5000)
    b = rng.standard_normal(65)
    got = _fftconv.oaconvolve(a, b, "full", [0])   # list axes must be accepted
    assert got is not None
    np.testing.assert_allclose(got, signal.fftconvolve(a, b, mode="full"),
                               rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize("mode", ["full", "same", "valid"])
@pytest.mark.parametrize("dtype", [np.complex128, np.complex64])
def test_fftconv_oaconvolve_complex_matches_reference(mode, dtype):
    from scipy.signal import _fftconv
    rng = np.random.default_rng(0)
    a = (rng.standard_normal(5000) + 1j*rng.standard_normal(5000)).astype(dtype)
    b = (rng.standard_normal(65) + 1j*rng.standard_normal(65)).astype(dtype)
    got = _fftconv.oaconvolve(a, b, mode, None)
    assert got is not None
    ref = signal.fftconvolve(a, b, mode=mode)
    rtol = 1e-5 if dtype == np.complex64 else 1e-10
    np.testing.assert_allclose(got, ref, rtol=rtol, atol=rtol)


def _random_of_dtype(rng, size, dtype):
    r = rng.standard_normal(size)
    if np.issubdtype(dtype, np.complexfloating):
        r = r + 1j*rng.standard_normal(size)
    return r.astype(dtype)


@pytest.mark.parametrize("mode", ["full", "same", "valid"])
@pytest.mark.parametrize("dtype",
                         [np.float64, np.complex128, np.float32, np.complex64])
@pytest.mark.parametrize("shape_axis",
                         [((8, 2000), 1), ((2000, 8), 0), ((4, 2000, 3), 1)])
def test_fftconv_oaconvolve_nd_single_axis(mode, dtype, shape_axis):
    from scipy.signal import _fftconv
    shape, axis = shape_axis
    rng = np.random.default_rng(0)
    kshape = list(shape)
    kshape[axis] = 65
    a = _random_of_dtype(rng, shape, dtype)
    b = _random_of_dtype(rng, tuple(kshape), dtype)
    got = _fftconv.oaconvolve(a, b, mode, axis)
    assert got is not None
    ref = signal.fftconvolve(a, b, mode=mode, axes=[axis])
    rtol = 1e-4 if dtype in (np.complex64, np.float32) else 1e-9
    np.testing.assert_allclose(got, ref, rtol=rtol, atol=rtol)


@pytest.mark.parametrize("mode", ["full", "same", "valid"])
@pytest.mark.parametrize("dtype", [np.float64, np.complex128])
def test_fftconv_oaconvolve_swapped_inputs(mode, dtype):
    # Kernel passed FIRST (in1 shorter than in2) exercises the internal
    # sig/ker reorientation; 'same' output length must follow the original
    # in1, not the reoriented one.
    from scipy.signal import _fftconv
    rng = np.random.default_rng(0)
    a = _random_of_dtype(rng, 65, dtype)
    b = _random_of_dtype(rng, 5000, dtype)
    got = _fftconv.oaconvolve(a, b, mode, None)
    assert got is not None
    ref = signal.fftconvolve(a, b, mode=mode)
    np.testing.assert_allclose(got, ref, rtol=1e-9, atol=1e-9)


def test_fftconv_oaconvolve_exact_multiple_block():
    # ns an exact multiple of the block step: with m=65 the block step is
    # 448 (block_size 512, overlap 64), so n=4480 fills the last block
    # exactly and must not read or add past the full-length output.
    from scipy.signal import _fftconv
    rng = np.random.default_rng(0)
    a = rng.standard_normal(4480)
    b = rng.standard_normal(65)
    got = _fftconv.oaconvolve(a, b, "full", None)
    assert got is not None
    np.testing.assert_allclose(got, signal.fftconvolve(a, b, mode="full"),
                               rtol=1e-9, atol=1e-9)


@pytest.mark.parametrize("mode", ["full", "same", "valid"])
@pytest.mark.parametrize("dtype", [np.float64, np.complex128])
def test_fftconv_oaconvolve_nd_kernel_longer(mode, dtype):
    # N-D single-axis with the KERNEL longer than the signal along the conv
    # axis: exercises the internal swap in the batched (outer, n, inner) core.
    from scipy.signal import _fftconv
    rng = np.random.default_rng(0)
    a = _random_of_dtype(rng, (4, 65, 3), dtype)
    b = _random_of_dtype(rng, (4, 2000, 3), dtype)
    got = _fftconv.oaconvolve(a, b, mode, 1)
    assert got is not None
    ref = signal.fftconvolve(a, b, mode=mode, axes=[1])
    np.testing.assert_allclose(got, ref, rtol=1e-9, atol=1e-9)


@pytest.mark.parametrize("da,db", [
    (np.complex128, np.complex64),
    (np.float64, np.complex128),
    (np.float64, np.float32),
])
def test_fftconv_oaconvolve_mixed_dtype_declines(da, db):
    from scipy.signal import _fftconv
    rng = np.random.default_rng(0)
    a = _random_of_dtype(rng, 2000, da)
    b = _random_of_dtype(rng, 65, db)
    assert _fftconv.oaconvolve(a, b, "full", None) is None   # mixed dtype -> decline
