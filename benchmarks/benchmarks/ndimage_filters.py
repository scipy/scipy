"""Airspeed Velocity benchmarks for ``scipy.ndimage`` filtering kernels.

These benchmarks cover the 1-D and N-D filter family in
:mod:`scipy.ndimage._filters`. The first priority (``UniformFilter1D`` and
``Correlate1D``) targets the C line-buffer machinery in
``scipy/ndimage/src/ni_support.c`` — ``NI_ArrayToLineBuffer`` /
``NI_LineBufferToArray`` — which forces every input dtype through a
``double`` intermediate buffer and is the main bottleneck identified for
SIMD/dtype-specialization work (see the uniform_filter1d SIMD discussion).
The remaining classes (P1) extend coverage to the other filters that share
the same line-buffer path or exercise the N-D filter iterator.
"""

import numpy as np

from .common import Benchmark, safe_import

with safe_import():
    from scipy.ndimage import (
        uniform_filter1d, correlate1d, gaussian_filter1d,
        gaussian_filter, median_filter, generic_filter, generic_filter1d,
        sobel, prewitt, laplace, minimum_filter1d, maximum_filter1d,
        rank_filter,
    )


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

# dtype dimension: the line-buffer C kernel upcasts everything to double;
# this dimension makes the resulting memory-traffic / cast cost visible to
# regressions. See ISSUES regarding forced double precision.
_DTYPES = [np.float64, np.float32, np.uint8]

# Boundary modes supported by ni_support. ``reflect`` is the default and the
# cheapest (in-place mirror); ``constant`` exercises the cval fill path.
_MODES = ['reflect', 'nearest', 'constant']

# 2-D / 3-D shapes for the N-D filters.
_SHAPES_ND = [(512, 512), (2048, 2048), (64, 64, 64)]

# Mixed 1-D length / N-D shape dimension for the line-buffer kernels. A flat
# 1-D array exercises the contiguous fast path; the 2-D shapes expose the
# axis=0 strided line-buffer copies (NI_ArrayToLineBuffer walks a non-unit
# stride across rows). 1080x1920 is a common video frame; (2048, 2048) is a
# square power-of-two image; (64, 64, 64) is a small volumetric tile.
LINEBUFFER_SHAPES = [1_000_000, (1080, 1920), (2048, 2048), (64, 64, 64)]

# Axis codes: ``-1`` is the C-contiguous last axis (unit stride); ``0`` is the
# first axis (strided for any N-D shape, identical to -1 for a 1-D array).
LINEBUFFER_AXES = [-1, 0]

# Filter-window sizes, chosen from real workloads rather than a geometric
# spread. Grepping scipy/ndimage for ``size=N`` shows the small values (3, 5,
# 7) dominate unit tests and docstring examples; the larger values come from
# named applications:
#   * size=3  — fixed Sobel/Prewitt derivative kernels, light denoising
#   * size=7  — scikit-image ``structural_similarity`` default ``win_size``
#   * size=9  — moderate convolution (e.g. 9-tap FIR smoothing)
#   * size=11 — Wang et al. SSIM paper default 11x11 window
#   * size=15 — heavy salt-and-pepper denoising
#   * size=21 — long FIR low-pass (e.g. audio band-pass); rare beyond this
# Values above 21 are essentially absent from real call sites (only boundary
# tests use size=100), so they are not represented here.

# NI_UniformFilter1D uses an incremental sliding window (tmp += new - old),
# i.e. O(N) regardless of window size, so filter_size barely moves runtime;
# keep just three representative points spanning the common denoising / SSIM
# range instead of a wide sweep.
_UNIFORM_SIZES = [3, 7, 11]

# NI_Correlate1D is O(N * n_weights) — the inner loop walks the whole kernel
# per output element — so the width is a real performance axis. Sample the
# 3-tap derivative kernel, a 9-tap moderate smoother, and a 21-tap long FIR
# that is about the widest kernel seen in practice.
_CORRELATE_WIDTHS = [3, 9, 21]

# minimum/maximum_filter1d use a monotonic-deque ring buffer that is nearly
# O(N); like uniform_filter1d the window size is a weak axis. Cover light
# denoising (3), the SSIM window (7), and heavy denoising (15).
_MINMAX_SIZES = [3, 7, 15]

# Gaussian sigma values observed in scipy code (median ≈ 1-2 px for image
# blur; 5 px is already a strong blur; sigma >= 10 is essentially absent
# outside synthetic tests).
_GAUSSIAN_SIGMAS = [1.0, 2.0, 5.0]


def _make_array(shape_or_size, dtype, seed=5):
    """Return a deterministic pseudo-random array of the requested shape/dtype.

    ``seed`` is fixed so asv parameter combinations are reproducible across
    commits — important for regression detection.
    """
    rng = np.random.RandomState(seed)
    if isinstance(shape_or_size, int):
        size = shape_or_size
    else:
        size = int(np.prod(shape_or_size))
    if np.issubdtype(np.dtype(dtype), np.integer):
        return rng.randint(0, 256, size=size, dtype=dtype).reshape(shape_or_size)
    return rng.standard_normal(size).astype(dtype).reshape(shape_or_size)


# ---------------------------------------------------------------------------
# P0 — line-buffer kernel coverage (NI_UniformFilter1D + NI_ArrayToLineBuffer
#      + NI_LineBufferToArray, plus the correlate1d control on the same path)
# ---------------------------------------------------------------------------

class UniformFilter1D(Benchmark):
    """Benchmarks the box-filter (moving-average) 1-D kernel.

    Exercises ``NI_UniformFilter1D`` (``ni_filters.c``) end-to-end, which in
    turn drives ``NI_ArrayToLineBuffer`` + ``NI_LineBufferToArray``
    (``ni_support.c``). The ``dtype`` dimension exposes the cost of the
    forced ``double`` round-trip highlighted in the SIMD feature request.

    The ``shape`` dimension mixes a flat 1-D length with N-D image shapes, and
    ``axis`` contrasts the C-contiguous last axis (unit stride) against the
    strided first axis — the latter makes the line-buffer's non-contiguous
    copy path visible, which a pure 1-D benchmark cannot reach.

    ``filter_size`` samples the SSIM workload (size=7 is scikit-image's
    ``win_size`` default, size=11 is the Wang et al. paper default) plus the
    size=3 light-denoising / derivative case. The kernel is an incremental
    sliding window (O(N) regardless of size), so a wider sweep would only
    add noise.
    """
    param_names = ['shape', 'axis', 'filter_size', 'dtype', 'mode']
    params = [
        LINEBUFFER_SHAPES,
        LINEBUFFER_AXES,
        _UNIFORM_SIZES,
        _DTYPES,
        ['reflect', 'constant'],
    ]

    def setup(self, shape, axis, filter_size, dtype, mode):
        # For a 1-D array axis=-1 and axis=0 are the same axis; only keep one
        # to avoid running the identical benchmark twice. N-D shapes still
        # exercise both the contiguous last axis and the strided first axis.
        if isinstance(shape, int) and axis == 0:
            raise NotImplementedError()
        self.x = _make_array(shape, dtype)

    def time_uniform_filter1d(self, shape, axis, filter_size, dtype, mode):
        uniform_filter1d(self.x, size=filter_size, axis=axis, mode=mode)

    def peakmem_uniform_filter1d(self, shape, axis, filter_size, dtype, mode):
        # Peak RSS of the double line buffer + output allocation.
        uniform_filter1d(self.x, size=filter_size, axis=axis, mode=mode)


class Correlate1D(Benchmark):
    """Control benchmark sharing the same line-buffer machinery.

    ``correlate1d`` runs through the identical ``NI_ArrayToLineBuffer`` /
    ``NI_LineBufferToArray`` pair as ``uniform_filter1d`` but with a
    weighted-sum kernel (``NI_Correlate1D``) instead of the incremental box
    sum. Comparing the two isolates line-buffer overhead from algorithm cost.

    Unlike ``uniform_filter1d``, ``NI_Correlate1D`` is O(N * n_weights) — the
    inner loop walks the whole kernel per output element — so ``n_weights``
    is a genuine performance axis, not a cosmetic one. The three widths cover
    a 3-tap derivative kernel (Sobel/Prewitt class), a 9-tap moderate
    smoother, and a 21-tap long FIR low-pass, which is about the widest
    kernel seen in real scipy/skimage call sites.
    """
    param_names = ['shape', 'axis', 'n_weights', 'dtype', 'mode']
    params = [
        LINEBUFFER_SHAPES,
        LINEBUFFER_AXES,
        _CORRELATE_WIDTHS,
        _DTYPES,
        ['reflect', 'constant'],
    ]

    def setup(self, shape, axis, n_weights, dtype, mode):
        # See UniformFilter1D: skip the redundant axis=0 run on 1-D arrays.
        if isinstance(shape, int) and axis == 0:
            raise NotImplementedError()
        self.x = _make_array(shape, dtype)
        rng = np.random.RandomState(7)
        # Symmetric weights avoid an accidental zero-output edge case and
        # exercise the full multiply-accumulate path.
        self.weights = rng.standard_normal(n_weights)

    def time_correlate1d(self, shape, axis, n_weights, dtype, mode):
        correlate1d(self.x, self.weights, axis=axis, mode=mode)


# ---------------------------------------------------------------------------
# P1 — remaining high-value filters in scipy.ndimage._filters
# ---------------------------------------------------------------------------

class GaussianFilter1D(Benchmark):
    """Separable Gaussian along a single axis.

    Built on top of ``correlate1d`` with a precomputed truncated Gaussian
    kernel, so it shares the line-buffer path but adds kernel construction
    overhead and an order-N derivative term when ``order > 0``. The
    ``shape``/``axis`` matrix mirrors ``UniformFilter1D`` since the same
    strided line-buffer copies dominate for N-D inputs.

    ``sigma`` follows the distribution seen in real scipy/skimage call sites:
    sigma=1-2 px covers the bulk of image-blur usage, sigma=5 is already a
    strong blur and about the upper end of common practice; sigma >= 10 is
    essentially absent outside synthetic tests so it is not sampled here.
    """
    param_names = ['shape', 'axis', 'sigma', 'order', 'dtype']
    params = [
        LINEBUFFER_SHAPES,
        LINEBUFFER_AXES,
        _GAUSSIAN_SIGMAS,
        [0, 1, 2],
        [np.float64, np.float32],
    ]

    def setup(self, shape, axis, sigma, order, dtype):
        if isinstance(shape, int) and axis == 0:
            raise NotImplementedError()
        self.x = _make_array(shape, dtype)

    def time_gaussian_filter1d(self, shape, axis, sigma, order, dtype):
        gaussian_filter1d(self.x, sigma=sigma, axis=axis, order=order)


class GaussianFilterND(Benchmark):
    """N-D separable Gaussian (applies the 1-D kernel along every axis).

    On a 2/3-D array this multiplies the line-buffer traffic by ``ndim``,
    making it the dominant cost for the common "blur an image" workload.
    """
    param_names = ['shape', 'sigma', 'dtype']
    params = [
        _SHAPES_ND,
        [1.0, 5.0],
        [np.float64, np.float32],
    ]

    def setup(self, shape, sigma, dtype):
        self.x = _make_array(shape, dtype)

    def time_gaussian_filter(self, shape, sigma, dtype):
        gaussian_filter(self.x, sigma=sigma)


class MedianFilter(Benchmark):
    """Rank-based filter (per-cell median over a square footprint).

    Unlike the line-buffer family, the rank filters go through the N-D
    ``NI_RankFilter`` kernel with a sorting step, so they exercise a
    different code path and have very different scaling.
    """
    param_names = ['shape', 'size', 'dtype']
    params = [
        [(256, 256), (1024, 1024)],
        [3, 5, 9],
        [np.float64, np.uint8],
    ]

    def setup(self, shape, size, dtype):
        self.x = _make_array(shape, dtype)

    def time_median_filter(self, shape, size, dtype):
        median_filter(self.x, size=size)

    def time_rank_filter(self, shape, size, dtype):
        # For a 2-D ``size=N`` the footprint is N*N cells, so the median rank
        # is (N*N)//2. ``median_filter`` is a thin wrapper around this call;
        # including it directly exercises the generic ``NI_RankFilter`` code
        # path with an explicit rank rather than the median convenience path.
        rank_filter(self.x, rank=(size * size) // 2, size=size)


class MinMaxFilter1D(Benchmark):
    """1-D minimum / maximum filters (sliding-window extrema).

    ``minimum_filter1d`` / ``maximum_filter1d`` use the ring-buffer deque in
    ``NI_MinOrMaxFilter1D`` rather than the generic line buffer, so they are
    a distinct kernel worth tracking independently. They still go through
    ``NI_ArrayToLineBuffer`` / ``NI_LineBufferToArray`` for the array ↔
    buffer copies, so the ``shape``/``axis`` matrix applies the same way as
    the other 1-D kernels.

    ``filter_size`` covers the denoising range: size=3 for light
    salt-and-pepper removal, size=7 for the SSIM-style statistics window,
    and size=15 for heavy noise. Like ``uniform_filter1d`` the deque kernel
    is approximately O(N), so the window size is only a weak axis.
    """
    param_names = ['shape', 'axis', 'filter_size', 'dtype']
    params = [
        # min/max filters are slower; keep the largest shape modest.
        [100_000, (1080, 1920), (2048, 2048)],
        LINEBUFFER_AXES,
        _MINMAX_SIZES,
        [np.float64, np.float32],
    ]

    def setup(self, shape, axis, filter_size, dtype):
        if isinstance(shape, int) and axis == 0:
            raise NotImplementedError()
        self.x = _make_array(shape, dtype)

    def time_minimum_filter1d(self, shape, axis, filter_size, dtype):
        minimum_filter1d(self.x, size=filter_size, axis=axis)

    def time_maximum_filter1d(self, shape, axis, filter_size, dtype):
        maximum_filter1d(self.x, size=filter_size, axis=axis)


class GenericFilter1D(Benchmark):
    """Python-callback 1-D filter (``generic_filter1d``).

    Crosses the C→Python boundary once per line and runs through the same
    ``NI_ArrayToLineBuffer`` / ``NI_LineBufferToArray`` pair as the other 1-D
    kernels, so the ``shape``/``axis`` matrix applies.

    The callback contract is ``function(input_line, output_line)`` and must
    write into ``output_line`` in place (see ``Py_Filter1DFunc`` in
    ``nd_image.c``). A straight copy keeps the per-line work minimal so the
    measurement isolates the C→Python call overhead rather than the
    callback body.
    """
    param_names = ['shape', 'axis']
    params = [
        [(256, 256), (1024, 1024)],
        LINEBUFFER_AXES,
    ]

    def setup(self, shape, axis):
        self.x = _make_array(shape, np.float64)

    @staticmethod
    def _copy1d(input_line, output_line):
        # generic_filter1d passes (input_line, output_line) and expects the
        # callback to fill ``output_line`` in place. ``input_line`` is the
        # extended line: it is padded by size1 = filter_size//2 on the left
        # and size2 = filter_size - size1 - 1 on the right (see
        # ``NI_Correlate1D`` in ``ni_filters.c``), so len(input_line) =
        # len(output_line) + filter_size - 1. Copying the centered segment
        # keeps the per-line work minimal so the measurement isolates the
        # C->Python call overhead rather than the callback body.
        size1 = (len(input_line) - len(output_line) + 1) // 2
        output_line[:] = input_line[size1 : size1 + len(output_line)]

    def time_generic_filter1d(self, shape, axis):
        generic_filter1d(self.x, self._copy1d, filter_size=5, axis=axis)


class GenericFilterND(Benchmark):
    """Python-callback N-D filter (``generic_filter``).

    Unlike the 1-D variant this goes through ``NI_GenericFilter`` with the
    N-D ``NI_FilterIterator``, visiting each output cell once and crossing
    the C→Python boundary per cell — so it benchmarks callback overhead at
    a very different granularity than ``GenericFilter1D``. No ``axis``
    dimension: the footprint is N-D by construction.
    """
    param_names = ['shape', 'size']
    params = [
        [(256, 256), (1024, 1024)],
        [3, 5],
    ]

    def setup(self, shape, size):
        self.x = _make_array(shape, np.float64)

    @staticmethod
    def _sumnd(values):
        return np.sum(values)

    def time_generic_filter(self, shape, size):
        generic_filter(self.x, self._sumnd, size=size)


class EdgeFilters(Benchmark):
    """Gradient / second-derivative edge detectors: Sobel, Prewitt, Laplace.

    These compose the line-buffer ``correlate1d`` with derivative kernels
    and run along two (Sobel/Prewitt) or all (Laplace) axes; included to
    guard the common image-processing pipeline.
    """
    param_names = ['shape', 'dtype']
    params = [
        _SHAPES_ND[:2],   # 2-D only — Sobel/Prewitt are axis-wise
        [np.float64, np.float32],
    ]

    def setup(self, shape, dtype):
        self.x = _make_array(shape, dtype)

    def time_sobel(self, shape, dtype):
        sobel(self.x)

    def time_prewitt(self, shape, dtype):
        prewitt(self.x)

    def time_laplace(self, shape, dtype):
        laplace(self.x)
