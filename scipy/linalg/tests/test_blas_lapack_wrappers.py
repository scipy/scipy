"""Contract tests for the hand-written C++ BLAS/LAPACK binding layer.

These do not test numerics -- ``test_blas.py``, ``test_fblas.py`` and ``test_lapack.py``
do that. They test the *wrapper machinery* in ``scipy/linalg/src/``: the parts that
``f2py`` used to generate from the ``.pyf.src`` declarations and that are now written by
hand, where a transcription slip produces a silently wrong answer rather than a build
failure.

Covered here:

* default ``n`` for the Level-1 BLAS routines with a stride and/or an offset
* ``overwrite_*`` semantics, including when a copy is forced regardless of the flag
* ``*_lwork`` workspace queries agreeing with the routines they size
* the ``gees``/``gges`` eigenvalue-sort callbacks, including failure unwinding
* the argument protocol (``parse_args``) and its CPython-compatible message order
* the permissive scalar coercions carried over from f2py's ``*_from_pyobj``
* f2py's rank reinterpretation for arrays whose rank does not match the declaration
* zero-sized inputs
* a smoke pass over every routine both modules expose
"""
import gc
import sys

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_equal
from pytest import raises as assert_raises

from scipy.linalg import get_lapack_funcs
from scipy._lib._testutils import IS_WASM

try:
    from scipy.linalg import _fblas as fblas
except ImportError:
    fblas = None

try:
    from scipy.linalg import _flapack as flapack
except ImportError:
    flapack = None

REAL_DTYPES = [np.float32, np.float64]
COMPLEX_DTYPES = [np.complex64, np.complex128]
DTYPES = REAL_DTYPES + COMPLEX_DTYPES

# float32/complex64 need a loose tolerance for the accumulating routines.
TOL = {np.float32: 1e-5, np.complex64: 1e-5, np.float64: 1e-12, np.complex128: 1e-12}


def _vector(n, dtype, rng):
    x = rng.uniform(size=n)
    if np.issubdtype(dtype, np.complexfloating):
        x = x + 1j * rng.uniform(size=n)
    return x.astype(dtype)


_PREFIX = {np.float32: 's', np.float64: 'd', np.complex64: 'c', np.complex128: 'z'}

# Level-1 families whose Python name is not simply <prefix><base>. `get_blas_funcs`
# cannot reach several of these (it would build `casum`, `diamax`, `zrot`), so the
# wrapper-layer tests below address the module directly.
_L1_NAMES = {
    'iamax': {np.float32: 'isamax', np.float64: 'idamax',
              np.complex64: 'icamax', np.complex128: 'izamax'},
    'asum': {np.float32: 'sasum', np.float64: 'dasum',
             np.complex64: 'scasum', np.complex128: 'dzasum'},
    'nrm2': {np.float32: 'snrm2', np.float64: 'dnrm2',
             np.complex64: 'scnrm2', np.complex128: 'dznrm2'},
    'rot': {np.float32: 'srot', np.float64: 'drot',
            np.complex64: 'csrot', np.complex128: 'zdrot'},
    'dot': {np.float32: 'sdot', np.float64: 'ddot',
            np.complex64: 'cdotu', np.complex128: 'zdotu'},
}


def _blas(base, dtype):
    """The `_fblas` routine for `base` at `dtype`, by its real Python name."""
    name = _L1_NAMES.get(base, {}).get(dtype, _PREFIX[dtype] + base)
    return getattr(fblas, name)


# --------------------------------------------------------------------------------------
# Default `n` for the Level-1 routines.
#
# The wrappers spell it `(len(x) - 1 - offx) / abs(incx) + 1`, which is exactly
# `len(x[offx::incx])`. An earlier revision used `(len(x) - offx) / abs(incx)`, which is
# one too few whenever `len(x) - offx` is not a multiple of `incx` -- a silently short
# read that the rest of the suite never noticed because it always passes `incx=1`.
# --------------------------------------------------------------------------------------
STRIDES = [(0, 1), (0, 2), (0, 3), (1, 2), (2, 2), (2, 3), (3, 1), (1, 4)]


@pytest.mark.skipif(fblas is None, reason='_fblas not built')
class TestDefaultN:

    @pytest.mark.parametrize('offx,incx', STRIDES)
    def test_iamax_visits_whole_strided_segment(self, offx, incx):
        # |x| strictly increasing, so iamax returns the index of the *last* element
        # visited: a direct read-out of the `n` the wrapper defaulted to.
        x = np.arange(1.0, 8.0)
        iamax = _blas('iamax', np.float64)
        assert_equal(iamax(x, offx=offx, incx=incx), len(x[offx::incx]) - 1)

    @pytest.mark.parametrize('dtype', DTYPES)
    @pytest.mark.parametrize('offx,incx', STRIDES)
    def test_asum_default_n(self, dtype, offx, incx):
        rng = np.random.default_rng(1234)
        x = _vector(7, dtype, rng)
        seg = x[offx::incx]
        asum = _blas('asum', dtype)
        # BLAS ?asum sums |re| + |im|; for a real flavor that is just |x|.
        expected = np.abs(seg.real).sum() + np.abs(seg.imag).sum()
        assert_allclose(asum(x, offx=offx, incx=incx), expected, rtol=TOL[dtype])

    @pytest.mark.parametrize('dtype', DTYPES)
    @pytest.mark.parametrize('offx,incx', STRIDES)
    def test_nrm2_default_n(self, dtype, offx, incx):
        rng = np.random.default_rng(5678)
        x = _vector(7, dtype, rng)
        nrm2 = _blas('nrm2', dtype)
        assert_allclose(nrm2(x, offx=offx, incx=incx),
                        np.linalg.norm(x[offx::incx]), rtol=TOL[dtype])

    @pytest.mark.parametrize('dtype', DTYPES)
    @pytest.mark.parametrize('offx,incx', STRIDES)
    def test_dot_default_n(self, dtype, offx, incx):
        rng = np.random.default_rng(99)
        x = _vector(7, dtype, rng)
        y = _vector(7, dtype, rng)
        dot = _blas('dot', dtype)
        # x and y are stepped independently, both defaulting off `len(x)`.
        seg_x = x[offx::incx]
        seg_y = y[offx::incx]
        assert_allclose(dot(x, y, offx=offx, incx=incx, offy=offx, incy=incx),
                        np.dot(seg_x, seg_y), rtol=TOL[dtype])

    @pytest.mark.parametrize('dtype', DTYPES)
    @pytest.mark.parametrize('offx,incx', STRIDES)
    def test_copy_default_n(self, dtype, offx, incx):
        rng = np.random.default_rng(4321)
        x = _vector(7, dtype, rng)
        y = np.zeros(7, dtype=dtype)
        copy = _blas('copy', dtype)
        out = copy(x, y, offx=offx, incx=incx, offy=offx, incy=incx)
        expected = np.zeros(7, dtype=dtype)
        expected[offx::incx] = x[offx::incx]
        assert_allclose(out, expected, rtol=TOL[dtype])

    @pytest.mark.parametrize('dtype', DTYPES)
    def test_explicit_n_matches_default(self, dtype):
        # The default is documented as `len(x[offx::incx])`; passing that value
        # explicitly must be indistinguishable from omitting it.
        rng = np.random.default_rng(7)
        x = _vector(9, dtype, rng)
        asum = _blas('asum', dtype)
        for offx, incx in STRIDES:
            n = len(x[offx::incx])
            assert_allclose(asum(x, offx=offx, incx=incx),
                            asum(x, n=n, offx=offx, incx=incx), rtol=TOL[dtype])

    @pytest.mark.parametrize('n_extra', [1, 2, 5])
    def test_n_beyond_the_segment_is_rejected(self, n_extra):
        # `CHECK(len(x) - offx > (n - 1) * abs(incx), n)`
        x = np.arange(1.0, 8.0)
        asum = _blas('asum', np.float64)
        n = len(x[::2]) + n_extra
        with assert_raises(ValueError):
            asum(x, n=n, incx=2)

    @pytest.mark.parametrize('offx', [-1, 7, 100])
    def test_offset_out_of_range_is_rejected(self, offx):
        x = np.arange(1.0, 8.0)
        asum = _blas('asum', np.float64)
        with assert_raises(ValueError):
            asum(x, offx=offx)

    def test_zero_increment_is_rejected(self):
        x = np.arange(1.0, 8.0)
        asum = _blas('asum', np.float64)
        with assert_raises(ValueError):
            asum(x, incx=0)

    def test_nrm2_rejects_negative_increment(self):
        # Unlike other Level-1 wrappers, nrm2 requires incx > 0.
        # Accelerate aborts when asum has a negative increment.
        x = np.arange(1.0, 8.0)
        nrm2 = _blas('nrm2', np.float64)
        with assert_raises(ValueError):
            nrm2(x, incx=-1)


# --------------------------------------------------------------------------------------
# overwrite_* semantics.
# --------------------------------------------------------------------------------------
@pytest.mark.skipif(fblas is None or flapack is None, reason='modules not built')
class TestOverwriteSemantics:

    @pytest.mark.parametrize('dtype', DTYPES)
    def test_gesv_respects_overwrite_flags(self, dtype):
        rng = np.random.default_rng(11)
        a0 = np.asfortranarray(np.eye(4, dtype=dtype) * 3 + 1)
        b0 = np.asfortranarray(_vector(4, dtype, rng).reshape(4, 1))
        gesv = get_lapack_funcs('gesv', dtype=dtype)

        a, b = a0.copy(order='F'), b0.copy(order='F')
        gesv(a, b, overwrite_a=0, overwrite_b=0)
        assert_equal(a, a0)
        assert_equal(b, b0)

        a, b = a0.copy(order='F'), b0.copy(order='F')
        lu, piv, x, info = gesv(a, b, overwrite_a=1, overwrite_b=1)
        assert info == 0
        # With the flag set and the buffers already suitable, LAPACK wrote through them.
        assert np.shares_memory(a, lu)
        assert np.shares_memory(b, x)
        assert not np.array_equal(a, a0)

    @pytest.mark.parametrize('dtype', DTYPES)
    def test_overwrite_ignored_when_layout_unsuitable(self, dtype):
        # A C-ordered 2-D array is not Fortran-contiguous, so the wrapper must copy even
        # with overwrite_a=1. Silently writing through would corrupt the caller's data.
        a0 = np.ascontiguousarray(np.eye(4, dtype=dtype) * 3 + 1)
        assert not a0.flags.f_contiguous
        b0 = np.asfortranarray(np.ones((4, 1), dtype=dtype))
        gesv = get_lapack_funcs('gesv', dtype=dtype)

        a = a0.copy(order='C')
        lu, piv, x, info = gesv(a, b0.copy(order='F'), overwrite_a=1, overwrite_b=1)
        assert info == 0
        assert_equal(a, a0)
        assert not np.shares_memory(a, lu)

    @pytest.mark.parametrize('dtype', DTYPES)
    def test_overwrite_ignored_on_dtype_mismatch(self, dtype):
        # An int array must be converted, so the caller's buffer cannot be written.
        a = np.asfortranarray(np.eye(4, dtype=np.int32) * 3 + 1)
        b = np.asfortranarray(np.ones((4, 1), dtype=np.int32))
        gesv = get_lapack_funcs('gesv', dtype=dtype)
        lu, piv, x, info = gesv(a, b, overwrite_a=1, overwrite_b=1)
        assert info == 0
        assert_equal(a, np.eye(4, dtype=np.int32) * 3 + 1)
        assert lu.dtype == dtype

    @pytest.mark.parametrize('dtype', DTYPES)
    def test_axpy_writes_y_in_place(self, dtype):
        # axpy's `y` is intent(in,out) with no `copy` and no overwrite_y flag, so it is
        # always written through when the layout allows. This mirrors f2py exactly.
        rng = np.random.default_rng(13)
        x = _vector(5, dtype, rng)
        y = _vector(5, dtype, rng)
        y_before = y.copy()
        axpy = _blas('axpy', dtype)
        out = axpy(x, y, a=2.0)
        assert np.shares_memory(y, out)
        assert_allclose(y, y_before + 2.0 * x, rtol=TOL[dtype])

    @pytest.mark.parametrize('dtype', DTYPES)
    def test_rot_honours_its_overwrite_flags(self, dtype):
        # rot, unlike axpy, does carry overwrite_x/overwrite_y.
        rng = np.random.default_rng(17)
        x = _vector(5, dtype, rng)
        y = _vector(5, dtype, rng)
        x0, y0 = x.copy(), y.copy()
        rot = _blas('rot', dtype)

        rot(x, y, 0.6, 0.8)
        assert_equal(x, x0)
        assert_equal(y, y0)

        xo, yo = rot(x, y, 0.6, 0.8, overwrite_x=1, overwrite_y=1)
        assert np.shares_memory(x, xo)
        assert np.shares_memory(y, yo)

    @pytest.mark.parametrize('dtype', REAL_DTYPES)
    def test_sysvx_writes_af_and_ipiv_in_place(self, dtype):
        # The expert drivers take `af` and `ipiv` as intent(in,out) with no `copy` and
        # no overwrite_ guard of their own, so a supplied buffer is written through.
        # That is what the .pyf declared; this pins it so it cannot drift silently.
        n = 4
        a = np.asfortranarray(np.eye(n, dtype=dtype) * 4 + 1)
        b = np.asfortranarray(np.ones((n, 1), dtype=dtype))
        sysvx = get_lapack_funcs('sysvx', dtype=dtype)

        af = np.zeros((n, n), dtype=dtype, order='F')
        ipiv = np.zeros(n, dtype=np.int32 if not _is_ilp64(sysvx) else np.int64)
        out = sysvx(a, b, af=af, ipiv=ipiv)
        af_out, ipiv_out = out[1], out[2]
        assert np.shares_memory(af, af_out)
        assert np.shares_memory(ipiv, ipiv_out)


def _is_ilp64(func):
    return func.int_dtype == np.dtype(np.int64)


# --------------------------------------------------------------------------------------
# Workspace queries.
# --------------------------------------------------------------------------------------
class TestLworkQueries:

    @pytest.mark.parametrize('dtype', DTYPES)
    @pytest.mark.parametrize('shape', [(6, 4), (4, 6), (5, 5)])
    def test_geqrf_lwork_is_accepted_by_geqrf(self, dtype, shape):
        m, n = shape
        rng = np.random.default_rng(21)
        a = np.asfortranarray(_vector(m * n, dtype, rng).reshape(m, n))
        geqrf, geqrf_lwork = get_lapack_funcs(('geqrf', 'geqrf_lwork'), dtype=dtype)

        work, info = geqrf_lwork(m, n)
        assert info == 0
        lwork = int(np.real(work))
        assert lwork >= 1

        qr_q, tau_q, w_q, info_q = geqrf(a, lwork=lwork)
        qr_d, tau_d, w_d, info_d = geqrf(a)
        assert info_q == 0 and info_d == 0
        # The factorization must not depend on how the workspace was sized.
        assert_allclose(qr_q, qr_d, rtol=TOL[dtype])
        assert_allclose(tau_q, tau_d, rtol=TOL[dtype])

    @pytest.mark.parametrize('dtype', DTYPES)
    def test_gesdd_lwork_is_accepted_by_gesdd(self, dtype):
        m, n = 6, 4
        rng = np.random.default_rng(23)
        a = np.asfortranarray(_vector(m * n, dtype, rng).reshape(m, n))
        gesdd, gesdd_lwork = get_lapack_funcs(('gesdd', 'gesdd_lwork'), dtype=dtype)

        work, info = gesdd_lwork(m, n)
        assert info == 0
        lwork = int(np.real(work))

        out_q = gesdd(a, lwork=lwork)
        out_d = gesdd(a)
        assert out_q[-1] == 0 and out_d[-1] == 0
        assert_allclose(out_q[1], out_d[1], rtol=TOL[dtype])   # singular values

    @pytest.mark.parametrize('dtype', DTYPES)
    def test_geqrf_inline_query_matches_lwork_helper(self, dtype):
        # `lwork=-1` on the routine itself must report the same size the *_lwork
        # helper does; they call the same LAPACK entry point in query mode.
        m, n = 7, 3
        rng = np.random.default_rng(29)
        a = np.asfortranarray(_vector(m * n, dtype, rng).reshape(m, n))
        geqrf, geqrf_lwork = get_lapack_funcs(('geqrf', 'geqrf_lwork'), dtype=dtype)

        helper, info_h = geqrf_lwork(m, n)
        _, _, work, info_i = geqrf(a, lwork=-1)
        assert info_h == 0 and info_i == 0
        assert int(np.real(helper)) == int(np.real(work[0]))


# --------------------------------------------------------------------------------------
# gees / gges eigenvalue-sort callbacks.
# --------------------------------------------------------------------------------------
@pytest.mark.skipif(flapack is None, reason='_flapack not built')
class TestSortCallbacks:

    @staticmethod
    def _matrix(dtype, n=4):
        rng = np.random.default_rng(31)
        a = rng.uniform(size=(n, n))
        if np.issubdtype(dtype, np.complexfloating):
            a = a + 1j * rng.uniform(size=(n, n))
        return np.asfortranarray(a.astype(dtype))

    @pytest.mark.parametrize('dtype', DTYPES)
    def test_callback_arity_variants(self, dtype):
        # LAPACK offers two scalars for a real gees and one for a complex one; f2py
        # passed only as many as the callable declared, and so must the C++ trampoline.
        gees = get_lapack_funcs('gees', dtype=dtype)
        a = self._matrix(dtype)

        seen = []

        def zero_arg():
            seen.append(0)
            return True

        def one_arg(w):
            seen.append(1)
            return True

        for fn in (zero_arg, one_arg):
            seen.clear()
            out = gees(fn, a, sort_t=1)
            assert out[-1] == 0
            assert seen, 'the callback was never invoked'

    @pytest.mark.parametrize('dtype', REAL_DTYPES)
    def test_real_callback_may_take_both_scalars(self, dtype):
        gees = get_lapack_funcs('gees', dtype=dtype)
        a = self._matrix(dtype)
        recorded = []

        def two_arg(wr, wi):
            recorded.append((wr, wi))
            return wr >= 0

        out = gees(two_arg, a, sort_t=1)
        assert out[-1] == 0
        assert recorded
        assert all(isinstance(p, float) and isinstance(q, float) for p, q in recorded)

    @pytest.mark.parametrize('dtype', DTYPES)
    def test_callback_needing_too_many_arguments_is_refused(self, dtype):
        # Checked up front, before the Fortran call starts -- as f2py did, and even
        # when sort_t=0 means the callback would never actually be invoked.
        gees = get_lapack_funcs('gees', dtype=dtype)
        a = self._matrix(dtype)

        def greedy(a1, a2, a3, a4, a5):
            return True

        with assert_raises(TypeError):
            gees(greedy, a, sort_t=1)
        with assert_raises(TypeError):
            gees(greedy, a, sort_t=0)

    @pytest.mark.parametrize('dtype', DTYPES)
    def test_non_callable_select_is_refused(self, dtype):
        gees = get_lapack_funcs('gees', dtype=dtype)
        a = self._matrix(dtype)
        for bad in (None, 42, 'nope'):
            with assert_raises(TypeError):
                gees(bad, a, sort_t=1)

    @pytest.mark.parametrize('dtype', DTYPES)
    def test_exception_in_callback_propagates(self, dtype):
        # The original exception must reach the caller rather than being swallowed
        # into a "do not select" result.
        gees = get_lapack_funcs('gees', dtype=dtype)
        a = self._matrix(dtype)

        class Boom(Exception):
            pass

        def raiser(*args):
            raise Boom('from the callback')

        with assert_raises(Boom):
            gees(raiser, a, sort_t=1)

    @pytest.mark.parametrize('dtype', DTYPES)
    def test_callback_result_with_raising_bool_propagates(self, dtype):
        gees = get_lapack_funcs('gees', dtype=dtype)
        a = self._matrix(dtype)

        class NoTruth:
            def __bool__(self):
                raise ValueError('no truth value')

        with assert_raises(ValueError):
            gees(lambda *args: NoTruth(), a, sort_t=1)

    @pytest.mark.parametrize('dtype', DTYPES)
    def test_extra_args_are_appended(self, dtype):
        gees = get_lapack_funcs('gees', dtype=dtype)
        a = self._matrix(dtype)
        extra_kw = gees.typecode + 'select_extra_args'
        captured = []

        def with_extra(w, threshold):
            captured.append(threshold)
            return True

        out = gees(with_extra, a, sort_t=1, **{extra_kw: (0.5,)})
        assert out[-1] == 0
        assert captured and all(t == 0.5 for t in captured)

    @pytest.mark.slow
    @pytest.mark.thread_unsafe(reason="sys.getrefcount sees other threads' references")
    @pytest.mark.parametrize('dtype', DTYPES)
    def test_callback_failure_does_not_leak_the_callable(self, dtype):
        gees = get_lapack_funcs('gees', dtype=dtype)
        a = self._matrix(dtype)

        def raiser(*args):
            raise RuntimeError

        gc.collect()
        before = sys.getrefcount(raiser)
        for _ in range(20):
            with assert_raises(RuntimeError):
                gees(raiser, a, sort_t=1)
        gc.collect()
        assert sys.getrefcount(raiser) == before


# --------------------------------------------------------------------------------------
# Argument protocol: parse_args reproduces CPython's checks, order and wording.
# --------------------------------------------------------------------------------------
@pytest.mark.skipif(fblas is None, reason='_fblas not built')
class TestArgumentProtocol:

    @staticmethod
    def _ddot():
        # kwlist: x, y, n, offx, incx, offy, incy   (2 required)
        return fblas.ddot, np.arange(3.0), np.arange(3.0)

    def test_too_many_arguments(self):
        ddot, x, y = self._ddot()
        with assert_raises(TypeError, match=r'takes at most 7 arguments \(8 given\)'):
            ddot(x, y, 3, 0, 1, 0, 1, 99)

    def test_missing_required_argument(self):
        ddot, x, _ = self._ddot()
        with assert_raises(TypeError, match=r"missing required argument 'y' \(pos 2\)"):
            ddot(x)

    def test_given_by_name_and_position(self):
        ddot, x, y = self._ddot()
        msg = r"given by name \('x'\) and position \(1\)"
        with assert_raises(TypeError, match=msg):
            ddot(x, y, x=x)

    def test_unexpected_keyword(self):
        ddot, x, y = self._ddot()
        with assert_raises(TypeError, match=r"unexpected keyword argument 'nope'"):
            ddot(x, y, nope=1)

    def test_missing_beats_clash(self):
        # CPython validates in separate passes, so a clash at an earlier index does
        # not pre-empt a missing-required at a later one: reports 'y', not the clash.
        ddot, x, _ = self._ddot()
        with assert_raises(TypeError, match=r"missing required argument 'y'"):
            ddot(x, x=x)

    def test_too_many_beats_everything(self):
        # The count is positional + keyword total, so an unexpected keyword that pushes
        # the total over the limit reports "takes at most", not "unexpected keyword".
        ddot, x, y = self._ddot()
        with assert_raises(TypeError, match=r'takes at most 7 arguments \(8 given\)'):
            ddot(x, y, 3, 0, 1, 0, 1, nope=1)


# --------------------------------------------------------------------------------------
# Scalar coercion: the deliberately permissive ports of f2py's *_from_pyobj.
# --------------------------------------------------------------------------------------
@pytest.mark.skipif(fblas is None, reason='_fblas not built')
class TestScalarCoercion:

    def setup_method(self):
        self.x = np.arange(1.0, 8.0)
        self.asum = fblas.dasum

    def test_float_truncates_to_int(self):
        assert_allclose(self.asum(self.x, n=3.7), self.asum(self.x, n=3))

    def test_complex_contributes_its_real_part(self):
        assert_allclose(self.asum(self.x, n=complex(3, 9)), self.asum(self.x, n=3))

    def test_sequence_contributes_element_zero(self):
        assert_allclose(self.asum(self.x, n=[3, 99]), self.asum(self.x, n=3))
        assert_allclose(self.asum(self.x, n=np.array([3, 99])), self.asum(self.x, n=3))

    def test_numeric_string_is_accepted(self):
        # PyNumber_Long accepts it, and so did f2py.
        assert_allclose(self.asum(self.x, n='3'), self.asum(self.x, n=3))

    @pytest.mark.skipif(IS_WASM, reason="overflows the smaller WASM stack before "
                                         "the recursion limit is reached")
    def test_self_referential_sequence_raises_recursion_error(self):
        # Pure-C descent into element 0 would otherwise overflow the C stack; the
        # converters opt into the interpreter's depth accounting instead.
        bad = []
        bad.append(bad)
        with assert_raises(RecursionError):
            self.asum(self.x, n=bad)

    def test_unconvertible_scalar_raises(self):
        with assert_raises((TypeError, ValueError)):
            self.asum(self.x, n=object())

    def test_overwrite_flag_requires_a_true_integer(self):
        # The overwrite_* flags use __index__, not the permissive path: a float would
        # otherwise be coerced and silently change the meaning of the call.
        rot = fblas.drot
        x, y = np.arange(3.0), np.arange(3.0)
        for bad in (1.0, '1', None):
            with assert_raises(TypeError):
                rot(x, y, 0.6, 0.8, overwrite_x=bad)
        rot(x, y, 0.6, 0.8, overwrite_x=True)     # bool is an integer
        rot(x, y, 0.6, 0.8, overwrite_x=np.int64(1))

    @pytest.mark.skipif(flapack is None, reason='_flapack not built')
    def test_option_letter_accepts_str_bytes_and_numpy_str(self):
        a = np.asfortranarray(np.eye(3) * 2 + 1)
        lange = flapack.dlange
        reference = lange('1', a)
        for spelling in (b'1', np.str_('1'), ['1']):
            assert_allclose(lange(spelling, a), reference)

    @pytest.mark.skipif(flapack is None, reason='_flapack not built')
    def test_option_letter_rejects_a_number(self):
        # Unlike the numeric converters, the character one takes no numbers.
        a = np.asfortranarray(np.eye(3) * 2 + 1)
        with assert_raises((TypeError, ValueError)):
            flapack.dlange(1, a)


# --------------------------------------------------------------------------------------
# Rank reinterpretation (fix_rank): f2py reinterpreted rank mismatches rather than
# rejecting them, and code in the wild relies on it.
# --------------------------------------------------------------------------------------
@pytest.mark.skipif(flapack is None, reason='_flapack not built')
class TestRankReinterpretation:

    def test_vector_is_treated_as_a_column(self):
        # A 1-D `b` where (n, nrhs) is declared becomes an (n, 1) column.
        a = np.asfortranarray(np.eye(3) * 2 + 1)
        b = np.array([1.0, 2.0, 3.0])
        lu, piv, x, info = flapack.dgesv(a, b)
        assert info == 0
        assert x.shape == b.shape          # the caller-shaped original is handed back
        assert_allclose(a @ x.reshape(3), b, rtol=1e-12)

    def test_unit_axes_are_squeezed(self):
        # A (1, 3, 1) array where a rank-2 is declared collapses to (3, 1).
        a = np.asfortranarray(np.eye(3) * 2 + 1)
        b = np.array([1.0, 2.0, 3.0]).reshape(1, 3, 1)
        lu, piv, x, info = flapack.dgesv(a, b)
        assert info == 0
        assert_allclose(a @ x.reshape(3), b.reshape(3), rtol=1e-12)

    def test_scalar_broadcasts_to_a_1x1(self):
        # The existing suite relies on `dgemm(3, [3], [-4])` working at all.
        out = fblas.dgemm(3, [3], [-4])
        assert_allclose(np.asarray(out).reshape(-1), [-36.0])


# --------------------------------------------------------------------------------------
# Zero-sized inputs.
# --------------------------------------------------------------------------------------
@pytest.mark.skipif(flapack is None, reason='_flapack not built')
class TestZeroSized:

    @pytest.mark.parametrize('dtype', DTYPES)
    def test_gesv_on_an_empty_system(self, dtype):
        a = np.zeros((0, 0), dtype=dtype, order='F')
        b = np.zeros((0, 1), dtype=dtype, order='F')
        gesv = get_lapack_funcs('gesv', dtype=dtype)
        lu, piv, x, info = gesv(a, b)
        assert info == 0
        assert lu.shape == (0, 0)
        assert x.shape == (0, 1)

    @pytest.mark.parametrize('dtype', DTYPES)
    def test_getrf_on_an_empty_matrix(self, dtype):
        a = np.zeros((0, 0), dtype=dtype, order='F')
        getrf = get_lapack_funcs('getrf', dtype=dtype)
        lu, piv, info = getrf(a)
        assert info == 0
        assert lu.shape == (0, 0)
        assert piv.shape == (0,)


# --------------------------------------------------------------------------------------
# Smoke pass over every exposed routine.
# --------------------------------------------------------------------------------------
def _exposed(module):
    return sorted(n for n in dir(module) if not n.startswith('_'))


@pytest.mark.parametrize('modname', ['_fblas', '_flapack'])
class TestEveryRoutine:

    @staticmethod
    def _module(modname):
        mod = {'_fblas': fblas, '_flapack': flapack}[modname]
        if mod is None:
            pytest.skip(f'{modname} not built')
        return mod

    def test_names_are_callable_and_self_describing(self, modname):
        mod = self._module(modname)
        names = _exposed(mod)
        assert len(names) > 100, 'suspiciously few routines exposed'
        for name in names:
            func = getattr(mod, name)
            assert callable(func), name
            assert func.__name__ == name

    def test_every_routine_has_a_docstring(self, modname):
        # f2py generated one per routine from the .pyf declaration; the replacements are
        # written by hand, so a new routine can easily be added without one.
        mod = self._module(modname)
        missing = [n for n in _exposed(mod)
                   if not isinstance(getattr(mod, n).__doc__, str)
                   or not getattr(mod, n).__doc__.strip()]
        assert missing == [], f'routines with no docstring: {missing}'

    def test_docstring_starts_with_the_signature(self, modname):
        mod = self._module(modname)
        bad = [n for n in _exposed(mod)
               if not getattr(mod, n).__doc__.startswith(n + '(')]
        assert bad == [], f'docstring does not open with its own signature: {bad}'

    def test_calling_with_a_bogus_keyword_raises_type_error(self, modname):
        mod = self._module(modname)
        for name in _exposed(mod):
            with assert_raises(TypeError):
                getattr(mod, name)(definitely_not_an_argument=1)


# --------------------------------------------------------------------------------------
# Reference counting.
# --------------------------------------------------------------------------------------
@pytest.mark.slow
@pytest.mark.thread_unsafe(reason="sys.getrefcount sees other threads' references")
@pytest.mark.skipif(flapack is None, reason='_flapack not built')
@pytest.mark.parametrize('dtype', [np.float64, np.complex128])
def test_repeated_calls_do_not_leak_inputs(dtype):
    # ~650 hand-written acquisition/release paths; a missing DECREF on a success path
    # shows up here as a steadily climbing refcount on the caller's arrays.
    a = np.asfortranarray(np.eye(5, dtype=dtype) * 3 + 1)
    b = np.asfortranarray(np.ones((5, 2), dtype=dtype))
    gesv = get_lapack_funcs('gesv', dtype=dtype)

    gesv(a, b)   # warm up any lazily built state
    gc.collect()
    before_a, before_b = sys.getrefcount(a), sys.getrefcount(b)
    for _ in range(50):
        gesv(a, b)
    gc.collect()
    assert sys.getrefcount(a) == before_a
    assert sys.getrefcount(b) == before_b


@pytest.mark.slow
@pytest.mark.thread_unsafe(reason="sys.getrefcount sees other threads' references")
@pytest.mark.skipif(flapack is None, reason='_flapack not built')
def test_failed_calls_do_not_leak_inputs():
    a = np.asfortranarray(np.eye(5) * 3 + 1)
    b = np.asfortranarray(np.ones((4, 2)))   # wrong shape: rejected by a CHECKARRAY
    gc.collect()
    before_a, before_b = sys.getrefcount(a), sys.getrefcount(b)
    for _ in range(50):
        with assert_raises(ValueError):
            flapack.dgesv(a, b)
    gc.collect()
    assert sys.getrefcount(a) == before_a
    assert sys.getrefcount(b) == before_b
