import queue
import threading
import multiprocessing
import numpy as np
import pytest
from numpy.random import random
from numpy.testing import (
    assert_array_almost_equal, assert_array_equal, assert_allclose
)
from pytest import raises as assert_raises
import scipy.fft as fft
from scipy.conftest import (
    array_api_compatible,
    skip_if_array_api,
    skip_if_array_api_gpu,
    skip_if_array_api_backend,
    set_assert_allclose
)
from scipy._lib._array_api import (
    _assert_matching_namespace,
    array_namespace,
    size
)


def fft1(x):
    L = len(x)
    phase = -2j*np.pi*(np.arange(L)/float(L))
    phase = np.arange(L).reshape(-1, 1) * phase
    return np.sum(x*np.exp(phase), axis=1)


class TestFFTShift:

    @array_api_compatible
    def test_fft_n(self, xp):
        x = xp.asarray([1, 2, 3])
        if xp.__name__ == 'torch':
            assert_raises(RuntimeError, fft.fft, x, 0)
        else:
            assert_raises(ValueError, fft.fft, x, 0)


class TestFFT1D:

    @array_api_compatible
    def test_identity(self, xp):
        maxlen = 512
        x = xp.asarray(random(maxlen) + 1j*random(maxlen))
        xr = xp.asarray(random(maxlen))
        _assert_allclose = set_assert_allclose(xp)
        for i in range(1, maxlen):
            _assert_allclose(fft.ifft(fft.fft(x[0:i])), x[0:i],
                             rtol=1e-9, atol=0)
            _assert_allclose(fft.irfft(fft.rfft(xr[0:i]), i), xr[0:i],
                             rtol=1e-9, atol=0)

    @array_api_compatible
    def test_fft(self, xp):
        x = random(30) + 1j*random(30)
        expect = xp.asarray(fft1(x))
        x = xp.asarray(x)
        _assert_allclose = set_assert_allclose(xp)
        _assert_allclose(expect, fft.fft(x))
        _assert_allclose(expect, fft.fft(x, norm="backward"))
        _assert_allclose(expect / xp.sqrt(xp.asarray(30, dtype=xp.float64)),
                         fft.fft(x, norm="ortho"))
        _assert_allclose(expect / 30, fft.fft(x, norm="forward"))

    @array_api_compatible
    def test_ifft(self, xp):
        x = xp.asarray(random(30) + 1j*random(30))
        _assert_allclose = set_assert_allclose(xp)
        _assert_allclose(x, fft.ifft(fft.fft(x)))
        for norm in ["backward", "ortho", "forward"]:
            _assert_allclose(x, fft.ifft(fft.fft(x, norm=norm), norm=norm))

    @skip_if_array_api_gpu
    @array_api_compatible
    def test_fft2(self, xp):
        x = xp.asarray(random((30, 20)) + 1j*random((30, 20)))
        expect = fft.fft(fft.fft(x, axis=1), axis=0)
        assert_allclose(expect, fft.fft2(x))
        assert_allclose(expect, fft.fft2(x, norm="backward"))
        assert_allclose(
            expect / xp.sqrt(xp.asarray(30 * 20, dtype=xp.float64)),
            fft.fft2(x, norm="ortho")
        )
        assert_allclose(expect / (30 * 20), fft.fft2(x, norm="forward"))

    @skip_if_array_api_gpu
    @array_api_compatible
    def test_ifft2(self, xp):
        x = xp.asarray(random((30, 20)) + 1j*random((30, 20)))
        expect = fft.ifft(fft.ifft(x, axis=1), axis=0)
        assert_allclose(expect, fft.ifft2(x))
        assert_allclose(expect, fft.ifft2(x, norm="backward"))
        assert_allclose(
            expect * xp.sqrt(xp.asarray(30 * 20, dtype=xp.float64)),
            fft.ifft2(x, norm="ortho")
        )
        assert_allclose(expect * (30 * 20), fft.ifft2(x, norm="forward"))

    @array_api_compatible
    @skip_if_array_api_backend('torch')
    def test_fftn(self, xp):
        x = xp.asarray(random((30, 20, 10)) + 1j*random((30, 20, 10)))
        _assert_allclose = set_assert_allclose(xp)
        expect = fft.fft(fft.fft(fft.fft(x, axis=2), axis=1), axis=0)
        _assert_allclose(expect, fft.fftn(x))
        _assert_allclose(expect, fft.fftn(x, norm="backward"))
        _assert_allclose(
            expect / xp.sqrt(xp.asarray(30 * 20 * 10, dtype=xp.float64)),
            fft.fftn(x, norm="ortho")
        )
        _assert_allclose(expect / (30 * 20 * 10), fft.fftn(x, norm="forward"))

    @array_api_compatible
    @skip_if_array_api_backend('torch')
    def test_ifftn(self, xp):
        x = xp.asarray(random((30, 20, 10)) + 1j*random((30, 20, 10)))
        _assert_allclose = set_assert_allclose(xp)
        expect = fft.ifft(fft.ifft(fft.ifft(x, axis=2), axis=1), axis=0)
        _assert_allclose(expect, fft.ifftn(x))
        _assert_allclose(expect, fft.ifftn(x, norm="backward"))
        _assert_allclose(
            fft.ifftn(x) * xp.sqrt(xp.asarray(30 * 20 * 10, dtype=xp.float64)),
            fft.ifftn(x, norm="ortho")
        )
        _assert_allclose(expect * (30 * 20 * 10), fft.ifftn(x, norm="forward"))

    @array_api_compatible
    def test_rfft(self, xp):
        x = xp.asarray(random(29))
        _assert_allclose = set_assert_allclose(xp)
        for n in [size(x), 2*size(x)]:
            for norm in [None, "backward", "ortho", "forward"]:
                _assert_allclose(fft.fft(x, n=n, norm=norm)[:(n//2 + 1)],
                                 fft.rfft(x, n=n, norm=norm))
            _assert_allclose(
                fft.rfft(x, n=n) / xp.sqrt(xp.asarray(n, dtype=xp.float64)),
                fft.rfft(x, n=n, norm="ortho")
            )

    @array_api_compatible
    def test_irfft(self, xp):
        x = xp.asarray(random(30))
        _assert_allclose = set_assert_allclose(xp)
        _assert_allclose(x, fft.irfft(fft.rfft(x)))
        for norm in ["backward", "ortho", "forward"]:
            _assert_allclose(x, fft.irfft(fft.rfft(x, norm=norm), norm=norm))

    @skip_if_array_api_gpu
    @array_api_compatible
    def test_rfft2(self, xp):
        x = xp.asarray(random((30, 20)))
        expect = fft.fft2(x)[:, :11]
        assert_allclose(expect, fft.rfft2(x))
        assert_allclose(expect, fft.rfft2(x, norm="backward"))
        assert_allclose(
            expect / xp.sqrt(xp.asarray(30 * 20, dtype=xp.float64)),
            fft.rfft2(x, norm="ortho")
        )
        assert_allclose(expect / (30 * 20), fft.rfft2(x, norm="forward"))

    @skip_if_array_api_gpu
    @array_api_compatible
    def test_irfft2(self, xp):
        x = xp.asarray(random((30, 20)))
        assert_allclose(x, fft.irfft2(fft.rfft2(x)))
        for norm in ["backward", "ortho", "forward"]:
            assert_allclose(x, fft.irfft2(fft.rfft2(x, norm=norm), norm=norm))

    @array_api_compatible
    @skip_if_array_api_backend('torch')
    def test_rfftn(self, xp):
        x = xp.asarray(random((30, 20, 10)))
        _assert_allclose = set_assert_allclose(xp)
        expect = fft.fftn(x)[:, :, :6]
        _assert_allclose(expect, fft.rfftn(x))
        _assert_allclose(expect, fft.rfftn(x, norm="backward"))
        _assert_allclose(
            expect / xp.sqrt(xp.asarray(30 * 20 * 10, dtype=xp.float64)),
            fft.rfftn(x, norm="ortho")
        )
        _assert_allclose(expect / (30 * 20 * 10), fft.rfftn(x, norm="forward"))

    @array_api_compatible
    @skip_if_array_api_backend('torch')
    def test_irfftn(self, xp):
        x = xp.asarray(random((30, 20, 10)))
        _assert_allclose = set_assert_allclose(xp)
        _assert_allclose(x, fft.irfftn(fft.rfftn(x)))
        for norm in ["backward", "ortho", "forward"]:
            _assert_allclose(x, fft.irfftn(fft.rfftn(x, norm=norm), norm=norm))

    @array_api_compatible
    def test_hfft(self, xp):
        x = random(14) + 1j*random(14)
        x_herm = np.concatenate((random(1), x, random(1)))
        x = np.concatenate((x_herm, x[::-1].conj()))
        x = xp.asarray(x)
        x_herm = xp.asarray(x_herm)
        _assert_allclose = set_assert_allclose(xp)
        expect = fft.fft(x)
        xp_test = array_namespace(x)
        _assert_allclose(expect,
                         xp_test.asarray(fft.hfft(x_herm),
                                         dtype=xp.complex128))
        _assert_allclose(expect,
                         xp_test.asarray(fft.hfft(x_herm, norm="backward"),
                                         dtype=xp.complex128))
        _assert_allclose(expect / xp.sqrt(xp.asarray(30, dtype=xp.float64)),
                         xp_test.asarray(fft.hfft(x_herm, norm="ortho"),
                                         dtype=xp.complex128))
        _assert_allclose(expect / 30,
                         xp_test.asarray(fft.hfft(x_herm, norm="forward"),
                                         dtype=xp.complex128))

    @array_api_compatible
    def test_ihfft(self, xp):
        x = random(14) + 1j*random(14)
        x_herm = np.concatenate((random(1), x, random(1)))
        x = np.concatenate((x_herm, x[::-1].conj()))
        x = xp.asarray(x)
        x_herm = xp.asarray(x_herm)
        _assert_allclose = set_assert_allclose(xp)
        _assert_allclose(x_herm, fft.ihfft(fft.hfft(x_herm)))
        for norm in ["backward", "ortho", "forward"]:
            _assert_allclose(x_herm,
                             fft.ihfft(fft.hfft(x_herm, norm=norm), norm=norm))

    @skip_if_array_api_gpu
    @array_api_compatible
    def test_hfft2(self, xp):
        x = xp.asarray(random((30, 20)))
        assert_allclose(x, fft.hfft2(fft.ihfft2(x)))
        for norm in ["backward", "ortho", "forward"]:
            assert_allclose(x, fft.hfft2(fft.ihfft2(x, norm=norm), norm=norm))

    @skip_if_array_api_gpu
    @array_api_compatible
    def test_ihfft2(self, xp):
        x = xp.asarray(random((30, 20)))
        expect = fft.ifft2(x)[:, :11]
        assert_allclose(expect, fft.ihfft2(x))
        assert_allclose(expect, fft.ihfft2(x, norm="backward"))
        assert_allclose(
            expect * xp.sqrt(xp.asarray(30 * 20, dtype=xp.float64)),
            fft.ihfft2(x, norm="ortho")
        )
        assert_allclose(expect * (30 * 20), fft.ihfft2(x, norm="forward"))

    @skip_if_array_api_gpu
    @array_api_compatible
    def test_hfftn(self, xp):
        x = xp.asarray(random((30, 20, 10)))
        assert_allclose(x, fft.hfftn(fft.ihfftn(x)))
        for norm in ["backward", "ortho", "forward"]:
            assert_allclose(x, fft.hfftn(fft.ihfftn(x, norm=norm), norm=norm))

    @skip_if_array_api_gpu
    @array_api_compatible
    @skip_if_array_api_backend('torch')
    def test_ihfftn(self, xp):
        x = xp.asarray(random((30, 20, 10)))
        expect = fft.ifftn(x)[:, :, :6]
        assert_allclose(expect, fft.ihfftn(x))
        assert_allclose(expect, fft.ihfftn(x, norm="backward"))
        assert_allclose(
            expect * xp.sqrt(xp.asarray(30 * 20 * 10, dtype=xp.float64)),
            fft.ihfftn(x, norm="ortho")
        )
        assert_allclose(expect * (30 * 20 * 10), fft.ihfftn(x, norm="forward"))

    @array_api_compatible
    @pytest.mark.parametrize("op", [fft.fftn, fft.ifftn,
                                    fft.rfftn, fft.irfftn])
    def test_axes_standard(self, op, xp):
        if xp.__name__ == 'torch':
            # @skip_if_array_api_backend doesn't seem to work
            # with @pytest.mark.parametrize
            pytest.xfail("torch.fft not yet implemented by array-api-compat")
        x = xp.asarray(random((30, 20, 10)))
        axes = [(0, 1, 2), (0, 2, 1), (1, 0, 2),
                (1, 2, 0), (2, 0, 1), (2, 1, 0)]
        _assert_allclose = set_assert_allclose(xp)
        xp_test = array_namespace(x)
        for a in axes:
            op_tr = op(xp_test.permute_dims(x, axes=a))
            tr_op = xp_test.permute_dims(op(x, axes=a), axes=a)
            _assert_allclose(op_tr, tr_op)

    @skip_if_array_api_gpu
    @array_api_compatible
    @pytest.mark.parametrize("op", [fft.hfftn, fft.ihfftn])
    def test_axes_non_standard(self, op, xp):
        x = xp.asarray(random((30, 20, 10)))
        axes = [(0, 1, 2), (0, 2, 1), (1, 0, 2),
                (1, 2, 0), (2, 0, 1), (2, 1, 0)]
        xp_test = array_namespace(x)
        for a in axes:
            op_tr = op(xp_test.permute_dims(x, axes=a))
            tr_op = xp_test.permute_dims(op(x, axes=a), axes=a)
            assert_allclose(op_tr, tr_op)

    @array_api_compatible
    @pytest.mark.parametrize("op", [fft.fftn, fft.ifftn,
                                    fft.rfftn, fft.irfftn])
    def test_axes_subset_with_shape_standard(self, op, xp):
        if xp.__name__ == 'torch':
            # @skip_if_array_api_backend doesn't seem to work
            # with @pytest.mark.parametrize
            pytest.xfail("torch.fft not yet implemented by array-api-compat")
        x = xp.asarray(random((16, 8, 4)))
        axes = [(0, 1, 2), (0, 2, 1), (1, 2, 0)]
        _assert_allclose = set_assert_allclose(xp)
        xp_test = array_namespace(x)
        for a in axes:
            # different shape on the first two axes
            shape = tuple([2*x.shape[ax] if ax in a[:2] else x.shape[ax]
                           for ax in range(x.ndim)])
            # transform only the first two axes
            op_tr = op(xp_test.permute_dims(x, axes=a),
                       s=shape[:2], axes=(0, 1))
            tr_op = xp_test.permute_dims(op(x, s=shape[:2], axes=a[:2]),
                                         axes=a)
            _assert_allclose(op_tr, tr_op)

    @skip_if_array_api_gpu
    @array_api_compatible
    @pytest.mark.parametrize("op", [fft.fft2, fft.ifft2,
                                    fft.rfft2, fft.irfft2,
                                    fft.hfft2, fft.ihfft2,
                                    fft.hfftn, fft.ihfftn])
    def test_axes_subset_with_shape_non_standard(self, op, xp):
        x = xp.asarray(random((16, 8, 4)))
        axes = [(0, 1, 2), (0, 2, 1), (1, 2, 0)]
        xp_test = array_namespace(x)
        for a in axes:
            # different shape on the first two axes
            shape = tuple([2*x.shape[ax] if ax in a[:2] else x.shape[ax]
                           for ax in range(x.ndim)])
            # transform only the first two axes
            op_tr = op(xp_test.permute_dims(x, axes=a),
                       s=shape[:2], axes=(0, 1))
            tr_op = xp_test.permute_dims(op(x, s=shape[:2], axes=a[:2]),
                                         axes=a)
            assert_allclose(op_tr, tr_op)

    @skip_if_array_api_gpu
    @array_api_compatible
    def test_all_1d_norm_preserving(self, xp):
        # verify that round-trip transforms are norm-preserving
        x = xp.asarray(random(30))
        xp_test = array_namespace(x)
        x_norm = xp_test.linalg.vector_norm(x)
        n = size(x) * 2
        func_pairs = [(fft.fft, fft.ifft),
                      (fft.rfft, fft.irfft),
                      # hfft: order so the first function takes x.size samples
                      #       (necessary for comparison to x_norm above)
                      (fft.ihfft, fft.hfft),
                      ]
        for forw, back in func_pairs:
            for n in [size(x), 2*size(x)]:
                for norm in ['backward', 'ortho', 'forward']:
                    tmp = forw(x, n=n, norm=norm)
                    tmp = back(tmp, n=n, norm=norm)
                    assert_allclose(x_norm, xp_test.linalg.vector_norm(tmp))

    @skip_if_array_api
    def test_float16(self):
        x = random(30).astype(np.float16)
        assert_array_almost_equal(fft.ifft(fft.fft(x)), x)
        assert_array_almost_equal(fft.irfft(fft.rfft(x)), x)
        assert_array_almost_equal(fft.hfft(fft.ihfft(x), len(x)), x)

    @array_api_compatible
    def test_float32(self, xp):
        x = xp.asarray(random(30), dtype=xp.float32)
        _assert_allclose = set_assert_allclose(xp)
        _assert_allclose(fft.ifft(fft.fft(x)),
                         xp.asarray(x, dtype=xp.complex64),
                         rtol=1e-4, atol=0)
        _assert_allclose(fft.irfft(fft.rfft(x)), x,
                         rtol=1e-4, atol=0)
        _assert_allclose(fft.hfft(fft.ihfft(x), x.shape[0]), x,
                         rtol=1e-4, atol=0)

    @array_api_compatible
    def test_float64(self, xp):
        x = xp.asarray(random(30), dtype=xp.float64)
        _assert_allclose = set_assert_allclose(xp)
        _assert_allclose(fft.ifft(fft.fft(x)),
                         xp.asarray(x, dtype=xp.complex128))
        _assert_allclose(fft.irfft(fft.rfft(x)), x)
        _assert_allclose(fft.hfft(fft.ihfft(x), x.shape[0]), x)


@skip_if_array_api
@pytest.mark.parametrize(
    "dtype",
    [np.float32, np.float64, np.longfloat,
     np.complex64, np.complex128, np.longcomplex])
@pytest.mark.parametrize("order", ["F", 'non-contiguous'])
@pytest.mark.parametrize(
    "fft",
    [fft.fft, fft.fft2, fft.fftn,
     fft.ifft, fft.ifft2, fft.ifftn])
def test_fft_with_order(dtype, order, fft):
    # Check that FFT/IFFT produces identical results for C, Fortran and
    # non contiguous arrays
    rng = np.random.RandomState(42)
    X = rng.rand(8, 7, 13).astype(dtype, copy=False)
    if order == 'F':
        Y = np.asfortranarray(X)
    else:
        # Make a non contiguous array
        Y = X[::-1]
        X = np.ascontiguousarray(X[::-1])

    if fft.__name__.endswith('fft'):
        for axis in range(3):
            X_res = fft(X, axis=axis)
            Y_res = fft(Y, axis=axis)
            assert_array_almost_equal(X_res, Y_res)
    elif fft.__name__.endswith(('fft2', 'fftn')):
        axes = [(0, 1), (1, 2), (0, 2)]
        if fft.__name__.endswith('fftn'):
            axes.extend([(0,), (1,), (2,), None])
        for ax in axes:
            X_res = fft(X, axes=ax)
            Y_res = fft(Y, axes=ax)
            assert_array_almost_equal(X_res, Y_res)
    else:
        raise ValueError


class TestFFTThreadSafe:
    threads = 16
    input_shape = (800, 200)

    def _test_mtsame(self, func, *args, xp=None):
        def worker(args, q):
            q.put(func(*args))

        q = queue.Queue()
        expected = func(*args)

        # Spin off a bunch of threads to call the same function simultaneously
        t = [threading.Thread(target=worker, args=(args, q))
             for i in range(self.threads)]
        [x.start() for x in t]

        [x.join() for x in t]

        if xp is None:
            _assert_array_equal = assert_array_equal
        elif xp.__name__ == 'cupy':
            import cupy as cp
            _assert_array_equal = cp.testing.assert_array_equal
        elif xp.__name__ == 'torch':
            import torch
            for i in range(self.threads):
                torch.testing.assert_close(
                    q.get(timeout=5), expected,
                    msg='Function returned wrong value in multithreaded context'
                )
            return
        else:
            _assert_array_equal = assert_array_equal

        # Make sure all threads returned the correct value
        for i in range(self.threads):
            _assert_array_equal(
                q.get(timeout=5), expected,
                err_msg='Function returned wrong value in multithreaded context'
            )

    @array_api_compatible
    def test_fft(self, xp):
        a = xp.ones(self.input_shape, dtype=xp.complex128)
        self._test_mtsame(fft.fft, a, xp=xp)

    @array_api_compatible
    def test_ifft(self, xp):
        a = xp.full(self.input_shape, 1+0j)
        self._test_mtsame(fft.ifft, a, xp=xp)

    @array_api_compatible
    def test_rfft(self, xp):
        a = xp.ones(self.input_shape)
        self._test_mtsame(fft.rfft, a, xp=xp)

    @array_api_compatible
    def test_irfft(self, xp):
        a = xp.full(self.input_shape, 1+0j)
        self._test_mtsame(fft.irfft, a, xp=xp)

    @array_api_compatible
    def test_hfft(self, xp):
        a = xp.ones(self.input_shape, dtype=xp.complex64)
        self._test_mtsame(fft.hfft, a, xp=xp)

    @array_api_compatible
    def test_ihfft(self, xp):
        a = xp.ones(self.input_shape)
        self._test_mtsame(fft.ihfft, a, xp=xp)


@skip_if_array_api
@pytest.mark.parametrize("func", [fft.fft, fft.ifft, fft.rfft, fft.irfft])
def test_multiprocess(func):
    # Test that fft still works after fork (gh-10422)

    with multiprocessing.Pool(2) as p:
        res = p.map(func, [np.ones(100) for _ in range(4)])

    expect = func(np.ones(100))
    for x in res:
        assert_allclose(x, expect)


class TestIRFFTN:

    @array_api_compatible
    @skip_if_array_api_backend('torch')
    def test_not_last_axis_success(self, xp):
        ar, ai = np.random.random((2, 16, 8, 32))
        a = ar + 1j*ai
        a = xp.asarray(a)

        axes = (-2,)

        # Should not raise error
        fft.irfftn(a, axes=axes)


class TestNamespaces:

    @array_api_compatible
    @pytest.mark.parametrize("func", [fft.fft, fft.ifft])
    def test_fft_ifft(self, func, xp):
        x = xp.asarray(random(30) + 1j*random(30))
        _assert_matching_namespace(func(x), x)

    @array_api_compatible
    @pytest.mark.parametrize("func", [fft.fftn, fft.ifftn])
    def test_fftn_ifftn(self, func, xp):
        if xp.__name__ == 'torch':
            # @skip_if_array_api_backend doesn't seem to work
            # with @pytest.mark.parametrize
            pytest.xfail("torch.fft not yet implemented by array-api-compat")
        x = xp.asarray(random((30, 20, 10)) + 1j*random((30, 20, 10)))
        _assert_matching_namespace(func(x), x)

    @array_api_compatible
    def test_rfft(self, xp):
        x = xp.asarray(random(29))
        _assert_matching_namespace(fft.rfft(x), x)

    @array_api_compatible
    def test_irfft(self, xp):
        x = xp.asarray(random(30))
        _assert_matching_namespace(fft.irfft(x), x)

    @array_api_compatible
    @pytest.mark.parametrize("func", [fft.rfftn, fft.irfftn])
    def test_rfftn_irfftn(self, func, xp):
        if xp.__name__ == 'torch':
            # @skip_if_array_api_backend doesn't seem to work
            # with @pytest.mark.parametrize
            pytest.xfail("torch.fft not yet implemented by array-api-compat")
        x = xp.asarray(random((30, 20, 10)))
        _assert_matching_namespace(func(x), x)

    @array_api_compatible
    def test_hfft_ihfft(self, xp):
        x = random(14) + 1j*random(14)
        x_herm = np.concatenate((random(1), x, random(1)))
        x_herm = xp.asarray(x_herm)
        y = fft.hfft(x_herm)
        _assert_matching_namespace(y, x_herm)
        _assert_matching_namespace(fft.ihfft(y), y)


@array_api_compatible
@pytest.mark.parametrize("func", [fft.fft, fft.ifft, fft.rfft, fft.irfft,
                                  fft.fftn, fft.ifftn,
                                  fft.rfftn, fft.irfftn, fft.hfft, fft.ihfft])
def test_non_standard_params(func, xp):
    if xp.__name__ == 'torch':
        # @skip_if_array_api_backend doesn't seem to work
        # with @pytest.mark.parametrize
        pytest.xfail("torch.fft not yet implemented by array-api-compat")
    if xp.__name__ != 'numpy':
        x = xp.asarray([1, 2, 3])
        # func(x) should not raise an exception
        func(x)
        assert_raises(ValueError, func, x, overwrite_x=True)
        assert_raises(ValueError, func, x, workers=2)
        # `plan` param is not tested since SciPy does not use it currently
        # but should be tested if it comes into use
