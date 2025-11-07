import pytest
import numpy as np

from scipy import stats
from scipy.stats._quantile import (_xp_searchsorted, _iquantile_methods,
    _iquantile_discontinuous_methods, _iquantile_continuous_methods)
from scipy._lib._array_api import (
    xp_default_dtype,
    is_numpy,
    is_torch,
    is_jax,
    is_cupy,
    make_xp_test_case,
    SCIPY_ARRAY_API,
)
from scipy._lib._array_api_no_0d import xp_assert_close, xp_assert_equal
from scipy._lib._util import _apply_over_batch
import scipy._lib.array_api_extra as xpx
from scipy.stats._axis_nan_policy import _broadcast_arrays

skip_xp_backends = pytest.mark.skip_xp_backends

lazy_xp_modules = [stats]

@_apply_over_batch(('x', 1), ('p', 1))
def quantile_reference_last_axis(x, p, nan_policy, method):
    if nan_policy == 'omit':
        x = x[~np.isnan(x)]
    p_mask = np.isnan(p)
    p = p.copy()
    p[p_mask] = 0.5
    if method == 'harrell-davis':
        # hdquantiles returns masked element if length along axis is 1 (bug)
        res = (np.full_like(p, x[0]) if x.size == 1
               else stats.mstats.hdquantiles(x, p).data)
        if nan_policy == 'propagate' and np.any(np.isnan(x)):
            res[:] = np.nan
    else:
        res = np.quantile(x, p, method=method)
    res[p_mask] = np.nan
    return res


def quantile_reference(x, p, *, axis, nan_policy, keepdims, method):
    x, p = np.moveaxis(x, axis, -1), np.moveaxis(p, axis, -1)
    res = quantile_reference_last_axis(x, p, nan_policy, method)
    res = np.moveaxis(res, -1, axis)
    if not keepdims:
        res = np.squeeze(res, axis=axis)
    return res


@make_xp_test_case(stats.quantile)
class TestQuantile:

    def test_input_validation(self, xp):
        x = xp.asarray([1, 2, 3])
        p = xp.asarray(0.5)

        message = "`x` must have real dtype."
        with pytest.raises(ValueError, match=message):
            stats.quantile(xp.asarray([True, False]), p)
        with pytest.raises(ValueError):
            stats.quantile(xp.asarray([1+1j, 2]), p)

        message = "`p` must have real floating dtype."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, xp.asarray([0, 1]))

        message = "`axis` must be an integer or None."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, axis=0.5)
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, axis=(0, -1))

        message = "`axis` is not compatible with the shapes of the inputs."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, axis=2)

        if not is_jax(xp):  # no data-dependent input validation for lazy arrays
            message = "The input contains nan values"
            with pytest.raises(ValueError, match=message):
                stats.quantile(xp.asarray([xp.nan, 1, 2]), p, nan_policy='raise')

        message = "method` must be one of..."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, method='a duck')

        message = "If specified, `keepdims` must be True or False."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, p, keepdims=42)

        message = "`keepdims` may be False only if the length of `p` along `axis` is 1."
        with pytest.raises(ValueError, match=message):
            stats.quantile(x, xp.asarray([0.5, 0.6]), keepdims=False)

    @pytest.mark.parametrize('method',
         ['inverted_cdf', 'averaged_inverted_cdf', 'closest_observation',
          'hazen', 'interpolated_inverted_cdf', 'linear',
          'median_unbiased', 'normal_unbiased', 'weibull',
          '_lower', '_higher', '_midpoint', '_nearest'])
    @pytest.mark.parametrize('shape_x, shape_p, axis',
         [(10, None, -1), (10, 10, -1), (10, (2, 3), -1),
          ((10, 2), None, 0), ((10, 2), None, 0),])
    def test_against_numpy(self, method, shape_x, shape_p, axis, xp):
        dtype = xp_default_dtype(xp)
        rng = np.random.default_rng(23458924568734956)
        x = rng.random(size=shape_x)
        p = rng.random(size=shape_p)
        ref = np.quantile(x, p, axis=axis,
                          method=method[1:] if method.startswith('_') else method)
        x, p = xp.asarray(x, dtype=dtype), xp.asarray(p, dtype=dtype)
        res = stats.quantile(x, p, method=method, axis=axis)

        xp_assert_close(res, xp.asarray(ref, dtype=dtype))

    @skip_xp_backends(cpu_only=True, reason="PyTorch doesn't have `betainc`.",
                      exceptions=['cupy', 'jax.numpy'])
    @pytest.mark.parametrize('axis', [0, 1])
    @pytest.mark.parametrize('keepdims', [False, True])
    @pytest.mark.parametrize('nan_policy', ['omit', 'propagate', 'marray'])
    @pytest.mark.parametrize('dtype', ['float32', 'float64'])
    @pytest.mark.parametrize('method', ['linear', 'harrell-davis'])
    def test_against_reference(self, axis, keepdims, nan_policy, dtype, method, xp):
        if is_jax(xp) and nan_policy == 'marray':  # mdhaber/marray#146
            pytest.skip("`marray` currently incompatible with JAX")
        rng = np.random.default_rng(23458924568734956)
        shape = (5, 6)
        x = rng.random(size=shape).astype(dtype)
        p = rng.random(size=shape).astype(dtype)
        mask = rng.random(size=shape) > 0.8
        assert np.any(mask)
        x[mask] = np.nan
        if not keepdims:
            p = np.mean(p, axis=axis, keepdims=True)

        # inject p = 0 and p = 1 to test edge cases
        # Currently would fail with CuPy/JAX (cupy/cupy#8934, jax-ml/jax#21900);
        # remove the `if` when those are resolved.
        if is_numpy(xp):
            p0 = p.ravel()
            p0[1] = 0.
            p0[-2] = 1.

        dtype = getattr(xp, dtype)

        if nan_policy == 'marray':
            if not SCIPY_ARRAY_API:
                pytest.skip("MArray is only available if SCIPY_ARRAY_API=1")
            marray = pytest.importorskip('marray')
            kwargs = dict(axis=axis, keepdims=keepdims, method=method)
            mxp = marray._get_namespace(xp)
            x_mp = mxp.asarray(x, mask=mask)
            res = stats.quantile(x_mp, mxp.asarray(p), **kwargs)
            ref = quantile_reference(x, p, nan_policy='omit', **kwargs)
            xp_assert_close(res.data, xp.asarray(ref, dtype=dtype))
            return

        kwargs = dict(axis=axis, keepdims=keepdims,
                      nan_policy=nan_policy, method=method)
        res = stats.quantile(xp.asarray(x), xp.asarray(p), **kwargs)
        ref = quantile_reference(x, p, **kwargs)
        xp_assert_close(res, xp.asarray(ref, dtype=dtype))

    def test_integer_input_output_dtype(self, xp):
        res = stats.quantile(xp.arange(10, dtype=xp.int64), 0.5)
        assert res.dtype == xp_default_dtype(xp)

    @pytest.mark.parametrize('x, p, ref, kwargs',
        [([], 0.5, np.nan, {}),
         ([1, 2, 3], [-1, 0, 1, 1.5, np.nan], [np.nan, 1, 3, np.nan, np.nan], {}),
         ([1, 2, 3], [], [], {}),
         ([[np.nan, 2]], 0.5, [np.nan, 2], {'nan_policy': 'omit'}),
         ([[], []], 0.5, np.full(2, np.nan), {'axis': -1}),
         ([[], []], 0.5, np.zeros((0,)), {'axis': 0, 'keepdims': False}),
         ([[], []], 0.5, np.zeros((1, 0)), {'axis': 0, 'keepdims': True}),
         ([], [0.5, 0.6], np.full(2, np.nan), {}),
         (np.arange(1, 28).reshape((3, 3, 3)), 0.5, [[[14.]]],
          {'axis': None, 'keepdims': True}),
         ([[1, 2], [3, 4]], [0.25, 0.5, 0.75], [[1.75, 2.5, 3.25]],
          {'axis': None, 'keepdims': True}),])
    def test_edge_cases(self, x, p, ref, kwargs, xp):
        default_dtype = xp_default_dtype(xp)
        x, p, ref = xp.asarray(x), xp.asarray(p), xp.asarray(ref, dtype=default_dtype)
        res = stats.quantile(x, p, **kwargs)
        xp_assert_equal(res, ref)

    @pytest.mark.parametrize('axis', [0, 1, 2])
    @pytest.mark.parametrize('keepdims', [False, True])
    def test_size_0(self, axis, keepdims, xp):
        shape = [3, 4, 0]
        out_shape = shape.copy()
        if keepdims:
            out_shape[axis] = 1
        else:
            out_shape.pop(axis)
        res = stats.quantile(xp.zeros(tuple(shape)), 0.5, axis=axis, keepdims=keepdims)
        assert res.shape == tuple(out_shape)

    @pytest.mark.parametrize('method',
        ['inverted_cdf', 'averaged_inverted_cdf', 'closest_observation',
         '_lower', '_higher', '_midpoint', '_nearest'])
    def test_transition(self, method, xp):
        # test that values of discontinuous estimators are correct when
        # p*n + m - 1 is integral.
        if method == 'closest_observation' and np.__version__ < '2.0.1':
            pytest.skip('Bug in np.quantile (numpy/numpy#26656) fixed in 2.0.1')
        x = np.arange(8., dtype=np.float64)
        p = np.arange(0, 1.03125, 0.03125)
        res = stats.quantile(xp.asarray(x), xp.asarray(p), method=method)
        ref = np.quantile(x, p, method=method[1:] if method.startswith('_') else method)
        xp_assert_equal(res, xp.asarray(ref, dtype=xp.float64))


@_apply_over_batch(('x', 1), ('y', 1))
def iquantile_reference_last_axis(x, y, nan_policy, method):
    i_nan = np.isnan(x)
    if nan_policy == 'propagate' and np.any(i_nan):
        return np.full_like(y, np.nan)
    elif nan_policy == 'omit':
        x = x[~i_nan]
    return stats.iquantile(x, y, keepdims=True, method=method)


def iquantile_reference(x, y, *, axis=0, nan_policy='propagate',
                        keepdims=None, method='linear'):
    x, y = _broadcast_arrays((x, y), axis=axis)
    x, y = np.moveaxis(x, axis, -1), np.moveaxis(y, axis, -1)
    res = iquantile_reference_last_axis(x, y, nan_policy, method)
    res = np.moveaxis(res, -1, axis)
    if not keepdims:
        res = np.squeeze(res, axis=axis)
    return res


@make_xp_test_case(stats.iquantile)
class TestIQuantile:
    def test_input_validation(self, xp):
        x = xp.asarray([1, 2, 3])
        y = xp.asarray(2)

        message = "`x` must have real dtype."
        with pytest.raises(ValueError, match=message):
            stats.iquantile(xp.asarray([True, False]), y)
        with pytest.raises(ValueError):
            stats.iquantile(xp.asarray([1+1j, 2]), y)

        message = "`y` must have real dtype."
        with pytest.raises(ValueError, match=message):
            stats.iquantile(x, xp.asarray([0+1j, 1]))

        message = "`axis` must be an integer or None."
        with pytest.raises(ValueError, match=message):
            stats.iquantile(x, y, axis=0.5)
        with pytest.raises(ValueError, match=message):
            stats.iquantile(x, y, axis=(0, -1))

        message = "`axis` is not compatible with the shapes of the inputs."
        with pytest.raises(ValueError, match=message):
            stats.iquantile(x, y, axis=2)

        if not is_jax(xp):  # no data-dependent input validation for lazy arrays
            message = "The input contains nan values"
            with pytest.raises(ValueError, match=message):
                stats.iquantile(xp.asarray([xp.nan, 1, 2]), y, nan_policy='raise')

        message = "method` must be one of..."
        with pytest.raises(ValueError, match=message):
            stats.iquantile(x, y, method='a duck')

        message = "If specified, `keepdims` must be True or False."
        with pytest.raises(ValueError, match=message):
            stats.iquantile(x, y, keepdims=42)

        message = "`keepdims` may be False only if the length of `y` along `axis` is 1."
        with pytest.raises(ValueError, match=message):
            stats.iquantile(x, xp.asarray([0.5, 0.6]), keepdims=False)

    @pytest.mark.parametrize('method', _iquantile_methods)
    @pytest.mark.parametrize('dtype', [None, 'float32', 'float64'])
    @pytest.mark.parametrize('x_shape', [2, 10, 11, 100, 1001, (2, 10), (2, 3, 11)])
    @pytest.mark.parametrize('y_shape', [None, 25])
    @pytest.mark.parametrize('ties', [False, True])
    def test_against_quantile(self, method, dtype, x_shape, y_shape, ties, xp):
        discontinuous = method in _iquantile_discontinuous_methods
        dtype = xp_default_dtype(xp) if dtype is None else getattr(xp, dtype)
        rng = np.random.default_rng(394529872549827485)
        y_shape = x_shape if y_shape is None else y_shape

        if ties:
            x = xp.asarray(rng.integers(9, size=x_shape), dtype=dtype)
        else:
            x = xp.asarray(rng.standard_normal(size=x_shape), dtype=dtype)

        p = xp.asarray(rng.random(size=y_shape), dtype=dtype)
        y = stats.quantile(x, p, method=method, axis=-1)
        res = stats.iquantile(x, y, method=method, axis=-1)
        ref = xp.broadcast_to(p, (*x.shape[:-1], y.shape[-1]))

        # check that `quantile` is the inverse of `iquantile`
        # note that for discontinuous methods, res is right on the cusp of a transition,
        # and there can be a tiny bit of error to the right or left. We shift it left
        # to ensure we're on the correct side of the transition, producing the same `y2`
        # as if the probability calculation were exact.
        res = res - 1e-6 if discontinuous else res
        y2 = stats.quantile(x, res, method=method, axis=-1)
        atol = 1e-6 if dtype == xp.float32 else 1e-12
        xp_assert_close(y2, y, atol=atol)

        # if there are ties or method is discontinuous, `quantile` is not invertible
        if ties or discontinuous:
            return

        # `quantile` is not invertible outside this domain
        a, b = _iquantile_continuous_methods[method]
        n = x.shape[-1]
        p_min = (1 - a) / (n + 1 - a - b)
        p_max = (n - a) / (n + 1 - a - b)
        i_very_low = y < xp.min(x, axis=-1, keepdims=True)
        i_very_high = y > xp.max(x, axis=-1, keepdims=True)
        i_low = (ref <= p_min) & ~i_very_low
        i_high = (ref >= p_max) & ~i_very_high
        i_ok = ~(i_low | i_high | i_very_low | i_very_high)

        # check for correct inversion within the domain
        xp_assert_close(res[i_ok], ref[i_ok])

        # check that all other values get mapped to bottom or top of range
        kwargs = dict(check_shape=False, check_dtype=False, check_0d=True)
        xp_assert_close(res[i_low], xp.asarray(p_min), **kwargs)
        xp_assert_close(res[i_high], xp.asarray(p_max), **kwargs)
        xp_assert_close(res[i_very_low], xp.asarray(0.0), **kwargs)
        xp_assert_close(res[i_very_high], xp.asarray(1.0), **kwargs)

    @pytest.mark.parametrize('axis', [0, 1])
    @pytest.mark.parametrize('keepdims', [False, True])
    @pytest.mark.parametrize('nan_policy', ['propagate', 'omit', 'marray'])
    @pytest.mark.parametrize('dtype', ['float32', 'float64'])
    @pytest.mark.parametrize('nans', [False, True])
    @pytest.mark.parametrize('meth', ['linear', 'inverted_cdf'])
    def test_against_reference(self, axis, keepdims, nan_policy, dtype, nans, meth, xp):
        if is_jax(xp) and nan_policy == 'marray':  # mdhaber/marray#146
            pytest.skip("`marray` currently incompatible with JAX")
        rng = np.random.default_rng(23458924568734956)
        shape = (5, 6)
        x = rng.standard_normal(size=shape).astype(dtype)
        y = rng.standard_normal(size=shape).astype(dtype)

        mask = None
        if nans:
            mask = rng.random(size=shape) > 0.8
            assert np.any(mask)
            x[mask] = np.nan

        if not keepdims:
            y = np.mean(y, axis=axis, keepdims=True)

        dtype = getattr(xp, dtype)

        if nan_policy == 'marray':
            if not SCIPY_ARRAY_API:
                pytest.skip("MArray is only available if SCIPY_ARRAY_API=1")
            marray = pytest.importorskip('marray')
            kwargs = dict(axis=axis, keepdims=keepdims, method=meth)
            mxp = marray._get_namespace(xp)
            x_mp = mxp.asarray(x, mask=mask)
            res = stats.iquantile(x_mp, mxp.asarray(y), **kwargs)
            ref = iquantile_reference(x, y, nan_policy='omit', **kwargs)
            xp_assert_close(res.data, xp.asarray(ref, dtype=dtype))
            return

        kwargs = dict(axis=axis, keepdims=keepdims,
                      nan_policy=nan_policy, method=meth)
        res = stats.iquantile(xp.asarray(x), xp.asarray(y), **kwargs)
        ref = iquantile_reference(x, y, **kwargs)
        xp_assert_close(res, xp.asarray(ref, dtype=dtype))

    @pytest.mark.parametrize('n', [50, 500])
    @pytest.mark.parametrize('method, ab', _iquantile_continuous_methods.items())
    def test_plotting_positions(self, n, method, ab, xp):
        a, b = ab
        rng = np.random.default_rng(539452987254982748)
        x = rng.standard_normal(n)

        mask = rng.random(n) < 0.1
        x[mask] = np.nan
        mask = xp.asarray(mask)

        x_masked = np.ma.masked_invalid(x)
        ref = stats.mstats.plotting_positions(x_masked, a, b)
        ref = xp.asarray(ref.data)

        x = xp.asarray(x)
        res = stats.iquantile(x, x, nan_policy='omit', method=method)

        xp_assert_close(res[~mask], ref[~mask])
        assert xp.all(xp.isnan(res[mask]))

    @pytest.mark.parametrize('ties', [False, True])
    def test_against_ecdf_percentileofscore(self, ties, xp):
        rng = np.random.default_rng(853945298725498274)
        n = 50
        dtype = xp_default_dtype(xp)
        x = rng.integers(10, size=n) if ties else rng.standard_normal(size=n)
        y = rng.integers(10, size=25) if ties else rng.standard_normal(size=25)
        ref = stats.ecdf(x).cdf.evaluate(y)
        ref2 = stats.percentileofscore(x, y, 'weak')
        x, y = xp.asarray(x, dtype=dtype), xp.asarray(y, dtype=dtype)
        res = stats.iquantile(x, y, method='inverted_cdf')
        ref, ref2 = xp.asarray(ref, dtype=dtype), xp.asarray(ref2, dtype=dtype)
        xp_assert_close(res, ref)
        xp_assert_close(res, ref2 / 100)

    def test_integer_input_output_dtype(self, xp):
        x = xp.arange(10, dtype=xp.int64)
        res = stats.iquantile(x, x)
        assert res.dtype == xp_default_dtype(xp)

    @pytest.mark.parametrize('nan_policy', ['propagate', 'omit', 'marray'])
    @pytest.mark.parametrize('method', _iquantile_methods)
    def test_size_one_sample(self, nan_policy, method, xp):
        discontinuous = method in _iquantile_discontinuous_methods
        x = xp.arange(10.)
        y = xp.asarray([0., -1., 1.])
        n = np.asarray(1.)
        with np.errstate(divide='ignore', invalid='ignore'):  # for method = 'linear'
            if discontinuous:
                ref = xp.asarray([1., 0., 1.])
            else:
                a, b = _iquantile_continuous_methods[method]
                ref = xp.asarray([(n - a) / (n + 1 - a - b), 0., 1.])

        if nan_policy == 'propagate':
            x = x[:1]
            kwargs = {'nan_policy': 'propagate'}
        elif nan_policy == 'omit':
            x = xpx.at(x)[1:].set(xp.nan)
            kwargs = {'nan_policy': 'omit'}
        elif nan_policy == 'marray':
            if is_jax(xp):
                pytest.skip("JAX currently incompatible with `marray`")
            if not SCIPY_ARRAY_API:
                pytest.skip("MArray is only available if SCIPY_ARRAY_API=1")
            marray = pytest.importorskip('marray')
            mxp = marray._get_namespace(xp)
            mask = (x > 0.)
            x = mxp.asarray(x, mask=mask)
            y = mxp.asarray(y)
            kwargs = {}

        with np.errstate(divide='ignore', invalid='ignore'):  # for method = 'linear'
            res = stats.iquantile(x, y, method=method, **kwargs)
        res = res.data if nan_policy == 'marray' else res
        xp_assert_close(res, ref)

    # skipping marray due to mdhaber/marray#24
    @pytest.mark.parametrize('nan_policy', ['propagate', 'omit'])
    @pytest.mark.parametrize('method', _iquantile_methods)
    def test_size_zero_sample(self, nan_policy, method, xp):
        x = xp.arange(10.)
        y = xp.asarray([0., -1., 1.])  # this should work
        ref = xp.full_like(y, xp.nan)

        if nan_policy == 'propagate':
            x = x[0:0]
            kwargs = {'nan_policy': 'propagate'}
        elif nan_policy == 'omit':
            x = xpx.at(x)[:].set(xp.nan)
            kwargs = {'nan_policy': 'omit'}
        elif nan_policy == 'marray':
            if not SCIPY_ARRAY_API:
                pytest.skip("MArray is only available if SCIPY_ARRAY_API=1")
            if is_jax(xp):
                pytest.skip("JAX currently incompatible with `marray`")
            marray = pytest.importorskip('marray')
            mxp = marray._get_namespace(xp)
            mask = (x >= 0.)
            x = mxp.asarray(x, mask=mask)
            y = mxp.asarray(y)
            kwargs = {}

        with np.errstate(divide='ignore', invalid='ignore'):  # for method = 'linear'
            res = stats.iquantile(x, y, method=method, **kwargs)

        if nan_policy == 'marray':
            assert xp.all(res.mask)
        else:
            xp_assert_close(res, ref)

    @pytest.mark.parametrize('x, y, ref, kwargs',
        [
         ([], 0.5, np.nan, {}),
         ([1, 2, 3], [0.999, 3.001, np.nan], [0., 1., np.nan], {}),
         ([1, 2, 3], [], [], {}),
         ([[np.nan, 2]], 2, [np.nan, 0.5], {'nan_policy': 'omit', 'method': 'weibull'}),
         ([[], []], 0.5, np.full(2, np.nan), {'axis': -1}),
         ([[], []], 0.5, np.zeros((0,)), {'axis': 0, 'keepdims': False}),
         ([[], []], 0.5, np.zeros((1, 0)), {'axis': 0, 'keepdims': True}),
         ([], [0.5, 0.6], np.full(2, np.nan), {}),
         (np.arange(1, 28).reshape((3, 3, 3)), 14., [[[0.5]]],
          {'axis': None, 'keepdims': True}),
         ([[1, 2], [3, 4]], [1.75, 2.5, 3.25], [[0.25, 0.5, 0.75]],
          {'axis': None, 'keepdims': True}),
         ([1, 2, 3], [-np.inf, np.inf], [0.0, 1.0], {}),
         # It is our choice how much effort and computational overhead we want to put
         # into adjusting for insane edge cases like when `x` contains infinite values,
         # especially when `y` does, too.
         # One practical argument would be that y = +/- inf should produce the same
         # results as an extremely large finite number, in which case the 0th element
         # of the 'linear' result should be `0`.
         # Another argument would be for producing NaN below wherever y = +/- inf.
         # Another would be that it's not appropriate to spend significant computation
         # correcting these edge cases; we should just document what we do.
         # In any case, this is the current status.
         ([-np.inf, -1, 0, 1, np.inf], [-np.inf, -2, -1, 0, 1, 2, np.inf],
          [0.25, 0.25, 0.25, 0.5 , 0.75, 0.75, np.nan], {'method':'linear'}),
         ([-np.inf, -1, 0, 1, np.inf], [-np.inf, -2, -1, 0, 1, 2, np.inf],
          [0.2, 0.2, 0.4, 0.6, 0.8, 0.8, 1.], {'method': 'inverted_cdf'}),
        ])
    def test_edge_cases(self, x, y, ref, kwargs, xp):
        if kwargs.get('method', None) == 'weibull' and is_torch(xp):
            pytest.skip('data-apis/array-api-compat#360')
        if kwargs.get('axis', None) == -1 and is_cupy(xp):
            pytest.skip('Fails; need to investigate.')
        default_dtype = xp_default_dtype(xp)
        x, y, ref = xp.asarray(x), xp.asarray(y), xp.asarray(ref, dtype=default_dtype)
        res = stats.iquantile(x, y, **kwargs)
        xp_assert_equal(res, ref)

    @pytest.mark.skip_xp_backends('jax.numpy', reason="arithmetic not exact for 2**k ?")
    @pytest.mark.parametrize('method', _iquantile_discontinuous_methods.keys())
    def test_transition(self, method, xp):
        # test that values of discontinuous estimators are as expected around
        # transition point
        x = np.arange(8., dtype=np.float64)
        xl, xr = np.nextafter(x, -np.inf), np.nextafter(x, np.inf)
        offset = 0.5 if method == 'closest_observation' else 0.0
        ref_r = (x + 1 + offset) / 8
        ref_r[-1] = 1.0  # value is greater than or equal to the maximum observation
        ref_l = (x + offset) / 8
        ref_l[0] = 0.0  # value is less than the minimum observation
        x, xl, xr = xp.asarray(x), xp.asarray(xl), xp.asarray(xr)
        ref_l, ref_r = xp.asarray(ref_l), xp.asarray(ref_r)
        xp_assert_equal(stats.iquantile(x, x, method=method), ref_r)
        xp_assert_equal(stats.iquantile(x, xr, method=method), ref_r)
        xp_assert_equal(stats.iquantile(x, xl, method=method), ref_l)
