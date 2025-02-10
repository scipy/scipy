import pytest
import numpy as np
from scipy import stats

from scipy._lib._array_api import xp_assert_close, xp_assert_equal
from scipy.stats._stats_py import _xp_mean, _xp_var

marray = pytest.importorskip('marray')
skip_backend = pytest.mark.skip_xp_backends


def get_arrays(n_arrays, *, dtype='float64', xp=np, shape=(7, 8), seed=84912165484321):
    mxp = marray._get_namespace(xp)
    rng = np.random.default_rng(seed)

    datas, masks = [], []
    for i in range(n_arrays):
        data = rng.random(size=shape)
        if dtype.startswith('complex'):
            data = 10*data * 10j*rng.standard_normal(size=shape)
        data = data.astype(dtype)
        datas.append(data)
        mask = rng.random(size=shape) > 0.75
        masks.append(mask)

    marrays = []
    nan_arrays = []
    for array, mask in zip(datas, masks):
        marrays.append(mxp.asarray(array, mask=mask))
        nan_array = array.copy()
        nan_array[mask] = xp.nan
        nan_arrays.append(nan_array)

    return mxp, marrays, nan_arrays


@skip_backend('dask.array', reason='Arrays need `device` attribute: dask/dask#11711')
@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@pytest.mark.parametrize('fun, kwargs', [(stats.gmean, {}),
                                         (stats.hmean, {}),
                                         (stats.pmean, {'p': 2})])
@pytest.mark.parametrize('axis', [0, 1])
def test_xmean(fun, kwargs, axis, xp):
    mxp, marrays, narrays = get_arrays(2, xp=xp)
    res = fun(marrays[0], weights=marrays[1], axis=axis, **kwargs)
    ref = fun(narrays[0], weights=narrays[1], nan_policy='omit', axis=axis, **kwargs)
    xp_assert_close(res.data, xp.asarray(ref))


@skip_backend('dask.array', reason='Arrays need `device` attribute: dask/dask#11711')
@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@pytest.mark.parametrize('axis', [0, 1, None])
@pytest.mark.parametrize('keepdims', [False, True])
def test_xp_mean(axis, keepdims, xp):
    mxp, marrays, narrays = get_arrays(2, xp=xp)
    kwargs = dict(axis=axis, keepdims=keepdims)
    res = _xp_mean(marrays[0], weights=marrays[1], **kwargs)
    ref = _xp_mean(narrays[0], weights=narrays[1], nan_policy='omit', **kwargs)
    xp_assert_close(res.data, xp.asarray(ref))


@skip_backend('dask.array', reason='Arrays need `device` attribute: dask/dask#11711')
@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@skip_backend('torch', reason="array-api-compat#242")
@pytest.mark.parametrize('fun, kwargs',
    [(stats.moment, {'order': 2}),
     (stats.skew, {}),
     (stats.skew, {'bias': False}),
     (stats.kurtosis, {}),
     (stats.kurtosis, {'bias': False}),
     (stats.sem, {}),
     (stats.kstat, {'n': 1}),
     (stats.kstat, {'n': 2}),
     (stats.kstat, {'n': 3}),
     (stats.kstat, {'n': 4}),
     (stats.kstatvar, {'n': 1}),
     (stats.kstatvar, {'n': 2}),
     (stats.circmean, {}),
     (stats.circvar, {}),
     (stats.circstd, {}),
     (_xp_var, {}),
     (stats.tmean, {'limits': (0.1, 0.9)}),
     (stats.tvar, {'limits': (0.1, 0.9)}),
     (stats.tmin, {'lowerlimit': 0.5}),
     (stats.tmax, {'upperlimit': 0.5}),
     (stats.tstd, {'limits': (0.1, 0.9)}),
     (stats.tsem, {'limits': (0.1, 0.9)}),
     ])
@pytest.mark.parametrize('axis', [0, 1, None])
def test_several(fun, kwargs, axis, xp):
    mxp, marrays, narrays = get_arrays(1, xp=xp)
    kwargs = dict(axis=axis) | kwargs
    res = fun(marrays[0], **kwargs)
    ref = fun(narrays[0], nan_policy='omit', **kwargs)
    xp_assert_close(res.data, xp.asarray(ref))


@skip_backend('dask.array', reason='Arrays need `device` attribute: dask/dask#11711')
@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@skip_backend('torch', reason="array-api-compat#242")
@pytest.mark.parametrize('axis', [0, 1])
@pytest.mark.parametrize('kwargs', [{}])
def test_describe(axis, kwargs, xp):
    mxp, marrays, narrays = get_arrays(1, xp=xp)
    kwargs = dict(axis=axis) | kwargs
    res = stats.describe(marrays[0], **kwargs)
    ref = stats.describe(narrays[0], nan_policy='omit', **kwargs)
    xp_assert_close(res.nobs.data, xp.asarray(ref.nobs))
    xp_assert_close(res.minmax[0].data, xp.asarray(ref.minmax[0].data))
    xp_assert_close(res.minmax[1].data, xp.asarray(ref.minmax[1].data))
    xp_assert_close(res.variance.data, xp.asarray(ref.variance.data))
    xp_assert_close(res.skewness.data, xp.asarray(ref.skewness.data))
    xp_assert_close(res.kurtosis.data, xp.asarray(ref.kurtosis.data))


@skip_backend('dask.array', reason='Arrays need `device` attribute: dask/dask#11711')
@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@skip_backend('torch', reason="array-api-compat#242")
@pytest.mark.parametrize('fun', [stats.zscore, stats.gzscore, stats.zmap])
@pytest.mark.parametrize('axis', [0, 1, None])
def test_zscore(fun, axis, xp):
    mxp, marrays, narrays = (get_arrays(2, xp=xp) if fun == stats.zmap
                             else get_arrays(1, xp=xp))
    res = fun(*marrays, axis=axis)
    ref = xp.asarray(fun(*narrays, nan_policy='omit', axis=axis))
    xp_assert_close(res.data[~res.mask], ref[~xp.isnan(ref)])
    xp_assert_equal(res.mask, marrays[0].mask)
