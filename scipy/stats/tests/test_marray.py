import numpy as np
import pytest
from scipy import stats

from scipy._lib._array_api import xp_assert_close, xp_assert_equal, _count_nonmasked
from scipy._lib._array_api import make_xp_pytest_param, make_xp_test_case
from scipy._lib._array_api import SCIPY_ARRAY_API, xp_default_dtype, is_torch
from scipy.stats._stats_py import _xp_mean, _xp_var


marray = pytest.importorskip('marray')
pytestmark = [
    pytest.mark.skipif(
        not SCIPY_ARRAY_API,
        reason=(
            "special function dispatch to marray required for these tests"
            " is hidden behind SCIPY_ARRAY_API flag."
        ),
    ),
]

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


@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@pytest.mark.parametrize('fun, kwargs', [
    make_xp_pytest_param(stats.gmean, {}),
    make_xp_pytest_param(stats.hmean, {}),
    make_xp_pytest_param(stats.pmean, {'p': 2}),
    make_xp_pytest_param(stats.expectile, {'alpha': 0.4})
])
@pytest.mark.parametrize('axis', [0, 1, None])
def test_xmean(fun, kwargs, axis, xp):
    mxp, marrays, narrays = get_arrays(2, xp=xp)
    res = fun(marrays[0], weights=marrays[1], axis=axis, **kwargs)
    ref = fun(narrays[0], weights=narrays[1], nan_policy='omit', axis=axis, **kwargs)
    xp_assert_close(res.data, xp.asarray(ref))


@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@pytest.mark.parametrize('axis', [0, 1, None])
@pytest.mark.parametrize('keepdims', [False, True])
@pytest.mark.uses_xp_capabilities(False, reason="private")
def test_xp_mean(axis, keepdims, xp):
    mxp, marrays, narrays = get_arrays(2, xp=xp)
    kwargs = dict(axis=axis, keepdims=keepdims)
    res = _xp_mean(marrays[0], weights=marrays[1], **kwargs)
    ref = _xp_mean(narrays[0], weights=narrays[1], nan_policy='omit', **kwargs)
    xp_assert_close(res.data, xp.asarray(ref))


@skip_backend('dask.array', reason='Arrays need `device` attribute: dask/dask#11711')
@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@pytest.mark.parametrize('fun, kwargs',
    [make_xp_pytest_param(stats.moment, {'order': 2}),
     make_xp_pytest_param(stats.skew, {}),
     make_xp_pytest_param(stats.skew, {'bias': False}),
     make_xp_pytest_param(stats.kurtosis, {}),
     make_xp_pytest_param(stats.kurtosis, {'bias': False}),
     make_xp_pytest_param(stats.sem, {}),
     make_xp_pytest_param(stats.kstat, {'n': 1}),
     make_xp_pytest_param(stats.kstat, {'n': 2}),
     make_xp_pytest_param(stats.kstat, {'n': 3}),
     make_xp_pytest_param(stats.kstat, {'n': 4}),
     make_xp_pytest_param(stats.kstatvar, {'n': 1}),
     make_xp_pytest_param(stats.kstatvar, {'n': 2}),
     make_xp_pytest_param(stats.circmean, {}),
     make_xp_pytest_param(stats.circvar, {}),
     make_xp_pytest_param(stats.circstd, {}),
     make_xp_pytest_param(stats.gstd, {}),
     make_xp_pytest_param(stats.variation, {}),
     pytest.param(
         _xp_var, {},
         marks=pytest.mark.uses_xp_capabilities(False, reason="private")
     ),
     make_xp_pytest_param(stats.tmean, {'limits': (0.1, 0.9)}),
     make_xp_pytest_param(stats.tvar, {'limits': (0.1, 0.9)}),
     make_xp_pytest_param(stats.tmin, {'lowerlimit': 0.5}),
     make_xp_pytest_param(stats.tmax, {'upperlimit': 0.5}),
     make_xp_pytest_param(stats.tstd, {'limits': (0.1, 0.9)}),
     make_xp_pytest_param(stats.tsem, {'limits': (0.1, 0.9)}),
     make_xp_pytest_param(stats.iqr, {}),
     make_xp_pytest_param(stats.median_abs_deviation, {}),
     ])
@pytest.mark.parametrize('axis', [0, 1, None])
def test_several(fun, kwargs, axis, xp):
    mxp, marrays, narrays = get_arrays(1, xp=xp)
    kwargs = dict(axis=axis) | kwargs
    res = fun(marrays[0], **kwargs)
    ref = fun(narrays[0], nan_policy='omit', **kwargs)
    xp_assert_close(res.data, xp.asarray(ref))


@make_xp_test_case(stats.describe)
@skip_backend('dask.array', reason='Arrays need `device` attribute: dask/dask#11711')
@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@pytest.mark.parametrize('axis', [0, 1, None])
@pytest.mark.parametrize('kwargs', [{}])
def test_describe(axis, kwargs, xp):
    mxp, marrays, narrays = get_arrays(1, xp=xp)
    kwargs = dict(axis=axis) | kwargs
    res = stats.describe(marrays[0], **kwargs)
    ref = stats.describe(narrays[0], nan_policy='omit', **kwargs)
    xp_assert_close(res.nobs.data, xp.asarray(ref.nobs))
    xp_assert_close(res.minmax[0].data, xp.asarray(ref.minmax[0].data))
    xp_assert_close(res.minmax[1].data, xp.asarray(ref.minmax[1].data))
    # copy reference arrays due to torch complaint about non-writeable buffer
    xp_assert_close(res.variance.data,
                    xp.asarray(np.asarray(ref.variance.data, copy=True)))
    xp_assert_close(res.skewness.data,
                    xp.asarray(np.asarray(ref.skewness.data, copy=True)))
    xp_assert_close(res.kurtosis.data,
                    xp.asarray(np.asarray(ref.kurtosis.data, copy=True)))


@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@pytest.mark.parametrize('fun', [make_xp_pytest_param(stats.zscore),
                                 make_xp_pytest_param(stats.gzscore),
                                 make_xp_pytest_param(stats.zmap)])
@pytest.mark.parametrize('axis', [0, 1, None])
def test_zscore(fun, axis, xp):
    mxp, marrays, narrays = (get_arrays(2, xp=xp) if fun == stats.zmap
                             else get_arrays(1, xp=xp))
    res = fun(*marrays, axis=axis)
    ref = xp.asarray(fun(*narrays, nan_policy='omit', axis=axis))
    xp_assert_close(res.data[~res.mask], ref[~xp.isnan(ref)])
    xp_assert_equal(res.mask, marrays[0].mask)


@skip_backend('dask.array', reason='Arrays need `device` attribute: dask/dask#11711')
@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@pytest.mark.parametrize('f', [make_xp_pytest_param(stats.ttest_1samp),
                               make_xp_pytest_param(stats.ttest_rel),
                               make_xp_pytest_param(stats.ttest_ind)])
@pytest.mark.parametrize('axis', [0, 1, None])
def test_ttest(f, axis, xp):
    f_name = f.__name__
    mxp, marrays, narrays = get_arrays(2, xp=xp)
    if f_name == 'ttest_1samp':
        marrays[1] = mxp.mean(marrays[1], axis=axis, keepdims=axis is not None)
        narrays[1] = np.nanmean(narrays[1], axis=axis, keepdims=axis is not None)
    res = f(*marrays, axis=axis)
    ref = f(*narrays, nan_policy='omit', axis=axis)
    xp_assert_close(res.statistic.data, xp.asarray(ref.statistic))
    xp_assert_close(res.pvalue.data, xp.asarray(ref.pvalue))
    res_ci = res.confidence_interval()
    ref_ci = ref.confidence_interval()
    xp_assert_close(res_ci.low.data, xp.asarray(ref_ci.low))
    xp_assert_close(res_ci.high.data, xp.asarray(ref_ci.high))


@skip_backend('dask.array', reason='Arrays need `device` attribute: dask/dask#11711')
@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@pytest.mark.filterwarnings("ignore::scipy.stats._axis_nan_policy.SmallSampleWarning")
@pytest.mark.parametrize('f', [make_xp_pytest_param(stats.skewtest),
                               make_xp_pytest_param(stats.kurtosistest),
                               make_xp_pytest_param(stats.normaltest),
                               make_xp_pytest_param(stats.jarque_bera)])
@pytest.mark.parametrize('axis', [0, 1, None])
def test_normality_tests(f, axis, xp):
    mxp, marrays, narrays = get_arrays(1, xp=xp, shape=(10, 11))

    res = f(*marrays, axis=axis)
    ref = f(*narrays, nan_policy='omit', axis=axis)

    xp_assert_close(res.statistic.data, xp.asarray(ref.statistic))
    xp_assert_close(res.pvalue.data, xp.asarray(ref.pvalue), atol=1e-16)


@skip_backend('dask.array', reason='Arrays need `device` attribute: dask/dask#11711')
@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@pytest.mark.parametrize('f, lambda_', [
    make_xp_pytest_param(stats.chisquare, {}),
    make_xp_pytest_param(stats.power_divergence, {'lambda_': 'pearson'}),
    make_xp_pytest_param(stats.power_divergence, {'lambda_': 'log-likelihood'}),
    make_xp_pytest_param(stats.power_divergence, {'lambda_': 'freeman-tukey'}),
    make_xp_pytest_param(stats.power_divergence,
                         {'lambda_': 'mod-log-likelihood'}),
    make_xp_pytest_param(stats.power_divergence, {'lambda_': 'neyman'}),
    make_xp_pytest_param(stats.power_divergence, {'lambda_': 'cressie-read'}),
])
@pytest.mark.parametrize('ddof', [0, 1])
@pytest.mark.parametrize('axis', [0, 1, None])
def test_power_divergence_chisquare(f, lambda_, ddof, axis, xp):
    mxp, marrays, narrays = get_arrays(2, xp=xp, shape=(5, 6))

    kwargs = dict(axis=axis, ddof=ddof)

    # test 1-arg
    res = f(marrays[0], **lambda_, **kwargs)
    ref = stats.power_divergence(narrays[0], nan_policy='omit', **lambda_, **kwargs)

    xp_assert_close(res.statistic.data, xp.asarray(ref[0]))
    xp_assert_close(res.pvalue.data, xp.asarray(ref[1]))

    # test 2-arg
    common_mask = np.isnan(narrays[0]) | np.isnan(narrays[1])
    normalize = (np.nansum(narrays[1] * ~common_mask, axis=axis, keepdims=True)
                 / np.nansum(narrays[0] * ~common_mask, axis=axis, keepdims=True))
    marrays[0] *= xp.asarray(normalize)
    narrays[0] *= normalize

    res = f(*marrays, **lambda_, **kwargs)
    ref = stats.power_divergence(*narrays, nan_policy='omit', **lambda_, **kwargs)

    xp_assert_close(res.statistic.data, xp.asarray(ref[0]))
    xp_assert_close(res.pvalue.data, xp.asarray(ref[1]))


@make_xp_test_case(stats.combine_pvalues)
@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@pytest.mark.parametrize('method', ['fisher', 'pearson', 'mudholkar_george',
                                    'tippett', 'stouffer'])
@pytest.mark.parametrize('axis', [0, 1, None])
def test_combine_pvalues(method, axis, xp):
    mxp, marrays, narrays = get_arrays(2, xp=xp, shape=(10, 11))

    kwargs = dict(method=method, axis=axis)
    res = stats.combine_pvalues(marrays[0], **kwargs)
    ref = stats.combine_pvalues(narrays[0], nan_policy='omit', **kwargs)

    xp_assert_close(res.statistic.data, xp.asarray(ref.statistic))
    xp_assert_close(res.pvalue.data, xp.asarray(ref.pvalue))

    if method != 'stouffer':
        return

    res = stats.combine_pvalues(marrays[0], weights=marrays[1], **kwargs)
    ref = stats.combine_pvalues(narrays[0], weights=narrays[1],
                                nan_policy='omit', **kwargs)

    xp_assert_close(res.statistic.data, xp.asarray(ref.statistic))
    xp_assert_close(res.pvalue.data, xp.asarray(ref.pvalue))


@make_xp_test_case(stats.ttest_ind_from_stats)
@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
def test_ttest_ind_from_stats(xp):
    shape = (10, 11)
    mxp, marrays, narrays = get_arrays(6, xp=xp, shape=shape)
    mask = np.sum(np.stack([np.isnan(arg) for arg in narrays]), axis=0).astype(bool)
    narrays = [arg[~mask] for arg in narrays]
    marrays[2], marrays[5] = marrays[2] * 100, marrays[5] * 100
    narrays[2], narrays[5] = narrays[2] * 100, narrays[5] * 100

    res = stats.ttest_ind_from_stats(*marrays)
    ref = stats.ttest_ind_from_stats(*narrays)

    mask = xp.asarray(mask)
    assert xp.any(mask) and xp.any(~mask)
    xp_assert_close(res.statistic.data[~mask], xp.asarray(ref.statistic))
    xp_assert_close(res.pvalue.data[~mask], xp.asarray(ref.pvalue))
    xp_assert_close(res.statistic.mask, mask)
    xp_assert_close(res.pvalue.mask, mask)
    assert res.statistic.shape == shape
    assert res.pvalue.shape == shape


def test_length_nonmasked_marray_iterable_axis_raises():
    xp = marray._get_namespace(np)

    data = [[1.0, 2.0], [3.0, 4.0]]
    mask = [[False, False], [True, False]]
    marr = xp.asarray(data, mask=mask)

    # Axis tuples are not currently supported for MArray input.
    # This test can be removed after support is added.
    with pytest.raises(NotImplementedError,
        match="`axis` must be an integer or None for use with `MArray`"):
        _count_nonmasked(marr, axis=(0, 1), xp=xp)


@make_xp_test_case(stats.directional_stats)
@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@pytest.mark.filterwarnings("ignore::RuntimeWarning")  # mdhaber/marray#120
@pytest.mark.parametrize('axis', [0, 1])
def test_directional_stats(xp, axis):
    mxp, marrays, narrays = get_arrays(1, shape=(19, 20, 3), xp=xp)
    res = stats.directional_stats(*marrays, axis=axis)

    x, = narrays
    if axis == 0:
        x = np.swapaxes(x, 0, 1)

    for i in range(x.shape[0]):
        xi = x[i]
        xi = xi[~np.any(np.isnan(xi), axis=-1), :]
        ref = stats.directional_stats(xi)
        xp_assert_close(res.mean_direction.data[i, ...],
                        xp.asarray(ref.mean_direction))
        xp_assert_close(res.mean_resultant_length.data[i],
                        xp.asarray(ref.mean_resultant_length)[()])
    assert not xp.any(res.mean_direction.mask)
    assert not xp.any(res.mean_resultant_length.mask)


@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@pytest.mark.parametrize('fun, kwargs', [
    make_xp_pytest_param(stats.wilcoxon,
                         {'method': 'asymptotic', 'zero_method': 'zsplit'}),
])
@pytest.mark.parametrize('axis', [0, 1, None])
def test_one_sample_tests(fun, kwargs, axis, xp):
    mxp, marrays, narrays = get_arrays(1, xp=xp, seed=84912165484322)
    res = fun(*marrays, axis=axis, **kwargs)
    ref = fun(*narrays, nan_policy='omit', axis=axis, **kwargs)
    xp_assert_close(res.statistic.data, xp.asarray(ref.statistic))
    xp_assert_close(res.pvalue.data, xp.asarray(ref.pvalue))


@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@pytest.mark.parametrize('fun, kwargs', [
    make_xp_pytest_param(stats.brunnermunzel, {}),
    make_xp_pytest_param(stats.mannwhitneyu, {'method': 'asymptotic'}),
    make_xp_pytest_param(stats.wilcoxon,
                         {'method': 'asymptotic', 'zero_method': 'zsplit'}),
    make_xp_pytest_param(stats.cramervonmises_2samp, {'method': 'exact'}),
])
@pytest.mark.parametrize('axis', [0, 1, None])
def test_two_sample_tests(fun, kwargs, axis, xp):
    if fun == stats.cramervonmises_2samp and axis is None:
        pytest.skip("Sample too large for exact method.")
    mxp, marrays, narrays = get_arrays(2, xp=xp, seed=84912165484322)
    res = fun(*marrays, axis=axis, **kwargs)
    ref = fun(*narrays, nan_policy='omit', axis=axis, **kwargs)
    xp_assert_close(res.statistic.data, xp.asarray(ref.statistic))
    xp_assert_close(res.pvalue.data, xp.asarray(ref.pvalue))


@skip_backend('dask.array', reason='Arrays need `device` attribute: dask/dask#11711')
@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@pytest.mark.parametrize('fun, kwargs', [
    make_xp_pytest_param(stats.bartlett, {}),
    make_xp_pytest_param(stats.alexandergovern, {}),
    make_xp_pytest_param(stats.levene, {'center': 'median'}),
    make_xp_pytest_param(stats.levene, {'center': 'mean'}),
    make_xp_pytest_param(stats.f_oneway, {'equal_var': True}),
    make_xp_pytest_param(stats.f_oneway, {'equal_var': False}),
    make_xp_pytest_param(stats.kruskal, {}),
    make_xp_pytest_param(stats.fligner, {'center': 'median'}),
    make_xp_pytest_param(stats.fligner, {'center': 'mean'}),
])
@pytest.mark.parametrize('axis', [0, 1, None])
def test_k_sample_tests(fun, kwargs, axis, xp):
    mxp, marrays, narrays = get_arrays(3, xp=xp)
    res = fun(*marrays, axis=axis, **kwargs)
    ref = fun(*narrays, nan_policy='omit', axis=axis, **kwargs)
    xp_assert_close(res.statistic.data, xp.asarray(ref.statistic))
    xp_assert_close(res.pvalue.data, xp.asarray(ref.pvalue))


@skip_backend('jax.numpy', reason="JAX currently incompatible with marray")
@pytest.mark.parametrize('fun, kwargs', [
    make_xp_pytest_param(stats.friedmanchisquare, {}),
])
@pytest.mark.parametrize('axis', [0, 1, None])
def test_k_sample_paired_tests(fun, kwargs, axis, xp):
    mxp, marrays, narrays = get_arrays(3, shape=(8, 9), xp=xp)
    res = fun(*marrays, axis=axis, **kwargs)
    ref = fun(*narrays, nan_policy='omit', axis=axis, **kwargs)
    xp_assert_close(res.statistic.data, xp.asarray(ref.statistic))
    xp_assert_close(res.pvalue.data, xp.asarray(ref.pvalue))


@skip_backend('dask.array', reason='Arrays need `device` attribute: dask/dask#11711')
@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@pytest.mark.parametrize('f', [
    make_xp_pytest_param(stats.pearsonr),
    make_xp_pytest_param(stats.pointbiserialr),
    make_xp_pytest_param(stats.spearmanrho),
])
@pytest.mark.parametrize('axis', [0, 1, None])
def test_pearsonr(f, axis, xp):
    mxp, marrays, narrays = get_arrays(2, xp=xp)
    res = f(*marrays, axis=axis)

    # `pearsonr` does not have `axis_nan_policy`, so do this manually
    x, y = narrays
    if axis == 0:
        x, y = x.T, y.T
    elif axis is None:
        x, y = x.ravel()[np.newaxis, :], y.ravel()[np.newaxis, :]

    for i in range(x.shape[0]):
        xi, yi = x[i, ...], y[i, ...]
        i = () if axis is None else i

        mask = np.isnan(xi) | np.isnan(yi)
        ref = f(xi[~mask], yi[~mask])

        atol = 1e-7 if (is_torch(xp) and xp_default_dtype(xp) == xp.float32
                        and f == stats.spearmanrho) else 0.
        xp_assert_close(res.statistic.data[i], xp.asarray(ref.statistic)[()])
        xp_assert_close(res.pvalue.data[i], xp.asarray(ref.pvalue)[()], atol=atol)

        if f == stats.pearsonr:
            res_ci_low, res_ci_high = res.confidence_interval()
            ref_ci_low, ref_ci_high = ref.confidence_interval()
            xp_assert_close(res_ci_low.data[i], xp.asarray(ref_ci_low)[()])
            xp_assert_close(res_ci_high.data[i], xp.asarray(ref_ci_high)[()])


@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@pytest.mark.parametrize('f', [make_xp_pytest_param(stats.linregress)])
@pytest.mark.parametrize('axis', [0, 1, None])
def test_linregress(f, axis, xp):
    mxp, marrays, narrays = get_arrays(2, seed=84912165484320, xp=xp)
    res = f(*marrays, axis=axis)
    ref = f(*narrays, nan_policy='omit', axis=axis)

    xp_assert_close(res.slope.data, xp.asarray(ref.slope))
    xp_assert_close(res.intercept.data, xp.asarray(ref.intercept))
    xp_assert_close(res.rvalue.data, xp.asarray(ref.rvalue))
    xp_assert_close(res.pvalue.data, xp.asarray(ref.pvalue))
    xp_assert_close(res.stderr.data, xp.asarray(ref.stderr))
    xp_assert_close(res.intercept_stderr.data, xp.asarray(ref.intercept_stderr))


@make_xp_test_case(stats.entropy)
@skip_backend('jax.numpy', reason="JAX doesn't allow item assignment.")
@pytest.mark.parametrize('qk', [False, True])
@pytest.mark.parametrize('axis', [0, 1, None])
def test_entropy(qk, axis, xp):
    mxp, marrays, narrays = get_arrays(2 if qk else 1, xp=xp)
    res = stats.entropy(*marrays, axis=axis)
    ref = stats.entropy(*narrays, nan_policy='omit', axis=axis)
    xp_assert_close(res.data, xp.asarray(ref))


@make_xp_test_case(stats.rankdata)
@pytest.mark.parametrize('axis', [0, 1, None])
@skip_backend('jax.numpy', reason="JAX currently incompatible with marray")
def test_rankdata(axis, xp):
    mxp, marrays, narrays = get_arrays(1, xp=xp)
    res = stats.rankdata(*marrays, axis=axis)
    ref = stats.rankdata(*narrays, nan_policy='omit', axis=axis)
    xp_assert_close(res.data[~res.mask], xp.asarray(ref[~np.isnan(ref)]))
    xp_assert_close(res.mask, xp.asarray(np.isnan(ref)))


@pytest.mark.parametrize('f', [
    make_xp_pytest_param(stats.levene),
    make_xp_pytest_param(stats.fligner),
])
def test_center_trimmed(f, xp):
    mxp, marrays, narrays = get_arrays(3, xp=xp)
    message = "`center='trimmed'` is incompatible with MArray."
    with pytest.raises(ValueError, match=message):
        f(*marrays, center='trimmed', axis=-1)
