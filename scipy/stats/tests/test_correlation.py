import pytest
import numpy as np

from scipy._lib._array_api import make_xp_test_case, xp_default_dtype, is_jax
from scipy._lib._array_api_no_0d import xp_assert_close
from scipy import stats
from scipy.stats._axis_nan_policy import SmallSampleWarning


@make_xp_test_case(stats.chatterjeexi)
class TestChatterjeeXi:
    @pytest.mark.parametrize('case', [
        dict(y_cont=True, statistic=-0.303030303030303, pvalue=0.9351329808526656),
        dict(y_cont=False, statistic=0.07407407407407396, pvalue=0.3709859367123997)])
    @pytest.mark.parametrize('dtype', ['float32', 'float64', None])
    def test_against_R_XICOR(self, case, dtype, xp):
        # Test against R package XICOR, e.g.
        # library(XICOR)
        # options(digits=16)
        # x = c(0.11027287231363914, 0.8154770102474279, 0.7073943466920335,
        #       0.6651317324378386, 0.6905752850115503, 0.06115250587536558,
        #       0.5209906494474178, 0.3155763519785274, 0.18405731803625924,
        #       0.8613557911541495)
        # y = c(0.8402081904493103, 0.5946972833914318, 0.23481606164114155,
        #       0.49754786197715384, 0.9146460831206026, 0.5848057749217579,
        #       0.7620801065573549, 0.31410063302647495, 0.7935620302236199,
        #       0.5423085761365468)
        # xicor(x, y, ties=FALSE, pvalue=TRUE)
        if dtype == 'float32' and np.__version__ < "2":
            pytest.skip("Scalar dtypes only respected after NEP 50.")
        dtype = xp_default_dtype(xp) if dtype is None else getattr(xp, dtype)
        rng = np.random.default_rng(25982435982346983)
        x = rng.random(size=10)
        y = (rng.random(size=10) if case['y_cont']
             else rng.integers(0, 5, size=10))

        x, y = xp.asarray(x, dtype=dtype), xp.asarray(y, dtype=dtype)
        res = stats.chatterjeexi(x, y, y_continuous=case['y_cont'])

        xp_assert_close(res.statistic, xp.asarray(case['statistic'], dtype=dtype))
        xp_assert_close(res.pvalue, xp.asarray(case['pvalue'], dtype=dtype))

    @pytest.mark.parametrize('y_continuous', (False, True))
    def test_permutation_asymptotic(self, y_continuous):
        # XICOR doesn't seem to perform the permutation test as advertised, so
        # compare the result of a permutation test against an asymptotic test.
        rng = np.random.default_rng(2524579827426)
        n = np.floor(rng.uniform(100, 150)).astype(int)
        shape = (2, n)
        x = rng.random(size=shape)
        y = (rng.random(size=shape) if y_continuous
             else rng.integers(0, 10, size=shape))
        method = stats.PermutationMethod(rng=rng)
        res = stats.chatterjeexi(x, y, method=method,
                                 y_continuous=y_continuous, axis=-1)
        ref = stats.chatterjeexi(x, y, y_continuous=y_continuous, axis=-1)
        np.testing.assert_allclose(res.statistic, ref.statistic, rtol=1e-15)
        np.testing.assert_allclose(res.pvalue, ref.pvalue, rtol=2e-2)

    def test_input_validation(self, xp):
        rng = np.random.default_rng(25932435798274926)
        x, y = rng.random(size=(2, 10))
        x, y = xp.asarray(x), xp.asarray(y)

        message = 'Array shapes are incompatible for broadcasting.|Incompatible shapes'
        with pytest.raises((ValueError, TypeError), match=message):
            stats.chatterjeexi(x, y[:-1])

        if not is_jax(xp):
            # jax misses out on some input validation from _axis_nan_policy decorator
            message = '...axis 10 is out of bounds for array...|out of range'
            with pytest.raises((ValueError, IndexError), match=message):
                stats.chatterjeexi(x, y, axis=10)

        message = '`y_continuous` must be boolean.'
        with pytest.raises(ValueError, match=message):
            stats.chatterjeexi(x, y, y_continuous='a herring')

        message = "`method` must be 'asymptotic' or"
        with pytest.raises(ValueError, match=message):
            stats.chatterjeexi(x, y, method='ekki ekii')

    @pytest.mark.skip_xp_backends('jax.numpy', reason='no SmallSampleWarning (lazy)')
    def test_special_cases(self, xp):
        message = 'One or more sample arguments is too small...'
        with pytest.warns(SmallSampleWarning, match=message):
            res = stats.chatterjeexi(xp.asarray([1]), xp.asarray([2]))

        assert xp.isnan(res.statistic)
        assert xp.isnan(res.pvalue)
