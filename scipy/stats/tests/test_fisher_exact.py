import pytest
import numpy as np
from numpy.testing import assert_equal, assert_allclose
from scipy.stats import fisher_exact
from scipy.stats._fisher_exact import _nc_hypergeom_mean
from .fisher_exact_results_from_r import data


class TestFisherExact:

    @pytest.mark.parametrize('table, expected_pvalue', [
        ([[14500, 20000], [30000, 40000]], 0.0110551177046),
        ([[100, 2], [1000, 5]], 0.13007593634300169),
        ([[2, 7], [8, 2]], 0.023014137565221155),
        ([[5, 1], [10, 10]], 0.19732441471571907),
        ([[5, 15], [20, 20]], 0.09580440012477634),
        ([[5, 16], [20, 25]], 0.1725864953812995),
        ([[10, 5], [10, 1]], 0.19732441471571907),
        ([[10, 5], [10, 0]], 0.06126482213438735),
        ([[5, 0], [1, 4]], 0.047619047619047616),
    ])
    def test_basic_pvalue_only(self, table, expected_pvalue):
        oddsratio, pvalue = fisher_exact(table)
        a, b, c, d = np.ravel(table)
        r = (a*d)/(b*c) if b*c > 0 else np.inf if a*d > 0 else np.nan
        assert_allclose(oddsratio, r, rtol=1e-15)
        assert_allclose(pvalue, expected_pvalue, rtol=1e-7)

    def test_basic_oddsratio_zero(self):
        oddsratio, pvalue = fisher_exact([[0, 1], [3, 2]])
        assert_equal(oddsratio, 0)
        assert_equal(pvalue, 1)
        oddsratio, pvalue = fisher_exact([[0, 2], [6, 4]])
        assert_equal(oddsratio, 0)
        assert_allclose(pvalue, 5/11, rtol=1e-11)

    @pytest.mark.parametrize('parameters, rresult', data)
    def test_results_from_r(self, parameters, rresult):
        alternative = parameters.alternative.replace('.', '-')
        result = fisher_exact(parameters.table,
                              alternative=alternative,
                              confidence_level=parameters.confidence_level)
        assert_allclose(result.pvalue, rresult.pvalue, rtol=1e-11)
        if result.conditional_odds_ratio < 400:
            or_rtol = 5e-4
            ci_rtol = 2e-2
        else:
            or_rtol = 5e-2
            ci_rtol = 1e-1
        assert_allclose(result.conditional_odds_ratio,
                        rresult.conditional_odds_ratio, rtol=or_rtol)
        assert_allclose(result.conditional_odds_ratio_ci,
                        rresult.conditional_odds_ratio_ci, rtol=ci_rtol)
        # Also do a self-check for the conditional odds ratio.
        # With the computed conditional odds ratio as the noncentrality
        # parameter of the noncentral hypergeometric distribution with
        # parameters table.sum(), table[0].sum(), and table[:,0].sum() as
        # total, ngood and nsample, respectively, the mean of the distribution
        # should equal table[0, 0].
        cor = result.conditional_odds_ratio
        table = np.array(parameters.table)
        total = table.sum()
        ngood = table[0].sum()
        nsample = table[:, 0].sum()
        # Avoid the warning from log(cor) when cor == 0.
        lognc = np.log(cor) if cor > 0 else -np.inf
        nchg_mean = _nc_hypergeom_mean(lognc, total, ngood, nsample)
        assert_allclose(nchg_mean, table[0, 0], rtol=1e-13)

    @pytest.mark.slow
    def test_large_numbers(self):
        # Test with some large numbers. Regression test for #1401
        pvals = [5.55976657259e-11, 2.66568479113e-11, 1.36320302448e-11]
        for pval, num in zip(pvals, [75, 76, 77]):
            res = fisher_exact([[17704, 496], [1065, num]])[1]
            assert_allclose(res, pval, rtol=1e-5)

        res = fisher_exact([[18000, 80000], [20000, 90000]])[1]
        assert_allclose(res, 0.275148170406, rtol=1e-5)

    def test_raises(self):
        # test we raise an error for wrong shape of input.
        with pytest.raises(ValueError, match='shape'):
            fisher_exact(np.arange(6).reshape(2, 3))

    @pytest.mark.parametrize('table', [
        [[0, 0], [5, 10]],
        [[5, 10], [0, 0]],
        [[0, 5], [0, 10]],
        [[5, 0], [10, 0]],
    ])
    def test_row_or_col_zero(self, table):
        result = fisher_exact(table)
        assert_equal(result.pvalue, 1.0)
        assert_equal(result.sample_odds_ratio, np.nan)
        assert_equal(result.conditional_odds_ratio, np.nan)
        assert_equal(result.conditional_odds_ratio_ci, (0, np.inf))

    @pytest.mark.parametrize(
        'table, exact_p_less, exact_p_greater',
        [([[0, 2], [3, 0]], 0.1, 1.0),
         ([[1, 1], [2, 1]], 0.7, 0.9),
         ([[2, 0], [1, 2]], 1.0, 0.3),
         ([[0, 1], [2, 3]], 2/3, 1.0),
         ([[1, 0], [1, 4]], 1.0, 1/3)]
    )
    def test_less_greater(self, table, exact_p_less, exact_p_greater):
        # Some tables with simple exact values
        # (includes regression test for ticket #1568).
        p_less = fisher_exact(table, alternative="less").pvalue
        assert_allclose(p_less, exact_p_less, atol=0, rtol=1e-11)
        p_greater = fisher_exact(table, alternative="greater").pvalue
        assert_allclose(p_greater, exact_p_greater, atol=0, rtol=1e-11)

    def test_gh3014(self):
        # Check if issue #3014 has been fixed.
        # Before, this would have raised a ValueError.
        odds, pvalue = fisher_exact([[1, 2], [9, 84419233]])
        assert_allclose(odds, 4689957.388888889, rtol=1e-15)
        assert_allclose(pvalue, 3.55369167322892e-07, rtol=1e-13)
