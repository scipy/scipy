import pytest
import numpy as np
from numpy.testing import assert_equal, assert_allclose
from scipy.stats import nchypergeom_fisher, hypergeom
from scipy.stats._odds_ratio import odds_ratio
from .fisher_exact_results_from_r import data


class TestOddsRatio:

    @pytest.mark.parametrize('parameters, rresult', data)
    def test_results_from_r(self, parameters, rresult):
        alternative = parameters.alternative.replace('.', '-')
        result = odds_ratio(parameters.table, alternative=alternative)
        # The results computed by R are not very precise.
        if result.odds_ratio < 400:
            or_rtol = 5e-4
            ci_rtol = 2e-2
        else:
            or_rtol = 5e-2
            ci_rtol = 1e-1
        assert_allclose(result.odds_ratio,
                        rresult.conditional_odds_ratio, rtol=or_rtol)
        assert_allclose(result.pvalue, rresult.pvalue, rtol=1e-6)
        ci = result.odds_ratio_ci(parameters.confidence_level)
        assert_allclose((ci.low, ci.high), rresult.conditional_odds_ratio_ci,
                        rtol=ci_rtol)

        # Also do a self-check for the conditional odds ratio.
        # With the computed conditional odds ratio as the noncentrality
        # parameter of the noncentral hypergeometric distribution with
        # parameters table.sum(), table[0].sum(), and table[:,0].sum() as
        # total, ngood and nsample, respectively, the mean of the distribution
        # should equal table[0, 0].
        cor = result.odds_ratio
        table = np.array(parameters.table)
        total = table.sum()
        ngood = table[0].sum()
        nsample = table[:, 0].sum()
        # nchypergeom_fisher does not allow the edge cases where the
        # noncentrality parameter is 0 or inf, so handle those values
        # separately here.
        if cor == 0:
            nchg_mean = hypergeom.support(total, ngood, nsample)[0]
        elif cor == np.inf:
            nchg_mean = hypergeom.support(total, ngood, nsample)[1]
        else:
            nchg_mean = nchypergeom_fisher.mean(total, ngood, nsample, cor)
        assert_allclose(nchg_mean, table[0, 0], rtol=1e-13)

    @pytest.mark.parametrize('table', [
        [[0, 0], [5, 10]],
        [[5, 10], [0, 0]],
        [[0, 5], [0, 10]],
        [[5, 0], [10, 0]],
    ])
    def test_row_or_col_zero(self, table):
        result = odds_ratio(table)
        assert_equal(result.odds_ratio, np.nan)
        ci = result.odds_ratio_ci()
        assert_equal((ci.low, ci.high), (0, np.inf))
