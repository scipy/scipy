import copy

import numpy as np
import pytest
from numpy.testing import assert_allclose

from scipy import stats
from scipy.stats._multicomp import _pvalue_dunnett, DunnettResult


class TestDunnett:

    @pytest.mark.parametrize(
        'rho, n_groups, df, statistic, pvalue, alternative',
        [
            # From Dunnett1995
            # Tables 1a and 1b pages 1117-1118
            (0.5, 1, 10, 1.81, 0.05, "greater"),  # different than two-sided
            (0.5, 3, 10, 2.34, 0.05, "greater"),
            (0.5, 2, 30, 1.99, 0.05, "greater"),
            (0.5, 5, 30, 2.33, 0.05, "greater"),
            (0.5, 4, 12, 3.32, 0.01, "greater"),
            (0.5, 7, 12, 3.56, 0.01, "greater"),
            (0.5, 2, 60, 2.64, 0.01, "greater"),
            (0.5, 4, 60, 2.87, 0.01, "greater"),
            (0.5, 4, 60, [2.87, 2.21], [0.01, 0.05], "greater"),
            # Tables 2a and 2b pages 1119-1120
            (0.5, 1, 10, 2.23, 0.05, "two-sided"),  # two-sided
            (0.5, 3, 10, 2.81, 0.05, "two-sided"),
            (0.5, 2, 30, 2.32, 0.05, "two-sided"),
            (0.5, 3, 20, 2.57, 0.05, "two-sided"),
            (0.5, 4, 12, 3.76, 0.01, "two-sided"),
            (0.5, 7, 12, 4.08, 0.01, "two-sided"),
            (0.5, 2, 60, 2.90, 0.01, "two-sided"),
            (0.5, 4, 60, 3.14, 0.01, "two-sided"),
            (0.5, 4, 60, [3.14, 2.55], [0.01, 0.05], "two-sided"),
        ],
    )
    def test_critical_values(
        self, rho, n_groups, df, statistic, pvalue, alternative
    ):
        rng = np.random.default_rng(165250594791731684851746311027739134893)
        rho = np.full((n_groups, n_groups), rho)
        np.fill_diagonal(rho, 1)

        statistic = np.array(statistic)
        res = _pvalue_dunnett(
            rho=rho, df=df, statistic=statistic,
            alternative=alternative,
            rng=rng
        )
        assert_allclose(res, pvalue, atol=5e-3)

    # For the following tests, p-values were computed using Matlab, e.g.
    #     sample = [18.  15.  18.  16.  17.  15.  14.  14.  14.  15.  15....
    #               14.  15.  14.  22.  18.  21.  21.  10.  10.  11.  9....
    #               25.  26.  17.5 16.  15.5 14.5 22.  22.  24.  22.5 29....
    #               24.5 20.  18.  18.5 17.5 26.5 13.  16.5 13.  13.  13....
    #               28.  27.  34.  31.  29.  27.  24.  23.  38.  36.  25....
    #               38. 26.  22.  36.  27.  27.  32.  28.  31....
    #               24.  27.  33.  32.  28.  19. 37.  31.  36.  36....
    #               34.  38.  32.  38.  32....
    #               26.  24.  26.  25.  29. 29.5 16.5 36.  44....
    #               25.  27.  19....
    #               25.  20....
    #               28.];
    #     j = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
    #          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
    #          0 0 0 0...
    #          1 1 1 1 1 1 1 1 1 1 1 1 1 1 1...
    #          2 2 2 2 2 2 2 2 2...
    #          3 3 3...
    #          4 4...
    #          5];
    #     [~, ~, stats] = anova1(sample, j, "off");
    #     [results, ~, ~, gnames] = multcompare(stats, ...
    #     "CriticalValueType", "dunnett", ...
    #     "Approximate", false);
    #     tbl = array2table(results, "VariableNames", ...
    #     ["Group", "Control Group", "Lower Limit", ...
    #     "Difference", "Upper Limit", "P-value"]);
    #     tbl.("Group") = gnames(tbl.("Group"));
    #     tbl.("Control Group") = gnames(tbl.("Control Group"))

    # Matlab doesn't report the statistic, so the statistics were
    # computed using R multcomp `glht`, e.g.:
    #     library(multcomp)
    #     options(digits=16)
    #     control < - c(18.0, 15.0, 18.0, 16.0, 17.0, 15.0, 14.0, 14.0, 14.0,
    #                   15.0, 15.0, 14.0, 15.0, 14.0, 22.0, 18.0, 21.0, 21.0,
    #                   10.0, 10.0, 11.0, 9.0, 25.0, 26.0, 17.5, 16.0, 15.5,
    #                   14.5, 22.0, 22.0, 24.0, 22.5, 29.0, 24.5, 20.0, 18.0,
    #                   18.5, 17.5, 26.5, 13.0, 16.5, 13.0, 13.0, 13.0, 28.0,
    #                   27.0, 34.0, 31.0, 29.0, 27.0, 24.0, 23.0, 38.0, 36.0,
    #                   25.0, 38.0, 26.0, 22.0, 36.0, 27.0, 27.0, 32.0, 28.0,
    #                   31.0)
    #     t < - c(24.0, 27.0, 33.0, 32.0, 28.0, 19.0, 37.0, 31.0, 36.0, 36.0,
    #             34.0, 38.0, 32.0, 38.0, 32.0)
    #     w < - c(26.0, 24.0, 26.0, 25.0, 29.0, 29.5, 16.5, 36.0, 44.0)
    #     x < - c(25.0, 27.0, 19.0)
    #     y < - c(25.0, 20.0)
    #     z < - c(28.0)
    #
    #     groups = factor(rep(c("control", "t", "w", "x", "y", "z"),
    #                         times=c(length(control), length(t), length(w),
    #                                 length(x), length(y), length(z))))
    #     df < - data.frame(response=c(control, t, w, x, y, z),
    #                       group=groups)
    #     model < - aov(response
    #     ~group, data = df)
    #     test < - glht(model=model,
    #                   linfct=mcp(group="Dunnett"),
    #                   alternative="g")
    #     summary(test)
    # p-values agreed with those produced by Matlab to at least atol=1e-3

    # From Matlab's documentation on multcompare
    samples_1 = [
        [
            24.0, 27.0, 33.0, 32.0, 28.0, 19.0, 37.0, 31.0, 36.0, 36.0,
            34.0, 38.0, 32.0, 38.0, 32.0
        ],
        [26.0, 24.0, 26.0, 25.0, 29.0, 29.5, 16.5, 36.0, 44.0],
        [25.0, 27.0, 19.0],
        [25.0, 20.0],
        [28.0]
    ]
    control_1 = [
        18.0, 15.0, 18.0, 16.0, 17.0, 15.0, 14.0, 14.0, 14.0, 15.0, 15.0,
        14.0, 15.0, 14.0, 22.0, 18.0, 21.0, 21.0, 10.0, 10.0, 11.0, 9.0,
        25.0, 26.0, 17.5, 16.0, 15.5, 14.5, 22.0, 22.0, 24.0, 22.5, 29.0,
        24.5, 20.0, 18.0, 18.5, 17.5, 26.5, 13.0, 16.5, 13.0, 13.0, 13.0,
        28.0, 27.0, 34.0, 31.0, 29.0, 27.0, 24.0, 23.0, 38.0, 36.0, 25.0,
        38.0, 26.0, 22.0, 36.0, 27.0, 27.0, 32.0, 28.0, 31.0
    ]
    pvalue_1 = [4.727e-06, 0.022346, 0.97912, 0.99953, 0.86579]
    statistic_1 = [5.27356, 2.91270, 0.60831, 0.27002, 0.96637]
    # These p-values were computeded using R multcomp `glht`
    pvalue_1_twosided = [1e-4, 0.02237, 0.97913, 0.99953, 0.86583]
    pvalue_1_greater = [1e-4, 0.011217, 0.768500, 0.896991, 0.577211]
    pvalue_1_less = [1, 1, 0.99660, 0.98398, .99953]

    # From Dunnett1995 comparing with R's DescTools: DunnettTest
    samples_2 = [
        [9.76, 8.80, 7.68, 9.36], [12.80, 9.68, 12.16, 9.20, 10.55]
    ]
    control_2 = [7.40, 8.50, 7.20, 8.24, 9.84, 8.32]
    pvalue_2 = [0.6201, 0.0058]
    statistic_2 = [0.85703, 3.69375]

    samples_3 = [
        [55, 64, 64], [55, 49, 52], [50, 44, 41]
    ]
    control_3 = [55, 47, 48]
    pvalue_3 = [0.0364, 0.8966, 0.4091]
    statistic_3 = [3.09073, 0.56195, -1.40488]

    # From Thomson and Short,
    # Mucociliary function in health, chronic obstructive airway disease,
    # and asbestosis, Journal of Applied Physiology, 1969. Table 1
    # Comparing with R's DescTools: DunnettTest
    samples_4 = [[3.8, 2.7, 4.0, 2.4], [2.8, 3.4, 3.7, 2.2, 2.0]]
    control_4 = [2.9, 3.0, 2.5, 2.6, 3.2]
    pvalue_4 = [0.5832, 0.9982]
    statistic_4 = [0.90875, -0.05007]

    @pytest.mark.parametrize(
        'samples, control, pvalue, statistic',
        [
            (samples_1, control_1, pvalue_1, statistic_1),
            (samples_2, control_2, pvalue_2, statistic_2),
            (samples_3, control_3, pvalue_3, statistic_3),
            (samples_4, control_4, pvalue_4, statistic_4),
        ]
    )
    def test_unbalanced(self, samples, control, pvalue, statistic):
        rng = np.random.default_rng(11681140010308601919115036826969764808)

        res = stats.dunnett(*samples, control=control, random_state=rng)

        assert isinstance(res, DunnettResult)
        assert_allclose(res.statistic, statistic, atol=1e-5)
        assert_allclose(res.pvalue, pvalue, atol=1e-3)

    @pytest.mark.parametrize("alternative, pvalue",
                             [["two-sided", pvalue_1_twosided],
                              ["greater", pvalue_1_greater],
                              ["less", pvalue_1_less]])
    def test_unbalanced_one_sided(self, alternative, pvalue):
        # Check p-value against R multcomp `glht` with one-sided alternatives
        rng = np.random.default_rng(19191150368269697648081168114001030860)
        res = stats.dunnett(*self.samples_1, control=self.control_1,
                            alternative=alternative, random_state=rng)
        assert_allclose(res.statistic, self.statistic_1, atol=1e-5)
        assert_allclose(res.pvalue, pvalue, atol=1e-4)

    @pytest.mark.parametrize(
        'alternative',
        ['two-sided', 'less', 'greater']
    )
    def test_ttest_ind(self, alternative):
        # check that `dunnett` agrees with `ttest_ind`
        # when there are only two groups
        rng = np.random.default_rng(114184017807316971636137493526995620351)

        for _ in range(10):
            sample = rng.integers(-100, 100, size=(10,))
            control = rng.integers(-100, 100, size=(10,))

            res = stats.dunnett(
                sample, control=control,
                alternative=alternative, random_state=rng
            )
            ref = stats.ttest_ind(
                sample, control,
                alternative=alternative, random_state=rng
            )

            assert_allclose(res.statistic, ref.statistic, atol=1e-3)
            assert_allclose(res.pvalue, ref.pvalue, atol=1e-3)

    @pytest.mark.parametrize(
        'alternative, statistic, pvalue',
        [
            ('less', 'low', 0),
            ('less', 'high', 1),
            ('greater', 'low', 1),
            ('greater', 'high', 0),
            ('two-sided', 'low', 0),
            ('two-sided', 'high', 0)
        ]
    )
    def test_alternatives(self, alternative, statistic, pvalue):
        rng = np.random.default_rng(114184017807316971636137493526995620351)

        for _ in range(10):
            sample_1 = rng.integers(0, 20, size=(10,))
            sample_2 = rng.integers(80, 100, size=(10,))

            if statistic == 'high':
                sample = sample_2
                control = sample_1
            else:
                sample = sample_1
                control = sample_2

            res = stats.dunnett(
                sample, control=control,
                alternative=alternative, random_state=rng
            )
            assert_allclose(res.pvalue, pvalue, atol=1e-7)

    @pytest.mark.parametrize(
        'alternative, allowance, ci_low, ci_high',
        [
            ('two-sided', 11, [0, -9, -16], [22, 13, 6]),
            ('less', 9, [2, -7, -14], [np.inf, np.inf, np.inf]),
            ('greater', 9, [-np.inf, -np.inf, -np.inf], [20, 11, 4])
        ]
    )
    def test_allowance(self, alternative, allowance, ci_low, ci_high):
        # Example (a) from Dunnett1995
        rng = np.random.default_rng(189117774084579816190295271136455278291)
        samples = [
            [55, 64, 64],
            [55, 49, 52],
            [50, 44, 41]
        ]
        control = [55, 47, 48]

        res = stats.dunnett(
            *samples, control=control, alternative=alternative,
            random_state=rng
        )
        allowance_ = res._allowance(confidence_level=0.95)
        assert_allclose(allowance_, allowance, atol=1)

        assert res._ci is None
        assert res._ci_cl is None
        ci = res.confidence_interval(confidence_level=0.95)
        assert_allclose(ci.low, ci_low, atol=1)
        assert_allclose(ci.high, ci_high, atol=1)

        # re-run to use the cached value "is" to check id as same object
        assert res._ci is ci
        assert res._ci_cl == 0.95
        ci_ = res.confidence_interval(confidence_level=0.95)
        assert ci_ is ci

        # check some str output
        res_str = str(res)
        assert '(Sample 2 - Control)' in res_str
        assert '95.0%' in res_str

        if alternative == 'less':
            assert 'inf' in res_str
            assert 'at least Lower' in res_str
            assert '-13.' in res_str
        elif alternative == 'greater':
            assert '-inf' in res_str
            assert 'at most Upper' in res_str
            assert '19.' in res_str
        else:
            assert 'inf' not in res_str
            assert 'between' in res_str
            assert '21.' in res_str

    def test_raises(self):
        samples = [
            [55, 64, 64],
            [55, 49, 52],
            [50, 44, 41]
        ]
        control = [55, 47, 48]

        # alternative
        with pytest.raises(ValueError, match="alternative must be"):
            stats.dunnett(*samples, control=control, alternative='bob')

        # 2D for a sample
        samples_ = copy.deepcopy(samples)
        samples_[0] = [samples_[0]]
        with pytest.raises(ValueError, match="must be 1D arrays"):
            stats.dunnett(*samples_, control=control)

        # 2D for control
        control_ = copy.deepcopy(control)
        control_ = [control_]
        with pytest.raises(ValueError, match="must be 1D arrays"):
            stats.dunnett(*samples, control=control_)

        # No obs in a sample
        samples_ = copy.deepcopy(samples)
        samples_[1] = []
        with pytest.raises(ValueError, match="at least 1 observation"):
            stats.dunnett(*samples_, control=control)

        # No obs in control
        control_ = []
        with pytest.raises(ValueError, match="at least 1 observation"):
            stats.dunnett(*samples, control=control_)

        res = stats.dunnett(*samples, control=control)
        with pytest.raises(ValueError, match="Confidence level must"):
            res.confidence_interval(confidence_level=3)
