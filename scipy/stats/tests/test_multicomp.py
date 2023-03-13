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

    # From Dunnett1995 comparing with R's DescTools: DunnettTest
    samples_2 = [
        [9.76, 8.80, 7.68, 9.36], [12.80, 9.68, 12.16, 9.20, 10.55]
    ]
    control_2 = [7.40, 8.50, 7.20, 8.24, 9.84, 8.32]
    pvalue_2 = [0.6201, 0.0058]

    samples_3 = [
        [55, 64, 64], [55, 49, 52], [50, 44, 41]
    ]
    control_3 = [55, 47, 48]
    pvalue_3 = [0.0364, 0.8966, 0.4091]

    # From Thomson and Short,
    # Mucociliary function in health, chronic obstructive airway disease,
    # and asbestosis, Journal of Applied Physiology, 1969. Table 1
    # Comparing with R's DescTools: DunnettTest
    samples_4 = [[3.8, 2.7, 4.0, 2.4], [2.8, 3.4, 3.7, 2.2, 2.0]]
    control_4 = [2.9, 3.0, 2.5, 2.6, 3.2]
    pvalue_4 = [0.5832, 0.9982]

    @pytest.mark.parametrize(
        'samples, control, pvalue',
        [
            (samples_1, control_1, pvalue_1),
            (samples_2, control_2, pvalue_2),
            (samples_3, control_3, pvalue_3),
            (samples_4, control_4, pvalue_4),
        ]
    )
    def test_unbalanced(self, samples, control, pvalue):
        rng = np.random.default_rng(11681140010308601919115036826969764808)

        res = stats.dunnett(*samples, control=control, random_state=rng)

        assert isinstance(res, DunnettResult)
        assert_allclose(res.pvalue, pvalue, atol=1e-3)

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
