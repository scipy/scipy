from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from scipy import stats
from scipy.stats._qmc import check_random_state
from scipy.stats._stats_py import _var

if TYPE_CHECKING:
    import numpy.typing as npt
    from scipy._lib._util import SeedType
    from typing import Literal, Sequence  # noqa: UP035


__all__ = [
    'dunnett'
]


@dataclass
class DunnettResult:
    statistic: np.ndarray
    pvalue: np.ndarray


def dunnett(
    *samples: npt.ArrayLike,
    control: npt.ArrayLike,
    alternative: Literal['two-sided', 'less', 'greater'] = "two-sided",
    random_state: SeedType = None
) -> DunnettResult:
    """Dunnett's test: multiple comparisons of means against a control group.

    This is an implementation of Dunnett's original, single-step test as
    described in [1]_.

    Parameters
    ----------
    sample1, sample2, ... : 1D array_like
        The sample measurements for each experimental group.
    control : 1D array_like
        The sample measurements for the control group.
    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis.

        The null hypothesis is that the means of the distributions underlying
        the samples and control are equal. The following alternative
        hypotheses are available (default is 'two-sided'):

        * 'two-sided': the means of the distributions underlying the samples
          and control are unequal.
        * 'less': the means of the distributions underlying the samples
          are less than the mean of the distribution underlying the control.
        * 'greater': the means of the distributions underlying the
          samples are greater than the mean of the distribution underlying
          the control.
    random_state : {None, int, `numpy.random.Generator`}, optional
        If `random_state` is an int or None, a new `numpy.random.Generator` is
        created using ``np.random.default_rng(random_state)``.
        If `random_state` is already a ``Generator`` instance, then the
        provided instance is used.

        The random number generator is used to control the randomized
        Quasi-Monte Carlo integration of the multivariate-t distribution.

    Returns
    -------
    res : DunnettResult
        An object containing attributes:

        statistic : float ndarray
            The computed statistic of the test for each comparison. The element
            at index ``i`` is the statistic for the comparison between
            groups ``i`` and the control.
        pvalue : float ndarray
            The computed p-value of the test for each comparison. The element
            at index ``i`` is the p-value for the comparison between
            group ``i`` and the control.

    See Also
    --------
    tukey_hsd : performs pairwise comparison of means.

    Notes
    -----
    Like the independent-sample t-test, Dunnett's test [1]_ is used to make
    inferences about the means of distributions from which samples were drawn.
    However, when multiple t-tests are performed at a fixed significance level,
    the "family-wise error rate" - the probability of incorrectly rejecting the
    null hypothesis in at least one test - will exceed the significance level.
    Dunnett's test is designed to perform multiple comparisons while
    controlling the family-wise error rate.

    Dunnett's test compares the means of multiple experimental groups
    against a single control group. Tukey's Honestly Significant Difference Test
    is another multiple-comparison test that controls the family-wise error
    rate, but `tukey_hsd` performs *all* pairwise comparisons between groups.
    When pairwise comparisons between experimental groups are not needed,
    Dunnett's test is preferable due to its higher power.


    The use of this test relies on several assumptions.

    1. The observations are independent within and among groups.
    2. The observations within each group are normally distributed.
    3. The distributions from which the samples are drawn have the same finite
       variance.

    References
    ----------
    .. [1] Charles W. Dunnett. "A Multiple Comparison Procedure for Comparing
       Several Treatments with a Control."
       Journal of the American Statistical Association, 50:272, 1096-1121,
       :doi:`10.1080/01621459.1955.10501294`, 1955.

    Examples
    --------
    In [1]_, the influence of drugs on blood count measurements on three groups
    of animal is investigated.

    The following table summarizes the results of the experiment in which
    two groups received different drugs, and one group acted as a control.
    Blood counts (in millions of cells per cubic millimeter) were recorded::

    >>> import numpy as np
    >>> control = np.array([7.40, 8.50, 7.20, 8.24, 9.84, 8.32])
    >>> drug_a = np.array([9.76, 8.80, 7.68, 9.36])
    >>> drug_b = np.array([12.80, 9.68, 12.16, 9.20, 10.55])

    We would like to see if the means between any of the groups are
    significantly different. First, visually examine a box and whisker plot.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(1, 1)
    >>> ax.boxplot([control, drug_a, drug_b])
    >>> ax.set_xticklabels(["Control", "Drug A", "Drug B"])  # doctest: +SKIP
    >>> ax.set_ylabel("mean")  # doctest: +SKIP
    >>> plt.show()

    Note the overlapping interquartile ranges of the drug A group and control
    group and the apparent separation between the drug B group and control
    group.

    Next, we will use Dunnett's test to assess whether the difference
    between group means is significant while controlling the family-wise error
    rate: the probability of making any false discoveries.
    Let the null hypothesis be that the experimental groups have the same
    mean as the control and the alternative be that an experimental group does
    not have the same mean as the control. We will consider a 5% family-wise
    error rate to be acceptable, and therefore we choose 0.05 as the threshold
    for significance.

    >>> from scipy.stats import dunnett
    >>> res = dunnett(drug_a, drug_b, control=control)
    >>> res.pvalue
    array([0.62004941, 0.0059035 ])  # may vary

    The p-value corresponding with the comparison between group A and control
    exceeds 0.05, so we do not reject the null hypothesis for that comparison.
    However, the p-value corresponding with the comparison between group B
    and control is less than 0.05, so we consider the experimental results
    to be evidence against the null hypothesis in favor of the alternative:
    group B has a different mean than the control group.

    """
    samples_, control_, rng = _iv_dunnett(
        samples=samples, control=control,
        alternative=alternative, random_state=random_state
    )

    rho, df, n_group = _params_dunnett(samples=samples_, control=control_)

    statistic = _statistic_dunnett(samples_, control_, df)

    pvalue = _pvalue_dunnett(
        rho=rho, df=df, statistic=statistic, alternative=alternative, rng=rng
    )

    return DunnettResult(
        statistic=statistic, pvalue=pvalue,
    )


def _iv_dunnett(
    samples: Sequence[npt.ArrayLike],
    control: npt.ArrayLike,
    alternative: Literal['two-sided', 'less', 'greater'],
    random_state: SeedType
) -> tuple[list[np.ndarray], np.ndarray, SeedType]:
    """Input validation for Dunnett's test."""
    rng = check_random_state(random_state)

    if alternative not in {'two-sided', 'less', 'greater'}:
        raise ValueError(
            "alternative must be 'less', 'greater' or 'two-sided'"
        )

    ndim_msg = "Control and samples groups must be 1D arrays"
    n_obs_msg = "Control and samples groups must have at least 1 observation"

    control = np.asarray(control)
    samples_ = [np.asarray(sample) for sample in samples]

    # samples checks
    samples_control: list[np.ndarray] = samples_ + [control]
    for sample in samples_control:
        if sample.ndim > 1:
            raise ValueError(ndim_msg)

        if sample.size < 1:
            raise ValueError(n_obs_msg)

    return samples_, control, rng


def _params_dunnett(
    samples: list[np.ndarray], control: np.ndarray
) -> tuple[np.ndarray, int, int]:
    """Specific parameters for Dunnett's test.

    Degree of freedom is the number of observations minus the number of groups
    including the control.
    """
    n_samples = np.array([sample.size for sample in samples])

    # From [1] p. 1100 d.f. = (sum N)-(p+1)
    n_sample = n_samples.sum()
    n_control = control.size
    n = n_sample + n_control
    n_groups = len(samples)
    df = n - n_groups - 1

    # From [1] p. 1103 rho_ij = 1/sqrt((N0/Ni+1)(N0/Nj+1))
    rho = n_control/n_samples + 1
    rho = 1/np.sqrt(rho[:, None] * rho[None, :])
    np.fill_diagonal(rho, 1)

    return rho, df, n_groups


def _statistic_dunnett(
    samples: list[np.ndarray], control: np.ndarray, df: int
) -> np.ndarray:
    """Statistic of Dunnett's test.

    Computation based on the original single-step test from [1].
    """
    n_control = control.size
    n_samples = np.array([sample.size for sample in samples])
    mean_control = np.mean(control)
    mean_samples = np.array([np.mean(sample) for sample in samples])
    all_samples = [control] + samples
    all_means = np.concatenate([[mean_control], mean_samples])

    # Variance estimate s^2 from [1] Eq. 1
    s2 = np.sum([_var(sample, mean=mean)*sample.size
                 for sample, mean in zip(all_samples, all_means)]) / df

    # z score inferred from [1] unlabeled equation after Eq. 1
    z = (mean_samples - mean_control) / np.sqrt(1/n_samples + 1/n_control)

    return z / np.sqrt(s2)


def _pvalue_dunnett(
    rho: np.ndarray, df: int, statistic: np.ndarray,
    alternative: Literal['two-sided', 'less', 'greater'],
    rng: SeedType = None
) -> np.ndarray:
    """pvalue from the multivariate t-distribution.

    Critical values come from the multivariate student-t distribution.
    """
    statistic = statistic.reshape(-1, 1)

    mvt = stats.multivariate_t(shape=rho, df=df, seed=rng)
    if alternative == "two-sided":
        statistic = abs(statistic)
        pvalue = 1 - mvt.cdf(statistic, lower_limit=-statistic)
    elif alternative == "greater":
        pvalue = 1 - mvt.cdf(statistic, lower_limit=-np.inf)
    else:
        pvalue = 1 - mvt.cdf(np.inf, lower_limit=statistic)

    return pvalue
