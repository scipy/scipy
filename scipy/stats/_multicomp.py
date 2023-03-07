from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from scipy import stats
from scipy.stats._qmc import check_random_state

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
    """Dunnett's test.

    Parameters
    ----------
    sample1, sample2, ... : 1D array_like
        The sample measurements for each experiment group.
    control : 1D array_like
        The sample measurements for the control group.
    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis.
        The following options are available (default is 'two-sided'):

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
            at index ``(i,)`` is the statistic for the comparison between
            groups ``i`` and the control.
        pvalue : float ndarray
            The computed p-value of the test for each comparison. The element
            at index ``(i,)`` is the p-value for the comparison between
            groups ``i`` and the control.

    See Also
    --------
    tukey_hsd : performs pairwise comparison of means.

    Notes
    -----
    Dunnett's test [1]_ compares the means of multiple experiment groups
    against a control group.
    `tukey_hsd` instead, performs pairwise comparison of means.
    It means Dunnett's test performs fewer tests, hence there is less p-value
    adjustment which makes the test more powerful.
    It should be preferred when there is control group.

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
    .. [2] K.S. Kwong, W. Liu. "Calculation of critical values for Dunnett
       and Tamhane's step-up multiple test procedure."
       Statistics and Probability Letters, 49, 411-416,
       :doi:`10.1016/S0167-7152(00)00076-6`, 2000.

    Examples
    --------
    In [1]_, the influence of drugs on blood count measurements on three groups
    of animal is investigated.

    The following table summarizes the results of the experiment in which
    two groups received different drug, and one group acted as a control.
    Blood counts (in millions of cells per cubic millimeter) were recorded::

         Control      Drug A      Drug B
           7.40        9.76        12.80
           8.50        8.80         9.68
           7.20        7.68        12.16
           8.24        9.36         9.20
           9.84                    10.55
           8.32

    >>> import numpy as np
    >>> control = np.array([7.40, 8.50, 7.20, 8.24, 9.84, 8.32])
    >>> drug_a = np.array([9.76, 8.80, 7.68, 9.36])
    >>> drug_b = np.array([12.80, 9.68, 12.16, 9.20, 10.55])

    The `dunnett` statistic is sensitive to the difference in means between
    the samples.

    We would like to see if the means between any of the groups are
    significantly different. First, visually examine a box and whisker plot.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(1, 1)
    >>> ax.boxplot([control, drug_a, drug_b])
    >>> ax.set_xticklabels(["Control", "Drug A", "Drug B"])  # doctest: +SKIP
    >>> ax.set_ylabel("mean")  # doctest: +SKIP
    >>> plt.show()

    From the box and whisker plot, we can see overlap in the interquartile
    ranges between the control group and the group from drug A.
    We can apply the `dunnett`
    test to determine if the difference between means is significant. We
    set a significance level of .05 to reject the null hypothesis.

    >>> from scipy.stats import dunnett
    >>> res = dunnett(drug_a, drug_b, control=control)
    >>> res.pvalue
    array([0.47773146, 0.00889328])  # random

    The null hypothesis is that each group has the same mean. The p-value for
    comparisons between ``control`` and ``drug_b`` do not exceed .05,
    so we reject the null hypothesis that they
    have the same means. The p-value of the comparison between ``control``
    and ``drug_a`` exceeds .05, so we accept the null hypothesis that there
    is not a significant difference between their means.

    """
    samples_, control, rng = iv_dunnett(
        samples=samples, control=control, random_state=random_state
    )

    rho, df, n_group = params_dunnett(samples=samples_, control=control)

    statistic = np.array([
        stats.ttest_ind(
            obs_, control, alternative=alternative, random_state=rng
        ).statistic
        for obs_ in samples_
    ])

    pvalue = pvalue_dunnett(
        rho=rho, df=df,
        statistic=statistic, alternative=alternative,
        rng=rng
    )

    return DunnettResult(
        statistic=statistic, pvalue=pvalue,
    )


def iv_dunnett(
    samples: Sequence[npt.ArrayLike],
    control: npt.ArrayLike,
    random_state: SeedType
) -> tuple[list[np.ndarray], np.ndarray, SeedType]:
    """Input validation for Dunnett's test."""
    rng = check_random_state(random_state)

    ndim_msg = "Control and samples groups must be 1D arrays"
    n_obs_msg = "Control and samples groups must have at least 1 observation"

    control = np.asarray(control)
    samples_ = [np.asarray(sample) for sample in samples]

    # samples checks
    samples_control: list[np.ndarray] = samples_ + [control]
    for sample in samples_control:
        if sample.ndim > 1:
            raise ValueError(ndim_msg)

        if len(sample) < 1:
            raise ValueError(n_obs_msg)

    return samples_, control, rng


def params_dunnett(
    samples: list[np.ndarray], control: np.ndarray
) -> tuple[np.ndarray, int, int]:
    """Specific parameters for Dunnett's test.

    Covariance matrix depends on the number of observations in each group:

    - All groups are equals (including the control), ``rho_ij=0.5`` except for
      the diagonal which is 1.
    - All groups but the control are equal, balanced design.
    - Groups are not equal, unbalanced design.

    Degree of freedom is the number of observations minus the number of groups
    including the control.
    """
    n_n_obs = np.array([len(obs_) for obs_ in samples])

    # From Dunnett1955 p. 1100 d.f. = (sum N)-(p+1)
    n_obs = n_n_obs.sum()
    n_control = len(control)
    n = n_obs + n_control
    n_groups = len(samples)
    df = n - n_groups - 1

    # rho_ij = 1/sqrt((N0/Ni+1)(N0/Nj+1))
    rho = n_control/n_n_obs + 1
    rho = 1/np.sqrt(rho[:, None] * rho[None, :])
    np.fill_diagonal(rho, 1)

    return rho, df, n_groups


def pvalue_dunnett(
    rho: np.ndarray, df: int, statistic: np.ndarray,
    alternative: Literal['two-sided', 'less', 'greater'],
    rng: SeedType = None
) -> np.ndarray:
    """pvalue from Dunnett critical value.

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
