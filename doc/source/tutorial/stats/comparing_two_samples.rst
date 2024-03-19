Comparing two samples
---------------------

In the following, we are given two samples, which can come either from the
same or from different distribution, and we want to test whether these
samples have the same statistical properties.


Comparing means
^^^^^^^^^^^^^^^

Test with sample with identical means:

    >>> rvs1 = stats.norm.rvs(loc=5, scale=10, size=500)
    >>> rvs2 = stats.norm.rvs(loc=5, scale=10, size=500)
    >>> stats.ttest_ind(rvs1, rvs2)
    Ttest_indResult(statistic=-0.5489036175088705, pvalue=0.5831943748663959)  # random

Test with sample with different means:

    >>> rvs3 = stats.norm.rvs(loc=8, scale=10, size=500)
    >>> stats.ttest_ind(rvs1, rvs3)
    Ttest_indResult(statistic=-4.533414290175026, pvalue=6.507128186389019e-06)  # random

Kolmogorov-Smirnov test for two samples ks_2samp
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the example, where both samples are drawn from the same distribution,
we cannot reject the null hypothesis, since the pvalue is high

    >>> stats.ks_2samp(rvs1, rvs2)
    KstestResult(statistic=0.026, pvalue=0.9959527565364388)  # random

In the second example, with different location, i.e., means, we can
reject the null hypothesis, since the pvalue is below 1%

    >>> stats.ks_2samp(rvs1, rvs3)
    KstestResult(statistic=0.114, pvalue=0.00299005061044668)  # random
