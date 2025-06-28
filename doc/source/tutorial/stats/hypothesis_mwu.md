---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.4
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

+++ {"tags": ["jupyterlite_sphinx_strip"]}

```{eval-rst}
.. notebooklite:: hypothesis_mwu.md
   :new_tab: True
```

(hypothesis_mwu)=

+++

# Mann-Whitney U Test

Consider the following data from [^1], which describes a study about the effect of diet protein choice (meat or fish) on blood plasma cholesterol concentrations. 

The `x` and `y` arrays below record measurements of cholesterol concentrations (mmol/L) from the two groups.

```{code-cell} ipython3
import numpy as np
# Fish eaters
x = [5.42, 5.86, 6.16, 6.55, 6.80, 7.00, 7.11]
# Meat eaters
y = [6.51, 7.56, 7.61, 7.84, 11.50]
```

We can analyze the data using the Mann-Whitney U test (`{class}`scipy.stats.mannwhitneyu`), a nonparametric statistic sensitive to differences between the distributions underlying the groups.

```{code-cell} ipython3
from scipy import stats
res = stats.mannwhitneyu(x, y)
res.statistic
```

The minimum value of the statistic, achieved when the maximum of `x` is less than the minimum of `y`, is 0. The maximum value of the statistic is $n_1 n_2$, the product of the number of observations in the two groups, and it is achieved when the minimum of `x` is greater than the maximum of `y`. The value of the statistic tends to be moderate (close to $n_1 n_2 / 2$) for samples drawn from the same statistical distribution, and low or high when the samples are drawn from different distributions.

The test is performed by comparing the observed value of the statistic against the null distribution: the distribution of statistic values derived under the null hypothesis that the samples were drawn from the same statistical distribution.

Although the set of all observations (`5.42, 5.86, 6.16, 6.55, 6.80, 7.00, 7.11, 6.51, 7.56, 7.61, 7.84, 11.50`) happened to be divided between the fish- and meat-eating groups as indicated, under the null hypothesis, these numerical values were equally likely to have been partitioned between the two groups in many other ways. Specifically, we can form the *null distribution* - the distribution of the test statistic when the null hypothesis is true - by partitioning the observations into groups of sizes seven and five in all possible ways, and computing the test statistic under each partition of the data.

```{code-cell} ipython3
from itertools import combinations
data = set(x + y)
n_1, n_2 = len(x), len(y)
null_distribution = []
# Compute the statistic for each partition of the data
for xi in combinations(data, n_1):  # `xi` is a group of size `n_1`
    yi = data - set(xi)  # The remaining data is in `yi`
    statistic = stats.mannwhitneyu(list(xi), list(yi)).statistic
    null_distribution.append(statistic)
```

```{code-cell} ipython3
# Plot the distribution of statistic values
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 5))
def plot(ax):  # we'll re-use this
    bins = np.arange(n_1*n_2 + 2) - 0.5  # center bins around values of statistic
    n, bins, patches = ax.hist(null_distribution, bins, density=True)
    ax.set_title(f"Mann-Whitney U Null Distribution (${n_1=}, {n_2=}$)")
    ax.set_xlabel("statistic")
    ax.set_ylabel("frequency density")
    return bins, patches
plot(ax)
plt.show()
```

The comparison is quantified by the p-value: the proportion of values in the null distribution as extreme or more extreme than the observed value of the statistic. In a two-sided test, the extremety of a statistic value is characterized by its absolute distance from the mean, $\frac{n_1 n_2}{2}$; more extreme values are further from the mean.

```{code-cell} ipython3
null_distribution = np.asarray(null_distribution)
mean = n_1 * n_2 / 2
as_extreme = abs(null_distribution - mean) >= abs(res.statistic - mean)
pvalue = np.count_nonzero(as_extreme) / null_distribution.size
assert pvalue == res.pvalue
pvalue
```

```{code-cell} ipython3
fig, ax = plt.subplots(figsize=(8, 5))
bins, patches = plot(ax)
null_unique = np.arange(2*mean + 1)
as_extreme = abs(null_unique - mean) >= abs(res.statistic - mean)
for patch, bin_edge, extreme in zip(patches, bins[:-1], as_extreme):  # bins[:-1] gives left edges of bars
    if extreme:
        patch.set_facecolor('red')
annotation = (f'p-value={pvalue:.4f}\n(shaded area)')
props = dict(facecolor='black', width=1, headwidth=5, headlength=8)
_ = ax.annotate(annotation, (31, 0.003), (20, 0.01), arrowprops=props)
_ = ax.annotate("", (4, 0.003), (19, 0.01), arrowprops=props)
```

If the p-value is "small" - that is, if there is low probability that chance alone would result in a statistic value at least as extreme as the one observed - this may be taken as evidence against the null hypothesis in favor of the alternative: the distributions of cholesterol concentrations for fish-eaters and meat-eaters are not the same. Note that:

- The inverse is not true; that is, the test is not used to provide evidence for the null hypothesis.
- The threshold for values that will be considered "small" is a choice that should be made before the data is analyzed [^3] with consideration of the risks of both false positives (incorrectly rejecting the null hypothesis) and false negatives (failure to reject a false null hypothesis).
- Small p-values are not evidence for a large effect; rather, they can only provide evidence for a "significant" effect, meaning that they are unlikely to have occurred under the null hypothesis.

+++

Suppose that before performing the experiment, the researchers had reason to believe that a fish diet would tend to lower plasma cholesterol concentrations. In that case, they might wish to consider only values of the statistic *less* than the observed value to be more extreme.

```{code-cell} ipython3
as_extreme = (null_distribution - mean) <= (res.statistic - mean)
pvalue = np.count_nonzero(as_extreme) / len(null_distribution)
res = stats.mannwhitneyu(x, y, alternative='less')
assert pvalue == res.pvalue
pvalue
```

In doing so, they would test the null hypothesis against a one-sided alternative: that plasma cholesterol concentrations from the population of fish eaters tend to be lower than plasma cholesterol concentrations from the population of meat eaters.

+++

The Mann-Whitney U test was invented at time when computational resources were much more limited than they are now. The statistic involves the use of ranks rather than the original data because it leads to a simplified method for computing the exact null distribution. This advantage is often irrelevant today, yet the use of ranks still has a few major downsides[^1]:

- It discards information, which may lead to lower test power, that is, ability to detect a real deviation from the null hypothesis.
- It complicates interpretation of the results. For instance, practitioners often wish to use the Mann-Whitney U test as a test for difference in population medians, but this is often not a valid interpretation.

In many cases, then, it may be preferable to use the same partitioning approach with a different statistic. If we are only interested in whether the population locations are distinct, for instance, a reasonable statistic might be the absolute value of the difference in sample means.

```{code-cell} ipython3
def statistic(x, y, axis):
    return np.abs(np.mean(x, axis=axis) - np.mean(y, axis=axis))
```

The function `scipy.stats.permutation_test` computes the statistic for all "permutations" (partitions) of the data. In this case, we consider statistic values *greater* than the observed statistic to be more extreme, so we pass argument `alternative='greater'`.

```{code-cell} ipython3
res = stats.permutation_test((x, y), statistic, alternative='greater')
res.pvalue
```

If we want to test the one-sided alternative that the mean of the distribution of plasma cholesterol concentrations of fish-eaters is less than the mean of the distribution of plasma cholesterol concentrations of meat-eaters, then an appropriate statistic would simply be the difference in means, without taking the absolute value). Because we would consider value of the statistic *less* than that observed value to be more extreme, we pass argument `alternative='less'`.

```{code-cell} ipython3
def statistic(x, y, axis):
    return np.mean(x, axis=axis) - np.mean(y, axis=axis)
    res = stats.permutation_test((x, y), statistic, alternative='less')
res.pvalue
```

Note that these p-values are substantially lower than the p-values provided by the Mann-Whitney U test, suggesting that this test may be more sensitive to differences in the underlying distributions. Also, the statistic is much easier to interpret than the one used by the Mann Whitney U test because it uses the data rather than the ranks. See [^1] for more information.

+++

## References

[^1]: Ludbrook, J. and Dudley, H. (1998). Why Permutation Tests Are Superior to t and F Tests in Biomedical Research. The American Statistician, 52(2), 127-132.

[^2]: Hollander, M., Wolfe, D. A., & Chicken, E. (2013). Nonparametric
statistical methods. John Wiley & Sons.

[^3]: Phipson, B. and Smyth, G. K. (2010). "Permutation P-values Should Never Be
Zero: Calculating Exact P-values When Permutations Are Randomly Drawn."
Statistical Applications in Genetics and Molecular Biology 9.1.
