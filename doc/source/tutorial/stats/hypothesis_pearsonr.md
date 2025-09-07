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
.. notebooklite:: hypothesis_pearsonr.md
   :new_tab: True
```

(hypothesis_pearsonr)=

+++

# Pearson's Correlation

Consider the following data from [^1], which studied the relationship between free proline (an amino acid) and total collagen (a protein often found in connective tissue) in unhealthy human livers.

The `x` and `y` arrays below record measurements of the two compounds. The observations are paired: each free proline measurement was taken from the same liver as the total collagen measurement at the same index.

```{code-cell} ipython3
import numpy as np
# total collagen (mg/g dry weight of liver)
x = np.array([7.1, 7.1, 7.2, 8.3, 9.4, 10.5, 11.4])
# free proline (μ mole/g dry weight of liver)
y = np.array([2.8, 2.9, 2.8, 2.6, 3.5, 4.6, 5.0])
```

```{code-cell} ipython3
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(x, y, '.')
ax.set_xlabel("total collagen (mg/g)")
ax.set_ylabel("free proline (μ mole/g)")
plt.show()
```

These data were analyzed in [^2] using Spearman’s correlation coefficient, a statistic sensitive to monotonic correlation between the samples. Here, we will analyze the data using Pearson's correlation coefficient (`{class}`scipy.stats.pearsonr`) which is sensitive to linear correlation.

```{code-cell} ipython3
from scipy import stats
res = stats.pearsonr(x, y)
res.statistic
```

The value of this statistic tends to be high (close to 1) for samples with a strongly positive linear correlation, low (close to -1) for samples with a strongly negative linear correlation, and small in magnitude (close to zero) for samples with weak linear correlation.

The test is performed by comparing the observed value of the statistic against the null distribution: the distribution of statistic values derived under the null hypothesis that total collagen and free proline measurements are drawn from independent normal distributions.

Under the null hypothesis, the population correlation coefficient is zero, and the sample correlation coefficent follows the beta distribution on the interval $(-1, 1)$ with shape parameters $a = b = \frac{n}{2}-1$, where $n$ is the number of observations in each sample.

```{code-cell} ipython3
n = len(x)  # len(x) == len(y)
a = b = n/2 - 1  # shape parameter
loc, scale = -1, 2  # support is (-1, 1)
dist = stats.beta(a=a, b=b, loc=loc, scale=scale)
r_vals = np.linspace(-1, 1, 1000)
pdf = dist.pdf(r_vals)
fig, ax = plt.subplots(figsize=(8, 5))
def plot(ax):  # we'll re-use this
    ax.plot(r_vals, pdf)
    ax.set_title("Pearson's R Test Null Distribution")
    ax.set_xlabel("statistic")
    ax.set_ylabel("probability density")
plot(ax)
plt.show()
```

The comparison is quantified by the p-value: the proportion of values in the null distribution as extreme or more extreme than the observed value of the statistic. In a two-sided test in which the statistic is positive, elements of the null distribution greater than the transformed statistic and elements of the null distribution less than the negative of the observed statistic are both considered "more extreme".

```{code-cell} ipython3
fig, ax = plt.subplots(figsize=(8, 5))
plot(ax)
rs = res.statistic  # original statistic
pvalue = dist.cdf(-rs) + dist.sf(rs)
annotation = (f'p-value={pvalue:.4f}\n(shaded area)')
props = dict(facecolor='black', width=1, headwidth=5, headlength=8)
_ = ax.annotate(annotation, (rs, 0.01), (0.25, 0.1), arrowprops=props)
i = r_vals >= rs
ax.fill_between(r_vals[i], y1=0, y2=pdf[i], color='C0')
i = r_vals <= -rs
ax.fill_between(r_vals[i], y1=0, y2=pdf[i], color='C0')
ax.set_xlim(-1.1, 1.1)
ax.set_ylim(0, 0.9)
plt.show()
```

```{code-cell} ipython3
res.pvalue  # two-sided p-value
```

If the p-value is "small" - that is, if there is a low probability of sampling data from independent normal distributions that produces such an extreme value of the statistic - this may be taken as evidence against the null hypothesis in favor of the alternative: the distribution of total collagen and free proline are not independent. Note that:

- The inverse is not true; that is, the test is not used to provide evidence for the null hypothesis.
- The threshold for values that will be considered "small" is a choice that should be made before the data is analyzed [^3] with consideration of the risks of both false positives (incorrectly rejecting the null hypothesis) and false negatives (failure to reject a false null hypothesis).
- Small p-values are not evidence for a large effect; rather, they can only provide evidence for a "significant" effect, meaning that they are unlikely to have occurred under the null hypothesis.

Suppose that before performing the experiment, the authors had reason to predict a positive linear correlation between the total collagen and free proline measurements, and that they had chosen to assess the plausibility of the null hypothesis against a one-sided alternative: free proline has a positive linear correlation with total collagen. In this case, only those values in the null distribution that are as great or greater than the observed statistic would be considered more extreme.

```{code-cell} ipython3
res = stats.pearsonr(x, y, alternative='greater')
res.statistic
fig, ax = plt.subplots(figsize=(8, 5))
plot(ax)
pvalue = dist.sf(rs)
annotation = (f'p-value={pvalue:.6f}\n(shaded area)')
props = dict(facecolor='black', width=1, headwidth=5, headlength=8)
_ = ax.annotate(annotation, (rs, 0.01), (0.86, 0.04), arrowprops=props)
i = r_vals >= rs
ax.fill_between(r_vals[i], y1=0, y2=pdf[i], color='C0')
ax.set_xlim(0.85, 1)
ax.set_ylim(0, 0.1)
plt.show()
```

```{code-cell} ipython3
res.pvalue  # one-sided p-value; half of the two-sided p-value
```

Note that the beta distribution is the exact null distribution for samples of any size under this null hypothesis. We can check this by computing a Monte Carlo null distribution: explicitly drawing samples from independent normal distributions and computing Pearson's statistic for each pair.

```{code-cell} ipython3
rng = np.random.default_rng(332520619051409741187796892627751113442)

def statistic(x, y, axis):
    return stats.pearsonr(x, y, axis=axis).statistic  # ignore pvalue
    
ref = stats.monte_carlo_test((x, y), rvs=(rng.standard_normal, rng.standard_normal),
                             statistic=statistic, alternative='greater', n_resamples=9999)

fig, ax = plt.subplots(figsize=(8, 5))
plot(ax)
ax.hist(ref.null_distribution, np.linspace(-1, 1, 26), density=True)
ax.legend(['exact null distribution (independent, normally-distributed observations)',
           f'Monte Carlo null distribution \n({len(ref.null_distribution)} permutations)'])
plt.show()
```

This is often a reasonable null hypothesis to test, but in other cases, it may be more appropriate to perform a permutation test: Under the null hypothesis that total collagen and free proline are independent (but not necessarily normally distributed), each of the free proline measurements are equally likely to be observed with any of the total collagen measurements. Therefore, we can form an exact null distribution by calculating the statistic under each possible pairing of elements between `x` and `y`. This is the null distribution used when we provide `pearsonr` with `method=stats.PermutationMethod()`.

```{code-cell} ipython3
res = stats.pearsonr(x, y, alternative='greater', method=stats.PermutationMethod())
res.pvalue
```

## References

[^1]: Kershenobich, D., Fierro, F. J., & Rojkind, M. (1970). The relationship
between the free pool of proline and collagen content in human liver cirrhosis.
The Journal of Clinical Investigation, 49(12), 2246-2249.

[^2]: Hollander, M., Wolfe, D. A., & Chicken, E. (2013). Nonparametric
statistical methods. John Wiley & Sons.

[^3]: Phipson, B. and Smyth, G. K. (2010). "Permutation P-values Should Never Be
Zero: Calculating Exact P-values When Permutations Are Randomly Drawn."
Statistical Applications in Genetics and Molecular Biology 9.1.
