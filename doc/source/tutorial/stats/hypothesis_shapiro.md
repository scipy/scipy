---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

+++ {"tags": ["jupyterlite_sphinx_strip"]}

```{eval-rst}
.. notebooklite:: hypothesis_shapiro.md
   :new_tab: True
```

(hypothesis_shapiro)=

+++

# Shapiro-Wilk test for normality

Suppose we wish to infer from measurements whether the weights of adult human
males in a medical study are not normally distributed [^1]. The weights (lbs)
are recorded in the array `x` below.

```{code-cell}
import numpy as np
x = np.array([148, 154, 158, 160, 161, 162, 166, 170, 182, 195, 236])
```

The normality test {func}`scipy.stats.shapiro` of [^1] and [^2] begins by
computing a statistic based on the relationship between the observations and the
expected order statistics of a normal distribution.

```{code-cell}
from scipy import stats
res = stats.shapiro(x)
res.statistic
```

The value of this statistic tends to be high (close to 1) for samples drawn from
a normal distribution.

The test is performed by comparing the observed value of the statistic against
the null distribution: the distribution of statistic values formed under the
null hypothesis that the weights were drawn from a normal distribution. For this
normality test, the null distribution is not easy to calculate exactly, so it is
usually approximated by Monte Carlo methods, that is, drawing many samples of
the same size as `x` from a normal distribution and computing the values of the
statistic for each.

```{code-cell}
def statistic(x):
    # Get only the `shapiro` statistic; ignore its p-value
    return stats.shapiro(x).statistic
ref = stats.monte_carlo_test(x, stats.norm.rvs, statistic,
                             alternative='less')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 5))
bins = np.linspace(0.65, 1, 50)

def plot(ax):  # we'll reuse this
    ax.hist(ref.null_distribution, density=True, bins=bins)
    ax.set_title("Shapiro-Wilk Test Null Distribution \n"
                 "(Monte Carlo Approximation, 11 Observations)")
    ax.set_xlabel("statistic")
    ax.set_ylabel("probability density")

plot(ax)
plt.show()
```

The comparison is quantified by the p-value: the proportion of values in the
null distribution less than or equal to the observed value of the statistic.

```{code-cell}
fig, ax = plt.subplots(figsize=(8, 5))
plot(ax)
annotation = (f'p-value={res.pvalue:.6f}\n(highlighted area)')
props = dict(facecolor='black', width=1, headwidth=5, headlength=8)
_ = ax.annotate(annotation, (0.75, 0.1), (0.68, 0.7), arrowprops=props)
i_extreme = np.where(bins <= res.statistic)[0]
for i in i_extreme:
    ax.patches[i].set_color('C1')
plt.xlim(0.65, 0.9)
plt.ylim(0, 4)
plt.show()
```

```{code-cell}
res.pvalue
```

If the p-value is "small" - that is, if there is a low probability of sampling
data from a normally distributed population that produces such an extreme value
of the statistic - this may be taken as evidence against the null hypothesis in
favor of the alternative: the weights were not drawn from a normal distribution.
Note that:

- The inverse is not true; that is, the test is not used to provide
  evidence *for* the null hypothesis.
- The threshold for values that will be considered "small" is a choice that
  should be made before the data is analyzed [^3] with consideration of the
  risks of both false positives (incorrectly rejecting the null hypothesis)
  and false negatives (failure to reject a false null hypothesis).

## References

[^1]: Shapiro, S. S. & Wilk, M.B (1965). An analysis of variance test for
normality (complete samples). Biometrika, Vol. 52, pp. 591-611.
{doi}`10.2307/2333709`
[^2]: NIST/SEMATECH e-Handbook of Statistical Methods, Section 7.2.1.3,
Anderson-Darling and Shapiro-Wilk tests,
https://www.itl.nist.gov/div898/handbook/prc/section2/prc213.htm
{doi}`10.18434/M32189`
[^3]: Phipson B., and Smyth, G. K. (2010). Permutation P-values Should Never Be
Zero: Calculating Exact P-values When Permutations Are Randomly Drawn.
Statistical Applications in Genetics and Molecular Biology, Vol.9.
{doi}`10.2202/1544-6115.1585`
