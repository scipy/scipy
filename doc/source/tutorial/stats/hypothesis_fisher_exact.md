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

```{eval-rst}
.. jupyterlite:: ../../_contents/hypothesis_fisher_exact.ipynb
   :new_tab: True
```

(hypothesis_fisher_exact)=
# Fisher's exact test

In [^1], the effective dose of acetazolamide for the prophylaxis of acute
mountain sickness was investigated. The study notably concluded:

> Acetazolamide 250 mg, 500 mg, and 750 mg daily were all efficacious for
> preventing acute mountain sickness. Acetazolamide 250 mg was the lowest
> effective dose with available evidence for this indication.

The following table summarizes the results of the experiment in which
some participants took a daily dose of acetazolamide 250 mg while others
took a placebo.

Cases of acute mountain sickness were recorded::

|                         | Acetazolamide | Control/Placebo V |
|-------------------------|---------------|-------------------|
| Acute mountain sickness |       7       |        17         |
| No                      |      15       |         5         |

Is there evidence that the acetazolamide 250 mg reduces the risk of
acute mountain sickness?

We begin by formulating a null hypothesis $H_0$:

> The odds of experiencing acute mountain sickness are the same with
> the acetazolamide treatment as they are with placebo.

Let's assess the plausibility of this hypothesis with
{class}`Fisher's test <scipy.stats.fisher_exact>`.

```{code-cell}
from scipy.stats import fisher_exact
res = fisher_exact([[7, 17], [15, 5]], alternative='less')
res.statistic
```

```{code-cell}
res.pvalue
```

Using a significance level of 5%, we would reject the null hypothesis in
favor of the alternative hypothesis: "The odds of experiencing acute
mountain sickness with acetazolamide treatment are less than the odds of
experiencing acute mountain sickness with placebo."

## Note

Because the null distribution of Fisher's exact test is formed under
the assumption that both row and column sums are fixed, the result of
the test are conservative when applied to an experiment in which the
row sums are not fixed.

In this case, the column sums are fixed; there are 22 subjects in each
group. But the number of cases of acute mountain sickness is not
(and cannot be) fixed before conducting the experiment. It is a
consequence.

{class}`Boschloo's test <scipy.stats.boschloo_exact>` does not depend on the
assumption that the row sums are fixed, and consequently, it provides a more
powerful test in this situation.

```{code-cell}
from scipy.stats import boschloo_exact
res = boschloo_exact([[7, 17], [15, 5]], alternative='less')
res.statistic
```

```{code-cell}
res.pvalue
```

We verify that the p-value is less than with ``fisher_exact``.

## References

[^1]: Low, Emma V. et al. (2012) "Identifying the lowest effective dose of
acetazolamide for the prophylaxis of acute mountain sickness: systematic review
and meta-analysis." BMJ, 345. {doi}`10.1136/bmj.e6779`
