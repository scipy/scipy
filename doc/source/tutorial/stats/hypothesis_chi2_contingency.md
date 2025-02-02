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
.. notebooklite:: hypothesis_chi2_contingency.md
   :new_tab: True
```

(hypothesis_chi2_contingency)=

+++

# Chi-square test of independence of variables in a contingency table

In [^1], the use of aspirin to prevent cardiovascular events in women and men
was investigated. The study notably concluded:

> ...aspirin therapy reduced the risk of a composite of
> cardiovascular events due to its effect on reducing the risk of
> ischemic stroke in women [...]

The article lists studies of various cardiovascular events. Let's focus on the
ischemic stoke in women.

The following table summarizes the results of the experiment in which
participants took aspirin or a placebo on a regular basis for several years.
Cases of ischemic stroke were recorded::

|                 | Aspirin | Control/Placebo |
|-----------------|---------|-----------------|
| Ischemic stroke |    176  |        230      |
| No stroke       |  21035  |      21018      |

Is there evidence that the aspirin reduces the risk of ischemic stroke? We begin
by formulating a null hypothesis $H_0$:

> The effect of aspirin is equivalent to that of placebo.

Let's assess the plausibility of this hypothesis with a
{class}`chi-square test <scipy.stats.contingency.chi2_contingency>` with the
observed [contingency table](https://en.wikipedia.org/wiki/Contingency_table)
as our input.

```{code-cell}
import numpy as np
from scipy.stats import chi2_contingency
table = np.array([[176, 230], [21035, 21018]])
res = chi2_contingency(table)
res.statistic
```

```{code-cell}
res.pvalue
```

Using a significance level of 5%, we would reject the null hypothesis in favor
of the alternative hypothesis: "the effect of aspirin is not equivalent to the
effect of placebo". Because {class}`scipy.stats.contingency.chi2_contingency`
performs a two-sided test, the alternative hypothesis does not indicate the
direction of the effect. We can use {class}`scipy.stats.contingency.odds_ratio`
to support the conclusion that aspirin *reduces* the risk of ischemic stroke.

## References

[^1]: Berger, Jeffrey S. et al. (2006) "Aspirin for the Primary Prevention of
      Cardiovascular Events in Women and Men: A Sex-Specific Meta-analysis of
      Randomized Controlled Trials." JAMA, 295(3):306-313.
      {doi}`10.1001/jama.295.3.306`
