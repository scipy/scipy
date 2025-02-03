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
.. notebooklite:: hypothesis_odds_ratio.md
   :new_tab: True
```

(hypothesis_odds_ratio)=

+++

# Odds ratio for a contingency table

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

|                 |    Aspirin    | Control/Placebo |
|-----------------|---------------|-----------------|
| Ischemic stroke |       176     |        230      |
| No stroke       |     21035     |      21018      |

The question we ask is "Is there evidence that the aspirin reduces the risk of
ischemic stroke?"

First, compute the {class}`odds ratio <scipy.stats.contingency.odds_ratio>` for
this contingency table:

```{code-cell}
from scipy.stats.contingency import odds_ratio
res = odds_ratio([[176, 230], [21035, 21018]])
res.statistic
```

For this sample, the odds of getting an ischemic stroke for those who have been
taking aspirin are 0.76 times that of those who have received the placebo.

To make statistical inferences about the population under study, we can compute
the 95% confidence interval for the odds ratio:

```{code-cell}
res.confidence_interval(confidence_level=0.95)
```

The 95% confidence interval for the conditional odds ratio is approximately
(0.62, 0.94).

The fact that the entire 95% confidence interval falls below 1 supports the
authors' conclusion that the aspirin was associated with a statistically
significant reduction in ischemic stroke.

## References

[^1]: Berger, Jeffrey S. et al. (2006). Aspirin for the Primary Prevention of
Cardiovascular Events in Women and Men: A Sex-Specific Meta-analysis of
Randomized Controlled Trials. JAMA, 295(3):306-313.
{doi}`10.1001/jama.295.3.306`
