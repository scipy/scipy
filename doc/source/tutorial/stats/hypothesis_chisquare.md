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
.. notebooklite:: hypothesis_chisquare.md
   :new_tab: True
```

(hypothesis_chisquare)=

+++

# Chi-square test

The {class}`chi-square test <scipy.stats.chisquare>` tests the null hypothesis
that a given set of categorical data has the given frequencies.

In [^1], bird foraging behavior was investigated in an old-growth forest of
Oregon. In the forest, 44% of the canopy volume was Douglas fir, 24% was
ponderosa pine, 29% was grand fir, and 3% was western larch. The authors
observed the behavior of several species of birds, one of which was the
red-breasted nuthatch. They made 189 observations of this species foraging,
recording 43 ("23%") of observations in Douglas fir, 52 ("28%") in ponderosa
pine, 54 ("29%") in grand fir, and 40 ("21%") in western larch.

Using a chi-square test, we can test the null hypothesis that the proportions of
foraging events are equal to the proportions of canopy volume. The authors of
the paper considered a p-value less than 1% to be significant.

Using the above proportions of canopy volume and observed events, we can infer
expected frequencies.

```{code-cell} ipython3
import numpy as np
f_exp = np.array([44, 24, 29, 3]) / 100 * 189
```

The observed frequencies of foraging were:

```{code-cell} ipython3
f_obs = np.array([43, 52, 54, 40])
```

We can now compare the observed frequencies with the expected frequencies.

```{code-cell} ipython3
from scipy.stats import chisquare
chisquare(f_obs=f_obs, f_exp=f_exp)
```

The p-value is well below the chosen significance level. Hence, the authors
considered the difference to be significant and concluded that the relative
proportions of foraging events were not the same as the relative proportions of
tree canopy volume.

## References

[^1]: Mannan, R. William and Meslow, E. Charles (1984) "Bird populations and
vegetation characteristics in managed and old-growth forests, northeastern
Oregon." Journal of Wildlife Management 48, 1219-1238. {doi}`10.2307/3801783`
