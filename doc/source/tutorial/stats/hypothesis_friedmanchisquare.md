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
.. jupyterlite:: ../../_contents/hypothesis_friedmanchisquare.ipynb
   :new_tab: True
```

(hypothesis_friedmanchisquare)=
# Friedman test for repeated samples

In [^1], the pulse rate (per minute) of a group of seven students was measured
before exercise, immediately after exercise and 5 minutes after exercise. Is
there evidence to suggest that the pulse rates on these three occasions are
similar?

We begin by formulating a null hypothesis $H_0$:

> The pulse rates are identical on these three occasions.

Let's assess the plausibility of this hypothesis with a Friedman test
({func}`scipy.stats.friedmanchisquare`.)

```{code-cell}
from scipy.stats import friedmanchisquare
before = [72, 96, 88, 92, 74, 76, 82]
immediately_after = [120, 120, 132, 120, 101, 96, 112]
five_min_after = [76, 95, 104, 96, 84, 72, 76]
res = friedmanchisquare(before, immediately_after, five_min_after)
res.statistic
```

```{code-cell}
res.pvalue
```

Using a significance level of 5%, we would reject the null hypothesis in favor
of the alternative hypothesis: "the pulse rates are different on these three
occasions".

## References

[^1]: Sprent, P. and Smeeton, N.C. (2000), "Applied Nonparametric Statistical
Methods, Third Edition". Chapter 6, Section 6.3.2.
