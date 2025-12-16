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
.. notebooklite:: hypothesis_dunnett.md
   :new_tab: True
```

(hypothesis_dunnett)=

+++

# Dunnett's test

Dunnett's test compares the means of multiple experimental groups against a
single control group. In [^1], the influence of drugs on blood count
measurements on three groups of animals is investigated.

The following table summarizes the results of the experiment in which two groups
received different drugs, and one group acted as a control. Blood counts (in
millions of cells per cubic millimeter) were recorded:

```{code-cell}
import numpy as np
control = np.array([7.40, 8.50, 7.20, 8.24, 9.84, 8.32])
drug_a = np.array([9.76, 8.80, 7.68, 9.36])
drug_b = np.array([12.80, 9.68, 12.16, 9.20, 10.55])
```

We would like to see if the means between any of the groups are
significantly different. First, visually examine a box and whisker plot.

```{code-cell}
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 1)
ax.boxplot([control, drug_a, drug_b])
ax.set_xticklabels(["Control", "Drug A", "Drug B"])
ax.set_ylabel("mean")
plt.show()
```

Note the overlapping interquartile ranges of the drug A group and control group
and the apparent separation between the drug B group and control group.

Next, we will use {func}`Dunnett's test <scipy.stats.dunnett>` to assess whether
the difference between group means is significant while controlling the
family-wise error rate: the probability of making any false discoveries.

Let the null hypothesis be that the experimental groups have the same mean as
the control and the alternative be that an experimental group does not have the
same mean as the control. We will consider a 5% family-wise error rate to be
acceptable, and therefore we choose 0.05 as the threshold for significance.

```{code-cell}
from scipy.stats import dunnett
res = dunnett(drug_a, drug_b, control=control)
res.pvalue
```

The p-value corresponding with the comparison between group A and control
exceeds 0.05, so we do not reject the null hypothesis for that comparison.
However, the p-value corresponding with the comparison between group B and
control is less than 0.05, so we consider the experimental results to be
evidence against the null hypothesis in favor of the alternative: group B has a
different mean than the control group.

## References

[^1]: Dunnett, Charles W. (1955) "A Multiple Comparison Procedure for Comparing
Several Treatments with a Control." Journal of the American Statistical
Association, 50:272, 1096-1121. {doi}`10.1080/01621459.1955.10501294`
