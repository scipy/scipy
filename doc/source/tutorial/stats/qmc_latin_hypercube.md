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
orphan: true
---

```{eval-rst}
.. jupyterlite:: ../../_contents/qmc_latin_hypercube.ipynb
   :new_tab: True
```

(qmc_latin_hypercube)=
# Latin Hypercube sampling

In [^1], a Latin Hypercube sampling (LHS) strategy was used to sample a
parameter space to study the importance of each parameter of an epidemic model.
Such analysis is also called a sensitivity analysis.

Since the dimensionality of the problem is high (6), it is computationally
expensive to cover the space. When numerical experiments are costly, QMC enables
analysis that may not be possible if using a grid.

The six parameters of the model represented the probability of illness, the
probability of withdrawal, and four contact probabilities. The authors assumed
uniform distributions for all parameters and generated 50 samples.

Using {class}`scipy.stats.qmc.LatinHypercube` to replicate the protocol, the
first step is to create a sample in the unit hypercube:

```{code-cell}
from scipy.stats import qmc
sampler = qmc.LatinHypercube(d=6)
sample = sampler.random(n=50)
```

Then the sample can be scaled to the appropriate bounds:

```{code-cell}
l_bounds = [0.000125, 0.01, 0.0025, 0.05, 0.47, 0.7]
u_bounds = [0.000375, 0.03, 0.0075, 0.15, 0.87, 0.9]
sample_scaled = qmc.scale(sample, l_bounds, u_bounds)
```

Such a sample was used to run the model 50 times, and a polynomial response
surface was constructed. This allowed the authors to study the relative
importance of each parameter across the range of possibilities of every other
parameter.

In this computer experiment, they showed a 14-fold reduction in the number of
samples required to maintain an error below 2% on their response surface when
compared to a grid sampling.

## References

[^1]: Seaholm, Susan K. et al. (1988). Latin hypercube sampling and the
sensitivity analysis of a Monte Carlo epidemic model. Int J Biomed Comput,
23(1-2), 97-112. {doi}`10.1016/0020-7101(88)90067-0`
