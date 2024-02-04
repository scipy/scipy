"""Integration convergence.

The function is a synthetic example specifically designed
to verify the correctness of the implementation [1]_.

References
----------

.. [1] Art B. Owen. On dropping the first Sobol' point. arXiv 2008.08051,
   2020.

"""
from collections import namedtuple

import numpy as np
import matplotlib.pyplot as plt


n_conv = 99
ns_gen = 2 ** np.arange(4, 13)  # 13


def art_2(sample):
    # dim 3, true value 5/3 + 5*(5 - 1)/4
    return np.sum(sample, axis=1) ** 2


functions = namedtuple('functions', ['name', 'func', 'dim', 'ref'])
case = functions('Art 2', art_2, 5, 5 / 3 + 5 * (5 - 1) / 4)


def conv_method(sampler, func, n_samples, n_conv, ref):
    samples = [sampler(n_samples) for _ in range(n_conv)]
    samples = np.array(samples)

    evals = [np.sum(func(sample)) / n_samples for sample in samples]
    squared_errors = (ref - np.array(evals)) ** 2
    rmse = (np.sum(squared_errors) / n_conv) ** 0.5

    return rmse


# Analysis
sample_mc_rmse = []
rng = np.random.default_rng()

def sampler_mc(x):
    return rng.random((x, case.dim))

for ns in ns_gen:
    # Monte Carlo
    conv_res = conv_method(sampler_mc, case.func, ns, n_conv, case.ref)
    sample_mc_rmse.append(conv_res)

sample_mc_rmse = np.array(sample_mc_rmse)

# Plot
fig, ax = plt.subplots(figsize=(5, 3))
ax.set_aspect('equal')

ratio = sample_mc_rmse[0] / ns_gen[0] ** (-1 / 2)
ax.plot(ns_gen, ns_gen ** (-1 / 2) * ratio, ls='-', c='k')

ax.scatter(ns_gen, sample_mc_rmse)

ax.set_xlabel(r'$N_s$')
ax.set_xscale('log')
ax.set_xticks(ns_gen)
ax.set_xticklabels([fr'$2^{{{ns}}}$' for ns in np.arange(4, 13)])

ax.set_ylabel(r'$\log (\epsilon)$')
ax.set_yscale('log')

fig.tight_layout()
plt.show()
