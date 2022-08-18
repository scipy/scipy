"""Integration convergence comparison: MC vs Sobol'.

The function is a synthetic example specifically designed
to verify the correctness of the implementation [2]_.

References
----------

.. [1] I. M. Sobol. The distribution of points in a cube and the accurate
   evaluation of integrals. Zh. Vychisl. Mat. i Mat. Phys., 7:784-802,
   1967.
.. [2] Art B. Owen. On dropping the first Sobol' point. arXiv 2008.08051,
   2020.

"""
from collections import namedtuple

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import qmc

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
sample_sobol_rmse = []
rng = np.random.default_rng()

for ns in ns_gen:
    # Monte Carlo
    sampler_mc = lambda x: rng.random((x, case.dim))
    conv_res = conv_method(sampler_mc, case.func, ns, n_conv, case.ref)
    sample_mc_rmse.append(conv_res)

    # Sobol'
    engine = qmc.Sobol(d=case.dim, scramble=False)
    conv_res = conv_method(engine.random, case.func, ns, 1, case.ref)
    sample_sobol_rmse.append(conv_res)

sample_mc_rmse = np.array(sample_mc_rmse)
sample_sobol_rmse = np.array(sample_sobol_rmse)

# Plot
fig, ax = plt.subplots(figsize=(4, 4))
ax.set_aspect('equal')

# MC
ratio = sample_mc_rmse[0] / ns_gen[0] ** (-1 / 2)
ax.plot(ns_gen, ns_gen ** (-1 / 2) * ratio, ls='-', c='k')

ax.scatter(ns_gen, sample_mc_rmse, label="MC")

# Sobol'
ratio = sample_sobol_rmse[0] / ns_gen[0] ** (-2/2)
ax.plot(ns_gen, ns_gen ** (-2/2) * ratio, ls='-.', c='k')

ax.scatter(ns_gen, sample_sobol_rmse, label="Sobol' unscrambled")

ax.set_xlabel(r'$N_s$')
ax.set_xscale('log')
ax.set_xticks(ns_gen)
ax.set_xticklabels([fr'$2^{{{ns}}}$' for ns in np.arange(4, 13)])

ax.set_ylabel(r'$\log (\epsilon)$')
ax.set_yscale('log')

ax.legend(loc='upper right')
fig.tight_layout()
plt.show()
