"""MC vs QMC in terms of space filling."""
from scipy.stats import qmc
import numpy as np

import matplotlib.pyplot as plt


rng = np.random.default_rng()

n_sample = 256
dim = 2

sample = {}

# MC
sample['MC'] = rng.random((n_sample, dim))

# Sobol'
engine = qmc.Sobol(d=dim, seed=rng)
sample["Sobol'"] = engine.random(n_sample)

fig, axs = plt.subplots(1, 2, figsize=(8, 4))

for i, kind in enumerate(sample):
    axs[i].scatter(sample[kind][:, 0], sample[kind][:, 1])

    axs[i].set_aspect('equal')
    axs[i].set_xlabel(r'$x_1$')
    axs[i].set_ylabel(r'$x_2$')
    axs[i].set_title(f'{kind}—$C^2 = ${qmc.discrepancy(sample[kind]):.2}')

plt.tight_layout()
plt.show()
