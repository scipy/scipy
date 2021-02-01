from scipy.stats import qmc
from scipy.stats._qmc import check_random_state
import numpy as np

import matplotlib.pyplot as plt


seed = np.random.RandomState(12345)
rng = check_random_state(seed)

n_sample = 256
dim = 2

sample = {}

# MC
sample['MC'] = rng.random_sample((n_sample, dim))

# Sobol
engine = qmc.Sobol(d=dim, seed=seed)
sample["Sobol'"] = engine.random(n_sample)

fig, axs = plt.subplots(1, 2)

for i, kind in enumerate(sample):
    axs[i].scatter(sample[kind][:, 0], sample[kind][:, 1])

    axs[i].set_aspect('equal')
    axs[i].set_xlabel(r'$x_1$')
    axs[i].set_ylabel(r'$x_2$')
    axs[i].set_title(f'{kind}—$C^2 = ${qmc.discrepancy(sample[kind]):.2}')

plt.show()
