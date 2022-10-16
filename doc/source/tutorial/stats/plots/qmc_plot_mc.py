"""Multiple MC to show how it can be bad."""
from scipy.stats import qmc
from scipy.stats._qmc import check_random_state
import numpy as np

import matplotlib.pyplot as plt


rng = np.random.default_rng()

n_sample = 256
dim = 2

sample = {}

# MC
sample['MC 1'] = rng.random((n_sample, dim))
sample["MC 2"] = rng.random((n_sample, dim))

fig, axs = plt.subplots(1, 2, figsize=(8, 4))

for i, kind in enumerate(sample):
    axs[i].scatter(sample[kind][:, 0], sample[kind][:, 1])

    axs[i].set_aspect('equal')
    axs[i].set_xlabel(r'$x_1$')
    axs[i].set_ylabel(r'$x_2$')

plt.tight_layout()
plt.show()
