"""Calculate the discrepancy of 2 designs and compare them."""
import numpy as np
from scipy.stats import qmc

import matplotlib.pyplot as plt

space_1 = np.array([[1, 3], [2, 6], [3, 2], [4, 5], [5, 1], [6, 4]])
space_2 = np.array([[1, 5], [2, 4], [3, 3], [4, 2], [5, 1], [6, 6]])

l_bounds = [0.5, 0.5]
u_bounds = [6.5, 6.5]
space_1 = qmc.scale(space_1, l_bounds, u_bounds, reverse=True)
space_2 = qmc.scale(space_2, l_bounds, u_bounds, reverse=True)

sample = {'space_1': space_1, 'space_2': space_2}

fig, axs = plt.subplots(1, 2, figsize=(8, 4))

for i, kind in enumerate(sample):
    axs[i].scatter(sample[kind][:, 0], sample[kind][:, 1])

    axs[i].set_aspect('equal')
    axs[i].set_xlabel(r'$x_1$')
    axs[i].set_ylabel(r'$x_2$')
    axs[i].set_title(f'{kind}—$C^2 = ${qmc.discrepancy(sample[kind]):.5}')

plt.tight_layout()
plt.show()
