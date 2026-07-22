"""Visualize the curse-of-dimensionality.

It presents a saturated design in 1, 2 and 3 dimensions for a
given discretization.
"""
import matplotlib.pyplot as plt
import numpy as np

disc = 10

x = np.linspace(0, 1, disc)
y = np.linspace(0, 1, disc)
z = np.linspace(0, 1, disc)

xx, yy, zz = np.meshgrid(x, y, z)

fig = plt.figure(figsize=(12, 4))
ax = fig.add_subplot(131)
ax.set_aspect('equal')
ax.scatter(xx, yy * 0)
ax.set_xlabel(r'$x_1$')
ax.get_yaxis().set_visible(False)

ax = fig.add_subplot(132)
ax.set_aspect('equal')
ax.scatter(xx, yy)
ax.set_xlabel(r'$x_1$')
ax.set_ylabel(r'$x_2$')

ax = fig.add_subplot(133, projection='3d')
ax.scatter(xx, yy, zz)
ax.set_xlabel(r'$x_1$')
ax.set_ylabel(r'$x_2$')
ax.set_zlabel(r'$x_3$')

plt.tight_layout(pad=2)
plt.show()
