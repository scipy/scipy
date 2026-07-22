import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter

grids = 2
boxs = 5

voxelarray = np.zeros((boxs * grids, boxs * grids, boxs * grids))

i = 1
for xi in range(0, 2):
    for yi in range(0, 2):
        for zi in range(0, 2):
            voxelarray[
                xi * boxs: xi * boxs + boxs,
                yi * boxs: yi * boxs + boxs,
                zi * boxs: zi * boxs + boxs,
            ] = i
            i += 1

voxelarray = np.uint8(voxelarray * 255 / 8)

cmap = plt.get_cmap("YlGnBu")


def plot_voxels(varray, ax, title):
    colors = cmap(varray)
    ax.view_init(30, 200)
    ax.axis("off")
    ax.voxels(varray, facecolors=colors, edgecolor="#000000", linewidth=0.1)
    ax.set_title(title, fontsize=30)


fig = plt.figure(figsize=(16, 9))
ax1 = fig.add_subplot(1, 3, 1, projection="3d")
ax2 = fig.add_subplot(1, 3, 2, projection="3d")
ax3 = fig.add_subplot(1, 3, 3, projection="3d")

plot_voxels(voxelarray, ax1, title="Original")
voxelarray2 = gaussian_filter(voxelarray, sigma=1)
plot_voxels(voxelarray2, ax2, title="gaussian_filter \n sigma=1")
voxelarray3 = gaussian_filter(voxelarray, sigma=3)
plot_voxels(voxelarray3, ax3, title="gaussian_filter \n sigma=3")

plt.tight_layout()
plt.show()
