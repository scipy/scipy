import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage

# code for ball taken from
# https://github.com/scikit-image/scikit-image/blob/main/skimage/morphology/footprints.py#L225-L252
# and therefore same as `from skimage.morphology import ball`


def ball(radius, dtype=np.uint8):
    n = 2 * radius + 1
    Z, Y, X = np.mgrid[
        -radius: radius: n * 1j,
        -radius: radius: n * 1j,
        -radius: radius: n * 1j
    ]
    s = X ** 2 + Y ** 2 + Z ** 2
    return np.array(s <= radius * radius, dtype=dtype)


def plot_voxels(varray, ax, title):
    ax.view_init(20, 200)
    ax.voxels(varray, edgecolor="k")
    ax.set_title(title, fontsize=30)


voxelarray = np.full((11, 11, 11), 0)
voxelarray[5, 3, 5] = 1
voxelarray[5, 7, 5] = 1
img_morphed = scipy.ndimage.binary_dilation(voxelarray, ball(3))
img_morphed2 = scipy.ndimage.binary_erosion(img_morphed, ball(2))

fig = plt.figure(figsize=(16, 9))
ax1 = fig.add_subplot(1, 3, 1, projection="3d")
ax2 = fig.add_subplot(1, 3, 2, projection="3d")
ax3 = fig.add_subplot(1, 3, 3, projection="3d")

plot_voxels(voxelarray, ax1, title="a) Original")
plot_voxels(img_morphed, ax2, title="b) binary_dilation \nwith ball, radius 3")
plot_voxels(img_morphed2, ax3,
            title="c) binary_erosion of b \nwith ball, radius 2")

plt.tight_layout()
plt.show()
