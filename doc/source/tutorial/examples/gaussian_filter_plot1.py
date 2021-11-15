import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter

grids = 2
boxs = 5
x, y, z = np.indices((boxs * grids, boxs * grids, boxs * grids))


def make_cube(pos=(0, 0, 0)):  # make a boolean mask
    cube = (
        (pos[0] <= x)
        & (x < pos[0] + boxs)
        & (pos[1] <= y)
        & (y < pos[1] + boxs)
        & (pos[2] <= z)
        & (z < pos[2] + boxs)
    )
    return cube


voxelarray = np.zeros((boxs * grids, boxs * grids, boxs * grids))

i = 1
for xi in range(0, 2):
    for yi in range(0, 2):
        for zi in range(0, 2):
            cube = make_cube((xi * boxs, yi * boxs, zi * boxs))
            voxelarray[cube] = i
            i += 1

voxelarray = np.uint8(voxelarray * 255 / 8)

cmap = plt.get_cmap("YlGnBu")


def plot_voxels(varray, fig, subp, title):
    colors = cmap(varray)
    ax = fig.add_subplot(*subp, projection="3d")
    ax.view_init(30, 200)
    plt.axis("off")
    ax.voxels(varray, facecolors=colors, edgecolor="#000000", linewidth=0.1)
    ax.set_title(title)


plt.rcParams["figure.dpi"] = 250
fig = plt.figure(figsize=plt.figaspect(0.5))
plot_voxels(voxelarray, fig, subp=(1, 3, 1), title="Original")
voxelarray2 = gaussian_filter(voxelarray, sigma=1)
plot_voxels(voxelarray2, fig, subp=(1, 3, 2), title="gaussian_filter \n sigma=1")
voxelarray2 = gaussian_filter(voxelarray, sigma=3)
plot_voxels(voxelarray2, fig, subp=(1, 3, 3), title="gaussian_filter \n sigma=3")
