import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from scipy.ndimage import gaussian_filter

grids = 2
boxs = 5
x, y, z = np.indices((boxs * grids, boxs * grids, boxs * grids))

def cube(pos=(0, 0, 0)):
    o1 = pos[0]
    o2 = pos[1]
    o3 = pos[2]
    cube = (
        (o1 <= x)
        & (x < o1 + boxs)
        & (o2 <= y)
        & (y < o2 + boxs)
        & (o3 <= z)
        & (z < o3 + boxs)
    )
    return cube


voxelarray = np.zeros((boxs * grids, boxs * grids, boxs * grids))

voxelarray[cube((0, 0, 0))] = 1 * 255 / 8
voxelarray[cube((boxs, 0, 0))] = 2 * 255 / 8
voxelarray[cube((0, boxs, 0))] = 3 * 255 / 8
voxelarray[cube((boxs, boxs, 0))] = 4 * 255 / 8
voxelarray[cube((0, 0, boxs))] = 5 * 255 / 8
voxelarray[cube((boxs, 0, boxs))] = 6 * 255 / 8
voxelarray[cube((0, boxs, boxs))] = 7 * 255 / 8
voxelarray[cube((boxs, boxs, boxs))] = 8 * 255 / 8
voxelarray = np.uint8(voxelarray)

cmap = plt.get_cmap("YlGnBu")
cols256 = []
for i in range(cmap.N):
    rgb = cmap(i)[:3]
    cols256.append(mpl.colors.rgb2hex(rgb))

def get_color(x):
    return cols256[x]

get_color = np.vectorize(get_color)

colors = get_color(voxelarray)

def plot_voxels(varray, fig, subp, title):
    colors = get_color(varray)
    ax = fig.add_subplot(*subp, projection="3d")
    plt.axis("off")
    ax.voxels(varray, facecolors=colors, edgecolor="#000000", linewidth=0.1)
    ax.set_title(title)

fig = plt.figure(figsize=plt.figaspect(0.5))
plot_voxels(voxelarray, fig, subp=(1, 3, 1), title="Original")
voxelarray2 = gaussian_filter(voxelarray, sigma=1)
plot_voxels(voxelarray2, fig, subp=(1, 3, 2), title="gaussian_filter \n sigma=1")
voxelarray2 = gaussian_filter(voxelarray, sigma=3)
plot_voxels(voxelarray2, fig, subp=(1, 3, 3), title="gaussian_filter \n sigma=3")