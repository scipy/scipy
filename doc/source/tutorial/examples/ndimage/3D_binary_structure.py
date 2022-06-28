import matplotlib.pyplot as plt
import scipy.ndimage


def plot_voxels(varray, ax, title):
    ax.view_init(20, 200)
    ax.voxels(varray, edgecolor="k")
    ax.set_title(title, fontsize=30)


fig = plt.figure(figsize=(16, 9))

for i in [1, 2, 3]:
    ax = fig.add_subplot(1, 3, i, projection="3d")
    arrray = scipy.ndimage.generate_binary_structure(3, i)
    plot_voxels(arrray, ax, title=f"rank=3 \n connectivity={i}")

plt.tight_layout()
plt.show()
