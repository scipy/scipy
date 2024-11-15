import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import multiscale_graphcorr


def mgc_plot(x, y, mgc_dict):
    """Plot sim and MGC-plot"""
    plt.figure(figsize=(8, 8))
    ax = plt.gca()

    # local correlation map
    mgc_map = mgc_dict["mgc_map"]

    # draw heatmap
    ax.set_title("Local Correlation Map", fontsize=20)
    im = ax.imshow(mgc_map, cmap='YlGnBu')

    # colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("", rotation=-90, va="bottom")
    ax.invert_yaxis()

    # Turn spines off and create white grid.
    for _, spine in ax.spines.items():
        spine.set_visible(False)

    # optimal scale
    opt_scale = mgc_dict["opt_scale"]
    ax.scatter(opt_scale[0], opt_scale[1],
               marker='X', s=200, color='red')

    # other formatting
    ax.tick_params(bottom="off", left="off")
    ax.set_xlabel('#Neighbors for X', fontsize=15)
    ax.set_ylabel('#Neighbors for Y', fontsize=15)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)


rng = np.random.default_rng()
unif = np.array(rng.uniform(0, 5, size=100))
x = unif * np.cos(np.pi * unif)
y = unif * np.sin(np.pi * unif) + 0.4 * rng.random(x.size)

_, _, mgc_dict = multiscale_graphcorr(x, y, random_state=rng)
mgc_plot(x, y, mgc_dict)
