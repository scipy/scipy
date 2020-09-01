import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker


orders = [2, 3]
fig, axes = plt.subplots(1, len(orders), figsize=(11, 5))
n_cells = 7  # grid will be size (n_cells, n_cells)

# desired interpolation coordinate (xi, yi)
xi, yi = 3.3, 3.7


def get_start(cc, order):
    if order % 1 == 0:
        start = math.floor(cc) - order // 2
    else:
        start = math.floor(cc + 0.5) - order // 2
    return start


for ax, order in zip(axes, orders):
    # draw open circles at the locations of pixel centers
    for n in range(n_cells):
        ax.plot(np.arange(n_cells), -np.full(n_cells, n), 'ko',
                fillstyle='none')

    # draw pixel borders
    for n in range(n_cells + 1):
        ax.plot([n - 0.5, n - 0.5], [0.5, -n_cells + .5], 'k-')
        ax.plot([-0.5, n_cells - .5], [-n + 0.5, -n + 0.5], 'k-')

    # plot an example coordinate location to interpolate
    ax.plot([xi], [-yi], 'rx')

    # plot filled circles for the points that will be involved in the
    # interpolation
    startx = get_start(xi, order)
    starty = get_start(yi, order)
    xc = np.tile(np.arange(startx, startx + order + 1)[:, np.newaxis],
                 (1, order + 1)).ravel()
    yc = np.tile(np.arange(starty, starty + order + 1)[np.newaxis, :],
                 (order + 1, 1)).ravel()
    ax.plot(xc, -yc, 'ko')
    ax.set_title("Interpolation (order = {})".format(order),
                 fontdict=dict(size=16, weight='bold'))

    # set limits and ticks for 0, 0 voxel at upper left
    ax.axis('square')
    ax.set_xticks(np.arange(n_cells + 1))
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    yticks = ticker.FixedLocator(-np.arange(n_cells, -1, -1))
    ax.yaxis.set_major_locator(yticks)
    yticklabels = ticker.FixedFormatter(np.arange(n_cells, -1, -1))
    ax.yaxis.set_major_formatter(yticklabels)
    ax.set_ylim([-n_cells + 0.5, 0.5])
    ax.set_xlim([-0.5, n_cells - 0.5])

plt.tight_layout()
plt.plot()
