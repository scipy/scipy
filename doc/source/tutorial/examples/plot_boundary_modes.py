import numpy as np
import matplotlib.pyplot as plt

from scipy import ndimage

img = np.array([-2, -1, 0, 1, 2], float)
x = np.linspace(-2, 6, num=1000)

modes = ['constant', 'grid-constant', 'nearest', 'reflect', 'mirror', 'wrap',
         'grid-wrap']
fig, axes = plt.subplots(len(modes), 3, figsize=(11, 8), sharex=True,
                         sharey=True)


for mode, (ax0, ax1, ax2) in zip(modes, axes):

    y = ndimage.map_coordinates(img, [x], order=0, mode=mode)
    ax0.scatter(np.arange(img.size), img)
    ax0.plot(x, y, '-')
    ax0.set_title(f'mode={mode}, order=0')

    y2 = ndimage.map_coordinates(img, [x], order=1, mode=mode)
    ax1.scatter(np.arange(img.size), img)
    ax1.plot(x, y2, '-')
    ax1.set_title(f'mode={mode}, order=1')

    y3 = ndimage.map_coordinates(img, [x], order=3, mode=mode)
    ax2.scatter(np.arange(img.size), img)
    ax2.plot(x, y3, '-')
    ax2.set_title(f'mode={mode}, order=3')

    sz = len(img)
    for ax in (ax0, ax1, ax2):
        if mode in ['grid-wrap', 'reflect']:
            ax.plot([-0.5, -0.5], [-2.5, 2.5], 'k--')
            ax.plot([sz - 0.5, sz - 0.5], [-2.5, 2.5], 'k--')
        elif mode in ['wrap', 'mirror']:
            ax.plot([0, 0], [-2.5, 2.5], 'k--')
            ax.plot([sz - 1, sz - 1], [-2.5, 2.5], 'k--')

    if mode != 'constant':
        for xx in range(int(x[0]), int(x[-1] + 1)):
            if (xx < 0) or (xx > img.size - 1):
                idx = np.argmin(np.abs(x - xx))

                for y_vals, ax in zip((y, y2, y3), (ax0, ax1, ax2)):
                    ax.scatter(
                        [x[idx]], [y_vals[idx]], facecolors='none',
                        edgecolor='#0343df', marker='o'
                    )

plt.tight_layout()
plt.show()
