import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import multiscale_graphcorr


def compute_mgc(x, y):
    """Compute MGC"""
    stat, pvalue, mgc_dict = multiscale_graphcorr(x, y)
    print("MGC test statistic: ", round(stat, 1))
    print("P-value: ", round(pvalue, 1))
    return mgc_dict


def mgc_plot(x, y, sim_name):
    """Plot sim and MGC-plot"""
    # run MGC
    mgc_dict = compute_mgc(x, y)

    # local correlation map
    plt.figure(figsize=(8, 8))
    mgc_map = mgc_dict["mgc_map"]

    ax = plt.gca()

    # draw heatmap
    plt.title("Local Correlation Map", fontsize=17)
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
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 60)
    plt.show()


x = np.array([-0.69198779, 0.18236784, -0.55349325, -0.29817661, -0.18634447,
              -0.87385644, -0.01675855, 0.54961639, 0.94876906, 0.84360142,
              -0.31630326, 0.26506057, -0.46969755, 0.17495724, 0.55487294,
              0.7197614, 0.6877137, -0.71202396, -0.8485548, -0.15804884,
              0.35631348, 0.69911379, 0.62262905, -0.35320649, -0.22129935,
              -0.40791716, -0.794978, -0.69907243, -0.05848131, -0.51222227,
              0.44514037, -0.23736801, 0.46391126, -0.49347417, 0.38685186,
              0.37943616, -0.60447922, -0.38036988, -0.53064025, 0.53489706,
              0.5080103, 0.56439519, 0.0677455, -0.49572507, 0.37072501,
              -0.29463516, -0.70315084, -0.67595525, 0.70258537, -0.58428249,
              0.14102167, -0.8983592, -0.78810629, -0.38004984, -0.23656869,
              0.65481118, -0.30425405, 0.11929743, -0.49289453, -0.67434041])


y = np.array([-0.6908471, 0.16075167, -0.67384637, -0.34711818, -0.26276349,
              -0.9370918, -0.06475708, 0.54963171, 0.92910889, 0.66925604,
              -0.42844048, 0.27750248, -0.38895249, 0.24754639, 0.69432075,
              0.8844475, 0.67773753, -0.72245891, -0.90505028, -0.15128523,
              0.6783898, 0.58426834, 0.52803583, -0.58496232, -0.08003326,
              -0.54283684, -0.72224492, -0.84101252, -0.13234972, -0.46774943,
              0.37698692, -0.20006398, 0.53392985, -0.56931348, 0.43570656,
              0.3181415, -0.87795924, -0.34648145, -0.6360669, 0.52274641,
              0.38759297, 0.74337648, 0.11119666, -0.68673302, 0.27804636,
              -0.38563257, -0.6312467, -0.6512263, 0.73456059, -0.58419793,
              0.18944069, -0.84614839, -0.81814372, -0.37322543, -0.27000976,
              0.65305384, -0.35459847, 0.30447249, -0.55376159, -0.82961632])


mgc_plot(x, y, "Linear")
